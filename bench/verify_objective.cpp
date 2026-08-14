// Does the vectorised fast path + fallback compute the SAME objective score as
// the robust scalar solver? The MCMC only ever sees this score, so if the two
// agree to rounding, the sampler cannot tell them apart.
#include <cstdio>
#include <cmath>
#include <random>
#include <vector>
#include <algorithm>
using namespace std;

static long long n_fallback = 0, n_line_evals = 0;

inline void ev(double x, double a, double b, double c, double d, double &f, double &fp) {
    f = d + x * (c + x * (b + x * a));
    fp = c + x * (2.0 * b + x * 3.0 * a);
}

// ---- Robust reference solver (as shipped in linefit.cpp) -------------------
inline double inverse_unit_poly(double y, double A, double B, double C, double D) {
    const double Ds = D - y;
    const double ftol = 1e-14 * (fabs(y) + 1.0);
    double u = -Ds / C;
    if (B != 0.0) {
        const double disc = C * C - 4.0 * B * Ds;
        if (disc >= 0.0) {
            const double q = -0.5 * (C + copysign(sqrt(disc), C));
            const double r1 = q / B, r2 = (q != 0.0) ? (Ds / q) : r1;
            u = (fabs(r1 - u) < fabs(r2 - u)) ? r1 : r2;
        }
    }
    if (!(u > -10.0 && u < 10.0)) u = 0.5;
    double best_u = u, best_f = HUGE_VAL;
    for (int it = 0; it < 12; ++it) {
        double f, fp;
        ev(u, A, B, C, Ds, f, fp);
        const double af = fabs(f);
        if (af < best_f) { best_f = af; best_u = u; }
        if (af < ftol) return u;
        if (fabs(fp) < 1e-30) break;
        double du = f / fp;
        if (du > 1.0) du = 1.0; else if (du < -1.0) du = -1.0;
        u -= du;
        if (fabs(du) < 1e-13) return u;
    }
    double lo = -2.0, hi = 2.0, flo, fhi, fpd;
    ev(lo, A, B, C, Ds, flo, fpd); ev(hi, A, B, C, Ds, fhi, fpd);
    if (flo * fhi > 0.0) return best_u;
    for (int it = 0; it < 60; ++it) {
        const double m = 0.5 * (lo + hi); double fm;
        ev(m, A, B, C, Ds, fm, fpd);
        if (fabs(fm) < ftol) return m;
        if (fm * flo < 0.0) { hi = m; fhi = fm; } else { lo = m; flo = fm; }
        if (hi - lo < 1e-15) return m;
    }
    return 0.5 * (lo + hi);
}

struct Ctx { const double* cy; const double* lines; int n_comp, n_lines, n_keep;
             double x_lo, span, t_scale; double* resid; double* uroot; };

static double score_from(const Ctx& ctx) {
    nth_element(ctx.resid, ctx.resid + ctx.n_keep, ctx.resid + ctx.n_lines);
    double sum = 0.0;
    for (int i = 0; i < ctx.n_keep; ++i) sum += ctx.resid[i];
    return sum * 1000.0 / ctx.n_keep;
}

static void abcd(double cu, double qu, double sp, double wl, const Ctx& c,
                 double& A, double& B, double& C, double& D) {
    const double S = c.span, x = c.x_lo;
    A = cu * S * S * S;
    B = S * S * (qu + 3.0 * cu * x);
    C = S * (sp + x * (2.0 * qu + 3.0 * cu * x));
    D = wl + x * (sp + x * (qu + x * cu));
}

// ---- Reference: robust solver for every line -------------------------------
double objective_ref(double cu, double qu, double sp, double wl, const Ctx& ctx) {
    double A, B, C, D; abcd(cu, qu, sp, wl, ctx, A, B, C, D);
    for (int i = 0; i < ctx.n_lines; ++i) {
        const double u = inverse_unit_poly(ctx.lines[i], A, B, C, D);
        if (!(u >= 0.0 && u <= 1.0)) { ctx.resid[i] = 1.0; continue; }
        const double t = u * ctx.t_scale;
        int j = (int)t; if (j >= ctx.n_comp - 1) j = ctx.n_comp - 2;
        const double fr = t - j;
        const double y = ctx.cy[j] + fr * (ctx.cy[j + 1] - ctx.cy[j]);
        ctx.resid[i] = (1.0 - y) * (1.0 - y);
    }
    return score_from(ctx);
}

// ---- Fast: monotonicity-gated vectorised Newton + verified fallback --------
static long long n_slowpath = 0, n_proposals = 0;

static inline bool dispersion_is_monotonic(double A3, double B2, double C) {
    const double d0 = C, d1 = C + B2 + A3;
    if (d0 == 0.0 || d1 == 0.0) return false;
    if ((d0 > 0.0) != (d1 > 0.0)) return false;
    if (A3 != 0.0) {
        const double uv = -B2 / (2.0 * A3);
        if (uv > 0.0 && uv < 1.0) {
            const double dv = C + uv * (B2 + uv * A3);
            if (dv == 0.0 || (dv > 0.0) != (d0 > 0.0)) return false;
        }
    }
    return true;
}

double objective_fast(double cu, double qu, double sp, double wl, const Ctx& ctx) {
    double A, B, C, D; abcd(cu, qu, sp, wl, ctx, A, B, C, D);
    const double B2 = 2.0 * B, A3 = 3.0 * A;
    const double* __restrict lines = ctx.lines;
    double* __restrict uroot = ctx.uroot;
    ++n_proposals;

    if (!dispersion_is_monotonic(A3, B2, C)) {
        ++n_slowpath;
        return objective_ref(cu, qu, sp, wl, ctx);
    }

    const double p0 = D, p1 = D + C + B + A;
    const double lam_lo = p0 < p1 ? p0 : p1, lam_hi = p0 < p1 ? p1 : p0;
    const double invC = 1.0 / C;

    #pragma omp simd
    for (int i = 0; i < ctx.n_lines; ++i) {
        double Ds = D - lines[i];
        double u = -Ds * invC;
        for (int it = 0; it < 6; ++it) {
            double fv = Ds + u * (C + u * (B + u * A));
            double fp = C + u * (B2 + u * A3);
            u -= fv / fp;
        }
        uroot[i] = u;
    }

    for (int i = 0; i < ctx.n_lines; ++i) {
        const double lam = lines[i];
        if (lam < lam_lo || lam > lam_hi) { ctx.resid[i] = 1.0; continue; }
        double u = uroot[i];
        ++n_line_evals;
        const double rf = (D - lam) + u * (C + u * (B + u * A));
        if (!(fabs(rf) < 1e-9 * (fabs(lam) + 1.0)) || !(u >= 0.0 && u <= 1.0)) {
            u = inverse_unit_poly(lam, A, B, C, D);
            ++n_fallback;
            if (!(u >= 0.0 && u <= 1.0)) { ctx.resid[i] = 1.0; continue; }
        }
        const double t = u * ctx.t_scale;
        int j = (int)t; if (j >= ctx.n_comp - 1) j = ctx.n_comp - 2;
        const double fr = t - j;
        const double y = ctx.cy[j] + fr * (ctx.cy[j + 1] - ctx.cy[j]);
        ctx.resid[i] = (1.0 - y) * (1.0 - y);
    }
    return score_from(ctx);
}

struct P { const char* n; const char* ll; double off, lin, quad, cub,
           ccov, scov, qcov, cubcov; int npix; };

int main() {
    vector<P> ps = {
        {"SOAR_930","linelists/SOAR_930_M1_FeAr.txt",
         3397.5,1725.0/1740.0,-7e-6,0.0,150,0.10,5e-5,3e-8,1740},
        {"ALFOSC_grism18","linelists/ALFOSC_grism18_ThAr.txt",
         3450.0,1900.0/2048.0,0.0,0.0,150,0.15,1e-4,3e-8,2048},
        {"CAFOS_b100","linelists/CAFOS_B100_arclamp.txt",
         2552.0,2547.0/2048.0,5.3429e-4,-1.116e-7,141,0.37,3.1e-4,8.17e-8,2048},
        {"EFOSC_grism19","linelists/EFOSC_grism19_HeAr.txt",
         4444.0,680.0/1024.0,5e-5,-1.3e-8,100,0.25,5e-4,3e-7,1024},
        {"EFOSC_grism7","linelists/EFOSC_grism7_HeAr.txt",
         3250.0,1800.0/1024.0,5e-5,-1.3e-8,100,0.50,1e-3,6e-7,1024},
    };

    for (const auto& p : ps) {
        vector<double> lines;
        FILE* fp = fopen(p.ll, "r");
        if (!fp) { printf("%-16s SKIP (no linelist)\n", p.n); continue; }
        double w; char rest[128];
        while (fscanf(fp, "%lf %127s", &w, rest) == 2) lines.push_back(w);
        fclose(fp);

        // Synthetic arc lamp at the nominal solution, same shape as the real one
        const int N = p.npix;
        vector<double> cy(N, 0.0);
        for (int i = 0; i < N; ++i) {
            double wl = p.off + i * (p.lin + i * (p.quad + i * p.cub));
            for (double l : lines) cy[i] += 0.8 * exp(-0.5 * pow((wl - l) / 1.6, 2));
        }
        double mx = 0; for (double v : cy) mx = max(mx, v);
        if (mx > 0) for (auto& v : cy) v /= mx;

        vector<double> resid(lines.size()), uroot(lines.size());
        Ctx ctx{cy.data(), lines.data(), N, (int)lines.size(),
                max(1, ((int)lines.size() * 4 + 4) / 5),
                0.0, (double)(N - 1), (double)(N - 1), resid.data(), uroot.data()};

        mt19937 g(2024);
        uniform_real_distribution<double> u01(0, 1);
        double max_abs = 0.0, max_rel = 0.0;
        long long n = 0, n_diff = 0;
        n_fallback = n_line_evals = 0; n_slowpath = n_proposals = 0;

        for (int tr = 0; tr < 200000; ++tr) {
            double off = p.off + (u01(g) - 0.5) * p.ccov;
            double lin = p.lin + (u01(g) - 0.5) * p.scov;
            double qu  = p.quad + (u01(g) - 0.5) * p.qcov;
            double cu  = p.cub + (u01(g) - 0.5) * p.cubcov;

            double s_ref = objective_ref(cu, qu, lin, off, ctx);
            double s_fst = objective_fast(cu, qu, lin, off, ctx);
            ++n;
            double d = fabs(s_ref - s_fst);
            if (d != 0.0) ++n_diff;
            if (d > max_abs) max_abs = d;
            double r = d / (fabs(s_ref) + 1e-300);
            if (r > max_rel) max_rel = r;
        }
        printf("%-16s proposals=%lld  differing=%lld  max|dScore|=%.3e  "
               "max rel=%.3e\n", p.n, n, n_diff, max_abs, max_rel);
        printf("%-16s fallback = %.4f%% of line solves | robust path = %.3f%% of proposals\n\n",
               "", 100.0 * (double)n_fallback / (double)n_line_evals,
               100.0 * (double)n_slowpath / (double)n_proposals);
    }
    return 0;
}
