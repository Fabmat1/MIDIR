// Direct equivalence check: does the new normalised-coordinate root solver
// return the same pixel positions as the old raw-pixel one, and how many
// polynomial evaluations does each burn?
#include <cstdio>
#include <cmath>
#include <random>
#include <vector>
#include <algorithm>
using namespace std;

static long long evals_old = 0, evals_new = 0;

inline void ev(double x, double a, double b, double c, double d, double &f, double &fp) {
    f = d + x * (c + x * (b + x * a));
    fp = c + x * (2.0 * b + x * 3.0 * a);
}

// ---- OLD (shipped) ---------------------------------------------------------
double inverse_cub_poly(double y, double a, double b, double c, double d) {
    double ds = d - y;
    double aa = fabs(a), ab = fabs(b), ac = fabs(c);
    if (aa < 1e-20 && ab < 1e-15) { if (ac < 1e-20) return NAN; return -ds / c; }
    if (aa < 1e-20) {
        double disc = c * c - 4.0 * b * ds;
        if (disc < 0) return NAN;
        double sq = sqrt(disc);
        double r1 = (-c + sq) / (2.0 * b), r2 = (-c - sq) / (2.0 * b), la = -ds / c;
        return (fabs(r1 - la) < fabs(r2 - la)) ? r1 : r2;
    }
    double x;
    if (ac > 1e-20) x = -ds / c;
    else if (ab > 1e-20) { double dd = -ds / b; x = (dd > 0) ? sqrt(dd) : 0.5; }
    else x = cbrt(-ds / a);
    x = fmax(-10.0, fmin(10.0, x));
    double f, fp, prev = x, best = x;
    ev(x, a, b, c, ds, f, fp); ++evals_old;
    double bf = fabs(f);
    for (int i = 0; i < 30; ++i) {
        ev(x, a, b, c, ds, f, fp); ++evals_old;
        if (fabs(f) < bf) { bf = fabs(f); best = x; }
        if (fabs(f) < 1e-14 * (fabs(y) + 1.0)) return x;
        if (fabs(fp) < 1e-30) { x += (f > 0 ? -1e-6 : 1e-6); continue; }
        double dx = f / fp;
        if (fabs(dx) > 1.0) dx = (dx > 0) ? 1.0 : -1.0;
        prev = x; x -= dx;
        if (fabs(dx) < 1e-14 * (fabs(x) + 1e-10)) return x;
        if (i > 5 && fabs(x - prev) < 1e-13) break;
    }
    double lo = -2.0, hi = 2.0, flo, fhi, fpd;
    ev(lo, a, b, c, ds, flo, fpd); ev(hi, a, b, c, ds, fhi, fpd); evals_old += 2;
    for (int i = 0; i < 10 && flo * fhi > 0; ++i) {
        lo *= 2.0; hi *= 2.0;
        ev(lo, a, b, c, ds, flo, fpd); ev(hi, a, b, c, ds, fhi, fpd); evals_old += 2;
    }
    if (flo * fhi > 0) return best;
    for (int i = 0; i < 60; ++i) {
        double m = 0.5 * (lo + hi), fm;
        ev(m, a, b, c, ds, fm, fpd); ++evals_old;
        if (fabs(fm) < 1e-14 * (fabs(y) + 1.0)) return m;
        if (fm * flo < 0) { hi = m; fhi = fm; } else { lo = m; flo = fm; }
        if (hi - lo < 1e-14 * (fabs(m) + 1e-10)) return m;
    }
    return 0.5 * (lo + hi);
}

// ---- NEW (normalised) ------------------------------------------------------
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
        ev(u, A, B, C, Ds, f, fp); ++evals_new;
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
    ev(lo, A, B, C, Ds, flo, fpd); ev(hi, A, B, C, Ds, fhi, fpd); evals_new += 2;
    if (flo * fhi > 0.0) return best_u;
    for (int it = 0; it < 60; ++it) {
        const double m = 0.5 * (lo + hi); double fm;
        ev(m, A, B, C, Ds, fm, fpd); ++evals_new;
        if (fabs(fm) < ftol) return m;
        if (fm * flo < 0.0) { hi = m; fhi = fm; } else { lo = m; flo = fm; }
        if (hi - lo < 1e-15) return m;
    }
    return 0.5 * (lo + hi);
}

struct Preset { const char *name; double off, lin, quad, cub, c_cov, s_cov, q_cov, cub_cov; int npix; };

int main() {
    // Search boxes taken from the shipped presets/*.json
    vector<Preset> presets = {
        {"SOAR_930",       3397.5, 1725.0/1740.0, -7e-6,     0.0,      150, 0.10, 5e-5,  3e-8,  1740},
        {"ALFOSC_grism18", 3450.0, 1900.0/2048.0,  0.0,      0.0,      150, 0.15, 1e-4,  3e-8,  2048},
        {"CAFOS_b100",     2552.0, 2547.0/2048.0,  5.3429e-4,-1.116e-7,141, 0.37, 3.1e-4,8.17e-8,2048},
        {"EFOSC_grism19",  4444.0,  680.0/1024.0,  5e-5,     -1.3e-8,  100, 0.25, 5e-4,  3e-7,  1024},
        {"EFOSC_grism7",   3250.0, 1800.0/1024.0,  5e-5,     -1.3e-8,  100, 0.50, 1e-3,  6e-7,  1024},
    };

    mt19937 gen(7);
    uniform_real_distribution<double> u01(0.0, 1.0);

    for (const auto &p : presets) {
        evals_old = evals_new = 0;
        long long n = 0, n_inrange = 0, n_disagree = 0;
        double worst_px = 0.0;
        // How often does each solver recover the pixel the probe wavelength was
        // actually generated from? (Only meaningful where the dispersion is
        // monotonic; where it is not, "the" root is genuinely ambiguous.)
        long long old_wrong = 0, new_wrong = 0, nonmono = 0;
        const double S = p.npix - 1;

        for (int trial = 0; trial < 40000; ++trial) {
            double off = p.off + (u01(gen) - 0.5) * p.c_cov;
            double lin = p.lin + (u01(gen) - 0.5) * p.s_cov;
            double qu  = p.quad + (u01(gen) - 0.5) * p.q_cov;
            double cu  = p.cub + (u01(gen) - 0.5) * p.cub_cov;

            // Is the dispersion monotonic over the chip? f'(x) = 3a x^2+2b x+c
            bool mono = true;
            for (int k = 0; k <= 64; ++k) {
                double xx = S * k / 64.0;
                if (lin + xx * (2.0 * qu + 3.0 * cu * xx) <= 0.0) { mono = false; break; }
            }
            if (!mono) ++nonmono;

            // Probe wavelengths spread across the nominal coverage
            for (int k = 0; k < 20; ++k) {
                double frac = (k + 0.5) / 20.0;
                double x_true = S * frac;
                double y = off + x_true * (lin + x_true * (qu + x_true * cu));

                double x_old = inverse_cub_poly(y, cu, qu, lin, off);

                double A = cu * S * S * S, B = S * S * qu, C = S * lin, D = off;
                double x_new = inverse_unit_poly(y, A, B, C, D) * S;

                ++n;
                bool in_old = (x_old == x_old && x_old >= 0.0 && x_old <= S);
                bool in_new = (x_new == x_new && x_new >= 0.0 && x_new <= S);
                if (in_old != in_new) { ++n_disagree; continue; }
                if (!in_old) continue;
                ++n_inrange;
                double d = fabs(x_old - x_new);
                if (d > worst_px) worst_px = d;
                if (mono) {
                    if (fabs(x_old - x_true) > 1e-4) ++old_wrong;
                    if (fabs(x_new - x_true) > 1e-4) ++new_wrong;
                }
            }
        }
        printf("%-16s in-range=%lld/%lld  range-disagreements=%lld  "
               "max |dx| = %.3e px\n", p.name, n_inrange, n, n_disagree, worst_px);
        printf("%-16s non-monotonic draws=%lld/40000 | wrong root (monotonic only): "
               "old=%lld  new=%lld\n", "", nonmono, old_wrong, new_wrong);
        printf("%-16s poly-evals/call: old=%.1f  new=%.1f  (%.1fx fewer)\n\n", "",
               (double)evals_old / n, (double)evals_new / n,
               (double)evals_old / evals_new);
    }
    return 0;
}
