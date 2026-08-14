// Instruments the current inverse_cub_poly(): how many polynomial evaluations
// does it burn per line, and does Newton ever actually converge?
#include <cstdio>
#include <cmath>
#include <random>
#include <vector>
using namespace std;

static long long n_newton_evals = 0, n_bisect_evals = 0, n_calls = 0,
                 n_newton_converged = 0, n_bisection_used = 0;

inline void eval_poly_and_deriv(double x, double a, double b, double c, double d,
                                double &f, double &fp) {
    f = d + x * (c + x * (b + x * a));
    fp = c + x * (2.0 * b + x * 3.0 * a);
}

// Verbatim copy of the shipped inverse_cub_poly(), plus counters.
double inverse_cub_poly(double y, double a, double b, double c, double d) {
    ++n_calls;
    double d_shifted = d - y;
    double abs_a = fabs(a), abs_b = fabs(b), abs_c = fabs(c);
    if (abs_a < 1e-20 && abs_b < 1e-15) { if (abs_c < 1e-20) return NAN; return -d_shifted / c; }
    if (abs_a < 1e-20) {
        double disc = c * c - 4.0 * b * d_shifted;
        if (disc < 0) return NAN;
        double sq = sqrt(disc);
        double r1 = (-c + sq) / (2.0 * b), r2 = (-c - sq) / (2.0 * b);
        double la = -d_shifted / c;
        return (fabs(r1 - la) < fabs(r2 - la)) ? r1 : r2;
    }
    double x;
    if (abs_c > 1e-20) x = -d_shifted / c;
    else if (abs_b > 1e-20) { double disc = -d_shifted / b; x = (disc > 0) ? sqrt(disc) : 0.5; }
    else x = cbrt(-d_shifted / a);
    x = fmax(-10.0, fmin(10.0, x));
    double f, fp, prev_x = x, best_x = x;
    eval_poly_and_deriv(x, a, b, c, d_shifted, f, fp); ++n_newton_evals;
    double best_f = fabs(f);
    for (int iter = 0; iter < 30; ++iter) {
        eval_poly_and_deriv(x, a, b, c, d_shifted, f, fp); ++n_newton_evals;
        if (fabs(f) < best_f) { best_f = fabs(f); best_x = x; }
        if (fabs(f) < 1e-14 * (fabs(y) + 1.0)) { ++n_newton_converged; return x; }
        if (fabs(fp) < 1e-30) { x += (f > 0 ? -1e-6 : 1e-6); continue; }
        double dx = f / fp;
        if (fabs(dx) > 1.0) dx = (dx > 0) ? 1.0 : -1.0;
        prev_x = x; x -= dx;
        if (fabs(dx) < 1e-14 * (fabs(x) + 1e-10)) { ++n_newton_converged; return x; }
        if (iter > 5 && fabs(x - prev_x) < 1e-13) break;
    }
    ++n_bisection_used;
    double x_lo = -2.0, x_hi = 2.0, f_lo, f_hi, fpd;
    eval_poly_and_deriv(x_lo, a, b, c, d_shifted, f_lo, fpd); ++n_bisect_evals;
    eval_poly_and_deriv(x_hi, a, b, c, d_shifted, f_hi, fpd); ++n_bisect_evals;
    for (int i = 0; i < 10 && f_lo * f_hi > 0; ++i) {
        x_lo *= 2.0; x_hi *= 2.0;
        eval_poly_and_deriv(x_lo, a, b, c, d_shifted, f_lo, fpd);
        eval_poly_and_deriv(x_hi, a, b, c, d_shifted, f_hi, fpd);
        n_bisect_evals += 2;
    }
    if (f_lo * f_hi > 0) return best_x;
    for (int iter = 0; iter < 60; ++iter) {
        double x_mid = 0.5 * (x_lo + x_hi), f_mid;
        eval_poly_and_deriv(x_mid, a, b, c, d_shifted, f_mid, fpd); ++n_bisect_evals;
        if (fabs(f_mid) < 1e-14 * (fabs(y) + 1.0)) return x_mid;
        if (f_mid * f_lo < 0) { x_hi = x_mid; f_hi = f_mid; }
        else { x_lo = x_mid; f_lo = f_mid; }
        if (x_hi - x_lo < 1e-14 * (fabs(x_mid) + 1e-10)) return x_mid;
    }
    return 0.5 * (x_lo + x_hi);
}

int main() {
    // SOAR_930-like search box, as in bench_mcmc.py
    double off0 = 3397.5, lin0 = 1725.0 / 1740.0, quad0 = -7e-6, cub0 = 0.0;
    double c_cov = 150, s_cov = 0.1, q_cov = 5e-5, cub_cov = 3e-8;
    vector<double> lines;
    FILE *fp = fopen("linelists/SOAR_930_M1_FeAr.txt", "r");
    double w; char rest[64];
    while (fscanf(fp, "%lf %63s", &w, rest) == 2) lines.push_back(w);
    fclose(fp);
    printf("n_lines = %zu\n", lines.size());

    mt19937 gen(42);
    uniform_real_distribution<double> u(0.0, 1.0);
    for (int trial = 0; trial < 20000; ++trial) {
        // Random parameter draw from the search box, like an MCMC proposal
        double off = off0 + (u(gen) - 0.5) * c_cov;
        double lin = lin0 + (u(gen) - 0.5) * s_cov;
        double qu = quad0 + (u(gen) - 0.5) * q_cov;
        double cu = cub0 + (u(gen) - 0.5) * cub_cov;
        for (double l : lines) inverse_cub_poly(l, cu, qu, lin, off);
    }
    printf("calls                 = %lld\n", n_calls);
    printf("newton poly-evals     = %lld  (%.1f per call)\n", n_newton_evals,
           (double)n_newton_evals / n_calls);
    printf("bisection poly-evals  = %lld  (%.1f per call)\n", n_bisect_evals,
           (double)n_bisect_evals / n_calls);
    printf("total poly-evals/call = %.1f\n",
           (double)(n_newton_evals + n_bisect_evals) / n_calls);
    printf("newton converged      = %.2f%% of calls\n",
           100.0 * n_newton_converged / n_calls);
    printf("fell back to bisection= %.2f%% of calls\n",
           100.0 * n_bisection_used / n_calls);
    return 0;
}
