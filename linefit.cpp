#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <random>
#include <mutex>
#include <atomic>
#include <thread>
#include <cstdio>
#include <numeric>
#include <limits>
#include <cstring>
#include <memory>

std::mutex cout_mutex;
std::atomic<long long> global_steps_completed(0);
std::atomic<long long> global_steps_accepted(0);
std::atomic<int> threads_running(0);
std::atomic<bool> stop_reporter(false);
// Set from Python (request_abort) when the user aborts the reduction. The
// sampling chains poll this and unwind at their next checkpoint.
std::atomic<bool> abort_requested(false);

using namespace std;

// Thread-safe random device initialization
static thread_local mt19937* thread_local_gen = nullptr;

mt19937& get_thread_gen() {
    if (!thread_local_gen) {
        random_device rd;
        thread_local_gen = new mt19937(rd() ^ (std::hash<std::thread::id>{}(std::this_thread::get_id())));
    }
    return *thread_local_gen;
}

// ============================================================================
// POLYNOMIAL FUNCTIONS
// ============================================================================

double cub_poly(double x, double a, double b, double c, double d) {
    return x * x * x * a + x * x * b + x * c + d;
}

// ============================================================================
// NUMERICALLY STABLE INVERSE CUBIC POLYNOMIAL
// ============================================================================

// Evaluate polynomial and its derivative simultaneously
inline void eval_poly_and_deriv(double x, double a, double b, double c, double d,
                                 double& f, double& fp) {
    // f(x) = ax³ + bx² + cx + d
    // f'(x) = 3ax² + 2bx + c
    // Use Horner's method for numerical stability
    f = d + x * (c + x * (b + x * a));
    fp = c + x * (2.0 * b + x * 3.0 * a);
}

// ---------------------------------------------------------------------------
// Inverse of the dispersion polynomial, solved in the NORMALISED pixel
// coordinate u = (x - x_lo) / (x_hi - x_lo), so that the detector spans
// u ∈ [0, 1].
//
// This normalisation is what makes the solver fast. The step damping below
// (|du| ≤ 1) and the bisection bracket ([-2, 2]) are only meaningful if one
// unit of the abscissa is comparable to the width of the domain. Feeding this
// routine raw pixel indices (0…2000) instead - as the caller used to - caps
// Newton at one *pixel* per iteration while the initial guess is ~100 px off,
// so it never converges and every single inversion falls through to a 60-step
// bisection: ~94 polynomial evaluations per line instead of ~4.
//
// Caller passes the coefficients already recast onto u:
//     A u³ + B u² + C u + D = y
// which for x = x_lo + u·S (S = x_hi - x_lo) means
//     A = a·S³, B = S²(b + 3a·x_lo), C = S(c + 2b·x_lo + 3a·x_lo²),
//     D = d + c·x_lo + b·x_lo² + a·x_lo³.
// In these units C is the full wavelength extent (~10³ Å) and B, A are the
// total quadratic/cubic excursions across the chip, so the problem is well
// conditioned regardless of instrument.
//
// Returns NaN if no usable root was found; the caller treats anything outside
// [0, 1] as a non-match.
// ---------------------------------------------------------------------------
inline double inverse_unit_poly(double y, double A, double B, double C, double D) {
    const double D_shifted = D - y;
    const double ftol = 1e-14 * (fabs(y) + 1.0);

    // --- Initial guess -----------------------------------------------------
    // Start from the linear solution, then upgrade it with the quadratic term
    // when that term is large enough to matter (it is, for grisms with strong
    // curvature). This puts us inside Newton's quadratic-convergence basin.
    double u = -D_shifted / C;
    if (B != 0.0) {
        const double disc = C * C - 4.0 * B * D_shifted;
        if (disc >= 0.0) {
            // Stable quadratic roots: avoids cancellation when 4BD ≪ C².
            const double q = -0.5 * (C + copysign(sqrt(disc), C));
            const double r1 = q / B;
            const double r2 = (q != 0.0) ? (D_shifted / q) : r1;
            u = (fabs(r1 - u) < fabs(r2 - u)) ? r1 : r2;
        }
    }
    // Also catches NaN/Inf from a degenerate C or B.
    if (!(u > -10.0 && u < 10.0)) u = 0.5;

    // --- Newton-Raphson ----------------------------------------------------
    double best_u = u;
    double best_f = HUGE_VAL;

    for (int iter = 0; iter < 12; ++iter) {
        double f, fp;
        eval_poly_and_deriv(u, A, B, C, D_shifted, f, fp);

        const double af = fabs(f);
        if (af < best_f) { best_f = af; best_u = u; }
        if (af < ftol) return u;

        if (fabs(fp) < 1e-30) break;  // stationary point - hand over to bisection

        double du = f / fp;
        if (du > 1.0) du = 1.0;
        else if (du < -1.0) du = -1.0;
        u -= du;

        // One unit of u is the whole detector, so this is ~1e-10 px: far below
        // any precision the interpolation downstream can make use of.
        if (fabs(du) < 1e-13) return u;
    }

    // --- Bisection fallback ------------------------------------------------
    // [-2, 2] brackets the detector with a wide margin. A root outside it maps
    // outside [0, 1] and is rejected by the caller anyway, so there is nothing
    // to gain from expanding the bracket.
    double u_lo = -2.0, u_hi = 2.0, f_lo, f_hi, fp_dummy;
    eval_poly_and_deriv(u_lo, A, B, C, D_shifted, f_lo, fp_dummy);
    eval_poly_and_deriv(u_hi, A, B, C, D_shifted, f_hi, fp_dummy);

    if (f_lo * f_hi > 0.0) {
        // No sign change: either no root here, or an even number of them.
        // Fall back on Newton's best effort, as before.
        return best_u;
    }

    for (int iter = 0; iter < 60; ++iter) {
        const double u_mid = 0.5 * (u_lo + u_hi);
        double f_mid;
        eval_poly_and_deriv(u_mid, A, B, C, D_shifted, f_mid, fp_dummy);

        if (fabs(f_mid) < ftol) return u_mid;

        if (f_mid * f_lo < 0.0) { u_hi = u_mid; f_hi = f_mid; }
        else                    { u_lo = u_mid; f_lo = f_mid; }

        if (u_hi - u_lo < 1e-15) return u_mid;
    }

    return 0.5 * (u_lo + u_hi);
}

// ============================================================================
// OBJECTIVE FUNCTION - Robust version using trimmed mean + match counting
// ============================================================================

struct ObjectiveResult {
    double score;           // Primary score (lower is better)
    int n_matched;          // Number of lines with value > threshold
    double trimmed_chisq;   // Trimmed chi-squared (ignores worst outliers)
};


// Everything about the comparison spectrum and the line list that does not
// change from one MCMC proposal to the next. Built once per chain so the inner
// loop touches raw pointers only.
struct ObjectiveCtx {
    const double* compspec_y;
    const double* lines;
    int n_comp;
    int n_lines;
    int n_keep;         // how many residuals survive the 80% trim
    double x_lo;
    double span;        // x_hi - x_lo
    double t_scale;     // n_comp - 1: normalised position -> grid index
    double* resid;      // scratch, n_lines doubles, owned by the chain
    double* uroot;      // scratch, n_lines doubles, owned by the chain
};

// Squared residual for a line that landed at normalised position u ∈ [0,1].
static inline double residual_at(double u, double t_scale, int n_comp,
                                 const double* __restrict cy) {
    // u ∈ [0,1] already, so the grid index needs no lower clamp.
    const double t = u * t_scale;
    int j = (int)t;
    if (j >= n_comp - 1) j = n_comp - 2;

    const double frac = t - j;
    const double y = cy[j] + frac * (cy[j + 1] - cy[j]);

    const double resid = 1.0 - y;
    return resid * resid;
}

// Is the dispersion strictly monotonic across the whole chip? p'(u) is a
// quadratic, so it suffices to check both ends plus the vertex if it falls
// inside [0,1]. When this holds, a wavelength has at most one pre-image on the
// detector, which is what lets the fast path below trust any root it finds.
static inline bool dispersion_is_monotonic(double A3, double B2, double C) {
    const double d0 = C;                    // p'(0)
    const double d1 = C + B2 + A3;          // p'(1)
    if (d0 == 0.0 || d1 == 0.0) return false;
    if ((d0 > 0.0) != (d1 > 0.0)) return false;

    if (A3 != 0.0) {
        const double u_vertex = -B2 / (2.0 * A3);
        if (u_vertex > 0.0 && u_vertex < 1.0) {
            const double dv = C + u_vertex * (B2 + u_vertex * A3);
            if (dv == 0.0 || (dv > 0.0) != (d0 > 0.0)) return false;
        }
    }
    return true;
}

double compute_objective(double cubic_fac, double quadratic_fac, double spacing,
                         double wl_start, const ObjectiveCtx& ctx) {
    const int n_lines = ctx.n_lines;
    const int n_comp = ctx.n_comp;
    const double* __restrict cy = ctx.compspec_y;
    const double* __restrict lines = ctx.lines;
    double* __restrict residuals = ctx.resid;

    // Recast the dispersion polynomial onto the normalised pixel coordinate
    // u = (x - x_lo)/span once per proposal, rather than fighting the raw
    // pixel scale inside every line's root solve. See inverse_unit_poly().
    const double S = ctx.span;
    const double x_lo = ctx.x_lo;
    const double A = cubic_fac * S * S * S;
    const double B = S * S * (quadratic_fac + 3.0 * cubic_fac * x_lo);
    const double C = S * (spacing + x_lo * (2.0 * quadratic_fac + 3.0 * cubic_fac * x_lo));
    const double D = wl_start + x_lo * (spacing + x_lo * (quadratic_fac + x_lo * cubic_fac));

    const double t_scale = ctx.t_scale;
    const double A3 = 3.0 * A, B2 = 2.0 * B;

    if (dispersion_is_monotonic(A3, B2, C)) {
        // ---- Fast path ----------------------------------------------------
        // The dispersion is invertible over the chip, so a line is observable
        // iff its wavelength lies between the two chip edges, and the root in
        // [0,1] - if there is one - is unique. That makes a plain fixed-count
        // Newton safe: any converged root inside [0,1] is necessarily *the*
        // root, so we can drop the branches and let this vectorise.
        const double p0 = D;                    // wavelength at u = 0
        const double p1 = D + C + B + A;        // wavelength at u = 1
        const double lam_lo = p0 < p1 ? p0 : p1;
        const double lam_hi = p0 < p1 ? p1 : p0;

        const double invC = 1.0 / C;
        double* __restrict uroot = ctx.uroot;

        #pragma omp simd
        for (int i = 0; i < n_lines; ++i) {
            double Ds = D - lines[i];
            double u = -Ds * invC;
            // Hand-unrolled: an inner loop counts as control flow and would
            // stop the vectoriser. Six undamped iterations from the linear
            // guess reach ~1e-12 px for every shipped preset; the rare
            // non-convergent line is caught below.
            for (int it = 0; it < 6; ++it) {
                double fv = Ds + u * (C + u * (B + u * A));
                double fp = C + u * (B2 + u * A3);
                u -= fv / fp;
            }
            uroot[i] = u;
        }

        for (int i = 0; i < n_lines; ++i) {
            const double lam = lines[i];

            // Off the chip: no root in [0,1] exists, no solve needed.
            if (lam < lam_lo || lam > lam_hi) {
                residuals[i] = 1.0;
                continue;
            }

            double u = uroot[i];
            // Did the fast path land on the root? A non-converged lane, a NaN,
            // or a root outside [0,1] fails this and is redone properly.
            const double fv = (D - lam) + u * (C + u * (B + u * A));
            if (!(fabs(fv) < 1e-9 * (fabs(lam) + 1.0)) || !(u >= 0.0 && u <= 1.0)) {
                u = inverse_unit_poly(lam, A, B, C, D);
                if (!(u >= 0.0 && u <= 1.0)) {
                    residuals[i] = 1.0;
                    continue;
                }
            }

            residuals[i] = residual_at(u, t_scale, n_comp, cy);
        }
    } else {
        // ---- Robust path --------------------------------------------------
        // Non-monotonic dispersion: several pixels can share a wavelength, so
        // "the" root is ambiguous and which one you land on matters. Use the
        // careful solver so the objective stays exactly what it always was.
        for (int i = 0; i < n_lines; ++i) {
            const double u = inverse_unit_poly(lines[i], A, B, C, D);

            // Rejects NaN as well: any comparison with NaN is false.
            if (!(u >= 0.0 && u <= 1.0)) {
                residuals[i] = 1.0;
                continue;
            }

            residuals[i] = residual_at(u, t_scale, n_comp, cy);
        }
    }

    if (n_lines == 0) return 1e6;

    const int n_keep = ctx.n_keep;

    // Partial sort: elements before nth are <= nth, elements after are >= nth
    nth_element(residuals, residuals + n_keep, residuals + n_lines);

    double sum = 0.0;
    for (int i = 0; i < n_keep; ++i) {
        sum += residuals[i];
    }

    return sum * 1000.0 / n_keep;
}

// ============================================================================
// PARALLEL TEMPERING MCMC
// ============================================================================

// Single parallel tempering chain group (multiple temperatures)
void run_parallel_tempering_chain(
    const vector<double>& compspec_x, const vector<double>& compspec_y,
    const vector<double>& lines,
    int n_samples, int n_burn_in,
    double wl_init, double spacing_init, double quad_init, double cub_init,
    double wl_lo, double wl_hi,
    double spacing_lo, double spacing_hi,
    double quad_lo, double quad_hi,
    double cub_lo, double cub_hi,
    double* out_wl, double* out_spacing, double* out_quad, double* out_cub, double* out_score,
    mt19937& gen) {

    threads_running.fetch_add(1);

    uniform_real_distribution<double> uniform(0.0, 1.0);
    normal_distribution<double> normal(0.0, 1.0);

    // Build the invariant part of the objective once for this chain.
    const int n_lines = (int)lines.size();
    const int n_comp = (int)compspec_x.size();
    vector<double> resid_scratch(n_lines > 0 ? n_lines : 1);
    vector<double> uroot_scratch(n_lines > 0 ? n_lines : 1);

    ObjectiveCtx ctx;
    ctx.compspec_y = compspec_y.data();
    ctx.lines      = lines.data();
    ctx.n_comp     = n_comp;
    ctx.n_lines    = n_lines;
    ctx.n_keep     = max(1, (n_lines * 4 + 4) / 5);  // keep 80%, round up
    ctx.x_lo       = compspec_x[0];
    ctx.span       = compspec_x[n_comp - 1] - compspec_x[0];
    ctx.t_scale    = (double)(n_comp - 1);
    ctx.resid      = resid_scratch.data();
    ctx.uroot      = uroot_scratch.data();

    const int n_temps = 6;
    const double temperatures[6] = {1.0, 3.0, 10.0, 30.0, 100.0, 300.0};
    const double inv_temperatures[6] = {1.0, 1.0/3.0, 0.1, 1.0/30.0, 0.01, 1.0/300.0};

    // Initialize chains
    struct Chain { double wl, sp, qu, cu, score; };
    Chain chains[6];
    
    for (int t = 0; t < n_temps; ++t) {
        if (t == 0) {
            chains[t] = {wl_init, spacing_init, quad_init, cub_init, 0.0};
        } else {
            chains[t].wl = wl_lo + uniform(gen) * (wl_hi - wl_lo);
            chains[t].sp = spacing_lo + uniform(gen) * (spacing_hi - spacing_lo);
            chains[t].qu = quad_lo + uniform(gen) * (quad_hi - quad_lo);
            chains[t].cu = cub_lo + uniform(gen) * (cub_hi - cub_lo);
        }
        chains[t].score = compute_objective(chains[t].cu, chains[t].qu,
                                            chains[t].sp, chains[t].wl, ctx);
    }

    // Step sizes
    double step_wl[6], step_sp[6], step_qu[6], step_cu[6];
    const double base_frac = 0.01;
    for (int t = 0; t < n_temps; ++t) {
        double scale = sqrt(temperatures[t]);
        step_wl[t] = (wl_hi - wl_lo) * base_frac * scale;
        step_sp[t] = (spacing_hi - spacing_lo) * base_frac * scale;
        step_qu[t] = (quad_hi - quad_lo) * base_frac * scale;
        step_cu[t] = (cub_hi - cub_lo) * base_frac * scale;
    }

    // Acceptance counters
    int accepts[6] = {0}, attempts[6] = {0};

    int local_steps = 0, local_accepted = 0;
    int sample_idx = 0;
    const int total_steps = n_samples + n_burn_in;

    // Pre-generate random numbers in batches for speed
    const int batch_size = 256;
    double rand_uniform[batch_size];
    double rand_normal[batch_size * 4];
    int rand_idx_u = batch_size, rand_idx_n = batch_size * 4;

    #define NEXT_UNIFORM() (rand_idx_u >= batch_size ? \
        (rand_idx_u = 0, [&]{ for(int i=0;i<batch_size;++i) rand_uniform[i] = uniform(gen); }(), rand_uniform[rand_idx_u++]) : \
        rand_uniform[rand_idx_u++])
    
    #define NEXT_NORMAL() (rand_idx_n >= batch_size*4 ? \
        (rand_idx_n = 0, [&]{ for(int i=0;i<batch_size*4;++i) rand_normal[i] = normal(gen); }(), rand_normal[rand_idx_n++]) : \
        rand_normal[rand_idx_n++])

    for (int step = 0; step < total_steps; ++step) {
        // Abort checkpoint. A relaxed load of a line that is only written on
        // abort costs nothing next to a step, and the collected samples are
        // discarded by the caller, so bailing out here is safe.
        if ((step & 63) == 0 && abort_requested.load(std::memory_order_relaxed)) break;

        ++local_steps;

        // Update all temperature chains
        for (int t = 0; t < n_temps; ++t) {
            ++attempts[t];

            // Propose new state
            double new_wl = chains[t].wl + NEXT_NORMAL() * step_wl[t];
            double new_sp = chains[t].sp + NEXT_NORMAL() * step_sp[t];
            double new_qu = chains[t].qu + NEXT_NORMAL() * step_qu[t];
            double new_cu = chains[t].cu + NEXT_NORMAL() * step_cu[t];

            // Reflect into bounds (inline)
            #define REFLECT(val, lo, hi) { \
                while (val < lo || val > hi) { \
                    if (val < lo) val = 2.0*lo - val; \
                    if (val > hi) val = 2.0*hi - val; \
                } \
            }
            REFLECT(new_wl, wl_lo, wl_hi);
            REFLECT(new_sp, spacing_lo, spacing_hi);
            REFLECT(new_qu, quad_lo, quad_hi);
            REFLECT(new_cu, cub_lo, cub_hi);
            #undef REFLECT

            double new_score = compute_objective(new_cu, new_qu, new_sp, new_wl, ctx);

            // Metropolis acceptance
            double delta = new_score - chains[t].score;
            
            // exp(-delta / T) = exp(-delta * inv_T)
            // Fast accept if better or probabilistic if worse
            bool accept = (delta <= 0.0) || 
                          (NEXT_UNIFORM() < exp(-delta * inv_temperatures[t]));

            if (accept) {
                chains[t].wl = new_wl;
                chains[t].sp = new_sp;
                chains[t].qu = new_qu;
                chains[t].cu = new_cu;
                chains[t].score = new_score;
                ++accepts[t];
                if (t == 0) ++local_accepted;
            }
        }

        // Replica exchange every 5 steps
        if ((step & 3) == 0) {  // step % 4 == 0, slightly more frequent
            int t = (int)(NEXT_UNIFORM() * (n_temps - 1));
            
            double delta_beta = inv_temperatures[t] - inv_temperatures[t + 1];
            double delta_E = chains[t].score - chains[t + 1].score;
            
            if (NEXT_UNIFORM() < exp(delta_beta * delta_E)) {
                Chain tmp = chains[t];
                chains[t] = chains[t + 1];
                chains[t + 1] = tmp;
            }
        }

        // Adapt step sizes during burn-in
        if (step < n_burn_in && step > 0 && (step & 1023) == 0) {  // step % 1024
            for (int t = 0; t < n_temps; ++t) {
                if (attempts[t] > 100) {
                    double rate = (double)accepts[t] / attempts[t];
                    double target = 0.234 + 0.1 * t / (n_temps - 1);  // Higher for hot chains
                    
                    double factor = rate / target;
                    if (factor < 0.5) factor = 0.5;
                    if (factor > 2.0) factor = 2.0;

                    step_wl[t] *= factor;
                    step_sp[t] *= factor;
                    step_qu[t] *= factor;
                    step_cu[t] *= factor;

                    // Clamp
                    double max_frac = 0.3, min_frac = 1e-6;
                    #define CLAMP_STEP(s, lo, hi) s = fmax(min_frac*(hi-lo), fmin(max_frac*(hi-lo), s))
                    CLAMP_STEP(step_wl[t], wl_lo, wl_hi);
                    CLAMP_STEP(step_sp[t], spacing_lo, spacing_hi);
                    CLAMP_STEP(step_qu[t], quad_lo, quad_hi);
                    CLAMP_STEP(step_cu[t], cub_lo, cub_hi);
                    #undef CLAMP_STEP

                    accepts[t] = 0;
                    attempts[t] = 0;
                }
            }
        }

        // Record sample from cold chain
        if (step >= n_burn_in) {
            out_wl[sample_idx] = chains[0].wl;
            out_spacing[sample_idx] = chains[0].sp;
            out_quad[sample_idx] = chains[0].qu;
            out_cub[sample_idx] = chains[0].cu;
            out_score[sample_idx] = chains[0].score;
            ++sample_idx;
        }

        // Progress update (less frequent)
        if (local_steps >= 2000) {
            global_steps_completed.fetch_add(local_steps, std::memory_order_relaxed);
            global_steps_accepted.fetch_add(local_accepted, std::memory_order_relaxed);
            local_steps = 0;
            local_accepted = 0;
        }
    }

    #undef NEXT_UNIFORM
    #undef NEXT_NORMAL

    global_steps_completed.fetch_add(local_steps, std::memory_order_relaxed);
    global_steps_accepted.fetch_add(local_accepted, std::memory_order_relaxed);
    threads_running.fetch_sub(1);
}

// ============================================================================
// PROGRESS REPORTER
// ============================================================================

void progress_reporter(int total_threads, long long total_steps) {
    const char* RESET = "\033[0m";
    const char* GREEN = "\033[32m";
    const char* YELLOW = "\033[33m";
    const char* RED = "\033[31m";

    while (!stop_reporter.load()) {
        long long steps = global_steps_completed.load();
        long long accepted = global_steps_accepted.load();
        int running = threads_running.load();

        double acc_rate = (steps > 0) ? 100.0 * (double)accepted / (double)steps : 0.0;

        const char* color;
        if (acc_rate >= 15.0 && acc_rate <= 35.0) color = GREEN;
        else if (acc_rate >= 8.0 && acc_rate <= 50.0) color = YELLOW;
        else color = RED;

        fprintf(stderr,
                "\r\033[K[Wavelength Solver] Threads %2d/%2d | Steps %10lld/%10lld | Accept %s%5.1f%%%s",
                running, total_threads, steps, total_steps, color, acc_rate, RESET);
        fflush(stderr);

        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }

    long long steps = global_steps_completed.load();
    long long accepted = global_steps_accepted.load();
    double acc_rate = (steps > 0) ? 100.0 * (double)accepted / (double)steps : 0.0;

    fprintf(stderr,
            "\r\033[K[Wavelength Solver] Threads %2d/%2d | Steps %10lld/%10lld | Accept %5.1f%%\n",
            threads_running.load(), total_threads, steps, total_steps, acc_rate);
    fflush(stderr);
}

// ============================================================================
// HISTOGRAM COMPUTATION
// ============================================================================

// Reduce the per-thread sample buffers straight into the output histograms.
//
// Consumes thread_samples in place: at production settings (16 threads x 2.5M
// samples) gathering them into one contiguous array first costs an extra 1.6 GB
// and a full copy, for nothing. Two parallel passes replace the old six serial
// ones, and the per-sample weights are formed on the fly instead of being
// materialised into another n_total-sized array.
//
// Layout: thread_samples[t] holds 5 columns of n_samples doubles - the four
// parameters followed by the score.
static void compute_weighted_histograms(
    const vector<unique_ptr<double[]>>& thread_samples,
    int n_threads, long long n_samples,
    int nbins,
    double* hist_output) {

    // ---- Pass 1: data ranges ----------------------------------------------
    double min_score = HUGE_VAL, max_score = -HUGE_VAL;
    double pmin[4], pmax[4];
    for (int p = 0; p < 4; ++p) { pmin[p] = HUGE_VAL; pmax[p] = -HUGE_VAL; }

    #pragma omp parallel
    {
        double l_min_s = HUGE_VAL, l_max_s = -HUGE_VAL;
        double l_min[4], l_max[4];
        for (int p = 0; p < 4; ++p) { l_min[p] = HUGE_VAL; l_max[p] = -HUGE_VAL; }

        #pragma omp for schedule(static)
        for (int t = 0; t < n_threads; ++t) {
            const double* base = thread_samples[t].get();
            const double* sc = base + 4 * n_samples;
            for (long long j = 0; j < n_samples; ++j) {
                const double s = sc[j];
                if (s < l_min_s) l_min_s = s;
                if (s > l_max_s) l_max_s = s;
            }
            for (int p = 0; p < 4; ++p) {
                const double* col = base + (long long)p * n_samples;
                for (long long j = 0; j < n_samples; ++j) {
                    const double v = col[j];
                    if (v < l_min[p]) l_min[p] = v;
                    if (v > l_max[p]) l_max[p] = v;
                }
            }
        }

        #pragma omp critical
        {
            if (l_min_s < min_score) min_score = l_min_s;
            if (l_max_s > max_score) max_score = l_max_s;
            for (int p = 0; p < 4; ++p) {
                if (l_min[p] < pmin[p]) pmin[p] = l_min[p];
                if (l_max[p] > pmax[p]) pmax[p] = l_max[p];
            }
        }
    }

    double score_range = max_score - min_score;
    if (score_range < 1e-10) score_range = 1.0;

    double data_min[4], range[4];
    for (int p = 0; p < 4; ++p) {
        double dmin = pmin[p], dmax = pmax[p];
        double r = dmax - dmin;
        if (r < 1e-15) {
            r = fabs(dmin) * 1e-6;
            if (r < 1e-15) r = 1e-15;
            dmin -= r / 2;
            dmax += r / 2;
            r = dmax - dmin;
        }
        data_min[p] = dmin;
        range[p] = r;
    }

    // ---- Pass 2: weighted binning -----------------------------------------
    // Weights: exp(-normalized_score * 5) to strongly favor low scores
    vector<double> hist((size_t)4 * nbins, 0.0);
    double weight_sum = 0.0;

    double bin_scale[4];
    for (int p = 0; p < 4; ++p) bin_scale[p] = (double)nbins / range[p];

    #pragma omp parallel
    {
        vector<double> l_hist((size_t)4 * nbins, 0.0);
        double l_wsum = 0.0;

        #pragma omp for schedule(static)
        for (int t = 0; t < n_threads; ++t) {
            const double* base = thread_samples[t].get();
            const double* sc = base + 4 * n_samples;
            for (long long j = 0; j < n_samples; ++j) {
                const double w = exp(-((sc[j] - min_score) / score_range) * 5.0);
                l_wsum += w;
                for (int p = 0; p < 4; ++p) {
                    const double* col = base + (long long)p * n_samples;
                    int bin_idx = (int)((col[j] - data_min[p]) * bin_scale[p]);
                    bin_idx = max(0, min(nbins - 1, bin_idx));
                    l_hist[(size_t)p * nbins + bin_idx] += w;
                }
            }
        }

        #pragma omp critical
        {
            weight_sum += l_wsum;
            for (size_t k = 0; k < l_hist.size(); ++k) hist[k] += l_hist[k];
        }
    }

    // Normalising here rather than per sample is equivalent and much cheaper.
    const double inv_wsum = (weight_sum > 0.0) ? 1.0 / weight_sum : 1.0;

    for (int p = 0; p < 4; ++p) {
        double* out_ptr = hist_output + p * 2 * nbins;
        for (int i = 0; i < nbins; ++i) {
            out_ptr[i] = data_min[p] + (i + 0.5) * range[p] / nbins;
            out_ptr[nbins + i] = hist[(size_t)p * nbins + i] * inv_wsum;
        }
    }
}

// ============================================================================
// EXTERNAL C API (signatures unchanged)
// ============================================================================

extern "C" {

int get_num_threads() {
    return omp_get_num_procs();
}

// --- Cooperative abort -------------------------------------------------------
// request_abort() may be called from any thread (typically the GUI thread while
// run_mcmc() is executing on a worker thread). run_mcmc() then returns early
// without touching hist_output; the caller is expected to discard the results.

void request_abort() {
    abort_requested.store(true, std::memory_order_relaxed);
}

void reset_abort() {
    abort_requested.store(false, std::memory_order_relaxed);
}

int is_aborted() {
    return abort_requested.load(std::memory_order_relaxed) ? 1 : 0;
}

int get_histogram_nbins(int n_samples) {
    int n_threads = omp_get_num_procs();
    long long n_total = (long long)n_threads * (long long)n_samples;
    return (int)ceil(2.0 * pow((double)n_total, 1.0 / 3.0));
}

void run_mcmc(const double* compspec_x, const double* compspec_y, int compspec_size,
              const double* lines, int lines_size,
              int n_samples, double wl_start, double spacing,
              double quadratic_fac, double cubic_fac,
              double wl_stepsize, double spacing_stepsize,
              double quad_stepsize, double cub_stepsize,
              double wl_cov, double spacing_cov,
              double quad_cov, double cub_cov,
              double acc_param,
              double* hist_output, int nbins, int* out_n_threads) {

    int n_threads = omp_get_num_procs();
    *out_n_threads = n_threads;

    // Abort may already have been requested before we got here.
    if (abort_requested.load(std::memory_order_relaxed)) {
        fprintf(stderr, "[Wavelength Solver] Aborted before start.\n");
        return;
    }

    global_steps_completed.store(0);
    global_steps_accepted.store(0);
    threads_running.store(0);
    stop_reporter.store(false);

    double wl_lo = wl_start - wl_cov / 2;
    double wl_hi = wl_start + wl_cov / 2;
    double spacing_lo = spacing - spacing_cov / 2;
    double spacing_hi = spacing + spacing_cov / 2;
    double quad_lo = quadratic_fac - quad_cov / 2;
    double quad_hi = quadratic_fac + quad_cov / 2;
    double cub_lo = cubic_fac - cub_cov / 2;
    double cub_hi = cubic_fac + cub_cov / 2;

    int n_burn_in = min(500000, max(50000, n_samples / 4));

    // One flat buffer per thread: 5 columns (wl, spacing, quad, cub, score) of
    // n_samples each. Left uninitialised on purpose - every element is written
    // before it is read, and this is ~1.6 GB at production settings.
    vector<unique_ptr<double[]>> thread_samples(n_threads);
    for (int i = 0; i < n_threads; ++i)
        thread_samples[i].reset(new double[(size_t)5 * (size_t)n_samples]);

    long long total_steps = (long long)n_threads * (long long)(n_samples + n_burn_in);
    std::thread reporter(progress_reporter, n_threads, total_steps);

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n_threads; ++i) {
        random_device rd;
        mt19937 gen(rd() ^ (i * 999983 + 12345));

        vector<double> compspec_x_vec(compspec_x, compspec_x + compspec_size);
        vector<double> compspec_y_vec(compspec_y, compspec_y + compspec_size);
        vector<double> lines_vec(lines, lines + lines_size);

        double* base = thread_samples[i].get();
        run_parallel_tempering_chain(
            compspec_x_vec, compspec_y_vec, lines_vec,
            n_samples, n_burn_in,
            wl_start, spacing, quadratic_fac, cubic_fac,
            wl_lo, wl_hi, spacing_lo, spacing_hi,
            quad_lo, quad_hi, cub_lo, cub_hi,
            base,
            base + (size_t)n_samples,
            base + (size_t)2 * n_samples,
            base + (size_t)3 * n_samples,
            base + (size_t)4 * n_samples,
            gen);
    }

    stop_reporter.store(true);
    reporter.join();

    if (abort_requested.load(std::memory_order_relaxed)) {
        // The chains stopped early, so the samples are incomplete: leave
        // hist_output untouched and let the caller raise.
        fprintf(stderr, "[Wavelength Solver] Aborted by user.\n");
        return;
    }

    // Report best
    double best_score = HUGE_VAL;
    int best_t = 0;
    long long best_j = 0;
    for (int t = 0; t < n_threads; ++t) {
        const double* sc = thread_samples[t].get() + (size_t)4 * n_samples;
        for (long long j = 0; j < n_samples; ++j) {
            if (sc[j] < best_score) { best_score = sc[j]; best_t = t; best_j = j; }
        }
    }

    {
        const double* b = thread_samples[best_t].get();
        fprintf(stderr, "[Wavelength Solver] Best: WL=%.4f Sp=%.8f Qu=%.4e Cu=%.4e Score=%.2f\n",
                b[best_j], b[(size_t)n_samples + best_j],
                b[(size_t)2 * n_samples + best_j], b[(size_t)3 * n_samples + best_j],
                best_score);
    }

    compute_weighted_histograms(thread_samples, n_threads, n_samples, nbins, hist_output);

    fprintf(stderr, "[Wavelength Solver] Done.\n");
}

} // extern "C"