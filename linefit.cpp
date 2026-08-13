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

std::mutex cout_mutex;
std::atomic<long long> global_steps_completed(0);
std::atomic<long long> global_steps_accepted(0);
std::atomic<int> threads_running(0);
std::atomic<bool> stop_reporter(false);

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

// Robust inverse using Newton-Raphson with bisection fallback
double inverse_cub_poly(double y, double a, double b, double c, double d) {
    // Solve: a*x³ + b*x² + c*x + d = y
    // Rearrange to: a*x³ + b*x² + c*x + (d - y) = 0
    
    double d_shifted = d - y;
    
    // Handle degenerate cases first
    double abs_a = fabs(a);
    double abs_b = fabs(b);
    double abs_c = fabs(c);
    
    // Essentially linear (most common case for small cubic/quadratic)
    if (abs_a < 1e-20 && abs_b < 1e-15) {
        if (abs_c < 1e-20) return NAN;
        return -d_shifted / c;
    }
    
    // Essentially quadratic
    if (abs_a < 1e-20) {
        double disc = c * c - 4.0 * b * d_shifted;
        if (disc < 0) return NAN;
        double sqrt_disc = sqrt(disc);
        // Return root closest to linear solution
        double r1 = (-c + sqrt_disc) / (2.0 * b);
        double r2 = (-c - sqrt_disc) / (2.0 * b);
        double linear_approx = -d_shifted / c;
        return (fabs(r1 - linear_approx) < fabs(r2 - linear_approx)) ? r1 : r2;
    }
    
    // Full cubic case: use Newton-Raphson with good initial guess and bracketing
    
    // For wavelength calibration, x is typically in [0, 1] (normalized pixel)
    // The polynomial should be monotonic, so we can bracket the root
    
    // Initial guess: use linear approximation
    double x;
    if (abs_c > 1e-20) {
        x = -d_shifted / c;  // Linear approximation
    } else if (abs_b > 1e-20) {
        // Quadratic approximation
        double disc = -d_shifted / b;
        x = (disc > 0) ? sqrt(disc) : 0.5;
    } else {
        x = cbrt(-d_shifted / a);
    }
    
    // Clamp initial guess to reasonable range
    x = fmax(-10.0, fmin(10.0, x));
    
    // Newton-Raphson with safety checks
    double f, fp;
    double prev_x = x;
    double best_x = x;
    eval_poly_and_deriv(x, a, b, c, d_shifted, f, fp);
    double best_f = fabs(f);
    
    for (int iter = 0; iter < 30; ++iter) {
        eval_poly_and_deriv(x, a, b, c, d_shifted, f, fp);
        
        // Track best solution
        if (fabs(f) < best_f) {
            best_f = fabs(f);
            best_x = x;
        }
        
        // Check convergence
        if (fabs(f) < 1e-14 * (fabs(y) + 1.0)) {
            return x;
        }
        
        // Newton step
        if (fabs(fp) < 1e-30) {
            // Derivative too small, perturb and continue
            x += (f > 0 ? -1e-6 : 1e-6);
            continue;
        }
        
        double dx = f / fp;
        
        // Damping for large steps (prevents oscillation)
        if (fabs(dx) > 1.0) {
            dx = (dx > 0) ? 1.0 : -1.0;
        }
        
        prev_x = x;
        x -= dx;
        
        // Check for convergence
        if (fabs(dx) < 1e-14 * (fabs(x) + 1e-10)) {
            return x;
        }
        
        // Check for oscillation (stuck between two values)
        if (iter > 5 && fabs(x - prev_x) < 1e-13) {
            break;
        }
    }
    
    // If Newton didn't converge well, try bisection in a reasonable range
    // Find brackets by scanning
    double x_lo = -2.0, x_hi = 2.0;
    double f_lo, f_hi, fp_dummy;
    eval_poly_and_deriv(x_lo, a, b, c, d_shifted, f_lo, fp_dummy);
    eval_poly_and_deriv(x_hi, a, b, c, d_shifted, f_hi, fp_dummy);
    
    // Expand range if needed to bracket root
    for (int i = 0; i < 10 && f_lo * f_hi > 0; ++i) {
        x_lo *= 2.0;
        x_hi *= 2.0;
        eval_poly_and_deriv(x_lo, a, b, c, d_shifted, f_lo, fp_dummy);
        eval_poly_and_deriv(x_hi, a, b, c, d_shifted, f_hi, fp_dummy);
    }
    
    if (f_lo * f_hi > 0) {
        // Couldn't bracket - return best Newton result
        return best_x;
    }
    
    // Bisection
    for (int iter = 0; iter < 60; ++iter) {
        double x_mid = 0.5 * (x_lo + x_hi);
        double f_mid;
        eval_poly_and_deriv(x_mid, a, b, c, d_shifted, f_mid, fp_dummy);
        
        if (fabs(f_mid) < 1e-14 * (fabs(y) + 1.0)) {
            return x_mid;
        }
        
        if (f_mid * f_lo < 0) {
            x_hi = x_mid;
            f_hi = f_mid;
        } else {
            x_lo = x_mid;
            f_lo = f_mid;
        }
        
        if (x_hi - x_lo < 1e-14 * (fabs(x_mid) + 1e-10)) {
            return x_mid;
        }
    }
    
    return 0.5 * (x_lo + x_hi);
}

// ============================================================================
// OBJECTIVE FUNCTION - Robust version using trimmed mean + match counting
// ============================================================================

struct ObjectiveResult {
    double score;           // Primary score (lower is better)
    int n_matched;          // Number of lines with value > threshold
    double trimmed_chisq;   // Trimmed chi-squared (ignores worst outliers)
};


double compute_objective(double cubic_fac, double quadratic_fac, double spacing,
                         double wl_start, const vector<double>& lines,
                         const vector<double>& compspec_x,
                         const vector<double>& compspec_y) {
    const int n_lines = (int)lines.size();
    const int n_comp = (int)compspec_x.size();
    const double x_lo = compspec_x[0];
    const double x_hi = compspec_x[n_comp - 1];
    const double x_range_inv = 1.0 / (x_hi - x_lo);
    
    // Thread-local static buffer to avoid repeated allocation
    thread_local vector<double> residuals;
    residuals.resize(n_lines);
    
    int valid_count = 0;
    
    for (int i = 0; i < n_lines; ++i) {
        double line_wl = lines[i];
        double tl = inverse_cub_poly(line_wl, cubic_fac, quadratic_fac, spacing, wl_start);
        
        if (tl != tl || tl < x_lo || tl > x_hi) {  // NaN check: tl != tl
            residuals[valid_count++] = 1.0;
            continue;
        }
        
        // Fast interpolation using normalized position
        double t = (tl - x_lo) * x_range_inv * (n_comp - 1);
        int j = (int)t;
        if (j < 0) j = 0;
        if (j >= n_comp - 1) j = n_comp - 2;
        
        double frac = t - j;
        double y = compspec_y[j] + frac * (compspec_y[j + 1] - compspec_y[j]);
        
        double resid = 1.0 - y;
        residuals[valid_count++] = resid * resid;
    }
    
    if (valid_count == 0) return 1e6;
    
    // Use nth_element for O(n) trimmed sum instead of O(n log n) sort
    int n_keep = (valid_count * 4 + 4) / 5;  // Keep 80%, round up
    if (n_keep < 1) n_keep = 1;
    
    // Partial sort: elements before nth are <= nth, elements after are >= nth
    nth_element(residuals.begin(), residuals.begin() + n_keep, residuals.begin() + valid_count);
    
    double sum = 0.0;
    for (int i = 0; i < n_keep; ++i) {
        sum += residuals[i];
    }
    
    return sum * 1000.0 / n_keep;
}


inline double interpolate_lines_chisq(double cubic_fac, double quadratic_fac, double spacing,
                                       double wl_start, const vector<double>& lines,
                                       const vector<double>& compspec_x,
                                       const vector<double>& compspec_y) {
    return compute_objective(cubic_fac, quadratic_fac, spacing, wl_start,
                            lines, compspec_x, compspec_y);
}

// ============================================================================
// PARALLEL TEMPERING MCMC
// ============================================================================

struct StepSizes {
    double wl;
    double spacing;
    double quad;
    double cub;
};

struct ChainState {
    double wl;
    double spacing;
    double quad;
    double cub;
    double score;
};

// Reflect value into bounds
inline double reflect_into_bounds(double val, double lo, double hi) {
    if (val >= lo && val <= hi) return val;
    double range = hi - lo;
    if (range <= 0) return (lo + hi) / 2;
    
    while (val < lo || val > hi) {
        if (val < lo) val = 2.0 * lo - val;
        if (val > hi) val = 2.0 * hi - val;
    }
    return val;
}

// Propose new state with Gaussian steps, occasionally large jumps
ChainState propose_state(const ChainState& current, const StepSizes& steps,
                         double wl_lo, double wl_hi,
                         double spacing_lo, double spacing_hi,
                         double quad_lo, double quad_hi,
                         double cub_lo, double cub_hi,
                         double temperature, mt19937& gen) {
    normal_distribution<double> normal(0.0, 1.0);
    uniform_real_distribution<double> uniform(0.0, 1.0);

    // Scale step sizes by sqrt(temperature) for better mixing at high T
    double temp_scale = sqrt(temperature);

    // Occasionally make a large jump (10% of the time)
    double jump_scale = (uniform(gen) < 0.1) ? 5.0 : 1.0;

    ChainState proposed;
    proposed.wl = current.wl + normal(gen) * steps.wl * temp_scale * jump_scale;
    proposed.spacing = current.spacing + normal(gen) * steps.spacing * temp_scale * jump_scale;
    proposed.quad = current.quad + normal(gen) * steps.quad * temp_scale * jump_scale;
    proposed.cub = current.cub + normal(gen) * steps.cub * temp_scale * jump_scale;

    // Reflect into bounds
    proposed.wl = reflect_into_bounds(proposed.wl, wl_lo, wl_hi);
    proposed.spacing = reflect_into_bounds(proposed.spacing, spacing_lo, spacing_hi);
    proposed.quad = reflect_into_bounds(proposed.quad, quad_lo, quad_hi);
    proposed.cub = reflect_into_bounds(proposed.cub, cub_lo, cub_hi);

    return proposed;
}

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
                                            chains[t].sp, chains[t].wl,
                                            lines, compspec_x, compspec_y);
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

            double new_score = compute_objective(new_cu, new_qu, new_sp, new_wl,
                                                 lines, compspec_x, compspec_y);

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

static void compute_weighted_histograms(
    const vector<vector<double>>& all_samples,
    int nbins,
    double* hist_output) {

    long long n_total = (long long)all_samples[0].size();
    const vector<double>& scores = all_samples[4];

    // Find score range for normalization
    double min_score = *min_element(scores.begin(), scores.end());
    double max_score = *max_element(scores.begin(), scores.end());
    double score_range = max_score - min_score;
    if (score_range < 1e-10) score_range = 1.0;

    // Weights: exp(-normalized_score * 5) to strongly favor low scores
    vector<double> weights(n_total);
    double weight_sum = 0.0;
    for (long long j = 0; j < n_total; ++j) {
        double norm_score = (scores[j] - min_score) / score_range;
        weights[j] = exp(-norm_score * 5.0);
        weight_sum += weights[j];
    }
    for (long long j = 0; j < n_total; ++j) {
        weights[j] /= weight_sum;
    }

    for (int p = 0; p < 4; ++p) {
        const vector<double>& col = all_samples[p];

        double data_min = *min_element(col.begin(), col.end());
        double data_max = *max_element(col.begin(), col.end());
        double range = data_max - data_min;
        if (range < 1e-15) {
            range = fabs(data_min) * 1e-6;
            if (range < 1e-15) range = 1e-15;
            data_min -= range / 2;
            data_max += range / 2;
            range = data_max - data_min;
        }

        vector<double> hist(nbins, 0.0);
        for (long long j = 0; j < n_total; ++j) {
            int bin_idx = (int)((col[j] - data_min) / range * (double)nbins);
            bin_idx = max(0, min(nbins - 1, bin_idx));
            hist[bin_idx] += weights[j];
        }

        double* out_ptr = hist_output + p * 2 * nbins;
        for (int i = 0; i < nbins; ++i) {
            out_ptr[i] = data_min + (i + 0.5) * range / nbins;
            out_ptr[nbins + i] = hist[i];
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

    vector<vector<vector<double>>> thread_samples(n_threads,
        vector<vector<double>>(5, vector<double>(n_samples, 0.0)));

    long long total_steps = (long long)n_threads * (long long)(n_samples + n_burn_in);
    std::thread reporter(progress_reporter, n_threads, total_steps);

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n_threads; ++i) {
        random_device rd;
        mt19937 gen(rd() ^ (i * 999983 + 12345));

        vector<double> compspec_x_vec(compspec_x, compspec_x + compspec_size);
        vector<double> compspec_y_vec(compspec_y, compspec_y + compspec_size);
        vector<double> lines_vec(lines, lines + lines_size);

        run_parallel_tempering_chain(
            compspec_x_vec, compspec_y_vec, lines_vec,
            n_samples, n_burn_in,
            wl_start, spacing, quadratic_fac, cubic_fac,
            wl_lo, wl_hi, spacing_lo, spacing_hi,
            quad_lo, quad_hi, cub_lo, cub_hi,
            thread_samples[i][0].data(),
            thread_samples[i][1].data(),
            thread_samples[i][2].data(),
            thread_samples[i][3].data(),
            thread_samples[i][4].data(),
            gen);
    }

    stop_reporter.store(true);
    reporter.join();

    long long n_total = (long long)n_threads * (long long)n_samples;
    vector<vector<double>> all_samples(5, vector<double>(n_total));

    for (int t = 0; t < n_threads; ++t) {
        long long offset = (long long)t * (long long)n_samples;
        for (int c = 0; c < 5; ++c) {
            memcpy(all_samples[c].data() + offset,
                   thread_samples[t][c].data(),
                   n_samples * sizeof(double));
        }
    }

    thread_samples.clear();
    thread_samples.shrink_to_fit();

    // Report best
    double best_score = all_samples[4][0];
    long long best_idx = 0;
    for (long long j = 1; j < n_total; ++j) {
        if (all_samples[4][j] < best_score) {
            best_score = all_samples[4][j];
            best_idx = j;
        }
    }

    fprintf(stderr, "[Wavelength Solver] Best: WL=%.4f Sp=%.8f Qu=%.4e Cu=%.4e Score=%.2f\n",
            all_samples[0][best_idx], all_samples[1][best_idx],
            all_samples[2][best_idx], all_samples[3][best_idx], best_score);

    compute_weighted_histograms(all_samples, nbins, hist_output);

    fprintf(stderr, "[Wavelength Solver] Done.\n");
}

} // extern "C"