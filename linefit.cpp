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

random_device rd;

double cub_poly(double x, double a, double b, double c, double d){
    return x*x*x*a + x*x*b + x*c + d;
}

double inverse_cub_poly(double x, double a, double b, double c, double d)
{
    b /= a;
    c /= a;
    d = (d - x) / a;
    double q = (3.0 * c - b * b) / 9.0;
    double r = (9.0 * b * c - 27.0 * d - 2.0 * b * b * b) / 54.0;
    double disc = q * q * q + r * r;
    double root_b = b / 3.0;

    double r0, r1, r2;

    if (disc > 0) {
        double s = r + sqrt(disc);
        s = (s < 0) ? -cbrt(-s) : cbrt(s);
        double t = r - sqrt(disc);
        t = (t < 0) ? -cbrt(-t) : cbrt(t);
        r0 = r1 = r2 = -root_b + s + t;
    } else if (disc == 0) {
        double r13 = (r < 0) ? -cbrt(-r) : cbrt(r);
        r0 = -root_b + 2.0 * r13;
        r1 = r2 = -root_b - r13;
    } else {
        double neg_q = -q;
        double dum1 = neg_q * neg_q * neg_q;
        dum1 = acos(r / sqrt(dum1));
        double r13 = 2.0 * sqrt(neg_q);
        r0 = -root_b + r13 * cos(dum1 / 3.0);
        r1 = -root_b + r13 * cos((dum1 + 2.0 * M_PI) / 3.0);
        r2 = -root_b + r13 * cos((dum1 + 4.0 * M_PI) / 3.0);
    }

    // Sort 3 elements — get the median with no branches
    if (r0 > r1) { double tmp = r0; r0 = r1; r1 = tmp; }
    if (r1 > r2) { double tmp = r1; r1 = r2; r2 = tmp; }
    if (r0 > r1) { double tmp = r0; r0 = r1; r1 = tmp; }

    return r1;
}

double interpolate_lines_chisq(double cubic_fac, double quadratic_fac, double spacing,
                               double wl_start, const vector<double>& lines,
                               const vector<double>& compspec_x,
                               const vector<double>& compspec_y) {
    double sum = 0.0;
    const int n_lines = (int)lines.size();
    const int n_comp = (int)compspec_x.size();
    const double x_lo = compspec_x[0];
    const double x_hi = compspec_x[n_comp - 1];

    for (int i = 0; i < n_lines; ++i) {
        double tl = inverse_cub_poly(lines[i], cubic_fac, quadratic_fac, spacing, wl_start);

        if (tl < x_lo || tl > x_hi) {
            sum += 1.0;
            continue;
        }

        auto it = lower_bound(compspec_x.begin(), compspec_x.end(), tl);
        int j = (int)(it - compspec_x.begin());

        double y;
        if (j == 0 || j >= n_comp - 1) {
            y = 0.0;
        } else {
            double x0 = compspec_x[j - 1];
            double x1 = compspec_x[j];
            double y0 = compspec_y[j - 1];
            double y1 = compspec_y[j];
            y = y0 + (y1 - y0) * ((tl - x0) / (x1 - x0));
        }
        sum += (y - 1.0) * (y - 1.0);
    }
    return (sum != 0.0) ? sum : 1000000.0;
}

// Optimized version — same distribution, fewer RNG calls, faster math
double levyRejectionSampling(double c, mt19937& local_gen) {
    // Pre-compute constant factor outside the loop
    const double c_over_2pi = c * (1.0 / (2.0 * M_PI));
    const double half_c = 0.5 * c;

    uniform_real_distribution<double> u_dist(0.0, 1.0);
    normal_distribution<double> n_dist(0.0, 1.0);

    while (true) {
        double v = n_dist(local_gen);
        double v2 = v * v;
        
        // x_candidate = c / v² (since mu is always 0 in your calls)
        double x_cand = c / v2;
        double x_cand_inv = v2 / c;  // 1/x_cand, avoids division later

        // p = sqrt(c/(2π)) * exp(-c/(2*x)) / x^(3/2)
        // Rewrite to minimize expensive operations:
        //   x^(3/2) = x * sqrt(x)
        //   p = sqrt(c/(2π)) * exp(-c/(2x)) * x_inv * x_inv_sqrt
        //   p = sqrt(c/(2π) * x_inv³) * exp(-c/(2x))
        //   p = sqrt(c/(2π)) * (v²/c)^(3/2) * exp(-half_c * v²/c)
        //     = sqrt(c/(2π)) * v³/(c^(3/2)) * exp(-v²/2)
        //     = v³ * exp(-v²/2) / (c * sqrt(2π))
        //
        // But v can be negative, and v³ would be negative, so use |v|:
        double abs_v = fabs(v);
        double abs_v2 = abs_v * abs_v;
        double p = abs_v * abs_v2 * exp(-0.5 * abs_v2) / (c * sqrt(2.0 * M_PI));

        // Single uniform for accept, single uniform for sign
        // We can combine: use one uniform [0,1) for accept,
        // and the sign of the original v for the flip (v is symmetric around 0)
        // BUT that changes the correlation structure. Keep two draws to be safe:
        double u2 = u_dist(local_gen);
        if (u2 <= p) {
            // Sign flip: use another uniform
            double u3 = u_dist(local_gen);
            return (u3 > 0.5) ? x_cand : -x_cand;
        }
    }
}
void progress_reporter(int total_threads, long long total_steps_per_thread, int n_burn_in) {
    long long total_steps = (long long)total_threads * ((long long)total_steps_per_thread + (long long)n_burn_in);

    const char* RESET  = "\033[0m";
    const char* GREEN  = "\033[32m";
    const char* YELLOW = "\033[33m";
    const char* RED    = "\033[31m";

    while (!stop_reporter.load()) {
        long long steps = global_steps_completed.load();
        long long accepted = global_steps_accepted.load();
        int running = threads_running.load();
        double acc_rate = (steps > 0) ? 100.0 * (double)accepted / (double)steps : 0.0;

        const char* color;
        if (acc_rate >= 10.0 && acc_rate <= 30.0) color = GREEN;
        else if (acc_rate >= 5.0 && acc_rate <= 40.0) color = YELLOW;
        else color = RED;

        fprintf(stderr,
            "\r\033[K[Wavelength Solver] Threads running %2d/%2d | Step %10lld/%10lld | Acceptance Rate %s%6.2f%%%s",
            running, total_threads, steps, total_steps, color, acc_rate, RESET);
        fflush(stderr);

        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }

    long long steps = global_steps_completed.load();
    long long accepted = global_steps_accepted.load();
    int running = threads_running.load();
    double acc_rate = (steps > 0) ? 100.0 * (double)accepted / (double)steps : 0.0;

    const char* color;
    if (acc_rate >= 10.0 && acc_rate <= 30.0) color = GREEN;
    else if (acc_rate >= 5.0 && acc_rate <= 40.0) color = YELLOW;
    else color = RED;

    fprintf(stderr,
        "\r\033[K[Wavelength Solver] Threads running %2d/%2d | Step %10lld/%10lld | Acceptance Rate %s%6.2f%%%s\n",
        running, total_threads, steps,
        (long long)total_threads * ((long long)total_steps_per_thread + (long long)n_burn_in),
        color, acc_rate, RESET);
    fflush(stderr);
}

// Per-thread MCMC sampler. Writes results into the provided per-thread output arrays.
void fitlines_mkcmk(vector<double>& compspec_x_vec, vector<double>& compspec_y_vec,
                     vector<double>& lines_vec,
                     int lines_size, int compspec_size, int n_samples, double wl_start,
                     double spacing, double quadratic_fac, double cubic_fac,
                     double wl_stepsize, double spacing_stepsize, double quad_stepsize,
                     double cub_stepsize, double wl_cov, double spacing_cov,
                     double quad_cov, double cub_cov, double acc_param,
                     double* out_wl, double* out_spacing, double* out_quad,
                     double* out_cub, double* out_corr){

    mt19937 local_gen(rd() + omp_get_thread_num() * 12345);
    threads_running.fetch_add(1);

    double this_correlation = interpolate_lines_chisq(cubic_fac, quadratic_fac, spacing, wl_start,
                                                lines_vec, compspec_x_vec, compspec_y_vec);

    const double wl_lo = wl_start - wl_cov / 2;
    const double wl_hi = wl_start + wl_cov / 2;
    const double spacing_lo = spacing - spacing_cov / 2;
    const double spacing_hi = spacing + spacing_cov / 2;
    const double quad_lo = quadratic_fac - quad_cov / 2;
    const double quad_hi = quadratic_fac + quad_cov / 2;
    const double cub_lo = cubic_fac - cub_cov / 2;
    const double cub_hi = cubic_fac + cub_cov / 2;

    double step_st, step_sp, step_quad, step_cub, step_num;
    double next_correlation;
    double nextwl, nextspacing, nextquad, nextcub;

    double m = 1 / (1 - acc_param);
    double n_val = acc_param / (acc_param - 1);

    uniform_real_distribution<> step_dist(0., 1.);

    int n_burn_in = 1000000;

    int local_steps = 0;
    int local_accepted = 0;

    for (int j = 0; j < n_samples + n_burn_in; ++j) {
        local_steps++;
        if (local_steps >= 1000) {
            global_steps_completed.fetch_add(local_steps, std::memory_order_relaxed);
            global_steps_accepted.fetch_add(local_accepted, std::memory_order_relaxed);
            local_steps = 0;
            local_accepted = 0;
        }

        step_num = step_dist(local_gen);

        step_st   = levyRejectionSampling(wl_stepsize, local_gen);
        step_sp   = levyRejectionSampling(spacing_stepsize, local_gen);
        step_quad = levyRejectionSampling(quad_stepsize, local_gen);
        step_cub  = levyRejectionSampling(cub_stepsize, local_gen);

        nextwl = wl_start + step_st;
        nextspacing = spacing + step_sp;
        nextquad = quadratic_fac + step_quad;
        nextcub = cubic_fac + step_cub;

        if (!(wl_lo < nextwl && nextwl < wl_hi && spacing_lo < nextspacing && nextspacing < spacing_hi &&
              quad_lo < nextquad && nextquad < quad_hi && cub_lo < nextcub && nextcub < cub_hi)) {
            if (j >= n_burn_in) {
                int idx = j - n_burn_in;
                out_wl[idx] = wl_start;
                out_spacing[idx] = spacing;
                out_quad[idx] = quadratic_fac;
                out_cub[idx] = cubic_fac;
                out_corr[idx] = this_correlation;
            }
            continue;
        }

        next_correlation = interpolate_lines_chisq(nextcub, nextquad, nextspacing, nextwl,
                                                   lines_vec, compspec_x_vec, compspec_y_vec);

        if (next_correlation < this_correlation) {
            wl_start = nextwl;
            spacing = nextspacing;
            quadratic_fac = nextquad;
            cubic_fac = nextcub;
            this_correlation = next_correlation;
            local_accepted++;
        } else if (step_num < (m * next_correlation / this_correlation) + n_val) {
            wl_start = nextwl;
            spacing = nextspacing;
            quadratic_fac = nextquad;
            cubic_fac = nextcub;
            this_correlation = next_correlation;
            local_accepted++;
        }

        if (j >= n_burn_in) {
            int idx = j - n_burn_in;
            out_wl[idx] = wl_start;
            out_spacing[idx] = spacing;
            out_quad[idx] = quadratic_fac;
            out_cub[idx] = cubic_fac;
            out_corr[idx] = this_correlation;
        }
    }

    // Flush remaining counts
    global_steps_completed.fetch_add(local_steps, std::memory_order_relaxed);
    global_steps_accepted.fetch_add(local_accepted, std::memory_order_relaxed);

    threads_running.fetch_sub(1);
}

// Compute weighted histograms identically to the Python code:
//   nbins = ceil(2 * (n_total ^ (1/3)))
//   weights = 1 / correlation
//   bin_edges = linspace(min, max, nbins+1)
//   bin_idx = int((val - min) / (max - min) * nbins), clamped to nbins-1
//   bin_centers = (edges[:-1] + edges[1:]) / 2
static void compute_weighted_histograms(
    const vector<vector<double>>& all_samples,  // [5][n_total]: wl, spacing, quad, cub, corr
    int nbins,
    double* hist_output)  // layout: [centers0(nbins), hist0(nbins), centers1, hist1, ..., centers3, hist3]
{
    long long n_total = (long long)all_samples[0].size();
    const vector<double>& corr = all_samples[4];

    // Precompute weights = 1/correlation
    vector<double> weights(n_total);
    for (long long j = 0; j < n_total; ++j) {
        weights[j] = 1.0 / corr[j];
    }

    for (int p = 0; p < 4; ++p) {
        const vector<double>& col = all_samples[p];

        // Find min/max
        double data_min = *min_element(col.begin(), col.end());
        double data_max = *max_element(col.begin(), col.end());

        // Compute bin edges: linspace(data_min, data_max, nbins+1)
        vector<double> bin_edges(nbins + 1);
        for (int i = 0; i <= nbins; ++i) {
            bin_edges[i] = data_min + (double)i * (data_max - data_min) / (double)nbins;
        }

        // Accumulate weighted histogram
        vector<double> hist(nbins, 0.0);
        double range = data_max - data_min;
        if (range == 0.0) range = 1.0;  // avoid division by zero

        for (long long j = 0; j < n_total; ++j) {
            int bin_idx = (int)((col[j] - data_min) / range * (double)nbins);
            if (bin_idx >= nbins) bin_idx = nbins - 1;
            if (bin_idx < 0) bin_idx = 0;
            hist[bin_idx] += weights[j];
        }

        // Compute bin centers
        vector<double> bin_centers(nbins);
        for (int i = 0; i < nbins; ++i) {
            bin_centers[i] = (bin_edges[i] + bin_edges[i + 1]) / 2.0;
        }

        // Write to output: [centers, hist] for this parameter
        double* out_ptr = hist_output + p * 2 * nbins;
        for (int i = 0; i < nbins; ++i) {
            out_ptr[i] = bin_centers[i];
        }
        for (int i = 0; i < nbins; ++i) {
            out_ptr[nbins + i] = hist[i];
        }
    }
}

extern "C" {

int get_num_threads() {
    return omp_get_num_procs();
}

// Compute nbins from total sample count: ceil(2 * (n_total^(1/3)))
// n_total = n_threads * n_samples
int get_histogram_nbins(int n_samples) {
    int n_threads = omp_get_num_procs();
    long long n_total = (long long)n_threads * (long long)n_samples;
    return (int)ceil(2.0 * pow((double)n_total, 1.0 / 3.0));
}

// Main entry point.
// hist_output must be pre-allocated with size: 4 * 2 * nbins doubles
// Layout: [param0_centers(nbins), param0_hist(nbins),
//          param1_centers(nbins), param1_hist(nbins),
//          param2_centers(nbins), param2_hist(nbins),
//          param3_centers(nbins), param3_hist(nbins)]
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

    int n_burn_in = 1000000;

    // Allocate per-thread sample storage internally
    // thread_samples[thread][column][sample_idx]
    // columns: 0=wl, 1=spacing, 2=quad, 3=cub, 4=corr
    vector<vector<vector<double>>> thread_samples(n_threads,
        vector<vector<double>>(5, vector<double>(n_samples, 0.0)));

    std::thread reporter(progress_reporter, n_threads, (long long)n_samples, n_burn_in);

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n_threads; ++i) {
        vector<double> compspec_x_vec(compspec_x, compspec_x + compspec_size);
        vector<double> compspec_y_vec(compspec_y, compspec_y + compspec_size);
        vector<double> lines_vec(lines, lines + lines_size);

        fitlines_mkcmk(compspec_x_vec, compspec_y_vec, lines_vec,
                        lines_size, compspec_size, n_samples,
                        wl_start, spacing, quadratic_fac, cubic_fac,
                        wl_stepsize, spacing_stepsize, quad_stepsize, cub_stepsize,
                        wl_cov, spacing_cov, quad_cov, cub_cov, acc_param,
                        thread_samples[i][0].data(),
                        thread_samples[i][1].data(),
                        thread_samples[i][2].data(),
                        thread_samples[i][3].data(),
                        thread_samples[i][4].data());
    }

    stop_reporter.store(true);
    reporter.join();

    // Merge all thread samples into one contiguous array per column
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

    // Free thread storage
    thread_samples.clear();
    thread_samples.shrink_to_fit();

    // Compute histograms
    fprintf(stderr, "[Wavelength Solver] Computing histograms (%d bins, %lld samples)...\n",
            nbins, n_total);
    fflush(stderr);

    compute_weighted_histograms(all_samples, nbins, hist_output);

    fprintf(stderr, "[Wavelength Solver] Done.\n");
    fflush(stderr);
}

} // extern "C"