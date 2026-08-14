"""Benchmark + correctness harness for the C++ wavelength solver.

Builds a synthetic comparison-lamp spectrum from a *known* cubic dispersion
solution, runs run_mcmc() through the same ctypes path MIDIR uses, and reports
both wall time and how well the true solution was recovered.

Usage:
    python bench/bench_mcmc.py <path-to-lib.so> [n_samples] [repeats] [preset|all]
"""
import ctypes
import sys
import time

import numpy as np
from numpy.ctypeslib import ndpointer
from scipy.ndimage import gaussian_filter, maximum_filter
from scipy.optimize import curve_fit

# Per-instrument test cases. Each is a plausible true dispersion solution that
# sits inside the search box the corresponding shipped preset gives the solver,
# so the benchmark measures a problem the solver can actually solve.
# "opts" mirror presets/<name>.json.
PRESETS = {
    "SOAR_930": dict(
        npix=1740, linelist="linelists/SOAR_930_M1_FeAr.txt",
        true=dict(offset=3397.5018518647, linear=0.99, quad=-1.5e-05, cub=5.0424132976e-10),
        opts=dict(extent_zero=1725, quad_zero=-7e-6, cube_zero=0.0,
                  offset_stepsize=0.1, linear_stepsize=0.0001,
                  quad_stepsize=1e-7, cube_stepsize=1e-10,
                  c_cov=150, s_cov=0.1, q_cov=5e-5, cub_cov=3e-8,
                  accept_param=1.05, lampfilterwindow=50)),
    "ALFOSC_grism18": dict(
        npix=2048, linelist="linelists/ALFOSC_grism18_ThAr.txt",
        true=dict(offset=3450.0, linear=0.9277, quad=1.2e-05, cub=-2.0e-12),
        opts=dict(extent_zero=1900.0, quad_zero=0.0, cube_zero=0.0,
                  offset_stepsize=0.01, linear_stepsize=1e-05,
                  quad_stepsize=1e-08, cube_stepsize=1e-10,
                  c_cov=150, s_cov=0.15, q_cov=1e-4, cub_cov=3e-8,
                  accept_param=1.05, lampfilterwindow=50)),
    "EFOSC_grism7": dict(
        npix=1024, linelist="linelists/EFOSC_grism7_HeAr.txt",
        true=dict(offset=3250.0, linear=1.758, quad=5e-05, cub=-1.3e-08),
        opts=dict(extent_zero=1800.0, quad_zero=5e-05, cube_zero=-1.3e-08,
                  offset_stepsize=0.1, linear_stepsize=0.0005,
                  quad_stepsize=1e-06, cube_stepsize=1e-09,
                  c_cov=100, s_cov=0.5, q_cov=1e-3, cub_cov=6e-7,
                  accept_param=1.015, lampfilterwindow=50)),
}

# Filled in by main() from the selected preset.
TRUE = PRESETS["SOAR_930"]["true"]
NPIX = PRESETS["SOAR_930"]["npix"]
LINELIST = PRESETS["SOAR_930"]["linelist"]
OPTS = PRESETS["SOAR_930"]["opts"]


def load(libpath):
    lib = ctypes.CDLL(libpath)
    lib.get_histogram_nbins.restype = ctypes.c_int
    lib.get_histogram_nbins.argtypes = [ctypes.c_int]
    lib.run_mcmc.restype = None
    lib.run_mcmc.argtypes = [
        ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
        ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
        ctypes.c_int,
        ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
        ctypes.c_int, ctypes.c_int,
    ] + [ctypes.c_double] * 13 + [
        ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
        ctypes.c_int, ctypes.POINTER(ctypes.c_int),
    ]
    try:
        lib.reset_abort.restype = None
        lib.reset_abort.argtypes = []
        lib.reset_abort()
    except AttributeError:
        pass
    return lib


def poly(x, off, lin, quad, cub):
    return off + x * (lin + x * (quad + x * cub))


def make_compspec(lines, seed=1234):
    """Synthetic arc lamp: Gaussian emission lines at the true solution."""
    rng = np.random.default_rng(seed)
    pixel = np.arange(NPIX, dtype=np.float64)
    wl = poly(pixel, TRUE["offset"], TRUE["linear"], TRUE["quad"], TRUE["cub"])
    flux = np.zeros(NPIX)
    amps = rng.uniform(0.3, 1.0, size=len(lines))
    for a, lam in zip(amps, lines):
        flux += a * np.exp(-0.5 * ((wl - lam) / 1.6) ** 2)
    flux += 0.02 * rng.standard_normal(NPIX) ** 2  # a little noise floor
    flux *= 1000.0
    # Same preprocessing as data_reduction.call_fitlines_markov
    flux /= np.maximum(maximum_filter(flux, OPTS["lampfilterwindow"]), 1e-10)
    flux = gaussian_filter(flux, 1)
    return pixel, flux


def markov_gaussian(x, amp, mean, std):
    return amp * np.exp(-(x - mean) ** 2 / (2 * std ** 2))


def run(lib, pixel, flux, lines, n_samples):
    nbins = lib.get_histogram_nbins(n_samples)
    hist_output = np.zeros(4 * 2 * nbins, dtype=np.float64)
    nthr = ctypes.c_int(0)
    offset_zero = TRUE["offset"] - 3.0  # deliberately imperfect starting point

    t0 = time.perf_counter()
    lib.run_mcmc(
        np.ascontiguousarray(pixel), np.ascontiguousarray(flux), len(pixel),
        np.ascontiguousarray(lines), len(lines), n_samples,
        offset_zero, OPTS["extent_zero"] / NPIX, OPTS["quad_zero"], OPTS["cube_zero"],
        OPTS["offset_stepsize"], OPTS["linear_stepsize"],
        OPTS["quad_stepsize"], OPTS["cube_stepsize"],
        OPTS["c_cov"], OPTS["s_cov"], OPTS["q_cov"], OPTS["cub_cov"],
        OPTS["accept_param"], hist_output, nbins, ctypes.byref(nthr))
    dt = time.perf_counter() - t0

    params = []
    for p in range(4):
        base = p * 2 * nbins
        centers = hist_output[base:base + nbins]
        hist = hist_output[base + nbins:base + 2 * nbins]
        try:
            popt, _ = curve_fit(markov_gaussian, centers, hist,
                                p0=[np.max(hist), centers[np.argmax(hist)],
                                    np.std(centers)], maxfev=1000000)
            params.append(popt[1])
        except RuntimeError:
            params.append(np.average(centers, weights=hist))
    return dt, params, nthr.value


def main():
    global TRUE, NPIX, LINELIST, OPTS
    libpath = sys.argv[1] if len(sys.argv) > 1 else "./liblinefit.so"
    n_samples = int(sys.argv[2]) if len(sys.argv) > 2 else 250000
    repeats = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    which = sys.argv[4] if len(sys.argv) > 4 else "SOAR_930"

    names = list(PRESETS) if which == "all" else [which]
    lib = load(libpath)

    for name in names:
        cfg = PRESETS[name]
        TRUE, NPIX, LINELIST, OPTS = (cfg["true"], cfg["npix"],
                                      cfg["linelist"], cfg["opts"])
        lines = np.genfromtxt(LINELIST)[:, 0]
        pixel, flux = make_compspec(lines)
        pixel_grid = np.arange(NPIX, dtype=np.float64)
        wl_true = poly(pixel_grid, TRUE["offset"], TRUE["linear"],
                       TRUE["quad"], TRUE["cub"])

        times, rmss, peaks = [], [], []
        for r in range(repeats):
            dt, params, nthr = run(lib, pixel, flux, lines, n_samples)
            times.append(dt)
            wl_fit = poly(pixel_grid, *params)
            rmss.append(np.sqrt(np.mean((wl_fit - wl_true) ** 2)))
            peaks.append(np.max(np.abs(wl_fit - wl_true)))

        print(f"\nSUMMARY {libpath}  preset={name}  n_samples={n_samples} "
              f"repeats={repeats} threads={nthr}")
        print(f"  time      : best={min(times):8.2f}s  "
              f"median={sorted(times)[len(times)//2]:8.2f}s")
        print(f"  RMS error : mean={np.mean(rmss):.4f} A  std={np.std(rmss):.4f} A  "
              f"worst={max(rmss):.4f} A")
        print(f"  max error : mean={np.mean(peaks):.4f} A  worst={max(peaks):.4f} A")


if __name__ == "__main__":
    main()
