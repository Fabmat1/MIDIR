# Wavelength-solver benchmark & verification harness

Tools for changing `linefit.cpp` without silently degrading the wavelength
solution. Run everything from the repository root.

## Speed + solution quality (end to end)

Builds a synthetic arc lamp from a *known* cubic dispersion solution, runs it
through the same `ctypes` path `data_reduction.py` uses, and reports both wall
time and how well the true solution was recovered.

```bash
make
python bench/bench_mcmc.py ./liblinefit.so 250000 3 all
```

Arguments: `<lib> [n_samples] [repeats] [preset|all]`. Presets cover a weakly
curved case (`SOAR_930`), a near-linear one (`ALFOSC_grism18`) and a strongly
curved one (`EFOSC_grism7`). Each preset's "true" solution is deliberately
placed *inside* the search box the corresponding `presets/*.json` gives the
solver — a truth outside the box only measures how fast the solver fails.

`RMS error` is the metric that matters: it is the RMS difference between the
recovered and true wavelength arrays, in Å. Compare it against a baseline build
of the previous `linefit.cpp`; run-to-run scatter is ~0.003 Å, so only
differences well above that mean anything.

To benchmark against the old code, build it from git and point the script at it:

```bash
git show <rev>:linefit.cpp > bench/linefit_baseline.cpp
g++ -O3 -march=native -std=c++17 -fPIC -fopenmp -shared \
    bench/linefit_baseline.cpp -o bench/liblinefit_baseline.so -lpthread
python bench/bench_mcmc.py ./bench/liblinefit_baseline.so 250000 3 all
```

## Root-solver equivalence

`verify_inverse.cpp` compares the shipped normalised-coordinate root solver
against the original raw-pixel one over the search box of every preset, and
counts polynomial evaluations per call.

```bash
g++ -O2 -std=c++17 bench/verify_inverse.cpp -o bench/verify_inverse && ./bench/verify_inverse
```

It reports how often each solver recovers the pixel a probe wavelength was
generated from. Large `max |dx|` values are expected for presets whose search
box admits *non-monotonic* dispersions: there, several pixels share a
wavelength and "the" root is genuinely ambiguous. The number to watch is
`wrong root (monotonic only)`, which must stay at 0.

## Objective equivalence

`verify_objective.cpp` checks that the vectorised fast path in
`compute_objective` returns the same score as the robust scalar solver — the
score is the only thing the MCMC ever sees.

```bash
g++ -O3 -march=native -std=c++17 -fopenmp bench/verify_objective.cpp \
    -o bench/verify_objective && ./bench/verify_objective
```

`max rel` should be ~1e-11 (pure floating-point rounding). It also reports how
often the monotonicity gate diverts a proposal to the robust path and how often
a line needs the per-line fallback.

## Profiling

`bench/bench_mcmc.py` goes through Python; to profile the solver alone, link a
small C driver against the library and run `perf` on that. `compute_objective`
should dominate — if the RNG (`normal_distribution`) climbs above ~10%, the
objective has got fast enough that the sampler itself is worth looking at.
