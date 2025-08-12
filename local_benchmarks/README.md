# Local Benchmarks

This directory contains benchmarks for profiling `fastp` internals. The build
uses [Google Benchmark](https://github.com/google/benchmark) and provides
helpers for comparing results against a saved baseline.

## Creating a baseline

1. Configure and build the benchmark:

   ```bash
   cmake -S local_benchmarks -B build
   cmake --build build --target run_benchmark
   ```

2. Copy the generated `build/benchmark_results.json` to
   `local_benchmarks/benchmark_baseline.json` to act as the baseline.

## Checking for regressions

After establishing a baseline, run:

```bash
cmake --build build --target check_benchmarks
```

This target reruns the benchmark and invokes
`scripts/compare_benchmarks.py` to compare the new results with the baseline.
If the regression in `cpu_time` or `items_per_second` exceeds the default 5%
threshold, the script exits with a non‑zero status.

Use the `--tolerance` flag in `compare_benchmarks.py` to adjust the threshold
as needed.
