#!/usr/bin/env python3
"""Compare Google Benchmark JSON results against a baseline.

Given a current benchmark result JSON and a baseline JSON, this script
computes percentage differences for selected metrics and exits with a
non‑zero status if any regression exceeds the configured tolerance.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def load_benchmarks(path: Path) -> dict[str, dict]:
    """Load benchmark results from *path* into a dict keyed by name."""
    with path.open() as fh:
        data = json.load(fh)
    benchmarks = {}
    for bench in data.get("benchmarks", []):
        name = bench.get("name")
        if name:
            benchmarks[name] = bench
    return benchmarks


def compare_metrics(current: dict, baseline: dict, tolerance: float) -> list[str]:
    """Compare metrics between current and baseline.

    Returns a list of regression messages exceeding *tolerance*.
    """
    regressions: list[str] = []
    cpu_cur = current.get("cpu_time")
    cpu_base = baseline.get("cpu_time")
    if cpu_cur is not None and cpu_base:
        diff = (cpu_cur - cpu_base) / cpu_base * 100
        if diff > tolerance:
            regressions.append(
                f"{current['name']}: cpu_time +{diff:.2f}% (baseline {cpu_base}, current {cpu_cur})"
            )
    ips_cur = current.get("items_per_second")
    ips_base = baseline.get("items_per_second")
    if ips_cur is not None and ips_base:
        diff = (ips_cur - ips_base) / ips_base * 100
        if diff < -tolerance:  # lower items/sec is a regression
            regressions.append(
                f"{current['name']}: items_per_second {diff:.2f}% (baseline {ips_base}, current {ips_cur})"
            )
    return regressions


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path, help="Baseline benchmark JSON file")
    parser.add_argument("current", type=Path, help="Current benchmark JSON file")
    parser.add_argument(
        "--tolerance",
        type=float,
        default=5.0,
        help="Allowed percentage regression before failing (default: 5)",
    )
    args = parser.parse_args()

    baseline_data = load_benchmarks(args.baseline)
    current_data = load_benchmarks(args.current)

    regressions: list[str] = []
    for name, current in current_data.items():
        base = baseline_data.get(name)
        if not base:
            continue
        regressions.extend(compare_metrics(current, base, args.tolerance))

    if regressions:
        print("Benchmark regressions detected:", file=sys.stderr)
        for msg in regressions:
            print(" -", msg, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
