#!/usr/bin/env python3
"""Compare a Mizoguchi run against the committed baseline record.

Runs Mizo_theta.py at each baseline time point and reports the delta, so a
stage of the reformulation can be accepted or rejected from one command.
"""
import argparse
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent

# (json path, label, tolerance, higher_is_better)
METRICS = [
    (("profile_metrics", "rmse"), "RMSE", 0.002, False),
    (("skill_vs_no_redistribution",), "Skill", 0.03, True),
    (("profile_metrics", "bias"), "bias", 0.005, None),
    (("profile_metrics", "max_abs_error"), "max|err|", 0.01, False),
    (("zero_crossing", "difference_m"), "front err [m]", 0.002, None),
    (("redistribution", "upper_gain_fraction_of_experiment"), "upper gain / exp", 0.15, None),
    (("solver", "nonlinear_rejection_fraction"), "NL reject frac", 0.05, False),
]


def dig(d, path):
    for k in path:
        if not isinstance(d, dict) or k not in d:
            return None
        d = d[k]
    return d


def compare_result(key, base, cur):
    """Compare one evaluated time against the baseline record.

    Returns the printable rows and the failures. A metric that is absent from
    either side is a failure and not a skip: a run that produced no metrics at
    all would otherwise compare nothing and be reported as agreement.
    """
    if base is None:
        return [], [f"{key}: no baseline record for this time"]
    if not cur:
        return [], [f"{key}: the run produced no result to compare"]

    rows, failures = [], []
    for path, label, tol, higher_better in METRICS:
        b, c = dig(base, path), dig(cur, path)
        if b is None and c is None:
            continue
        if b is None or c is None:
            missing = "baseline" if b is None else "run"
            failures.append(f"{key} {label}: not available in the {missing}")
            continue
        if not math.isfinite(b) or not math.isfinite(c):
            failures.append(f"{key} {label}: nonfinite metric")
            continue
        d = c - b
        flag = ""
        if abs(d) > tol:
            # a move in the good direction is reported, not failed
            good = (higher_better is True and d > 0) or (higher_better is False and d < 0)
            flag = "  better" if good else "  *** OVER TOL"
            if not good:
                failures.append(f"{key} {label}: {b:.6g} -> {c:.6g} (tol {tol})")
        rows.append((label, b, c, d, tol, flag))

    if not rows and not failures:
        failures.append(f"{key}: no metric could be compared")
    return rows, failures


def initial_condition_differences(baseline, current):
    """Identify incompatible recorded initial states before regression checks."""
    differences = []
    for key in ("initial_pressure_Pa", "initial_temperature_C"):
        if key in baseline and key in current:
            if not math.isclose(baseline[key], current[key], rel_tol=1e-10, abs_tol=1e-8):
                differences.append(f"{key}: baseline={baseline[key]}, run={current[key]}")
    return differences


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vtu_dir", type=Path, required=True)
    ap.add_argument("--snapshot_format", choices=("auto", "vtu", "hdf5"), default="auto",
                    help="Field-output layout to read; passed through to Mizo_theta.py")
    ap.add_argument("--baseline", type=Path, default=HERE / "mizo_baseline.json")
    ap.add_argument("--times", type=str, default="12,24,50")
    ap.add_argument("--label", type=str, default="run")
    args = ap.parse_args()

    baseline = json.loads(args.baseline.read_text())
    times = [t.strip() for t in args.times.split(",") if t.strip()]

    print(f"baseline: {args.baseline}")
    for field in ("source_run", "initial_condition"):
        if baseline.get(field) is not None:
            print(f"baseline {field}: {baseline[field]}")
    print(f"run: {args.vtu_dir} ({args.snapshot_format})")

    from Mizo_theta import _initial_condition
    differences = initial_condition_differences(
        baseline.get("initial_condition", {}), _initial_condition(args.vtu_dir)
    )
    if differences:
        print("INCOMPARABLE INITIAL CONDITIONS: this is not a code regression")
        for difference in differences:
            print("  -", difference)
        return 2

    compared = 0
    failures = []
    for t in times:
        key = f"{t}h"
        base = baseline.get("times", {}).get(key)
        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as tmp:
            out = Path(tmp.name)
        cmd = [sys.executable, str(HERE / "Mizo_theta.py"),
               "--vtu_dir", str(args.vtu_dir), "--time", t,
               "--snapshot_format", args.snapshot_format,
               "--json_output", str(out)]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            print(f"[FAIL] Mizo_theta.py failed at {key}\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}")
            failures.append(f"{key}: evaluation failed")
            out.unlink(missing_ok=True)
            continue
        cur = json.loads(out.read_text())
        out.unlink(missing_ok=True)

        print(f"\n=== t = {key}  ({args.label} vs baseline) ===")
        initial = (cur.get("initial_condition") or {}).get("initial_pressure_Pa")
        if initial is not None:
            print(f"run initial pressure: {initial:.4f} Pa")
        rows, time_failures = compare_result(key, base, cur)
        failures.extend(time_failures)
        if rows:
            compared += 1
            print(f"{'metric':22s} {'baseline':>12s} {'run':>12s} {'delta':>12s}  tol")
            for label, b, c, d, tol, flag in rows:
                print(f"{label:22s} {b:12.6g} {c:12.6g} {d:+12.6g}  {tol:g}{flag}")

    print("\n" + "=" * 60)
    if not compared:
        failures.append("no time could be compared against the baseline")
    if failures:
        print("REGRESSION:")
        for f in failures:
            print("  -", f)
        return 1
    print(f"OK: every metric of {compared} time(s) within tolerance of the baseline")
    return 0


if __name__ == "__main__":
    sys.exit(main())
