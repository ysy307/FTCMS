"""Validate FTCMS Mizoguchi water-content profiles and solver robustness."""

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import pyvista as pv


EXPERIMENTAL_DATA = {
    0.01: [0.345, 0.420, 0.375, 0.410],
    0.02: [0.345, 0.420, 0.400, 0.410],
    0.03: [0.345, 0.400, 0.395, 0.400],
    0.04: [0.345, 0.390, 0.395, 0.395],
    0.05: [0.345, 0.345, 0.405, 0.395],
    0.06: [0.345, 0.270, 0.400, 0.390],
    0.07: [0.345, 0.250, 0.345, 0.390],
    0.08: [0.345, 0.300, 0.260, 0.390],
    0.09: [0.340, 0.310, 0.280, 0.385],
    0.10: [0.345, 0.320, 0.290, 0.375],
    0.11: [0.345, 0.325, 0.300, 0.350],
    0.12: [0.345, 0.330, 0.305, 0.310],
    0.13: [0.350, 0.330, 0.310, 0.290],
    0.14: [0.345, 0.335, 0.315, 0.270],
    0.15: [0.345, 0.340, 0.320, 0.260],
    0.16: [0.355, 0.340, 0.320, 0.270],
    0.17: [0.350, 0.340, 0.325, 0.275],
    0.18: [0.350, 0.340, 0.325, 0.275],
    0.19: [0.350, 0.340, 0.325, 0.280],
    0.20: [0.350, 0.340, 0.330, 0.280],
}
TIME_INDEX = {0: 0, 12: 1, 24: 2, 50: 3}
ICE_TO_WATER_VOLUME_RATIO = 917.0 / 1000.0
DZ_DEFAULT = 0.01


def _output_interval_minutes(vtu_dir):
    """Read the field-output interval from the project's Output.json.

    The output index for a target time is time/interval, so a hard-coded interval
    silently reads the wrong snapshot whenever the project changes its output
    cadence. Search upward from the VTU directory for Input/Output.json.
    """
    unit_seconds = {"second": 1.0, "minute": 60.0, "hour": 3600.0, "day": 86400.0}
    for parent in [Path(vtu_dir).resolve(), *Path(vtu_dir).resolve().parents]:
        candidate = parent / "Input" / "Output.json"
        if not candidate.is_file():
            continue
        interval = json.loads(candidate.read_text())["field_output"]["output_interval"]
        unit = str(interval.get("unit", "second")).lower()
        if unit not in unit_seconds:
            raise ValueError(f"Unsupported output_interval unit '{unit}' in {candidate}")
        return float(interval["value"]) * unit_seconds[unit] / 60.0
    raise FileNotFoundError(
        f"Could not locate Input/Output.json above {vtu_dir}; "
        "pass --output_interval_minutes explicitly"
    )


# Cell types by the number of corners a cell has, for the two spatial
# dimensions the solver writes. The file records only the corner indices, so the
# element is identified by how many of them there are.
_CELL_TYPES = {
    (2, 3): pv.CellType.TRIANGLE,
    (2, 4): pv.CellType.QUAD,
    (3, 4): pv.CellType.TETRA,
    (3, 8): pv.CellType.HEXAHEDRON,
}


def _h5_path(output_dir):
    """The single HDF5 a run writes, searched for under the output directory."""
    output_dir = Path(output_dir)
    matches = sorted(set(output_dir.glob("*.h5")) | set(output_dir.glob("Files/*.h5")))
    if len(matches) != 1:
        raise FileNotFoundError(
            f"Expected exactly one HDF5 output under {output_dir}, found {len(matches)}"
        )
    return matches[0]


def _timestep_index(time_hours, output_interval_minutes):
    return round(time_hours * 60.0 / output_interval_minutes)


def _read_grid(path, index):
    """The mesh and the fields of one output time, as an unstructured grid.

    A run writes one self-contained HDF5: the vertices, the corner indices of
    every cell, and one dataset per quantity indexed by output time. The grid is
    assembled here so that everything downstream -- the resampling, the depth
    bins, the metrics -- is the same code whatever wrote the file.
    """
    with h5py.File(path, "r") as handle:
        points = np.asarray(handle["geometry/vertices"])
        cells = np.asarray(handle["viz/topology/cells"])
        fields = {}
        for name in ("WaterContent", "IceContent"):
            if f"fields/{name}" not in handle:
                raise ValueError(f"{name} not found in {path}")
            data = handle[f"fields/{name}"]
            if index >= data.shape[0]:
                raise FileNotFoundError(
                    f"{path} holds {data.shape[0]} output times, none at index {index}"
                )
            fields[name] = np.asarray(data[index]).reshape(data.shape[1], -1)[:, 0]

    if points.shape[1] == 2:
        points = np.column_stack([points, np.zeros(len(points))])
    # A mesh embedded in three coordinates is still two-dimensional when one of
    # them never varies: an x-z column is meshed that way so that gravity acts
    # along a coordinate the solve carries.
    varying = int(sum(np.ptp(points[:, axis]) > 0.0 for axis in range(3)))
    corners = cells.shape[1]
    cell_type = _CELL_TYPES.get((varying, corners))
    if cell_type is None:
        raise ValueError(
            f"No cell type for {corners} corners in {varying} dimensions in {path}"
        )

    counts = np.full((len(cells), 1), corners, dtype=np.int64)
    connectivity = np.hstack([counts, cells.astype(np.int64)]).ravel()
    grid = pv.UnstructuredGrid(
        connectivity, np.full(len(cells), cell_type, dtype=np.uint8), points
    )
    for name, values in fields.items():
        grid.point_data[name] = values
    return grid


def _sample_total_water(grid, resolution):
    water = np.asarray(grid.point_data["WaterContent"])
    ice = np.asarray(grid.point_data["IceContent"])
    grid.point_data["theta_tot"] = water + ICE_TO_WATER_VOLUME_RATIO * ice

    x_min, x_max, y_min, y_max, z_min, z_max = grid.bounds
    nx = max(2, int(np.ceil((x_max - x_min) / resolution)) + 1)
    nz = max(2, int(np.ceil((z_max - z_min) / resolution)) + 1)
    dy = 1.0 if y_min == y_max else y_max - y_min
    image = pv.ImageData(
        dimensions=(nx, 1, nz),
        spacing=(resolution, dy, resolution),
        origin=(x_min, y_min, z_min),
    )
    sampled = image.sample(grid)
    theta = np.asarray(sampled.point_data["theta_tot"])
    valid = np.isfinite(theta)
    if "vtkValidPointMask" in sampled.point_data:
        valid &= np.asarray(sampled.point_data["vtkValidPointMask"]).astype(bool)
    points = sampled.points[valid]
    if points.size == 0:
        raise ValueError("No valid sample points found in the output grid")
    return z_max - points[:, 2], theta[valid]


def _depth_profile(depths, theta, time_hours, dz):
    experiment_index = TIME_INDEX[time_hours]
    rows = []
    for depth_end, values in EXPERIMENTAL_DATA.items():
        depth_start = depth_end - dz
        mask = (depths >= depth_start) & (depths < depth_end)
        simulated = float(np.mean(theta[mask])) if np.any(mask) else None
        experimental = values[experiment_index]
        rows.append(
            {
                "depth_start_m": depth_start,
                "depth_end_m": depth_end,
                "simulated": simulated,
                "experimental": experimental,
                "error": None if simulated is None else simulated - experimental,
            }
        )
    return rows


def _print_profile(rows):
    print(f"{'Depth (m)':<20} | {'Simulated':<12} | {'Exp':<12} | {'Error':<12}")
    print("-" * 65)
    for row in rows:
        interval = f"{row['depth_start_m']:4.2f} - {row['depth_end_m']:4.2f}"
        if row["simulated"] is None:
            print(f"{interval:<20} | {'No Data':<12} | {row['experimental']:12.6f} | {'-':<12}")
        else:
            print(
                f"{interval:<20} | {row['simulated']:12.6f} | "
                f"{row['experimental']:12.6f} | {row['error']:+12.6f}"
            )


def _profile_metrics(rows):
    errors = np.asarray([row["error"] for row in rows if row["error"] is not None])
    if errors.size == 0:
        raise ValueError("No depth bins contain simulation data")
    return {
        "rmse": float(np.sqrt(np.mean(errors**2))),
        "mae": float(np.mean(np.abs(errors))),
        "bias": float(np.mean(errors)),
        "max_abs_error": float(np.max(np.abs(errors))),
    }


def _redistribution_metrics(initial_rows, target_rows, dz):
    simulated_changes = []
    experimental_changes = []
    for initial, target in zip(initial_rows, target_rows):
        if initial["simulated"] is None or target["simulated"] is None:
            simulated_changes.append(np.nan)
        else:
            simulated_changes.append(target["simulated"] - initial["simulated"])
        experimental_changes.append(target["experimental"] - initial["experimental"])
    simulated_changes = np.asarray(simulated_changes)
    experimental_changes = np.asarray(experimental_changes)
    depth_ends = np.asarray(list(EXPERIMENTAL_DATA))
    upper = depth_ends <= 0.05
    lower = (depth_ends > 0.05) & (depth_ends <= 0.12)
    valid = np.isfinite(simulated_changes)
    simulated_upper = float(np.sum(simulated_changes[upper & valid]) * dz)
    experimental_upper = float(np.sum(experimental_changes[upper]) * dz)
    return {
        "upper_0_005m_gain_m": simulated_upper,
        "lower_005_012m_change_m": float(np.sum(simulated_changes[lower & valid]) * dz),
        "column_0_020m_change_m": float(np.sum(simulated_changes[valid]) * dz),
        "experimental_upper_0_005m_gain_m": experimental_upper,
        "experimental_lower_005_012m_change_m": float(
            np.sum(experimental_changes[lower]) * dz
        ),
        "experimental_column_0_020m_change_m": float(np.sum(experimental_changes) * dz),
        "upper_gain_fraction_of_experiment": (
            simulated_upper / experimental_upper if experimental_upper != 0.0 else None
        ),
    }


def _bin_centres():
    return np.asarray([end - 0.5 * DZ_DEFAULT for end in EXPERIMENTAL_DATA])


def _zero_crossing_depth(changes):
    """Depth where the water-content change switches from gain to loss.

    Mizoguchi's columns wet near the cold end and dry below it, so the crossing
    depth is the front position the redistribution has reached. Returns None when
    the profile never changes sign (no redistribution, or all one way).
    """
    centres = _bin_centres()
    values = np.asarray(changes, dtype=float)
    finite = np.isfinite(values)
    if finite.sum() < 2:
        return None
    centres, values = centres[finite], values[finite]
    order = np.argsort(centres)
    centres, values = centres[order], values[order]
    for i in range(len(values) - 1):
        a, b = values[i], values[i + 1]
        if a > 0.0 >= b:
            # Linear interpolation of the sign change between bin centres.
            span = a - b
            if span == 0.0:
                return float(centres[i])
            return float(centres[i] + (centres[i + 1] - centres[i]) * a / span)
    return None


def _skill_score(rows, baseline_key="experimental_initial"):
    """Skill of the simulation against a no-redistribution reference.

    The reference prediction is "the profile never changed", taken from the
    measured initial profile so the reference is independent of the model.
    score = 1 - MSE(model)/MSE(reference); positive means the simulation beats
    assuming no redistribution at all, which is the minimum a freezing model has
    to clear to be saying anything.
    """
    initial = EXPERIMENTAL_DATA
    model_sq, ref_sq = [], []
    for row, values in zip(rows, initial.values()):
        if row["simulated"] is None:
            continue
        model_sq.append((row["simulated"] - row["experimental"]) ** 2)
        ref_sq.append((values[TIME_INDEX[0]] - row["experimental"]) ** 2)
    if not model_sq:
        return None
    mse_model = float(np.mean(model_sq))
    mse_ref = float(np.mean(ref_sq))
    if mse_ref == 0.0:
        return None
    return 1.0 - mse_model / mse_ref


def _profile_changes(initial_rows, target_rows):
    changes = []
    for initial, target in zip(initial_rows, target_rows):
        if initial["simulated"] is None or target["simulated"] is None:
            changes.append(np.nan)
        else:
            changes.append(target["simulated"] - initial["simulated"])
    return changes


def _experimental_changes(time_hours):
    idx0, idx = TIME_INDEX[0], TIME_INDEX[time_hours]
    return [values[idx] - values[idx0] for values in EXPERIMENTAL_DATA.values()]


def _sampling_sensitivity(vtu_dir, time_hours, dz, res, interval_minutes):
    """Largest bin-mean difference between res and res/2 sampling.

    The depth bins are 10 mm and the sampling grid is 1 mm, so a bin mean should
    not depend on the grid. If it does, the reported profile is a sampling
    artefact rather than a result.
    """
    grid = _read_grid(_h5_path(vtu_dir), _timestep_index(time_hours, interval_minutes))
    coarse = _depth_profile(*_sample_total_water(grid, res), time_hours, dz)
    fine = _depth_profile(*_sample_total_water(grid, 0.5 * res), time_hours, dz)
    deltas = [
        abs(c["simulated"] - f["simulated"])
        for c, f in zip(coarse, fine)
        if c["simulated"] is not None and f["simulated"] is not None
    ]
    return {
        "resolution_m": res,
        "refined_resolution_m": 0.5 * res,
        "max_abs_bin_difference": float(max(deltas)) if deltas else None,
    }


def _write_csv(path, results):
    import csv

    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            ["time_hours", "depth_start_m", "depth_end_m", "simulated", "experimental", "error"]
        )
        for result in results:
            for row in result["profile"]:
                writer.writerow(
                    [
                        result["time_hours"],
                        f"{row['depth_start_m']:.4f}",
                        f"{row['depth_end_m']:.4f}",
                        "" if row["simulated"] is None else f"{row['simulated']:.6f}",
                        f"{row['experimental']:.6f}",
                        "" if row["error"] is None else f"{row['error']:.6f}",
                    ]
                )


def _read_solver_history(path):
    names = None
    records = []
    with Path(path).open(encoding="utf-8") as stream:
        for line in stream:
            stripped = line.strip()
            if stripped.startswith("# attempt "):
                names = stripped[2:].split()
            elif stripped and not stripped.startswith("#"):
                if names is None:
                    raise ValueError(f"Missing solver-history schema in {path}")
                values = stripped.split()
                if len(values) != len(names):
                    raise ValueError(f"Malformed solver-history row in {path}: {stripped}")
                records.append(dict(zip(names, values)))
    if not records:
        raise ValueError(f"No solver attempts found in {path}")
    return records


def _solver_metrics(path):
    records = _read_solver_history(path)
    accepted = [record for record in records if int(record["accepted"]) == 1]
    rejected = [record for record in records if int(record["accepted"]) == 0]
    nonlinear_rejected = [record for record in rejected if record["status"] != "lte_rejected"]
    lte_rejected = [record for record in rejected if record["status"] == "lte_rejected"]
    accepted_lte = [
        float(record["lte_rel"])
        for record in accepted
        if float(record["lte_rel"]) >= 0.0
    ]
    return {
        "attempts": len(records),
        "accepted_steps": len(accepted),
        "rejected_steps": len(rejected),
        "nonlinear_rejections": len(nonlinear_rejected),
        "lte_rejections": len(lte_rejected),
        "nonlinear_rejection_fraction": len(nonlinear_rejected) / len(records),
        "final_accepted_time_s": max(float(record["time_accepted_s"]) for record in accepted),
        "max_accepted_lte": max(accepted_lte, default=-1.0),
    }


def calculate_theta_error_high_precision(
    vtu_base_dir,
    target_time_hours,
    dz=0.01,
    res=0.001,
    output_interval_minutes=5.0,
):
    if target_time_hours not in TIME_INDEX:
        raise ValueError(f"Invalid time: {target_time_hours} h")
    target_grid = _read_grid(
        _h5_path(vtu_base_dir),
        _timestep_index(target_time_hours, output_interval_minutes),
    )
    depths, theta = _sample_total_water(target_grid, res)
    rows = _depth_profile(depths, theta, target_time_hours, dz)

    print(f"Target Time: {target_time_hours} h (sampling resolution={res} m)")
    _print_profile(rows)
    metrics = _profile_metrics(rows)
    print(
        "Profile metrics: "
        f"RMSE={metrics['rmse']:.6f}, MAE={metrics['mae']:.6f}, "
        f"bias={metrics['bias']:+.6f}, max|error|={metrics['max_abs_error']:.6f}"
    )

    result = {"time_hours": target_time_hours, "profile": rows, "profile_metrics": metrics}

    skill = _skill_score(rows)
    result["skill_vs_no_redistribution"] = skill
    if skill is not None:
        print(f"Skill vs no-redistribution reference: {skill:+.4f}")

    if target_time_hours > 0:
        initial_grid = _read_grid(
            _h5_path(vtu_base_dir), _timestep_index(0, output_interval_minutes)
        )
        initial_depths, initial_theta = _sample_total_water(initial_grid, res)
        initial_rows = _depth_profile(initial_depths, initial_theta, 0, dz)
        redistribution = _redistribution_metrics(initial_rows, rows, dz)
        result["redistribution"] = redistribution

        changes = _profile_changes(initial_rows, rows)
        zc_sim = _zero_crossing_depth(changes)
        zc_exp = _zero_crossing_depth(_experimental_changes(target_time_hours))
        result["zero_crossing"] = {
            "simulated_m": zc_sim,
            "experimental_m": zc_exp,
            "difference_m": (
                None if zc_sim is None or zc_exp is None else zc_sim - zc_exp
            ),
        }
        print(
            "Zero-crossing depth: "
            f"simulated={'n/a' if zc_sim is None else f'{zc_sim:.4f} m'}, "
            f"experimental={'n/a' if zc_exp is None else f'{zc_exp:.4f} m'}"
            + (
                ""
                if zc_sim is None or zc_exp is None
                else f", difference={zc_sim - zc_exp:+.4f} m"
            )
        )
        print(
            "Redistribution: "
            f"upper gain={redistribution['upper_0_005m_gain_m']:+.6e} m, "
            f"lower change={redistribution['lower_005_012m_change_m']:+.6e} m, "
            f"column change={redistribution['column_0_020m_change_m']:+.6e} m"
        )
        print(
            "Experimental redistribution: "
            f"upper gain={redistribution['experimental_upper_0_005m_gain_m']:+.6e} m, "
            f"lower change={redistribution['experimental_lower_005_012m_change_m']:+.6e} m, "
            "simulated/experimental upper gain="
            f"{redistribution['upper_gain_fraction_of_experiment']:.3f}"
        )
    return result


def _gate_failures(result, time_hours, solver, args):
    """Gate checks for one evaluated time, returned as human-readable strings."""
    failures = []
    label = f"{time_hours} h"
    metrics = result["profile_metrics"]
    if args.rmse_limit is not None and metrics["rmse"] > args.rmse_limit:
        failures.append(f"{label}: RMSE {metrics['rmse']:.6f} exceeds {args.rmse_limit:.6f}")

    skill = result.get("skill_vs_no_redistribution")
    if time_hours > 0 and skill is not None and skill <= 0.0:
        failures.append(
            f"{label}: skill {skill:+.4f} is not positive, so the simulation is no "
            "better than assuming no redistribution"
        )

    if time_hours > 0:
        redistribution = result["redistribution"]
        if redistribution["upper_0_005m_gain_m"] <= 0.0:
            failures.append(f"{label}: no upward water gain in the upper 0.05 m")
        if args.min_upper_gain_fraction is not None:
            fraction = redistribution["upper_gain_fraction_of_experiment"]
            if fraction is None or fraction < args.min_upper_gain_fraction:
                failures.append(
                    f"{label}: upper water gain is "
                    f"{'n/a' if fraction is None else f'{fraction:.3f}'} of the experiment, "
                    f"below {args.min_upper_gain_fraction:.3f}"
                )
        if args.max_abs_column_change_m is not None:
            change = abs(redistribution["column_0_020m_change_m"])
            if change > args.max_abs_column_change_m:
                failures.append(
                    f"{label}: absolute column-water change {change:.6e} m exceeds "
                    f"{args.max_abs_column_change_m:.6e} m"
                )
        if args.max_zero_crossing_error_m is not None:
            difference = (result.get("zero_crossing") or {}).get("difference_m")
            if difference is None or abs(difference) > args.max_zero_crossing_error_m:
                failures.append(
                    f"{label}: zero-crossing error "
                    f"{'n/a' if difference is None else f'{difference:+.4f} m'} exceeds "
                    f"{args.max_zero_crossing_error_m:.4f} m"
                )

    sensitivity = result.get("sampling_sensitivity")
    if sensitivity is not None and args.max_sampling_difference is not None:
        difference = sensitivity["max_abs_bin_difference"]
        if difference is None or difference > args.max_sampling_difference:
            failures.append(
                f"{label}: sampling difference "
                f"{'n/a' if difference is None else f'{difference:.6f}'} exceeds "
                f"{args.max_sampling_difference:.6f}"
            )

    if solver is None:
        failures.append(f"{label}: solver history is required for --check")
    else:
        if solver["final_accepted_time_s"] + 1.0e-8 < time_hours * 3600.0:
            failures.append(f"{label}: simulation did not reach this validation time")
        if solver["nonlinear_rejection_fraction"] > args.max_nonlinear_rejection_fraction:
            failures.append(
                f"{label}: nonlinear rejection fraction "
                f"{solver['nonlinear_rejection_fraction']:.3f} exceeds "
                f"{args.max_nonlinear_rejection_fraction:.3f}"
            )
    return failures


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vtu_dir", type=Path, required=True)
    parser.add_argument(
        "--time",
        required=True,
        help="Target time in hours (%s) or 'all' to evaluate every gate time in order"
        % ", ".join(str(k) for k in sorted(TIME_INDEX)),
    )
    parser.add_argument("--res", type=float, default=0.001)
    parser.add_argument("--output_interval_minutes", type=float,
                        help="Overrides the interval read from Input/Output.json")
    parser.add_argument("--solver_history", type=Path)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--rmse_limit", type=float)
    parser.add_argument("--min_upper_gain_fraction", type=float)
    parser.add_argument("--max_abs_column_change_m", type=float)
    parser.add_argument(
        "--max_zero_crossing_error_m",
        type=float,
        help="Gate on |simulated - experimental| zero-crossing depth, e.g. 0.02 for two bins",
    )
    parser.add_argument("--max_sampling_difference", type=float)
    parser.add_argument("--max_nonlinear_rejection_fraction", type=float, default=0.2)
    parser.add_argument("--json_output", type=Path)
    parser.add_argument("--csv_output", type=Path, help="Per-bin profiles for every evaluated time")
    parser.add_argument(
        "--sampling_sensitivity",
        action="store_true",
        help="Also sample at res/2 and report the largest bin-mean difference",
    )
    args = parser.parse_args()

    if str(args.time).lower() == "all":
        target_times = sorted(TIME_INDEX)
    else:
        try:
            requested = int(args.time)
        except ValueError:
            parser.error(f"--time must be an integer or 'all', got {args.time!r}")
        if requested not in TIME_INDEX:
            parser.error(f"--time must be one of {sorted(TIME_INDEX)} or 'all'")
        target_times = [requested]

    interval_minutes = args.output_interval_minutes
    if interval_minutes is None:
        interval_minutes = _output_interval_minutes(args.vtu_dir)
        print(f"Output interval from Input/Output.json: {interval_minutes:g} min")

    history_path = args.solver_history
    if history_path is None:
        candidate = args.vtu_dir.parent / "solver_history.log"
        if candidate.is_file():
            history_path = candidate
    solver = _solver_metrics(history_path) if history_path is not None else None
    if solver is not None:
        print(
            "Solver metrics: "
            f"attempts={solver['attempts']}, accepted={solver['accepted_steps']}, "
            f"nonlinear rejects={solver['nonlinear_rejections']}, "
            f"LTE rejects={solver['lte_rejections']}, "
            f"nonlinear reject fraction={solver['nonlinear_rejection_fraction']:.3f}, "
            f"final accepted time={solver['final_accepted_time_s'] / 3600.0:.3f} h"
        )
        print()

    results = []
    failures = []
    for target_time in target_times:
        # A time the run never reached has no output file; report it as a failed
        # gate rather than dying on FileNotFoundError, so 'all' still summarises
        # everything the run did produce.
        try:
            result = calculate_theta_error_high_precision(
                args.vtu_dir,
                target_time,
                res=args.res,
                output_interval_minutes=interval_minutes,
            )
        except (FileNotFoundError, ValueError) as error:
            print(f"[{target_time} h] not evaluated: {error}")
            failures.append(f"{target_time} h could not be evaluated: {error}")
            print()
            continue

        if args.sampling_sensitivity:
            sensitivity = _sampling_sensitivity(
                args.vtu_dir, target_time, DZ_DEFAULT, args.res, interval_minutes
            )
            result["sampling_sensitivity"] = sensitivity
            difference = sensitivity["max_abs_bin_difference"]
            print(
                "Sampling sensitivity "
                f"({sensitivity['resolution_m']:g} m vs {sensitivity['refined_resolution_m']:g} m): "
                f"max|bin difference|={'n/a' if difference is None else f'{difference:.6f}'}"
            )

        if solver is not None:
            result["solver"] = solver
        results.append(result)

        if args.check:
            failures.extend(_gate_failures(result, target_time, solver, args))
        print()

    if args.csv_output is not None and results:
        _write_csv(args.csv_output, results)
        print(f"Wrote per-bin profiles to {args.csv_output}")

    if len(results) > 1:
        print("Summary")
        header = f"{'t [h]':>6} {'RMSE':>10} {'bias':>10} {'skill':>9} {'zc_sim':>9} {'zc_exp':>9}"
        print(header)
        for result in results:
            metrics = result["profile_metrics"]
            skill = result.get("skill_vs_no_redistribution")
            crossing = result.get("zero_crossing") or {}
            fields = [
                f"{result['time_hours']:6d}",
                f"{metrics['rmse']:10.6f}",
                f"{metrics['bias']:+10.6f}",
                "      n/a" if skill is None else f"{skill:+9.4f}",
            ]
            for key in ("simulated_m", "experimental_m"):
                value = crossing.get(key)
                fields.append("      n/a" if value is None else f"{value:9.4f}")
            print(" ".join(fields))
        print()

    if args.json_output is not None:
        payload = results[0] if len(results) == 1 else {"times": results}
        args.json_output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    if failures:
        for failure in failures:
            print(f"CHECK FAILED: {failure}")
        raise SystemExit(1)
    if args.check:
        print("All requested V&V checks passed.")


if __name__ == "__main__":
    main()
