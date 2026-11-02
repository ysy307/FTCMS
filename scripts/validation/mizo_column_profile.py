#!/usr/bin/env python3
"""Compare a mixed-column cell profile against the Mizoguchi water-content table.

The profile this reads is the one ``type_mixed_column%write_profile`` writes:
one row per cell, carrying the cell's centroid, its volume and the state the
closure evaluated there. The comparison stays in cells. A depth band's water
content is the volume-weighted mean of the cells whose centroid falls in it,
and a cell that straddles a band edge is not split -- the P0 unknown has no
sub-cell structure to split, and inventing one is how a comparison acquires
water the run never moved.

The measured table is ``Mizo_theta.EXPERIMENTAL_DATA``. Its digitisation is not
recorded anywhere in this repository, so every result here reports that
provenance alongside the numbers; an agreement measured against it is an
agreement with that table and not, by itself, with the experiment.
"""
import argparse
import csv
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from Mizo_theta import (  # noqa: E402
    EXPERIMENTAL_DATA,
    EXPERIMENTAL_REFERENCE,
    ICE_TO_WATER_VOLUME_RATIO,
    TIME_INDEX,
)

# The table is stated at 1 cm intervals measured down from the cooled surface,
# and its key is the lower edge of the band: 0.01 is the first centimetre.
BAND_HEIGHT = 0.01


class ProfileError(RuntimeError):
    """A profile that cannot be read, or cannot be compared to the table."""


def read_profile(path):
    """One written profile, as a list of per-cell dictionaries."""
    rows = []
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            try:
                cell = {
                    "time_s": float(row["time_s"]),
                    "z_m": float(row["z_m"]),
                    "volume_m2": float(row["volume_m2"]),
                    "temperature_C": float(row["temperature_C"]),
                    "pressure_water_Pa": float(row["pressure_water_Pa"]),
                    "pressure_liquid_Pa": float(row["pressure_liquid_Pa"]),
                    "theta_w": float(row["theta_w"]),
                    "theta_i": float(row["theta_i"]),
                    "stored_water": float(row["stored_water"]),
                    "flux_z_m_per_s": float(row["flux_z_m_per_s"]),
                }
            except (KeyError, TypeError, ValueError) as exc:
                raise ProfileError(f"{path}: unreadable row: {exc}") from exc
            rows.append(cell)
    if not rows:
        raise ProfileError(f"{path}: holds no cells")

    times = {round(c["time_s"], 6) for c in rows}
    if len(times) != 1:
        raise ProfileError(f"{path}: rows carry {len(times)} different times")

    for cell in rows:
        for key, value in cell.items():
            if not math.isfinite(value):
                raise ProfileError(f"{path}: nonfinite {key}")
        if cell["volume_m2"] <= 0.0:
            raise ProfileError(f"{path}: a cell has no volume")
    return rows


def total_water(cell):
    """The quantity the measured table is: liquid plus ice as its own volume."""
    return cell["theta_w"] + ICE_TO_WATER_VOLUME_RATIO * cell["theta_i"]


def band_profile(rows, height):
    """Volume-weighted cell means per 1 cm band, keyed by depth below the top.

    The key is the band's lower edge in metres, matching the table: 0.01 is the
    first centimetre below the cooled surface.
    """
    bands = {}
    for cell in rows:
        depth = height - cell["z_m"]
        index = int(math.floor(depth / BAND_HEIGHT))
        # A centroid exactly on the far face belongs to the last band.
        index = min(max(index, 0), int(round(height / BAND_HEIGHT)) - 1)
        key = round((index + 1) * BAND_HEIGHT, 6)
        entry = bands.setdefault(key, {"volume": 0.0, "cells": 0, "sums": {}})
        entry["volume"] += cell["volume_m2"]
        entry["cells"] += 1
        for name, value in (
            ("total_water", total_water(cell)),
            ("theta_w", cell["theta_w"]),
            ("theta_i", cell["theta_i"]),
            ("temperature_C", cell["temperature_C"]),
            ("pressure_liquid_Pa", cell["pressure_liquid_Pa"]),
            ("flux_z_m_per_s", cell["flux_z_m_per_s"]),
        ):
            entry["sums"][name] = entry["sums"].get(name, 0.0) + value * cell["volume_m2"]

    return {
        key: {name: total / entry["volume"] for name, total in entry["sums"].items()}
        | {"volume_m2": entry["volume"], "cells": entry["cells"]}
        for key, entry in sorted(bands.items())
    }


def column_inventory(rows):
    """What the profile itself says the column holds, per unit thickness."""
    water = sum(c["stored_water"] * c["volume_m2"] for c in rows)
    ice = sum(c["theta_i"] * c["volume_m2"] for c in rows)
    area = sum(c["volume_m2"] for c in rows)
    return {"stored_water_m2": water, "ice_m2": ice, "area_m2": area}


def compare_to_table(bands, hours):
    """Band means against the measured column for this time, where both exist."""
    if hours not in TIME_INDEX:
        raise ProfileError(
            f"the table carries no column for {hours} h; it has {sorted(TIME_INDEX)}"
        )
    index = TIME_INDEX[hours]

    pairs = []
    for depth, measured in sorted(EXPERIMENTAL_DATA.items()):
        if depth not in bands:
            continue
        pairs.append((depth, bands[depth]["total_water"], measured[index]))
    if not pairs:
        raise ProfileError("no depth band of the profile lines up with the table")

    errors = [computed - observed for _, computed, observed in pairs]
    return {
        "bands_compared": len(pairs),
        "rmse": math.sqrt(sum(e * e for e in errors) / len(errors)),
        "bias": sum(errors) / len(errors),
        "max_abs_error": max(abs(e) for e in errors),
        "per_band": [
            {"depth_m": d, "computed": c, "observed": o, "error": c - o}
            for d, c, o in pairs
        ],
    }


def freezing_front(bands):
    """The shallowest band that carries ice, and the 0 C crossing depth."""
    frozen = [d for d, b in sorted(bands.items()) if b["theta_i"] > 0.0]
    deepest_frozen = max(frozen) if frozen else None

    crossing = None
    ordered = sorted(bands.items())
    for (d0, b0), (d1, b1) in zip(ordered, ordered[1:]):
        if b0["temperature_C"] <= 0.0 < b1["temperature_C"]:
            span = b1["temperature_C"] - b0["temperature_C"]
            crossing = d0 + (d1 - d0) * (0.0 - b0["temperature_C"]) / span
            break
    return {"deepest_frozen_band_m": deepest_frozen, "zero_crossing_m": crossing}


def evaluate(path, hours, height):
    rows = read_profile(path)
    bands = band_profile(rows, height)
    result = {
        "profile": str(path),
        "time_hours": hours,
        "time_s": rows[0]["time_s"],
        "cells": len(rows),
        "column_height_m": height,
        "inventory": column_inventory(rows),
        "front": freezing_front(bands),
        "reference": EXPERIMENTAL_REFERENCE,
    }
    try:
        result["vs_table"] = compare_to_table(bands, hours)
    except ProfileError as exc:
        result["vs_table"] = {"unavailable": str(exc)}
    result["bands"] = {f"{d:.2f}": b for d, b in bands.items()}
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profile", type=Path, help="a CSV written by write_profile")
    parser.add_argument(
        "--hours",
        type=int,
        required=True,
        help="which column of the measured table this profile is to be read against",
    )
    parser.add_argument(
        "--height", type=float, default=0.20, help="column height in metres"
    )
    parser.add_argument("--json", type=Path, help="write the full result here")
    args = parser.parse_args(argv)

    try:
        result = evaluate(args.profile, args.hours, args.height)
    except (ProfileError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    if args.json:
        args.json.write_text(json.dumps(result, indent=2))

    inv = result["inventory"]
    print(f"{args.profile.name}: {result['cells']} cells, t = {result['time_s']:.1f} s")
    print(
        f"  stored water {inv['stored_water_m2']:.6e} m2   ice {inv['ice_m2']:.6e} m2"
    )
    front = result["front"]
    print(
        f"  deepest frozen band {front['deepest_frozen_band_m']}   "
        f"0 C crossing {front['zero_crossing_m']}"
    )

    table = result["vs_table"]
    if "unavailable" in table:
        print(f"  against the table: {table['unavailable']}")
    else:
        print(
            f"  against the table ({table['bands_compared']} bands): "
            f"RMSE {table['rmse']:.4f}  bias {table['bias']:+.4f}  "
            f"max|err| {table['max_abs_error']:.4f}"
        )
        print("  depth   computed  observed     error")
        for band in table["per_band"]:
            print(
                f"  {band['depth_m']:.2f}    {band['computed']:.4f}    "
                f"{band['observed']:.4f}   {band['error']:+.4f}"
            )
    print(f"  reference provenance: {result['reference']['provenance']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
