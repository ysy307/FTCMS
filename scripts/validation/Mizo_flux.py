"""Report the vertical liquid-water flux around the freezing front.

Gate 2 of the Dall'Amico driver work asks a single question: does the liquid
pressure p_l = -s_eff actually pull water upward at the front, and at the
order of magnitude the Mizoguchi column needs? The 1 cm-bin water-content
profile that Mizo_theta.py reports only answers that after hours of simulated
time; the flux answers it in the first minutes, which is what makes a short
smoke run a usable gate.

Reference magnitude: redistributing 2.85e-3 m of water into the upper layer
over 50 h needs a mean upward flux of about 1.6e-8 m/s, and the transient peak
is larger. A run whose front flux is orders below that cannot reproduce the
experiment no matter how long it runs.
"""

import argparse
from pathlib import Path

import numpy as np
import pyvista as pv

from ftcms_fields import read_fields

ICE_TO_WATER_VOLUME_RATIO = 917.0 / 1000.0
REFERENCE_FLUX_M_S = 1.74e-8


def _latest_vtu(vtu_dir):
    files = sorted(Path(vtu_dir).glob("Out_*.vtu"))
    if not files:
        raise FileNotFoundError(f"no Out_*.vtu under {vtu_dir}")
    return files[-1]


def _front_depth(z, temperature):
    """Depth of the 0 C isotherm, by linear interpolation on the T(z) profile."""
    order = np.argsort(z)
    z_sorted, t_sorted = z[order], temperature[order]
    sign_change = np.where(np.diff(np.signbit(t_sorted)))[0]
    if sign_change.size == 0:
        return None
    i = sign_change[-1]
    t0, t1 = t_sorted[i], t_sorted[i + 1]
    if t1 == t0:
        return z_sorted[i]
    return z_sorted[i] + (z_sorted[i + 1] - z_sorted[i]) * (0.0 - t0) / (t1 - t0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vtu", type=Path, help="Explicit .vtu file")
    parser.add_argument("--vtu_dir", type=Path, help="Directory; the last Out_*.vtu is used")
    parser.add_argument("--band", type=float, default=0.01,
                        help="Half-width [m] of the band around the front to report")
    args = parser.parse_args()

    if args.vtu is None and args.vtu_dir is None:
        parser.error("pass --vtu or --vtu_dir")
    path = args.vtu if args.vtu is not None else _latest_vtu(args.vtu_dir)

    # Cell fields (the engine's derived quantities) are averaged onto the nodes.
    grid = read_fields(path).cell_data_to_point_data(pass_cell_data=True)
    pd = grid.point_data
    z = grid.points[:, 2]
    temperature = np.asarray(pd["Temperature"])
    flux = np.asarray(pd["waterFlux"])
    qz = flux[:, 2]
    water = np.asarray(pd["WaterContent"])
    ice = np.asarray(pd["IceContent"])
    theta_tot = water + ICE_TO_WATER_VOLUME_RATIO * ice

    print(f"file: {path}")
    print(f"nodes: {z.size}   z range: {z.min():.4f} .. {z.max():.4f} m")
    print(f"T range: {temperature.min():+.4f} .. {temperature.max():+.4f} C")

    front = _front_depth(z, temperature)
    if front is None:
        print("freezing front: none (no 0 C crossing in the column)")
    else:
        print(f"freezing front (0 C isotherm): z = {front:.4f} m")
        band = np.abs(z - front) <= args.band
        if band.any():
            qz_band = qz[band]
            print(f"band |z-front| <= {args.band:.3f} m: {band.sum()} nodes")
            print(f"  q_z  max upward = {qz_band.max():+.3e} m/s")
            print(f"  q_z  mean       = {qz_band.mean():+.3e} m/s")
            print(f"  q_z  min        = {qz_band.min():+.3e} m/s")
            print(f"  ratio to reference {REFERENCE_FLUX_M_S:.2e} m/s: "
                  f"{qz_band.max() / REFERENCE_FLUX_M_S:+.2f}")

    print(f"column-wide q_z max upward = {qz.max():+.3e} m/s "
          f"at z = {z[np.argmax(qz)]:.4f} m")
    print(f"theta_tot: min {theta_tot.min():.4f}  max {theta_tot.max():.4f}  "
          f"mean {theta_tot.mean():.4f}")

    top = z >= z.max() - 1.0e-6
    if top.any():
        print(f"top node row: theta_tot max {theta_tot[top].max():.4f}  "
              f"T mean {temperature[top].mean():+.4f} C")


if __name__ == "__main__":
    main()
