#!/usr/bin/env python3
"""Write a structured rectangle as a Gmsh 4.1 mesh in the x-z plane.

Physical groups follow the Mizoguchi projects: surface 1 is the soil, curves
2, 3, 4 and 5 are the left (x = 0), bottom (z = 0), right (x = width) and top
(z = height) boundaries. The vertical spacing may be graded towards the top.
With `--interface Z` the soil above z = Z is a second layer, surface 6, and
each layer gets rows in proportion to its height.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import gmsh


def build(width: float, height: float, nx: int, nz: int, element: str, grading: float,
          interface: float | None) -> None:
    gmsh.model.add("rectangle")
    levels = [0.0, height] if interface is None else [0.0, interface, height]
    lower_left = gmsh.model.geo.addPoint(0.0, 0.0, 0.0)
    lower_right = gmsh.model.geo.addPoint(width, 0.0, 0.0)
    bottom = gmsh.model.geo.addLine(lower_left, lower_right)
    lower, across, lefts, rights, surfaces = bottom, [bottom], [], [], []
    for below, level in zip(levels[:-1], levels[1:]):
        upper_right = gmsh.model.geo.addPoint(width, 0.0, level)
        upper_left = gmsh.model.geo.addPoint(0.0, 0.0, level)
        right = gmsh.model.geo.addLine(lower_right, upper_right)
        top = gmsh.model.geo.addLine(upper_right, upper_left)
        left = gmsh.model.geo.addLine(upper_left, lower_left)
        surface = gmsh.model.geo.addPlaneSurface([gmsh.model.geo.addCurveLoop([lower, right, top, left])])
        rows = max(1, round(nz * (level - below) / height))
        # Progression > 1 on the right edge (bottom to top) and its inverse on the
        # left edge (top to bottom) both refine towards the top.
        gmsh.model.geo.mesh.setTransfiniteCurve(right, rows + 1, "Progression", 1.0 / grading)
        gmsh.model.geo.mesh.setTransfiniteCurve(left, rows + 1, "Progression", grading)
        across.append(top)
        lefts.append(left)
        rights.append(right)
        surfaces.append(surface)
        lower, lower_left, lower_right = -top, upper_left, upper_right
    for line in across:
        gmsh.model.geo.mesh.setTransfiniteCurve(line, nx + 1)
    for surface in surfaces:
        gmsh.model.geo.mesh.setTransfiniteSurface(surface)
        if element == "quad":
            gmsh.model.geo.mesh.setRecombine(2, surface)
    gmsh.model.geo.synchronize()

    gmsh.model.addPhysicalGroup(2, surfaces[:1], 1)
    if len(surfaces) > 1:
        gmsh.model.addPhysicalGroup(2, surfaces[1:], 6)
    gmsh.model.addPhysicalGroup(1, lefts, 2)
    gmsh.model.addPhysicalGroup(1, [bottom], 3)
    gmsh.model.addPhysicalGroup(1, rights, 4)
    gmsh.model.addPhysicalGroup(1, across[-1:], 5)
    gmsh.model.mesh.generate(2)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--width", type=float, default=0.08)
    parser.add_argument("--height", type=float, default=0.2)
    parser.add_argument("--nx", type=int, default=4)
    parser.add_argument("--nz", type=int, default=20)
    parser.add_argument("--element", choices=("tri", "quad"), default="tri")
    parser.add_argument("--grading", type=float, default=1.0,
                        help="ratio of successive cell heights towards the top (< 1 refines the top)")
    parser.add_argument("--interface", type=float, default=None,
                        help="height of a layer interface; the soil above it is surface 6")
    args = parser.parse_args()

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    try:
        build(args.width, args.height, args.nx, args.nz, args.element, args.grading, args.interface)
        gmsh.write(str(args.output))
    finally:
        gmsh.finalize()


if __name__ == "__main__":
    main()
