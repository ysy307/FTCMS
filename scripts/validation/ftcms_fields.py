"""Read an FTCMS field file in the physical frame and with the analysis names.

FTCMS writes its fields through PETSc's VTU viewer: on the intrinsic 2D mesh of
the analysis plane, with the `Output.json` variable names, a Lagrange potential
as point data and everything else as cell data. The analysis scripts work in
the project's (x, y, z) and with the names the in-house writer used; this puts
a grid in that form once, so every script reads either kind of file alike.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pyvista as pv

#: Output.json name -> the name the analysis scripts use.
ALIASES = {
    "temperature": "Temperature",
    "pressure": "Pressure",
    "water_content": "WaterContent",
    "ice_content": "IceContent",
    "water_flux": "waterFlux",
}

#: `calculate_type` of Basic.json for an analysis in the x-z plane.
XZ_PLANE = 2


def computation_type(path: Path) -> int | None:
    """`simulation_settings.calculate_type` of the project above `path`, if any."""
    for parent in [Path(path).resolve(), *Path(path).resolve().parents]:
        basic = parent / "Input" / "Basic.json"
        if basic.is_file():
            return int(json.loads(basic.read_text())["simulation_settings"]["calculate_type"])
    return None


#: VTK cell type of a cell with this many corners in the plane.
_PLANE_CELLS = {3: pv.CellType.TRIANGLE, 4: pv.CellType.QUAD}


def _read_hdf5(path: Path) -> pv.UnstructuredGrid:
    """One HDF5 field file of PETSc's visualisation layout, as a grid.

    The file holds /geometry/vertices, /viz/topology/cells, and one dataset per
    field under /vertex_fields or /cell_fields named `fields_<name>`.
    """
    import h5py

    with h5py.File(path, "r") as handle:
        points = np.asarray(handle["geometry/vertices"], dtype=float)
        cells = np.asarray(handle["viz/topology/cells"], dtype=np.int64)
        fields = {}
        for group, kind in (("vertex_fields", "point"), ("cell_fields", "cell")):
            if group not in handle:
                continue
            for key, data in handle[group].items():
                fields[(kind, key.removeprefix("fields_"))] = np.asarray(data)
    if points.shape[1] == 2:
        points = np.column_stack([points, np.zeros(len(points))])
    corners = cells.shape[1]
    connectivity = np.hstack([np.full((len(cells), 1), corners, dtype=np.int64), cells]).ravel()
    grid = pv.UnstructuredGrid(connectivity, np.full(len(cells), _PLANE_CELLS[corners], dtype=np.uint8), points)
    for (kind, name), values in fields.items():
        (grid.point_data if kind == "point" else grid.cell_data)[name] = values
    return grid


def read_fields(path: Path) -> pv.DataSet:
    """The grid of one field file, in physical coordinates and with both names of every field."""
    grid = _read_hdf5(path) if Path(path).suffix == ".h5" else pv.read(path)
    planar = grid.bounds[5] - grid.bounds[4] == 0.0
    in_xz = planar and computation_type(path) == XZ_PLANE
    if in_xz:
        points = np.asarray(grid.points)
        grid.points = np.column_stack([points[:, 0], np.zeros(len(points)), points[:, 1]])
    for data in (grid.point_data, grid.cell_data):
        # PETSc names a one-component point field "<name>.0".
        for written in [key for key in data.keys() if key.endswith(".0") and key[:-2] in ALIASES]:
            data[written[:-2]] = np.asarray(data[written])
        for name, alias in ALIASES.items():
            if name not in data or alias in data:
                continue
            values = np.asarray(data[name])
            if values.ndim == 2 and values.shape[1] in (2, 3):
                # A vector of the intrinsic plane: PETSc pads its two components with a zero.
                zeros = np.zeros(len(values))
                values = (np.column_stack([values[:, 0], zeros, values[:, 1]]) if in_xz
                          else np.column_stack([values[:, 0], values[:, 1], zeros]))
            data[alias] = values
    return grid
