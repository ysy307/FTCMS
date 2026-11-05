#!/usr/bin/env python3
"""Write a Mizoguchi freezing-column project from project/Mizo-xz-Convex.

Every variant the column studies ran with `-mizo_*` options is a project input
here, so FTCMS alone runs it:

* the top condition: `convective` (coolant behind a surface resistance),
  `plate` (the measured plate temperatures imposed on the soil) or `contact`
  (the plate behind a contact conductance, optionally switching when the
  surface freezes);
* the base condition: `leak` (the room behind a resistance) or `measured`
  (the measured 20 cm temperatures);
* the discretisation of the water equation, the ice impedance factor and the
  mesh (the project's own, a given file, or a structured rectangle).

The temperatures are those of docs/data/soil_temperature_profiles.csv; the
0 h plate value is the intercept of the exponential fitted to the measured ones.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "project" / "Mizo-xz-Convex" / "Input"

PLATE_HOURS = [0.0, 1.0, 3.0, 6.0, 12.0, 24.0, 45.0]
PLATE_TEMPERATURES = [-2.7, -4.2, -5.7, -6.5, -6.6, -6.7, -6.8]
BASE_TEMPERATURES = [6.7, 7.1, 5.9, 4.4, 2.4, 1.4, 1.1]
TOP_ID, BASE_ID = 5, 3

DISCRETIZATIONS = {
    "mixed": {"space": "mixed", "degree": 0, "coefficient_evaluation": "quadrature"},
    "lagrange": {"space": "lagrange", "degree": 1, "coefficient_evaluation": "quadrature"},
    "nodal": {"space": "lagrange", "degree": 1, "coefficient_evaluation": "nodal"},
}


def table(values: list) -> list[dict]:
    return [{"time": 3600.0 * hour, "value": value} for hour, value in zip(PLATE_HOURS, values)]


def top_condition(arguments: argparse.Namespace) -> dict:
    if arguments.top == "convective":
        return {"type": "robin", "value_type": "constant", "values": [arguments.top_coefficient, -6.0]}
    if arguments.top == "plate":
        return {"type": "dirichlet", "value_type": "table", "values": table(PLATE_TEMPERATURES)}
    switch = []
    if arguments.contact_frozen_coefficient is not None:
        switch = [arguments.contact_frozen_coefficient, arguments.switch_temperature, arguments.switch_width]
    rows = [[arguments.contact_coefficient, temperature] + switch for temperature in PLATE_TEMPERATURES]
    return {"type": "robin", "value_type": "table", "values": table(rows)}


def base_condition(arguments: argparse.Namespace) -> dict:
    if arguments.base == "leak":
        return {"type": "robin", "value_type": "constant", "values": [3.0, 6.7]}
    return {"type": "dirichlet", "value_type": "table", "values": table(BASE_TEMPERATURES)}


def write_project(arguments: argparse.Namespace) -> Path:
    target = arguments.output / "Input"
    target.mkdir(parents=True, exist_ok=True)

    basic = json.loads((BASE / "Basic.json").read_text())
    basic["discretization_settings"] = {
        "thermal": DISCRETIZATIONS["lagrange"],
        "hydraulic": DISCRETIZATIONS[arguments.hydraulic],
    }
    for material in basic["materials"]:
        material["hydraulic_conductivity_model"]["impedance_factor"] = arguments.impedance_factor
    if arguments.mesh is not None:
        shutil.copy(arguments.mesh, target / arguments.mesh.name)
        basic["geometry_settings"] = {"mesh_file_name": arguments.mesh.name}
    elif arguments.nz is not None:
        subprocess.run([sys.executable, str(ROOT / "scripts" / "mesh" / "rectangle_msh.py"), str(target / "column.msh"),
                        "--width", str(arguments.width), "--nx", str(arguments.nx), "--nz", str(arguments.nz),
                        "--element", arguments.element, "--grading", str(arguments.grading)], check=True)
        basic["geometry_settings"] = {"mesh_file_name": "column.msh"}
    else:
        shutil.copy(BASE / basic["geometry_settings"]["mesh_file_name"], target)
    (target / "Basic.json").write_text(json.dumps(basic, indent=4))

    conditions = json.loads((BASE / "Conditions.json").read_text())
    conditions["time_controls"]["simulation_period"] = {"unit": "hour", "start": 0.0, "end": arguments.hours}
    for condition in conditions["boundary_conditions"]:
        if condition["id"] == TOP_ID:
            condition["thermal"] = top_condition(arguments)
        elif condition["id"] == BASE_ID:
            condition["thermal"] = base_condition(arguments)
    (target / "Conditions.json").write_text(json.dumps(conditions, indent=4))

    output = json.loads((BASE / "Output.json").read_text())
    if arguments.mesh is None and arguments.nz is not None:
        # The project's observation points sit at x = 4 cm; keep them inside a narrower rectangle.
        for point in output["history_output"]["coordinates"]:
            point["x"] = 0.5 * arguments.width
    (target / "Output.json").write_text(json.dumps(output, indent=4))
    return arguments.output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("output", type=Path, help="project directory to write (its Input/ is created)")
    parser.add_argument("--top", choices=("convective", "plate", "contact"), default="plate")
    parser.add_argument("--top-coefficient", type=float, default=28.0, help="convective top, W/(m2 K)")
    parser.add_argument("--contact-coefficient", type=float, default=200.0, help="contact top, W/(m2 K)")
    parser.add_argument("--contact-frozen-coefficient", type=float, default=None,
                        help="contact coefficient once the surface froze, W/(m2 K)")
    parser.add_argument("--switch-temperature", type=float, default=0.0, help="C")
    parser.add_argument("--switch-width", type=float, default=0.5, help="K")
    parser.add_argument("--base", choices=("leak", "measured"), default="leak")
    parser.add_argument("--hydraulic", choices=sorted(DISCRETIZATIONS), default="mixed")
    parser.add_argument("--impedance-factor", type=float, default=7.0)
    parser.add_argument("--hours", type=float, default=50.0)
    parser.add_argument("--mesh", type=Path, default=None, help="a Gmsh file to use instead of the project's")
    parser.add_argument("--nx", type=int, default=2)
    parser.add_argument("--nz", type=int, default=None, help="build a structured rectangle with this many rows")
    parser.add_argument("--width", type=float, default=0.02)
    parser.add_argument("--element", choices=("tri", "quad"), default="quad")
    parser.add_argument("--grading", type=float, default=1.0)
    print(write_project(parser.parse_args()))


if __name__ == "__main__":
    main()
