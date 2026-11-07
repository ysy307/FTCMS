#!/usr/bin/env python3
"""Verification cases of the FTCMS engine, one per ctest entry.

Each case is a project, the changes it makes to it, and one statement about
the result: a quantity, its expected value and a tolerance, or the refusal the
input must meet. A case copies its project into a fresh directory, applies its
changes, runs FTCMS there and checks the observation files. A failure prints
the case, the expected value and the actual value.

    verification.py CASE --ftcms bin/FTCMS --work DIR [--mpirun MPIRUN] [--ranks N]
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

PROJECTS = Path(__file__).resolve().parent / "data" / "projects"
DENSITY = 1000.0
GRAVITY = 9.80665
ZERO_FLUX = {"type": "flux", "value_type": "constant", "values": [0.0]}


def condition(boundary: int, thermal: dict, hydraulic: dict) -> dict:
    return {"id": boundary, "calculate_thermal": True, "calculate_hydraulic": True,
            "calculate_mechanical": False, "thermal": thermal, "hydraulic": hydraulic}


def discretization(space: str, degree: int, evaluation: str) -> dict:
    return {"space": space, "degree": degree, "coefficient_evaluation": evaluation}


def isothermal(space: str, degree: int, evaluation: str, top: dict, days: float) -> dict:
    """A closed column of water at 6.7 C, with `top` as the condition at z = 0.2 m."""
    return {
        "Basic.json": {"analysis_controls.calculate_thermal": False,
                       "discretization_settings.hydraulic": discretization(space, degree, evaluation)},
        "Conditions.json": {"time_controls.simulation_period": {"unit": "day", "start": 0.0, "end": days},
                            "time_controls.time_stepping.max_step": 1.0e5,
                            "boundary_conditions": [condition(2, ZERO_FLUX, ZERO_FLUX),
                                                    condition(3, ZERO_FLUX, ZERO_FLUX),
                                                    condition(4, ZERO_FLUX, ZERO_FLUX),
                                                    condition(5, ZERO_FLUX, top)]},
        "Output.json": {"history_output.output_interval": {"unit": "day", "value": days},
                        "field_output.output_interval": {"unit": "day", "value": days}},
    }


BASE_PRESSURE = -20000.0
BASE_DIRICHLET = {"type": "dirichlet", "value_type": "constant", "values": [BASE_PRESSURE]}


def water_table(space: str, degree: int, evaluation: str) -> dict:
    """A water column held at a pore pressure at its base, closed elsewhere."""
    case = isothermal(space, degree, evaluation, ZERO_FLUX, 60.0)
    case["Conditions.json"]["boundary_conditions"][1]["hydraulic"] = BASE_DIRICHLET
    return case


LAYER_HEIGHTS = (0.19, 0.12, 0.08, 0.01)
LAYER_RATIO = 4.0


def layers(space: str, degree: int, evaluation: str) -> dict:
    """Steady saturated downward flow through two layers, the upper one LAYER_RATIO times less conductive."""
    inflow = {"type": "flux", "value_type": "constant", "values": [2.0e-7]}
    case = isothermal(space, degree, evaluation, inflow, 1.0)
    case["Conditions.json"]["boundary_conditions"][1]["hydraulic"] = {
        "type": "dirichlet", "value_type": "constant", "values": [10000.0]}
    case["Conditions.json"]["initial_conditions.hydraulic.value"] = 10000.0
    case["Output.json"]["history_output.coordinates"] = [{"x": 0.04, "y": 0.0, "z": z} for z in LAYER_HEIGHTS]
    return case


def conduction() -> dict:
    """Heat alone between 10 C at the top and 2 C at the bottom, the water held still."""
    return {
        "Basic.json": {"analysis_controls.calculate_hydraulic": False},
        "Conditions.json": {"time_controls.simulation_period": {"unit": "day", "start": 0.0, "end": 10.0},
                            "time_controls.time_stepping.max_step": 1.0e5,
                            "boundary_conditions": [
                                condition(2, ZERO_FLUX, ZERO_FLUX),
                                condition(3, {"type": "dirichlet", "value_type": "constant", "values": [2.0]},
                                          ZERO_FLUX),
                                condition(4, ZERO_FLUX, ZERO_FLUX),
                                condition(5, {"type": "dirichlet", "value_type": "constant", "values": [10.0]},
                                          ZERO_FLUX)]},
        "Output.json": {"history_output.output_interval": {"unit": "day", "value": 10.0},
                        "field_output.output_interval": {"unit": "day", "value": 10.0}},
    }


def freezing(hydraulic: dict) -> dict:
    """The column cooled through its top until ice forms there."""
    return {"Basic.json": {"discretization_settings.hydraulic": hydraulic},
            "Conditions.json": {"time_controls.simulation_period": {"unit": "hour", "start": 0.0, "end": 8.0}},
            "Output.json": {"history_output.output_interval": {"unit": "hour", "value": 8.0},
                            "field_output.output_interval": {"unit": "hour", "value": 8.0}}}


def jacobian(hydraulic: dict) -> dict:
    """Two steps of the column started frozen, each assembled Jacobian checked against finite differences."""
    case = freezing(hydraulic)
    case["Basic.json"]["solver_settings.linear_solver.petsc_options"] = "-snes_test_jacobian"
    case["Conditions.json"]["time_controls.simulation_period"] = {"unit": "second", "start": 0.0, "end": 600.0}
    case["Conditions.json"]["initial_conditions.thermal.value"] = -1.0
    return case


def no_advection(case: dict) -> dict:
    """The same case with the heat equation carrying no liquid sensible heat."""
    case["Basic.json"]["analysis_controls.enable_heat_advection"] = False
    return case


def radiation() -> dict:
    """Heat alone: 10 C at the bottom, long-wave exchange with a 0 C sky (the soil stays unfrozen) at the top."""
    case = conduction()
    case["Conditions.json"]["boundary_conditions"][3]["thermal"] = {
        "type": "radiation", "value_type": "constant", "values": [0.95, 0.0]}
    case["Conditions.json"]["boundary_conditions"][1]["thermal"] = {
        "type": "dirichlet", "value_type": "constant", "values": [10.0]}
    case["Output.json"]["history_output.variables"] = ["temperature", "thermal_conductivity"]
    return case


INFLOW = {"type": "flux", "value_type": "constant", "values": [1.0e-8]}

# name: (project, changes, check)
CASES: dict[str, tuple[str, list[dict], str]] = {
    "conduction": ("soil-column", [conduction()], "linear_temperature"),
    "radiation": ("soil-column", [radiation()], "radiation_balance"),
    "hydrostatic_mixed": ("soil-column", [isothermal("mixed", 0, "quadrature", ZERO_FLUX, 60.0)], "hydrostatic"),
    "hydrostatic_lagrange": ("soil-column", [isothermal("lagrange", 1, "quadrature", ZERO_FLUX, 60.0)],
                             "hydrostatic"),
    "water_table_mixed": ("soil-column", [water_table("mixed", 0, "quadrature")], "water_table"),
    "water_table_lagrange": ("soil-column", [water_table("lagrange", 1, "quadrature")], "water_table"),
    "layers_mixed": ("soil-layers", [layers("mixed", 0, "quadrature")], "layers"),
    "layers_lagrange": ("soil-layers", [layers("lagrange", 1, "quadrature")], "layers"),
    "inflow_mixed_lagrange": ("soil-column", [isothermal("mixed", 0, "quadrature", INFLOW, 1.0),
                                              isothermal("lagrange", 1, "quadrature", INFLOW, 1.0)],
                              "same_pressure"),
    "freezing_mixed": ("soil-column", [freezing(discretization("mixed", 0, "quadrature"))], "freezing"),
    "freezing_nodal": ("soil-column", [freezing(discretization("lagrange", 1, "nodal"))], "freezing"),
    "jacobian_mixed": ("soil-column", [jacobian(discretization("mixed", 0, "quadrature"))], "jacobian"),
    "jacobian_lagrange": ("soil-column", [jacobian(discretization("lagrange", 1, "quadrature"))], "jacobian"),
    "jacobian_nodal": ("soil-column", [jacobian(discretization("lagrange", 1, "nodal"))], "jacobian"),
    "jacobian_no_advection": ("soil-column", [no_advection(jacobian(discretization("mixed", 0, "quadrature")))],
                              "jacobian"),
    "atmosphere_surface": ("atmosphere-surface", [{}], "finite"),
    "refuse_mixed_heat": ("soil-column", [{"Basic.json": {"discretization_settings.thermal":
                                                          discretization("mixed", 0, "quadrature")}}],
                          "refused"),
    "refuse_nodal_isothermal": ("soil-column", [isothermal("lagrange", 1, "nodal", ZERO_FLUX, 1.0)], "refused"),
    "refuse_vapour": ("soil-column", [{"Basic.json": {"analysis_controls.enable_vapor_transport": True}}],
                      "refused"),
    "refuse_segregation": ("soil-column", [{"Basic.json": {"materials": None}}], "refused_segregation"),
}


def set_path(document: dict, path: str, value) -> None:
    if path == "materials" and value is None:
        # Switch on segregation in every material of the project.
        for material in document["materials"]:
            material["phase_change"]["equilibrium_model"]["segregation"] = True
        return
    keys = path.split(".")
    for key in keys[:-1]:
        document = document.setdefault(key, {})
    document[keys[-1]] = value


def prepare(project: str, changes: dict, directory: Path) -> Path:
    if directory.exists():
        shutil.rmtree(directory)
    shutil.copytree(PROJECTS / project / "Input", directory / "Input")
    for file_name, edits in changes.items():
        path = directory / "Input" / file_name
        document = json.loads(path.read_text())
        for key, value in edits.items():
            set_path(document, key, value)
        path.write_text(json.dumps(document, indent=4))
    return directory


def run(directory: Path, arguments: argparse.Namespace) -> subprocess.CompletedProcess:
    command = [arguments.mpirun, "-n", str(arguments.ranks), str(arguments.ftcms.resolve())]
    environment = dict(os.environ, FTCMS_PROJECT_PATH=str(directory))
    with open(directory / "run.log", "w") as log:
        return subprocess.run(command, cwd=directory, env=environment, stdout=log, stderr=subprocess.STDOUT)


def last_row(directory: Path, code: str) -> list[float]:
    """The final line of Output/obsf_<code>.dat: time, then one value per observation."""
    rows = [line.split() for line in (directory / "Output" / f"obsf_{code}.dat").read_text().splitlines()
            if line and not line.startswith(("#", "Time"))]
    return [float(value) for value in rows[-1]]


def expect(name: str, actual: float, expected: float, tolerance: float) -> bool:
    ok = math.isfinite(actual) and abs(actual - expected) <= tolerance
    print(f"{name}: expected {expected:.6g} +- {tolerance:.3g}, actual {actual:.6g} -> {'pass' if ok else 'FAIL'}")
    return ok


def check(kind: str, directories: list[Path], results: list[subprocess.CompletedProcess]) -> bool:
    if kind == "refused_segregation":
        kind = "refused"
    if kind == "refused":
        log = (directories[0] / "run.log").read_text()
        ok = results[0].returncode != 0 and ("Not implemented" in log or "not implemented" in log
                                             or "not solved" in log)
        print(f"refusal: exit {results[0].returncode}, message {'found' if ok else 'missing'} -> "
              f"{'pass' if ok else 'FAIL'}")
        return ok
    if any(result.returncode != 0 for result in results):
        print(f"run failed with exit {[result.returncode for result in results]}; see run.log")
        return False
    if kind == "linear_temperature":
        temperature = last_row(directories[0], "T")
        # Observations at z = 0.19, 0.10, 0.01 m between 10 C at 0.2 m and 2 C at 0 m.
        return all(expect(f"T(z={z})", temperature[k + 1], 2.0 + 8.0 * z / 0.2, 2.0e-3)
                   for k, z in enumerate((0.19, 0.10, 0.01)))
    if kind == "radiation_balance":
        temperature = last_row(directories[0], "T")
        conductivity = last_row(directories[0], "TC")
        # Observations at z = 0.19 and 0.10 m; the profile is linear at steady state.
        gradient = (temperature[1] - temperature[2]) / (0.19 - 0.10)
        surface = temperature[1] + gradient * (0.20 - 0.19)
        conducted = conductivity[1] * gradient
        radiated = -0.95 * 5.670374419e-8 * ((surface + 273.15) ** 4 - (0.0 + 273.15) ** 4)
        return expect("conducted heat at the surface [W/m2]", conducted, radiated, 0.02 * abs(radiated))
    if kind == "hydrostatic":
        pressure = last_row(directories[0], "P")
        head_difference = -DENSITY * GRAVITY * (0.19 - 0.01)
        return expect("p(0.19) - p(0.01) [Pa]", pressure[1] - pressure[3], head_difference,
                      0.02 * abs(head_difference))
    if kind == "water_table":
        pressure = last_row(directories[0], "P")
        # Hydrostatic above the base: p(z) = p_base - rho g z at z = 0.19, 0.10, 0.01 m.
        return all(expect(f"p(z={z}) [Pa]", pressure[k + 1], BASE_PRESSURE - DENSITY * GRAVITY * z, 0.03 * DENSITY
                          * GRAVITY * 0.2)
                   for k, z in enumerate((0.19, 0.10, 0.01)))
    if kind == "layers":
        pressure = last_row(directories[0], "P")
        head = [p / (DENSITY * GRAVITY) + z for p, z in zip(pressure[1:], LAYER_HEIGHTS)]
        # One flux through both layers, so the head gradients stand in the inverse ratio of the conductivities.
        upper = (head[0] - head[1]) / (LAYER_HEIGHTS[0] - LAYER_HEIGHTS[1])
        lower = (head[2] - head[3]) / (LAYER_HEIGHTS[2] - LAYER_HEIGHTS[3])
        return expect("upper / lower head gradient", upper / lower, LAYER_RATIO, 0.02 * LAYER_RATIO)
    if kind == "same_pressure":
        first, second = (last_row(directory, "P") for directory in directories)
        rise = abs(first[1] - (-44683.6901))
        return all(expect(f"pressure at observation {k}", first[k], second[k], 0.05 * rise + 1.0)
                   for k in range(1, len(first)))
    if kind == "freezing":
        temperature = last_row(directories[0], "T")
        ice = last_row(directories[0], "Qi")
        ok = all(math.isfinite(value) for value in temperature + ice)
        ok = ok and temperature[1] < temperature[3] and ice[1] > 0.0
        print(f"freezing: top T {temperature[1]:.4g} C below base {temperature[3]:.4g} C, top ice {ice[1]:.4g}"
              f" -> {'pass' if ok else 'FAIL'}")
        return ok
    if kind == "jacobian":
        # PETSc prints "||J - Jfd||_F/||J||_F = <ratio>, ||J - Jfd||_F = <norm>" per assembled Jacobian.
        log = (directories[0] / "run.log").read_text().splitlines()
        ratios = [float(line.split("=")[1].split(",")[0]) for line in log if "||J - Jfd||_F/||J||_F =" in line]
        # The closure tangents are themselves differences, so agreement is to their truncation, not to rounding.
        return bool(ratios) and expect(f"largest ||J - Jfd||/||J|| of {len(ratios)}", max(ratios), 0.0, 1.0e-3)
    if kind == "finite":
        temperature = last_row(directories[0], "T")
        ok = all(math.isfinite(value) for value in temperature)
        print(f"finite observations {temperature} -> {'pass' if ok else 'FAIL'}")
        return ok
    raise ValueError(kind)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("case", choices=sorted(CASES))
    parser.add_argument("--ftcms", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--mpirun", default="mpirun")
    parser.add_argument("--ranks", type=int, default=1)
    arguments = parser.parse_args()

    project, changes, kind = CASES[arguments.case]
    directories, results = [], []
    for index, change in enumerate(changes):
        directory = prepare(project, change,
                            arguments.work.resolve() / f"{arguments.case}_np{arguments.ranks}_{index}")
        directories.append(directory)
        results.append(run(directory, arguments))
    return 0 if check(kind, directories, results) else 1


if __name__ == "__main__":
    sys.exit(main())
