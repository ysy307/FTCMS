"""Reader equivalence and failure tests for the Mizoguchi validation scripts.

Two output layouts have to be measured with one instrument: one ``Out_NNN.vtu``
per output step, or one ``Out_NNN.h5`` per step in PETSc's visualisation layout
(`file_format: hdf5`). Every accuracy number in this project is only comparable
across those two if the *reader* is the only thing that differs and every
downstream step -- resampling, depth binning, metrics -- is shared.

These tests therefore assert two things:

1. equivalence: the same field on the same mesh, written once as VTU and once
   as HDF5, produces identical depth profiles and identical profile metrics;
2. failure: every way of asking for a snapshot that does not exist, or of
   asking ambiguously, raises rather than returning an empty-but-successful
   result.

Run with:
    .venv/bin/python -m unittest discover -s scripts/validation/tests -v
"""

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np
import pyvista as pv

SCRIPT_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = SCRIPT_DIR.parents[1]
sys.path.insert(0, str(SCRIPT_DIR))

import Mizo_theta as mt  # noqa: E402
import mizo_compare as mc  # noqa: E402

QGATE_ARCHIVE = (
    REPO_ROOT / "project/Mizo-xz-Convex/Output_archive_20260818_qgate_theta0330"
)


def _column_mesh(nx=5, nz=41, width=0.08, height=0.20):
    """Points and triangles of a rectangular x-z column, y identically zero."""
    xs = np.linspace(0.0, width, nx)
    zs = np.linspace(0.0, height, nz)
    points = np.array([[x, 0.0, z] for z in zs for x in xs], dtype=float)
    cells = []
    for k in range(nz - 1):
        for i in range(nx - 1):
            a = k * nx + i
            b = a + 1
            c = a + nx
            d = c + 1
            cells.append([a, b, d])
            cells.append([a, d, c])
    return points, np.asarray(cells, dtype=np.int64)


def _synthetic_fields(points, height=0.20):
    """A depth-dependent water/ice split with a sign change, as the column has."""
    depth = height - points[:, 2]
    water = 0.30 + 0.25 * depth
    ice = np.where(depth < 0.06, 0.10 - 1.2 * depth, 0.0)
    return water, np.maximum(ice, 0.0)


def _write_vtu(directory, points, cells, snapshots, *, stem="Out", fields=True):
    """Write one Out_NNN.vtu per snapshot, in the archived layout."""
    directory.mkdir(parents=True, exist_ok=True)
    connectivity = np.hstack(
        [np.full((len(cells), 1), cells.shape[1], dtype=np.int64), cells]
    ).ravel()
    written = []
    for index, (water, ice) in snapshots.items():
        grid = pv.UnstructuredGrid(
            connectivity,
            np.full(len(cells), pv.CellType.TRIANGLE, dtype=np.uint8),
            points,
        )
        if fields:
            grid.point_data["WaterContent"] = water
            grid.point_data["IceContent"] = ice
        path = directory / f"{stem}_{index:03d}.vtu"
        grid.save(path)
        written.append(path)
    return written


def _write_hdf5(directory, points, cells, snapshots, *, fields=("water_content", "ice_content")):
    """Write one Out_NNN.h5 per snapshot, in the layout FTCMS writes through PETSc."""
    directory.mkdir(parents=True, exist_ok=True)
    written = []
    for index, (water, ice) in snapshots.items():
        path = directory / f"Out_{index:03d}.h5"
        with h5py.File(path, "w") as handle:
            handle.create_dataset("geometry/vertices", data=points)
            handle.create_dataset("viz/topology/cells", data=cells.astype(np.int32))
            if "water_content" in fields:
                handle.create_dataset("vertex_fields/fields_water_content", data=water)
            if "ice_content" in fields:
                handle.create_dataset("vertex_fields/fields_ice_content", data=ice)
        written.append(path)
    return written


class ReaderEquivalenceTest(unittest.TestCase):
    """The VTU and the HDF5 reader must be interchangeable."""

    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        root = Path(cls._tmp.name)
        points, cells = _column_mesh()
        water0, ice0 = _synthetic_fields(points)
        snapshots = {0: (water0, ice0), 24: (water0 + 0.01, ice0)}
        cls.points, cls.cells, cls.snapshots = points, cells, snapshots
        cls.vtu_dir = root / "vtu" / "Files"
        _write_vtu(cls.vtu_dir, points, cells, snapshots)
        cls.h5_dir = root / "h5" / "Files"
        _write_hdf5(cls.h5_dir, points, cells, snapshots)

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def _profile(self, directory, snapshot_format, index, time_hours):
        source = mt.resolve_snapshot_source(directory, snapshot_format)
        grid = source.read_grid(index)
        depths, theta = mt._sample_total_water(grid, 0.001)
        return mt._depth_profile(depths, theta, time_hours, mt.DZ_DEFAULT)

    def test_formats_are_detected(self):
        self.assertEqual(mt.resolve_snapshot_source(self.vtu_dir, "auto").format, "vtu")
        self.assertEqual(mt.resolve_snapshot_source(self.h5_dir, "auto").format, "hdf5")

    def test_depth_profiles_are_identical(self):
        vtu_rows = self._profile(self.vtu_dir, "vtu", 24, 24)
        h5_rows = self._profile(self.h5_dir, "hdf5", 24, 24)
        self.assertEqual(len(vtu_rows), len(h5_rows))
        for a, b in zip(vtu_rows, h5_rows):
            self.assertIsNotNone(a["simulated"])
            self.assertIsNotNone(b["simulated"])
            self.assertAlmostEqual(a["simulated"], b["simulated"], places=12)

    def test_profile_metrics_are_identical(self):
        vtu = mt._profile_metrics(self._profile(self.vtu_dir, "vtu", 24, 24))
        h5 = mt._profile_metrics(self._profile(self.h5_dir, "hdf5", 24, 24))
        for key in ("rmse", "mae", "bias", "max_abs_error"):
            self.assertAlmostEqual(vtu[key], h5[key], places=12, msg=key)

    def test_snapshot_counts_and_indices(self):
        vtu = mt.resolve_snapshot_source(self.vtu_dir, "vtu")
        h5 = mt.resolve_snapshot_source(self.h5_dir, "hdf5")
        self.assertEqual(vtu.available_indices(), [0, 24])
        self.assertEqual(h5.snapshot_count(), 2)

    def test_source_description_records_provenance(self):
        source = mt.resolve_snapshot_source(self.vtu_dir, "vtu")
        described = source.describe(24)
        self.assertEqual(described["format"], "vtu")
        self.assertEqual(described["snapshot_index"], 24)
        self.assertIn("Out_024.vtu", described["snapshot_path"])


class SelectionFailureTest(unittest.TestCase):
    """Ambiguous or absent output must raise, never resolve silently."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.root = Path(self._tmp.name)
        self.points, self.cells = _column_mesh(nx=3, nz=5)
        water, ice = _synthetic_fields(self.points)
        self.snapshots = {0: (water, ice)}

    def tearDown(self):
        self._tmp.cleanup()

    def test_empty_directory_is_a_failure(self):
        empty = self.root / "empty"
        empty.mkdir()
        with self.assertRaises(mt.SnapshotError):
            mt.resolve_snapshot_source(empty, "auto")

    def test_missing_directory_is_a_failure(self):
        with self.assertRaises(mt.SnapshotError):
            mt.resolve_snapshot_source(self.root / "absent", "auto")

    def test_duplicate_vtu_indices_are_ambiguous(self):
        directory = self.root / "duplicates"
        _write_vtu(directory, self.points, self.cells, self.snapshots)
        _write_vtu(directory / "Files", self.points, self.cells, self.snapshots)
        with self.assertRaises(mt.SnapshotError):
            mt.resolve_snapshot_source(directory, "vtu")

    def test_both_formats_present_is_ambiguous(self):
        both = self.root / "both" / "Files"
        _write_vtu(both, self.points, self.cells, self.snapshots)
        _write_hdf5(both, self.points, self.cells, self.snapshots)
        with self.assertRaises(mt.SnapshotError) as raised:
            mt.resolve_snapshot_source(both, "auto")
        self.assertIn("--snapshot_format", str(raised.exception))
        # Naming the format resolves the ambiguity rather than guessing.
        self.assertEqual(mt.resolve_snapshot_source(both, "vtu").format, "vtu")
        self.assertEqual(mt.resolve_snapshot_source(both, "hdf5").format, "hdf5")

    def test_requesting_an_absent_format_is_a_failure(self):
        vtu_only = self.root / "vtu_only" / "Files"
        _write_vtu(vtu_only, self.points, self.cells, self.snapshots)
        with self.assertRaises(mt.SnapshotError):
            mt.resolve_snapshot_source(vtu_only, "hdf5")

    def test_duplicate_hdf5_indices_are_ambiguous(self):
        directory = self.root / "duplicates_h5"
        _write_hdf5(directory, self.points, self.cells, self.snapshots)
        _write_hdf5(directory / "Files", self.points, self.cells, self.snapshots)
        with self.assertRaises(mt.SnapshotError):
            mt.resolve_snapshot_source(directory, "hdf5")


class SnapshotAvailabilityTest(unittest.TestCase):
    """Asking for a time the run never wrote must fail with what is available."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.root = Path(self._tmp.name)
        self.points, self.cells = _column_mesh(nx=3, nz=5)
        water, ice = _synthetic_fields(self.points)
        self.snapshots = {0: (water, ice), 1: (water, ice)}

    def tearDown(self):
        self._tmp.cleanup()

    def test_vtu_missing_index_reports_what_exists(self):
        directory = self.root / "vtu" / "Files"
        _write_vtu(directory, self.points, self.cells, self.snapshots)
        source = mt.resolve_snapshot_source(directory, "vtu")
        with self.assertRaises(mt.SnapshotError) as raised:
            source.read_grid(24)
        self.assertIn("24", str(raised.exception))

    def test_hdf5_missing_index_reports_what_exists(self):
        directory = self.root / "h5" / "Files"
        _write_hdf5(directory, self.points, self.cells, self.snapshots)
        source = mt.resolve_snapshot_source(directory, "hdf5")
        with self.assertRaises(mt.SnapshotError) as raised:
            source.read_grid(24)
        self.assertIn("2", str(raised.exception))

    def test_vtu_without_water_content_is_a_failure(self):
        directory = self.root / "bare" / "Files"
        _write_vtu(directory, self.points, self.cells, self.snapshots, fields=False)
        source = mt.resolve_snapshot_source(directory, "vtu")
        with self.assertRaises((mt.SnapshotError, ValueError)):
            source.read_grid(0)

    def test_hdf5_without_ice_content_is_a_failure(self):
        directory = self.root / "partial" / "Files"
        _write_hdf5(directory, self.points, self.cells, self.snapshots, fields=("water_content",))
        source = mt.resolve_snapshot_source(directory, "hdf5")
        with self.assertRaises((mt.SnapshotError, ValueError)):
            source.read_grid(0)


class CommandLineFailureTest(unittest.TestCase):
    """A run that cannot be evaluated must exit non-zero and emit no clean JSON."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.root = Path(self._tmp.name)
        points, cells = _column_mesh(nx=3, nz=9)
        water, ice = _synthetic_fields(points)
        self.directory = self.root / "run" / "Output" / "Files"
        _write_vtu(self.directory, points, cells, {0: (water, ice)})

    def tearDown(self):
        self._tmp.cleanup()

    def _run(self, *extra):
        return subprocess.run(
            [sys.executable, str(SCRIPT_DIR / "Mizo_theta.py"),
             "--vtu_dir", str(self.directory),
             "--output_interval_minutes", "30", *extra],
            capture_output=True, text=True,
        )

    def test_unreachable_time_exits_non_zero(self):
        result = self._run("--time", "12")
        self.assertNotEqual(result.returncode, 0, msg=result.stdout)
        self.assertIn("not evaluated", result.stdout)

    def test_no_successful_json_is_written_when_nothing_evaluates(self):
        target = self.root / "out.json"
        result = self._run("--time", "12", "--json_output", str(target))
        self.assertNotEqual(result.returncode, 0, msg=result.stdout)
        if target.exists():
            payload = json.loads(target.read_text())
            self.assertTrue(payload.get("failures"), msg=payload)
            self.assertFalse(payload.get("times"), msg=payload)

    def test_evaluated_json_records_source_and_initial_condition(self):
        target = self.root / "ok.json"
        result = self._run("--time", "0", "--json_output", str(target))
        self.assertEqual(result.returncode, 0, msg=result.stdout + result.stderr)
        payload = json.loads(target.read_text())
        self.assertIn("source", payload)
        self.assertEqual(payload["source"]["format"], "vtu")
        self.assertIn("experimental_reference", payload)


class CompareGuardTest(unittest.TestCase):
    """mizo_compare must not report OK when it compared nothing."""

    def test_missing_metrics_are_reported_as_failures(self):
        base = {
            "profile_metrics": {"rmse": 0.02, "bias": 0.0, "max_abs_error": 0.05},
            "skill_vs_no_redistribution": 0.7,
        }
        rows, failures = mc.compare_result("12h", base, {})
        self.assertFalse(rows)
        self.assertTrue(failures)

    def test_matching_metrics_produce_no_failure(self):
        base = {"profile_metrics": {"rmse": 0.02}}
        rows, failures = mc.compare_result("12h", base, base)
        self.assertTrue(rows)
        self.assertFalse(failures)

    def test_different_initial_conditions_are_not_comparable(self):
        baseline = {"initial_pressure_Pa": -54065.7412, "initial_temperature_C": 6.7}
        current = {"initial_pressure_Pa": -44683.6901, "initial_temperature_C": 6.7}
        self.assertTrue(mc.initial_condition_differences(baseline, current))
        self.assertFalse(mc.initial_condition_differences(baseline, baseline))

    def test_nonfinite_metrics_are_failures(self):
        for value in (float("nan"), float("inf"), -float("inf")):
            with self.subTest(value=value):
                base = {"profile_metrics": {"rmse": 0.02}}
                cur = {"profile_metrics": {"rmse": value}}
                self.assertTrue(mc.compare_result("12h", base, cur)[1])

    def test_absent_baseline_time_is_a_failure(self):
        rows, failures = mc.compare_result("99h", None, {"profile_metrics": {"rmse": 0.02}})
        self.assertTrue(failures)


@unittest.skipUnless(QGATE_ARCHIVE.is_dir(), "archived qgate run not present")
class ArchiveReadableTest(unittest.TestCase):
    """The archived VTU run must be readable and carry its stated initial state."""

    def test_initial_snapshot_is_uniform_theta0330(self):
        source = mt.resolve_snapshot_source(QGATE_ARCHIVE / "Files", "auto")
        self.assertEqual(source.format, "vtu")
        self.assertEqual(source.snapshot_count(), 101)
        grid = source.read_grid(0)
        depths, theta = mt._sample_total_water(grid, 0.002)
        rows = mt._depth_profile(depths, theta, 0, mt.DZ_DEFAULT)
        values = [row["simulated"] for row in rows]
        self.assertTrue(all(value is not None for value in values))
        self.assertAlmostEqual(min(values), 0.330, places=3)
        self.assertAlmostEqual(max(values), 0.330, places=3)


if __name__ == "__main__":
    unittest.main()
