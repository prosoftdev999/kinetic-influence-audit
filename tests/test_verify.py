import json
import math
import os
import resource
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

EXPECTED = json.loads((Path(__file__).parent / "expected.json").read_text())
RESULT = Path(os.environ.get("TASK_RESULT", "/app/result.json"))
SPECIES = ["A0", "A1", "B0", "B1", "C0", "C1", "D0", "D1"]


def _require_keys(obj, keys):
    assert isinstance(obj, dict)
    assert set(keys).issubset(obj)


def assert_sequence_audit(result, expected):
    _require_keys(result, ["vial_order", "carryover_fraction_polar",
                           "carryover_fraction_sticky", "deep_wash_after_acquisition",
                           "chi_square"])
    assert result["vial_order"] == expected["vial_order"]
    assert result["carryover_fraction_polar"] == expected["carryover_fraction_polar"]
    assert result["carryover_fraction_sticky"] == expected["carryover_fraction_sticky"]
    assert result["deep_wash_after_acquisition"] == expected["deep_wash_after_acquisition"]
    score = result["chi_square"]
    assert type(score) in (int, float) and math.isfinite(score)
    assert abs(score - expected["chi_square"]) <= 0.05


def assert_core_campaign(result, expected):
    _require_keys(result, ["audit", "design_run_ids", "sequence_audit",
                           "stability", "parameters", "predictions"])
    assert_sequence_audit(result["sequence_audit"], expected["sequence_audit"])
    audit = result["audit"]
    _require_keys(audit, ["swapped_run_pairs", "time_scale_by_run", "exchange_site"])
    assert audit["swapped_run_pairs"] == expected["audit"]["swapped_run_pairs"]
    assert audit["time_scale_by_run"] == expected["audit"]["time_scale_by_run"]
    assert audit["exchange_site"] == expected["audit"]["exchange_site"]
    assert result["design_run_ids"] == expected["design_run_ids"]
    assert set(result["parameters"]) == set(expected["parameters"])
    for name, wanted in expected["parameters"].items():
        got = result["parameters"][name]
        assert type(got) in (int, float) and math.isfinite(got)
        tolerance = 0.025 * abs(wanted) if name.startswith("ea") else 0.035 * abs(wanted)
        assert abs(got - wanted) <= tolerance
    got_rows = result["predictions"]
    wanted_rows = expected["predictions"]
    assert len(got_rows) == len(wanted_rows) == 72
    wanted_keys = [(r["run_id"], r["time_min"], r["isotopologue"]) for r in wanted_rows]
    assert [(r.get("run_id"), r.get("time_min"), r.get("isotopologue"))
            for r in got_rows] == wanted_keys
    for got, wanted in zip(got_rows, wanted_rows):
        _require_keys(got, ["run_id", "time_min", "isotopologue", "concentration_mmol_L"])
        value = got["concentration_mmol_L"]
        assert type(value) in (int, float) and math.isfinite(value) and value >= -1e-7
        assert abs(value - wanted["concentration_mmol_L"]) <= 0.003
    for run in ("H01", "H02", "H03"):
        for time in (9, 27, 55):
            total = sum(r["concentration_mmol_L"] for r in got_rows
                        if r["run_id"] == run and r["time_min"] == time)
            charge = {"H01": 0.92, "H02": 1.08, "H03": 1.04}[run]
            assert abs(total - charge) <= 0.003


def _swap_map_from_expected(expected):
    mapping = {}
    for left, right in expected["audit"]["swapped_run_pairs"]:
        mapping[left] = right
        mapping[right] = left
    return mapping


def _check_stability(stability, expected_stability, expected, corrected_labels):
    _require_keys(stability, ["most_influential_run_id", "most_influential_pair",
                              "leave_one_trace_out", "leave_two_traces_out"])
    swap_map = _swap_map_from_expected(expected)

    def norm_id(run_id):
        return swap_map.get(run_id, run_id) if corrected_labels else run_id

    submitted_single_ids = [row.get("omitted_run_id") for row in stability["leave_one_trace_out"]]
    assert submitted_single_ids == sorted(submitted_single_ids)
    got_deletions = stability["leave_one_trace_out"]
    wanted_deletions = expected_stability["leave_one_trace_out"]
    assert len(got_deletions) == len(wanted_deletions) == 12
    got_by_id = {}
    for got in got_deletions:
        _require_keys(got, ["omitted_run_id", "parameter_log_shift_l2",
                            "max_prediction_shift_mmol_L", "parameters"])
        key = norm_id(got["omitted_run_id"])
        assert key not in got_by_id
        got_by_id[key] = got
    wanted_by_id = {row["omitted_run_id"]: row for row in wanted_deletions}
    assert set(got_by_id) == set(wanted_by_id)
    for run_id, wanted in wanted_by_id.items():
        got = got_by_id[run_id]
        assert abs(got["parameter_log_shift_l2"] - wanted["parameter_log_shift_l2"]) <= 2e-4
        assert abs(got["max_prediction_shift_mmol_L"] -
                   wanted["max_prediction_shift_mmol_L"]) <= 2e-5
        assert set(got["parameters"]) == set(wanted["parameters"])
        for name, target in wanted["parameters"].items():
            assert abs(got["parameters"][name] - target) <= 0.004 * abs(target)

    submitted_pairs = []
    for row in stability["leave_two_traces_out"]:
        pair = row.get("omitted_run_ids")
        assert isinstance(pair, list) and len(pair) == 2 and pair == sorted(pair)
        submitted_pairs.append(pair)
    assert submitted_pairs == sorted(submitted_pairs)
    got_pairs = stability["leave_two_traces_out"]
    wanted_pairs = expected_stability["leave_two_traces_out"]
    assert len(got_pairs) == len(wanted_pairs) == 66
    got_by_pair = {}
    for got in got_pairs:
        _require_keys(got, ["omitted_run_ids", "parameter_log_shift_l2",
                            "max_prediction_shift_mmol_L", "interaction_log_shift_excess",
                            "interaction_prediction_shift_excess_mmol_L", "parameters"])
        key = tuple(sorted(norm_id(x) for x in got["omitted_run_ids"]))
        assert key not in got_by_pair
        got_by_pair[key] = got
    wanted_by_pair = {tuple(row["omitted_run_ids"]): row for row in wanted_pairs}
    assert set(got_by_pair) == set(wanted_by_pair)
    for pair, wanted in wanted_by_pair.items():
        got = got_by_pair[pair]
        for field, tolerance in [
            ("parameter_log_shift_l2", 2e-4),
            ("max_prediction_shift_mmol_L", 2e-5),
            ("interaction_log_shift_excess", 3e-4),
            ("interaction_prediction_shift_excess_mmol_L", 3e-5),
        ]:
            assert abs(got[field] - wanted[field]) <= tolerance
        assert set(got["parameters"]) == set(wanted["parameters"])
        for name, target in wanted["parameters"].items():
            assert abs(got["parameters"][name] - target) <= 0.004 * abs(target)

    assert norm_id(stability["most_influential_run_id"]) == expected_stability["most_influential_run_id"]
    got_most_pair = sorted(norm_id(x) for x in stability["most_influential_pair"])
    assert got_most_pair == expected_stability["most_influential_pair"]


def test_result_matches_campaign():
    assert RESULT.is_file()
    result = json.loads(RESULT.read_text())
    assert_core_campaign(result, EXPECTED)
    stability = result["stability"]
    expected_stability = EXPECTED["stability"]

    # The campaign is scientifically identical whether a deletion row is labeled by
    # the exported trace ID or, after reconciliation, by that trace's true reactor ID.
    # Accept either convention when it is used consistently across the full table.
    accepted = False
    for corrected_labels in (False, True):
        try:
            _check_stability(stability, expected_stability, EXPECTED, corrected_labels)
            accepted = True
            break
        except (AssertionError, KeyError, TypeError, ValueError):
            pass
    assert accepted

def _drop_solver_privileges():
    resource.setrlimit(resource.RLIMIT_CPU, (300, 300))
    resource.setrlimit(resource.RLIMIT_AS, (3_000_000_000, 3_000_000_000))
    resource.setrlimit(resource.RLIMIT_FSIZE, (2_000_000, 2_000_000))
    resource.setrlimit(resource.RLIMIT_NOFILE, (64, 64))
    os.setgroups([])
    os.setgid(65534)
    os.setuid(65534)


def test_solver_generalizes_to_sealed_export_incident():
    solver = Path(os.environ.get("TASK_SOLVER", "/app/solver.py"))
    assert solver.is_file() and solver.stat().st_size < 100_000
    case_root = Path(tempfile.mkdtemp(prefix="kinetic-case-"))
    data_dir = case_root / "data"
    output = case_root / "result.json"
    shutil.copytree(Path(__file__).parent / "hidden_data", data_dir)
    local_verifier = os.environ.get("LOCAL_VERIFIER") == "1"
    if not local_verifier:
        for path in [case_root, data_dir, *data_dir.iterdir()]:
            os.chown(path, 65534, 65534)
        case_root.chmod(0o700)
        data_dir.chmod(0o500)
        for path in data_dir.iterdir():
            path.chmod(0o400)
    environment = {
        "PATH": "/usr/local/bin:/usr/bin:/bin",
        "HOME": str(case_root),
        "PYTHONNOUSERSITE": "1",
        "TASK_DATA_DIR": str(data_dir),
        "TASK_OUTPUT": str(output),
        "TASK_REPLAY": "1",
    }
    completed = subprocess.run(
        [sys.executable if local_verifier else "python", str(solver)],
        cwd=case_root, env=environment,
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        timeout=330, check=False,
        preexec_fn=None if local_verifier else _drop_solver_privileges,
    )
    assert completed.returncode == 0 and output.is_file()
    result = json.loads(output.read_text())
    hidden_expected = json.loads((Path(__file__).parent / "hidden_expected.json").read_text())
    assert_core_campaign(result, hidden_expected)
