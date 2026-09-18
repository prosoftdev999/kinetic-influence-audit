from __future__ import annotations

import csv
import itertools
import json
import os
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.optimize import least_squares
from scipy.optimize import linear_sum_assignment

DATA = Path(os.environ.get("TASK_DATA_DIR", "/app/data"))
OUTPUT = Path(os.environ.get("TASK_OUTPUT", "/app/result.json"))
REPLAY = os.environ.get("TASK_REPLAY") == "1"
SPECIES = ["A0", "A1", "B0", "B1", "C0", "C1", "D0", "D1"]


def rows(name):
    with (DATA / name).open() as handle:
        return list(csv.DictReader(handle))


cal = rows("calibration.csv")
rf = {}
for row in cal:
    key = row["session"], row["compound"]
    rf.setdefault(key, []).append(float(row["peak_area"]) / float(row["known_mmol_L"]))
rf = {key: float(np.mean(value)) for key, value in rf.items()}
manifest = {row["run_id"]: row for row in rows("run_manifest.csv")}
obs_rows = rows("chromatograms.csv")
observed = {}
for row in obs_rows:
    run, time, species = row["run_id"], float(row["time_min"]), row["isotopologue"]
    observed.setdefault(run, {}).setdefault(time, {})[species] = float(row["peak_area"])


def integrate(theta, run, times, exchange_site, time_scale=1.0):
    p = np.exp(theta)
    record = manifest[run]
    temp, catalyst = float(record["temperature_K"]), float(record["catalyst_mmol_L"])
    initial_frac = float(record["initial_A_label_fraction"])
    solvent_frac = float(record["solvent_label_fraction"])
    charge = float(record["initial_A_mmol_L"])
    gas = 0.008314462618
    rate = p[:5] * catalyst * np.exp(-p[5:] / gas * (1/temp - 1/328.15))
    k1, k2, k3, k4, kx = rate
    matrix = np.zeros((8, 8))
    matrix[0,0] = matrix[1,1] = -(k1+k2)
    matrix[2,0], matrix[3,1] = k1, k1
    matrix[2,2] = matrix[3,3] = -k3
    matrix[4,0], matrix[5,1] = k2, k2
    matrix[4,4] = matrix[5,5] = -k4
    matrix[6,2], matrix[6,4] = k3, k4
    matrix[7,3], matrix[7,5] = k3, k4
    offset = 2 if exchange_site == "B" else 4
    matrix[offset,offset] -= kx*solvent_frac
    matrix[offset,offset+1] += kx*(1-solvent_frac)
    matrix[offset+1,offset] += kx*solvent_frac
    matrix[offset+1,offset+1] -= kx*(1-solvent_frac)
    y0 = [charge*(1-initial_frac), charge*initial_frac, 0, 0, 0, 0, 0, 0]
    return np.array([expm(matrix*time*time_scale) @ y0 for time in times])


def residual(theta, assignment, time_scales, exchange_site, omitted_runs=()):
    if isinstance(omitted_runs, str):
        omitted_runs = (omitted_runs,)
    result = []
    for run in sorted(observed):
        if run in omitted_runs:
            continue
        times = sorted(observed[run])
        condition_run = assignment[run]
        predicted = integrate(theta, condition_run, times, exchange_site, time_scales[run])
        session = manifest[condition_run]["session"]
        for time, values in zip(times, predicted):
            result.extend((values[i] - observed[run][time][s] / rf[(session, s[0])]) / 0.001
                          for i, s in enumerate(SPECIES))
    return result


def assign_runs(theta, exchange_site):
    runs = sorted(observed)
    costs = np.zeros((len(runs), len(runs)))
    scale_choice = {}
    for row_index, exported_run in enumerate(runs):
        times = sorted(observed[exported_run])
        for col_index, condition_run in enumerate(runs):
            session = manifest[condition_run]["session"]
            measured = np.array([[observed[exported_run][time][s] / rf[(session, s[0])]
                                  for s in SPECIES] for time in times])
            choices = []
            for scale in (0.80, 1.00, 1.25):
                predicted = integrate(theta, condition_run, times, exchange_site, scale)
                choices.append((float(np.sum(((predicted - measured) / 0.001) ** 2)), scale))
            costs[row_index, col_index], scale_choice[exported_run, condition_run] = min(choices)
    rows_found, cols_found = linear_sum_assignment(costs)
    assignment = {runs[row]: runs[col] for row, col in zip(rows_found, cols_found)}
    scales = {run: scale_choice[run, condition] for run, condition in assignment.items()}
    return assignment, scales


start = np.log([0.04, 0.035, 0.018, 0.03, 0.015, 40, 50, 35, 45, 25])
lower = np.log([0.002]*5 + [10]*5)
upper = np.log([0.20]*5 + [100]*5)
candidates = []
for site in ("B", "C"):
    theta = start.copy()
    assignment, time_scales = assign_runs(theta, site)
    for _iteration in range(6):
        fit_try = least_squares(residual, theta, bounds=(lower, upper),
                                args=(assignment, time_scales, site), max_nfev=800,
                                xtol=2e-10, ftol=2e-10, gtol=2e-10)
        updated, updated_scales = assign_runs(fit_try.x, site)
        theta = fit_try.x
        if updated == assignment and updated_scales == time_scales:
            break
        assignment, time_scales = updated, updated_scales
    final_residual = residual(theta, assignment, time_scales, site)
    candidates.append((float(np.dot(final_residual, final_residual)), site,
                       assignment, time_scales, fit_try))
candidates.sort(key=lambda item: item[0])
if os.environ.get("KINETIC_DIAGNOSTICS"):
    for score, site, assignment, time_scales, _fit in candidates:
        print(f"site={site} scaled_sse={score:.6f} assignment={assignment} scales={time_scales}")
_, exchange_site, assignment, time_scales, fit = candidates[0]
swapped_pairs = []
for exported_run, condition_run in sorted(assignment.items()):
    if exported_run != condition_run and exported_run < condition_run:
        if assignment.get(condition_run) != exported_run:
            raise RuntimeError("best assignment is not a set of reciprocal swaps")
        swapped_pairs.append([exported_run, condition_run])
if len(swapped_pairs) != 3:
    raise RuntimeError(f"expected three swapped pairs, recovered {swapped_pairs}")
nonunit_scales = {run: scale for run, scale in sorted(time_scales.items()) if scale != 1.0}
if sorted(nonunit_scales.values()) != [0.8, 1.25]:
    raise RuntimeError(f"expected one clock fault of each kind, recovered {nonunit_scales}")
values = np.exp(fit.x)
names = ["k1_ref", "k2_ref", "k3_ref", "k4_ref", "kx_ref", "ea1_kJ_mol",
         "ea2_kJ_mol", "ea3_kJ_mol", "ea4_kJ_mol", "eax_kJ_mol"]

# A deletion diagnostic is reported because a global fit can look convincing while
# being carried by one chromatogram.  The discrete campaign repair is held fixed;
# each exported trace is removed in turn and only the shared kinetic quantities
# are refitted.  Prediction displacement is evaluated on the complete withheld
# schedule, not on the training rows used by the refit.
schedule = rows("prediction_times.csv")
prediction_grid = []
for run in sorted({row["run_id"] for row in schedule}):
    times = sorted(float(row["time_min"]) for row in schedule if row["run_id"] == run)
    prediction_grid.extend(integrate(fit.x, run, times, exchange_site).ravel())
prediction_grid = np.asarray(prediction_grid)


def withheld_grid(theta):
    grid = []
    for run in sorted({row["run_id"] for row in schedule}):
        times = sorted(float(row["time_min"]) for row in schedule if row["run_id"] == run)
        grid.extend(integrate(theta, run, times, exchange_site).ravel())
    return np.asarray(grid)

stability_rows = []
pair_stability_rows = []
most_influential_run_id = None
most_influential_pair = None
if not REPLAY:
    single_fits = {}
    for omitted_run in sorted(observed):
        deletion_fit = least_squares(
            residual, fit.x, bounds=(lower, upper),
            args=(assignment, time_scales, exchange_site, omitted_run), max_nfev=800,
            xtol=2e-10, ftol=2e-10, gtol=2e-10,
        )
        deletion_values = np.exp(deletion_fit.x)
        deletion_grid = withheld_grid(deletion_fit.x)
        single_fits[omitted_run] = deletion_fit.x
        stability_rows.append({
            "omitted_run_id": omitted_run,
            "parameter_log_shift_l2": float(np.linalg.norm(deletion_fit.x - fit.x)),
            "max_prediction_shift_mmol_L": float(
                np.max(np.abs(np.asarray(deletion_grid) - prediction_grid))
            ),
            "parameters": dict(zip(names, map(float, deletion_values))),
        })

    most_influential_run_id = max(
        stability_rows,
        key=lambda row: (row["parameter_log_shift_l2"] +
                         row["max_prediction_shift_mmol_L"] / 0.001,
                         row["omitted_run_id"]),
    )["omitted_run_id"]

    single_rows = {row["omitted_run_id"]: row for row in stability_rows}
    for first_run, second_run in itertools.combinations(sorted(observed), 2):
        pair_start = (single_fits[first_run] + single_fits[second_run]) / 2
        pair_fit = least_squares(
            residual, pair_start, bounds=(lower, upper),
            args=(assignment, time_scales, exchange_site, (first_run, second_run)),
            max_nfev=800, xtol=2e-10, ftol=2e-10, gtol=2e-10,
        )
        pair_log_shift = float(np.linalg.norm(pair_fit.x - fit.x))
        pair_prediction_shift = float(
            np.max(np.abs(withheld_grid(pair_fit.x) - prediction_grid))
        )
        pair_stability_rows.append({
            "omitted_run_ids": [first_run, second_run],
            "parameter_log_shift_l2": pair_log_shift,
            "max_prediction_shift_mmol_L": pair_prediction_shift,
            "interaction_log_shift_excess": float(
                pair_log_shift
                - single_rows[first_run]["parameter_log_shift_l2"]
                - single_rows[second_run]["parameter_log_shift_l2"]
            ),
            "interaction_prediction_shift_excess_mmol_L": float(
                pair_prediction_shift
                - single_rows[first_run]["max_prediction_shift_mmol_L"]
                - single_rows[second_run]["max_prediction_shift_mmol_L"]
            ),
            "parameters": dict(zip(names, map(float, np.exp(pair_fit.x)))),
        })

    most_influential_pair = max(
        pair_stability_rows,
        key=lambda row: (row["parameter_log_shift_l2"] +
                         row["max_prediction_shift_mmol_L"] / 0.001,
                         row["omitted_run_ids"]),
    )["omitted_run_ids"]

design_candidates = rows("design_candidates.csv")
design_scenarios = rows("design_scenarios.csv")
for candidate in design_candidates:
    manifest[candidate["candidate_id"]] = candidate
design_times = [6.0, 18.0, 42.0]
step = 0.0001
information = {}
for scenario in design_scenarios:
    factors = np.array([float(scenario[key]) for key in
                        ["k1_ref", "k2_ref", "k3_ref", "k4_ref", "kx_ref",
                         "ea1", "ea2", "ea3", "ea4", "eax"]])
    scenario_theta = np.log(values * factors)
    for candidate in design_candidates:
        candidate_id = candidate["candidate_id"]
        columns = []
        for index in range(10):
            upper_theta = scenario_theta.copy()
            lower_theta = scenario_theta.copy()
            upper_theta[index] += step
            lower_theta[index] -= step
            upper_state = integrate(upper_theta, candidate_id, design_times, exchange_site).ravel()
            lower_state = integrate(lower_theta, candidate_id, design_times, exchange_site).ravel()
            columns.append((upper_state - lower_state) / (2 * step))
        jacobian = np.column_stack(columns)
        information[candidate_id, scenario["scenario_id"]] = jacobian.T @ jacobian / 0.002**2

best_score = -np.inf
design_run_ids = None
design_rankings = []
for choice in itertools.combinations(design_candidates, 4):
    catalyst_total = sum(float(row["catalyst_mmol_L"]) for row in choice)
    temperatures = [float(row["temperature_K"]) for row in choice]
    solvent_fractions = {float(row["solvent_label_fraction"]) for row in choice}
    if catalyst_total > 3.80 + 1e-12 or max(temperatures) - min(temperatures) < 15.0:
        continue
    if len(solvent_fractions) < 2:
        continue
    ids = [row["candidate_id"] for row in choice]
    scenario_scores = []
    for scenario in design_scenarios:
        matrix = np.eye(10)
        for candidate_id in ids:
            matrix += information[candidate_id, scenario["scenario_id"]]
        sign, logdet = np.linalg.slogdet(matrix)
        if sign <= 0:
            raise RuntimeError("non-positive design information matrix")
        scenario_scores.append(float(logdet))
    score = min(scenario_scores)
    design_rankings.append((score, ids))
    if score > best_score + 1e-10 or (abs(score-best_score) <= 1e-10 and
                                     (design_run_ids is None or ids < design_run_ids)):
        best_score, design_run_ids = score, ids
if os.environ.get("KINETIC_DIAGNOSTICS"):
    for score, ids in sorted(design_rankings, key=lambda item: (-item[0], item[1]))[:5]:
        print(f"design={ids} worst_case_logdet={score:.12f}")



def solve_autosampler_sequence():
    nominal_rows = rows("autosampler_nominals.csv")
    sequence_rows = rows("autosampler_sequence.csv")
    channel_names = ["polar_1", "polar_2", "polar_3",
                     "sticky_1", "sticky_2", "sticky_3"]
    vial_ids = [row["vial_id"] for row in nominal_rows]
    nominal = np.asarray([[float(row[name]) for name in channel_names]
                          for row in nominal_rows], dtype=float)
    measured = np.asarray([[float(row[name]) for name in channel_names]
                           for row in sequence_rows], dtype=float)
    if len(vial_ids) != 10 or measured.shape != (10, 6):
        raise RuntimeError("unexpected autosampler qualification dimensions")

    sigma = 0.003
    grid = (0.18, 0.24, 0.30, 0.36, 0.42)
    n = len(vial_ids)
    full_mask = (1 << n) - 1
    infinity = float("inf")

    def best_for_fractions(rho_polar, rho_sticky, keep_path=False):
        no_wash = np.full((1 << n, n), infinity)
        used_wash = np.full((1 << n, n), infinity)
        predecessor_no_wash = {} if keep_path else None
        predecessor_used_wash = {} if keep_path else None

        for vial in range(n):
            no_wash[1 << vial, vial] = float(
                np.sum(((measured[0] - nominal[vial]) / sigma) ** 2)
            )

        for mask in range(1, 1 << n):
            position = mask.bit_count()
            if position >= n:
                continue
            for previous in range(n):
                cost_no_wash = no_wash[mask, previous]
                cost_used_wash = used_wash[mask, previous]
                if not np.isfinite(cost_no_wash) and not np.isfinite(cost_used_wash):
                    continue
                carry = nominal[previous].copy()
                carry[:3] *= rho_polar
                carry[3:] *= rho_sticky
                for current in range(n):
                    if mask & (1 << current):
                        continue
                    next_mask = mask | (1 << current)
                    ordinary = float(
                        np.sum(((measured[position] - nominal[current] - carry) / sigma) ** 2)
                    )
                    if cost_no_wash + ordinary < no_wash[next_mask, current]:
                        no_wash[next_mask, current] = cost_no_wash + ordinary
                        if keep_path:
                            predecessor_no_wash[next_mask, current] = (
                                mask, previous, False, False
                            )
                    if cost_used_wash + ordinary < used_wash[next_mask, current]:
                        used_wash[next_mask, current] = cost_used_wash + ordinary
                        if keep_path:
                            predecessor_used_wash[next_mask, current] = (
                                mask, previous, True, False
                            )

                    # Exactly one deep wash is required.  Using it on this boundary
                    # removes only the carryover into the current acquisition.
                    washed = float(
                        np.sum(((measured[position] - nominal[current]) / sigma) ** 2)
                    )
                    if cost_no_wash + washed < used_wash[next_mask, current]:
                        used_wash[next_mask, current] = cost_no_wash + washed
                        if keep_path:
                            predecessor_used_wash[next_mask, current] = (
                                mask, previous, False, True
                            )

        final_vial = int(np.argmin(used_wash[full_mask]))
        score = float(used_wash[full_mask, final_vial])
        if not keep_path:
            return score

        order = [final_vial]
        mask = full_mask
        previous = final_vial
        wash_already_used = True
        wash_after = None
        while mask.bit_count() > 1:
            table = predecessor_used_wash if wash_already_used else predecessor_no_wash
            parent_mask, parent_vial, parent_used_wash, crossed_wash = table[mask, previous]
            if crossed_wash:
                wash_after = parent_mask.bit_count()
            order.append(parent_vial)
            mask = parent_mask
            previous = parent_vial
            wash_already_used = parent_used_wash
        order.reverse()
        return score, order, wash_after

    candidates = []
    for rho_polar in grid:
        for rho_sticky in grid:
            candidates.append((best_for_fractions(rho_polar, rho_sticky),
                               rho_polar, rho_sticky))
    candidates.sort(key=lambda item: (item[0], item[1], item[2]))
    _, rho_polar, rho_sticky = candidates[0]
    score, order, wash_after = best_for_fractions(rho_polar, rho_sticky, keep_path=True)
    return {
        "vial_order": [vial_ids[index] for index in order],
        "carryover_fraction_polar": rho_polar,
        "carryover_fraction_sticky": rho_sticky,
        "deep_wash_after_acquisition": wash_after,
        "chi_square": score,
    }


sequence_audit = solve_autosampler_sequence()

predictions = []
for run in sorted({row["run_id"] for row in schedule}):
    times = sorted(float(row["time_min"]) for row in schedule if row["run_id"] == run)
    for time, state in zip(times, integrate(fit.x, run, times, exchange_site)):
        for species, concentration in zip(SPECIES, state):
            predictions.append({"run_id": run, "time_min": int(time), "isotopologue": species,
                                "concentration_mmol_L": float(concentration)})
OUTPUT.write_text(json.dumps({"audit": {"swapped_run_pairs": swapped_pairs,
                                         "time_scale_by_run": nonunit_scales,
                                         "exchange_site": exchange_site},
                              "design_run_ids": design_run_ids,
                              "sequence_audit": sequence_audit,
                              "stability": {"most_influential_run_id": most_influential_run_id,
                                            "most_influential_pair": most_influential_pair,
                                            "leave_one_trace_out": stability_rows,
                                            "leave_two_traces_out": pair_stability_rows},
                              "parameters": dict(zip(names, map(float, values))),
                              "predictions": predictions}, indent=2) + "\n")
