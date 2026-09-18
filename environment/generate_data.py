from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
DATA.mkdir(exist_ok=True)

RUNS = [
    ("R01", "S1", 318.15, 0.80, 0.25, 0.05, 1.00, True),
    ("R02", "S2", 318.15, 1.15, 0.25, 0.18, 1.00, True),
    ("R03", "S3", 323.15, 0.65, 0.25, 0.32, 1.00, True),
    ("R04", "S1", 323.15, 1.30, 0.25, 0.47, 1.00, True),
    ("R05", "S2", 328.15, 0.75, 0.25, 0.62, 1.00, True),
    ("R06", "S3", 328.15, 1.20, 0.25, 0.78, 1.00, True),
    ("R07", "S1", 333.15, 0.60, 0.25, 0.12, 1.00, True),
    ("R08", "S2", 333.15, 1.05, 0.25, 0.27, 1.00, True),
    ("R09", "S3", 338.15, 0.70, 0.25, 0.43, 1.00, True),
    ("R10", "S1", 338.15, 1.25, 0.25, 0.58, 1.00, True),
    ("R11", "S2", 343.15, 0.55, 0.25, 0.73, 1.00, True),
    ("R12", "S3", 343.15, 0.95, 0.25, 0.88, 1.00, True),
    ("H01", "S2", 323.15, 0.95, 0.18, 0.62, 0.92, False),
    ("H02", "S3", 333.15, 0.72, 0.30, 0.08, 1.08, False),
    ("H03", "S1", 343.15, 0.52, 0.47, 0.81, 1.04, False),
]
TIMES = np.array([0, 3, 7, 12, 20, 32, 48, 70], dtype=float)
SPECIES = ["A0", "A1", "B0", "B1", "C0", "C1", "D0", "D1"]
RF = {
    "S1": {"A": 1.00, "B": 0.77, "C": 1.22, "D": 0.91},
    "S2": {"A": 1.08, "B": 0.72, "C": 1.29, "D": 0.86},
    "S3": {"A": 0.94, "B": 0.83, "C": 1.15, "D": 0.97},
}
PARAM = {"k1_ref": 0.047, "k2_ref": 0.031, "k3_ref": 0.021,
         "k4_ref": 0.026, "kx_ref": 0.018, "ea1": 43.0, "ea2": 57.0,
         "ea3": 38.0, "ea4": 51.0, "eax": 29.0}


def rates(temp: float, catalyst: float):
    gas = 0.008314462618
    result = []
    for name, ea in [("k1_ref", "ea1"), ("k2_ref", "ea2"), ("k3_ref", "ea3"),
                     ("k4_ref", "ea4"), ("kx_ref", "eax")]:
        result.append(PARAM[name] * catalyst * np.exp(-PARAM[ea] / gas * (1 / temp - 1 / 328.15)))
    return result


def simulate(temp: float, catalyst: float, initial_fraction: float, solvent_fraction: float,
             charge: float, time_scale: float = 1.0):
    k1, k2, k3, k4, kx = rates(temp, catalyst)

    def rhs(_t, y):
        a0, a1, b0, b1, c0, c1, d0, d1 = y
        exchange = kx * (b0 * solvent_fraction - b1 * (1 - solvent_fraction))
        return [-(k1 + k2) * a0, -(k1 + k2) * a1,
                k1*a0 - k3*b0 - exchange, k1*a1 - k3*b1 + exchange,
                k2*a0 - k4*c0, k2*a1 - k4*c1,
                k3*b0 + k4*c0, k3*b1 + k4*c1]

    y0 = [charge*(1-initial_fraction), charge*initial_fraction, 0, 0, 0, 0, 0, 0]
    sample_times = TIMES * time_scale
    return solve_ivp(rhs, (0, sample_times[-1]), y0, t_eval=sample_times,
                     rtol=1e-11, atol=1e-13).y.T


with (DATA / "calibration.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["session", "compound", "known_mmol_L", "peak_area"])
    for session in RF:
        for compound in "ABCD":
            for concentration in (0.20, 0.55, 1.00):
                writer.writerow([session, compound, concentration,
                                 f"{RF[session][compound] * concentration * 100000:.3f}"])

with (DATA / "run_manifest.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["run_id", "session", "temperature_K", "catalyst_mmol_L",
                     "initial_A_label_fraction", "solvent_label_fraction", "initial_A_mmol_L", "role"])
    for run, session, temp, catalyst, initial_frac, solvent_frac, charge, training in RUNS:
        writer.writerow([run, session, temp, catalyst, initial_frac, solvent_frac, charge,
                         "training" if training else "prediction"])

rng = np.random.default_rng(70419)
with (DATA / "chromatograms.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["run_id", "time_min", "isotopologue", "peak_area"])
    exported_id = {"R02": "R09", "R09": "R02", "R04": "R11", "R11": "R04",
                   "R06": "R07", "R07": "R06"}
    for run, session, temp, catalyst, initial_frac, solvent_frac, charge, training in RUNS:
        if not training:
            continue
        clock_scale = {"R03": 0.80, "R10": 1.25}.get(run, 1.0)
        values = simulate(temp, catalyst, initial_frac, solvent_frac, charge, clock_scale)
        for time, row in zip(TIMES, values):
            for species, concentration in zip(SPECIES, row):
                compound = species[0]
                noise = rng.normal(0, 45.0) if time else 0.0
                area = max(0.0, RF[session][compound] * concentration * 100000 + noise)
                writer.writerow([exported_id.get(run, run), int(time), species, f"{area:.3f}"])

with (DATA / "prediction_times.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["run_id", "time_min"])
    for run, *_rest, training in RUNS:
        if not training:
            for time in (9, 27, 55):
                writer.writerow([run, time])

with (DATA / "design_candidates.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["candidate_id", "temperature_K", "catalyst_mmol_L",
                     "initial_A_label_fraction", "solvent_label_fraction", "initial_A_mmol_L"])
    candidate_number = 1
    for temperature in (318.15, 328.15, 338.15, 348.15):
        for catalyst in (0.50, 1.00, 1.40):
            for solvent_fraction in (0.15, 0.75):
                writer.writerow([f"E{candidate_number:02d}", temperature, catalyst, 0.25,
                                 solvent_fraction, 1.00])
                candidate_number += 1

with (DATA / "design_scenarios.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    names = ["k1_ref", "k2_ref", "k3_ref", "k4_ref", "kx_ref", "ea1", "ea2", "ea3", "ea4", "eax"]
    writer.writerow(["scenario_id", *names])
    writer.writerow(["central", *([1.0] * 10)])
    writer.writerow(["fast_parallel", 1.06, 1.05, 0.96, 0.97, 1.03, 0.98, 1.02, 1.01, 0.99, 1.03])
    writer.writerow(["slow_parallel", 0.95, 0.94, 1.05, 1.04, 0.97, 1.02, 0.98, 0.99, 1.01, 0.97])
    writer.writerow(["hot_sensitive", 1.02, 0.98, 1.03, 0.97, 1.04, 1.04, 1.05, 0.96, 0.95, 1.03])
    writer.writerow(["exchange_weak", 1.01, 1.00, 0.98, 1.02, 0.91, 0.99, 1.01, 1.02, 0.98, 0.94])

(DATA / "lab_record.md").write_text("""# Reactor campaign record

The campaign follows A through two observed intermediates, B and C, to a common product D. Each compound is integrated as its M+0 and M+1 isotopologue. Forward chemical steps preserve the label. The notebook does not settle whether reversible solvent exchange occurs in B or in C, but it occurs in exactly one of them and nowhere else. Runs use the same catalyst lot and the kinetic constants are shared across the campaign. Catalyst loading scales every elementary rate linearly.

The five elementary rates follow Arrhenius temperature dependence relative to 328.15 K. Report activation energies in kJ/mol. Concentrations in the kinetic model are mmol/L and time is minutes. The calibration injections have zero intercept; M+0 and M+1 share their parent compound's response factor within a session.

For exchange, the solvent-pool labelled fraction is listed separately in the manifest. The forward exchange flux is proportional to the unlabelled exchanging intermediate times that fraction; the reverse flux is proportional to its labelled isotopologue times the complementary fraction, with the same exchange constant. The reactor is closed, so A+B+C+D remains equal to the initial A charge within measurement noise.

The sample log records one export incident: three disjoint pairs of complete training chromatograms received each other's run IDs. Exactly six training IDs are affected. Within a chromatogram the time points and isotopologue channels stayed together, and the manifest itself is correct. All training reactors began with the same A concentration and isotope fraction, so the time-zero rows do not identify the pairs. The prediction runs were entered later and are not involved.

The pump-controller audit found two additional clock faults among the exported training chromatograms. For one trace, actual elapsed time equals recorded time times 0.80; for one other trace it equals recorded time times 1.25. Every remaining trace has scale 1.00. The two affected exported IDs are not among the six IDs involved in the transpositions.

After reconstruction, choose four follow-up runs from `design_candidates.csv`. For each uncertainty row in `design_scenarios.csv`, multiply the fitted parameters by that row's factors and compute central-difference sensitivities of all eight concentrations at 6, 18, and 42 minutes with respect to the natural logarithm of each of the ten parameters; use a log-step of 0.0001. For candidate c and scenario s, let J be the resulting 24 by 10 sensitivity matrix and set I(c,s) = J'J / 0.002^2. Score a four-run set by the minimum, over the five scenarios, of log(det(Identity + sum I(c,s))). The selected set must have total catalyst loading no greater than 3.80 mmol/L, span at least 15 K, and include both available solvent label fractions. Choose the feasible set with the largest score; if scores agree within 1e-10, choose the lexicographically smaller sorted candidate-ID list.
""")

(DATA / "stability_protocol.md").write_text("""# Campaign deletion analysis

The group needs to know whether the recovered model is being carried by one exported chromatogram. Start from the accepted full-data discrete repair and fitted parameter vector. For each exported training run in sorted ID order, remove that complete observed trace and refit the ten continuous log-parameters. Hold the exchange site, trace-to-condition assignment, and clock scales at their full-data values. Use the same calibrated concentration-space residuals and parameter bounds as the full fit.

For every deletion, report `omitted_run_id`, the Euclidean norm of the change in the ten natural-log parameters as `parameter_log_shift_l2`, the largest absolute concentration change anywhere on the complete withheld prediction schedule as `max_prediction_shift_mmol_L`, and the ten refitted values under `parameters`. Preserve the parameter names from `schema.json`.

Place the twelve rows in `leave_one_trace_out`, sorted by `omitted_run_id`. Set `most_influential_run_id` to the row maximizing `parameter_log_shift_l2 + max_prediction_shift_mmol_L / 0.001`; if scores tie, use the lexicographically larger run ID.

Single deletions can hide interactions between two individually modest traces. Next evaluate every unordered pair of exported training runs. Start each pair refit at the elementwise midpoint of its two single-deletion log-parameter vectors, remove both complete traces, and otherwise use the same fixed discrete repair, residuals, bounds, and withheld prediction grid. Report all 66 rows in `leave_two_traces_out`, ordered lexicographically by `omitted_run_ids`.

Each pair row must contain `omitted_run_ids`, `parameter_log_shift_l2`, `max_prediction_shift_mmol_L`, and the ten refitted `parameters`. Also report `interaction_log_shift_excess` as the pair's log-shift norm minus the two corresponding single-deletion norms. Report `interaction_prediction_shift_excess_mmol_L` in the same way from the maximum prediction displacements. Do not clamp negative interaction excesses. Set `most_influential_pair` to the pair maximizing `parameter_log_shift_l2 + max_prediction_shift_mmol_L / 0.001`; if scores tie, use the lexicographically larger pair.
""")

(DATA / "schema.json").write_text(json.dumps({
    "output": "/app/result.json",
    "required": {
        "audit": ["swapped_run_pairs", "time_scale_by_run", "exchange_site"],
        "design_run_ids": "four sorted IDs from design_candidates.csv",
        "stability": "leave-one-trace-out analysis from stability_protocol.md",
        "parameters": ["k1_ref", "k2_ref", "k3_ref", "k4_ref", "kx_ref",
                       "ea1_kJ_mol", "ea2_kJ_mol", "ea3_kJ_mol", "ea4_kJ_mol", "eax_kJ_mol"],
        "predictions": "one row per run_id/time_min/isotopologue in prediction_times.csv",
        "sequence_audit": {
            "vial_order": "all ten QC vial IDs in acquisition order",
            "carryover_fraction_polar": "one value from the qualification grid",
            "carryover_fraction_sticky": "one value from the qualification grid",
            "deep_wash_after_acquisition": "1-based acquisition number before the wash break",
            "chi_square": "minimum normalized residual sum from autosampler_protocol.md"
        }
    }
}, indent=2) + "\n")


# Autosampler qualification.  This is generated after the reactor files so the
# original campaign RNG stream and every pre-existing data file remain unchanged.
qc_rng = np.random.default_rng(45177)
qc_channels = ["polar_1", "polar_2", "polar_3", "sticky_1", "sticky_2", "sticky_3"]
qc_base = np.array([0.42, 0.31, 0.24, 0.18, 0.13, 0.09])
qc_nominal = []
for index in range(10):
    response = (qc_base.copy() + qc_rng.normal(0, 0.055, size=6)
                + 0.035 * np.sin((index + 1) * np.arange(1, 7) * 1.7))
    qc_nominal.append(np.clip(response, 0.025, None))
qc_nominal = np.asarray(qc_nominal)
qc_order = [7, 1, 9, 3, 0, 6, 2, 8, 4, 5]
qc_rho_polar, qc_rho_sticky = 0.36, 0.24
qc_wash_after = 6
qc_sigma = 0.003
qc_measured = []
for position, vial in enumerate(qc_order):
    response = qc_nominal[vial].copy()
    if position > 0 and position != qc_wash_after:
        carry = qc_nominal[qc_order[position - 1]].copy()
        carry[:3] *= qc_rho_polar
        carry[3:] *= qc_rho_sticky
        response += carry
    response += qc_rng.normal(0, qc_sigma, size=6)
    qc_measured.append(response)

with (DATA / "autosampler_nominals.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["vial_id", *qc_channels])
    for index, response in enumerate(qc_nominal, 1):
        writer.writerow([f"Q{index:02d}", *[f"{value:.9f}" for value in response]])

with (DATA / "autosampler_sequence.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["acquisition_index", *qc_channels])
    for index, response in enumerate(qc_measured, 1):
        writer.writerow([index, *[f"{value:.9f}" for value in response]])

(DATA / "autosampler_protocol.md").write_text("""# Autosampler carryover qualification

The autosampler qualification immediately preceded the reactor campaign. Its sequence export kept the acquisition positions but lost the vial IDs. `autosampler_nominals.csv` contains the six-channel normalized response expected from each qualification vial when injected after a clean needle. `autosampler_sequence.csv` contains the measured six-channel responses in acquisition order. Every listed vial was injected exactly once.

The qualification model includes one-step needle memory. Channels `polar_1` through `polar_3` share one carryover fraction, and channels `sticky_1` through `sticky_3` share a second fraction. Each fraction is one of the five values 0.18, 0.24, 0.30, 0.36, or 0.42. Except across a deep wash, the predicted response at acquisition i is the current vial's nominal response plus the previous vial's nominal response multiplied by the appropriate carryover fraction. The first acquisition has no previous-vial contribution.

Exactly one programmed deep wash occurred between two adjacent acquisitions. For the acquisition immediately after that wash, the previous-vial contribution is zero. That acquisition then becomes the previous vial for the following acquisition, so ordinary carryover resumes. The independent standard deviation of each normalized channel is 0.003.

Recover the vial order, the two carryover fractions, and the wash location by minimizing

`chi_square = sum(((measured - predicted) / 0.003)^2)`

over all vial permutations, both carryover-grid choices, and every possible wash boundary. Use each vial exactly once. Report the 1-based acquisition number after which the wash occurred. If hypotheses agree within 1e-9 in chi-square, choose the lexicographically smaller `vial_order`; if that is also tied, choose the smaller wash index, then the smaller polar carryover, then the smaller sticky carryover.
""")
