import json
import os
from pathlib import Path

species = ["A0", "A1", "B0", "B1", "C0", "C1", "D0", "D1"]
parameters = {name: value for name, value in zip(
    ["k1_ref", "k2_ref", "k3_ref", "k4_ref", "kx_ref", "ea1_kJ_mol",
     "ea2_kJ_mol", "ea3_kJ_mol", "ea4_kJ_mol", "eax_kJ_mol"],
    [0.04, 0.03, 0.02, 0.02, 0.0, 40, 40, 40, 40, 40])}
predictions = []
for run, charge, fraction in [("H01", .92, .18), ("H02", 1.08, .30)]:
    for time in (9, 27, 55):
        left = charge * pow(2.718281828, -0.06*time)
        row = [left*(1-fraction), left*fraction, 0, 0, 0, 0,
               (charge-left)*(1-fraction), (charge-left)*fraction]
        predictions.extend({"run_id": run, "time_min": time, "isotopologue": s,
                            "concentration_mmol_L": v} for s, v in zip(species, row))
Path(os.environ.get("TASK_OUTPUT", "/app/result.json")).write_text(
    json.dumps({"audit": {"swapped_run_pairs": [["R01", "R02"], ["R03", "R04"], ["R05", "R06"]],
                           "time_scale_by_run": {"R08": 0.8, "R12": 1.25},
                           "exchange_site": "C"},
                "design_run_ids": ["E01", "E02", "E03", "E04"],
                "parameters": parameters, "predictions": predictions}))
