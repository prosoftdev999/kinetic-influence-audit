# Campaign deletion analysis

The group needs to know whether the recovered model is being carried by one exported chromatogram. Start from the accepted full-data discrete repair and fitted parameter vector. For each exported training run in sorted ID order, remove that complete observed trace and refit the ten continuous log-parameters. Hold the exchange site, trace-to-condition assignment, and clock scales at their full-data values. Use the same calibrated concentration-space residuals and parameter bounds as the full fit.

For every deletion, report `omitted_run_id`, the Euclidean norm of the change in the ten natural-log parameters as `parameter_log_shift_l2`, the largest absolute concentration change anywhere on the complete withheld prediction schedule as `max_prediction_shift_mmol_L`, and the ten refitted values under `parameters`. Preserve the parameter names from `schema.json`.

Place the twelve rows in `leave_one_trace_out`, sorted by `omitted_run_id`. Set `most_influential_run_id` to the row maximizing `parameter_log_shift_l2 + max_prediction_shift_mmol_L / 0.001`; if scores tie, use the lexicographically larger run ID.

Single deletions can hide interactions between two individually modest traces. Next evaluate every unordered pair of exported training runs. Start each pair refit at the elementwise midpoint of its two single-deletion log-parameter vectors, remove both complete traces, and otherwise use the same fixed discrete repair, residuals, bounds, and withheld prediction grid. Report all 66 rows in `leave_two_traces_out`, ordered lexicographically by `omitted_run_ids`.

Each pair row must contain `omitted_run_ids`, `parameter_log_shift_l2`, `max_prediction_shift_mmol_L`, and the ten refitted `parameters`. Also report `interaction_log_shift_excess` as the pair's log-shift norm minus the two corresponding single-deletion norms. Report `interaction_prediction_shift_excess_mmol_L` in the same way from the maximum prediction displacements. Do not clamp negative interaction excesses. Set `most_influential_pair` to the pair maximizing `parameter_log_shift_l2 + max_prediction_shift_mmol_L / 0.001`; if scores tie, use the lexicographically larger pair.
