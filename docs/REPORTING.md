# Task reports and current status

These rules implement the user's accepted reporting format. They do not change
the model, apparatus inputs or acceptance tolerances.

## Task completion report

- At most 30 lines: what changed (2–4 lines), at most five key numbers,
  tests in one line, changed limitations, questions requiring a decision,
  and the next step.
- Always include **Для проверки**, at most ten lines: unmet/qualified criteria
  with observed value and tolerance; decisions without a source, one per line;
  assumption/limitation changes since the previous report.
- When none apply, say so explicitly in one line.
- Put detailed tables, all numbers, source lists, reasoning and path mappings
  in `results/<task>/details.md` and JSON. Do not repeat UNKNOWNS in task reports.
- The reorganization report additionally includes the top-level directory tree.

## STATUS.md

Overwrite the root file after each task, at most40 lines. Include5–7 capabilities,
3–5 verification lines, at mostfive main limitations, open work/next steps,
and2–3 launch commands. Do not append historical progress logs.

Historical stage generators remain reproducibility utilities for archived
reports. They are not the template for subsequent short task reports.
