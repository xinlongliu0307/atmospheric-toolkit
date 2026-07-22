# Local Coding Agent: Capability Observations

A record of what the local agent (Ollama-served qwen2.5-coder, 7B and 32B)
can and cannot reliably do, measured across a graduated series of tasks.
The purpose is to guide delegation: which real tasks to hand to which model,
and which to keep for manual work.

## Model-selection rule

- **7B**: scaffolding new standalone code, and diagnosing bare, isolated
  bugs (a single wrong operator in a short function). Fast and reliable
  on its strong ground.
- **32B**: implementing a formula against a precise test contract, and
  most modification work. Reliable where the 7B is not, for
  implementation and coordinated edits.
- **Manual**: diagnosing a bug embedded in a real formula. Unreliable on
  both models (see Tier 3). A human fix is faster and safer.

## Toy levels (metcalc / sandbox)

- **Level 1 (scaffold, 7B)**: correct code, but self-written tests had
  wrong reference values. Lesson: verify test *content*, not just green.
- **Level 2 (modify, 7B then 32B)**: coordinated code+test change defeated
  both single-shot; 7B added a signature without logic, 32B wrote tests
  without the code. Surfaced a harness bug (single-quoted tool calls not
  recovered), since fixed.
- **Level 3 (implement to spec, 7B fail, 32B pass)**: 7B produced a wrong
  algorithm with dataset-artifact markers; 32B implemented cleanly.
  Implement-to-spec is past the 7B edge.
- **Level 4 (diagnostic, 7B pass)**: bare +/- sign in a two-line function;
  7B read it correctly and fixed it, running its own check.

## Atmospheric tiers (atmospheric-toolkit)

- **Tier 1 (scaffold wind_speed, 7B)**: clean implementation; needed a
  pytest pythonpath fix for the package to import. Green on CI.
- **Tier 2 (potential temperature vs MetPy, 32B)**: implemented the
  Poisson relation correctly; verified against MetPy as oracle within
  rel=1e-3 across 300-1000 hPa. The most valuable pattern: test a
  hand-rolled formula against an authoritative library. Green on CI.
- **Tier 3 (dewpoint sign bug)**: BOTH models failed the diagnosis.
  7B froze re-reading the file and was halted by the progress detector.
  32B, denied a MetPy oracle by the system-shell environment, abandoned
  diagnosis and thrashed through three whole-formula rewrites, each
  worse. Progress detector did not fire because each result differed.
  Fixed manually with a one-character sign change. Green on CI.

## Cross-cutting lessons

- **Passing is not correct.** A test with a wrong reference passes
  falsely. Anchor to an authoritative source; use a library as oracle
  where one exists.
- **Agreement is bounded by constants.** Physical implementations differ
  at the 4th-5th significant figure; use relative tolerances, not
  machine precision.
- **The agent runs in the system shell.** It cannot import MetPy or
  pytest; those live only in the project .venv. Do not ask it to run
  tests that need them; verify manually.
- **Safety machinery held throughout.** Confinement, per-action approval,
  and snapshots protected every run. The progress detector caught the
  7B dead-end; careful mode let bad edits be declined.
- **Close the loop with a push.** Work is only tracked when committed
  AND pushed; CI is the backstop that proves tests pass on a clean
  machine.
