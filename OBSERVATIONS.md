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

## Tier 4 (multi-file coordination, 32B): partial, not complete

The 32B made the two harder edits correctly — added wind_direction to
speed.py with its import and formula, and updated the diagnostics.py
import to reference it — but left the DIAGNOSTICS dictionary itself
un-updated, an internally inconsistent state: it imported a function it
never registered. Its third edit (the tests) was emitted as an embedded
tool call in the final-answer slot and did not execute, so it never
landed. The suite reported 11 passed, but that green was hollow: it ran
only the pre-existing tests. The test that would have caught the missing
registration was the very edit that failed to run. Completed manually by
registering wind_direction and adding its tests (then 13 passed).

Lesson: coordinated multi-file change sits at the reliable edge even for
the 32B. It gets the mechanically obvious edits but can miss a step, and
a passing suite can conceal the omission when the guarding test is itself
the un-landed edit. Verify content and file count, not the exit code.

Harness note: an embedded tool call in the final-answer slot was dropped
for the third time. The single-quote recovery fix handles calls the
extractor processes, but a call that lands in the final-answer branch
bypasses it. Bounded, known limitation; careful manual verification
catches its consequences.

## Workflow note: environment errors can silently skip verification

During the SAM extension, a stale UV_CACHE_DIR pointing at an unplugged
external drive killed the pytest step in a chained sequence, and the
commit and push proceeded with unverified code. CI caught it (green, as
it happened), and the local run afterwards confirmed 4 passed. Guard:
when a verify command errors for environmental reasons, stop the chain
and re-verify before committing. Backstop: CI. The cache setting is now
fixed permanently in .zshrc.

## ASL extension: emission failure confirmed as reliable, not occasional

The 32B read the ASL contract and produced a correct sector-masked
minimum-finder — but emitted the entire write_file call into the
final-answer slot, unrecoverable because the docstring's unescaped
quotes make the blob unparseable as JSON or a Python literal. Fourth
occurrence of the final-answer emission failure. Completed by hand;
20 passed with the correct changed-file set. Conclusion: single-shot
write-heavy tasks on the 32B need the human fallback ready, and the
harness limitation is bounded but not fixable by better parsing alone.
