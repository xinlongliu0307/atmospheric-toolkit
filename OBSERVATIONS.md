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

## ZW3 extension: analytic contract passed at machine precision

Fourier ZW3 index landed with all seven contract tests green at near
machine precision, including wave-2 orthogonality (leakage < 1e-9) and
zonal-mean invariance. [Record here whether the 32B produced the file
or the hand fallback was used, and update the emission-failure count.]
Suite at 27; commit 913e1e1.

## ZW3 extension: analytic contract passed at machine precision

Fourier ZW3 index landed with all seven contract tests green at near
machine precision, including wave-2 orthogonality (leakage < 1e-9) and
zonal-mean invariance. The 32B run hit the final-answer emission failure
again (fifth occurrence); the hand fallback was used, and the registry
was updated manually as planned. Suite at 27; commit 913e1e1.

## ZW3 real-data validation against ERA5 (1979-2025, 564 months, 49S)

Amplitude seasonality: JJA 54.6 m, DJF 40.7 m, July maximum 55.7 m,
December minimum 32.1 m. Winter maximum confirmed, consistent with the
documented winter enhancement of SH stationary wave 3.

Phase: circular mean 50.7 deg on the 120-deg cycle, placing ridges near
51E, 171E, and 291E. These fall within a few degrees of the reference
longitudes used by Raphael (2004), namely 49S at 50E, 166E, and 76W
(doi:10.1029/2004GL020365; locations as stated in Raphael 2007,
doi:10.1029/2006JD007852). Phase clustering R = 0.51.

Significance: the analytic contract could not test the phase sign
convention, because the test and the implementation shared it. Recovering
independently chosen published ridge longitudes provides the external
check that the synthetic oracle structurally could not.

## SAM and ASL from ERA5 MSLP (1979-2025, 564 months)

Climatology reproduces three published features independently:
  - Absolute central pressure shows the semiannual oscillation (minima
    April 975.1 and October 970.6 hPa; maxima January 981.0 and June
    977.4 hPa).
  - Relative central pressure shows instead a winter minimum (August
    -11.7 hPa) and summer maximum (December/January -7.6 hPa), the
    contrast documented by Hosking et al. 2013 (doi:10.1175/JCLI-D-12-00813.1).
  - Longitude migrates from ~224E in JJA to ~242E in DJF, consistent with
    the documented ~220E winter to ~250E summer shift and with Turner
    et al. 2013 (doi:10.1002/joc.3558).

Correlations with deseasonalised anomalies, autocorrelation-adjusted:
  absolute central pressure r = -0.770 (n_eff 477, p 7e-95)
  relative central pressure r = -0.031 (p 0.48)
  ASL longitude             r = +0.090 (p 0.035, fails Bonferroni; treat as null)

Interpretation: the ASL's absolute depth is largely the local expression
of SAM; its depth relative to the sector background is independent of it.
This independently recovers the motivation for the relative-pressure
metric.

CAVEAT PENDING TEST: the ASL sector is 36% of the 65S circle, so the ASL
contributes to the southern node of a zonal-mean SAM index. The -0.770 is
partly circular. Sensitivity test with the sector excluded is in
scripts/test_sam_asl_circularity.py.

## Circularity test and the ZW3-ASL relationship

Sector-excluded SAM (231 of 360 longitudes) correlates 0.871 with the
full-circle index. SAM vs absolute central pressure drops from -0.770 to
-0.585 (34% of variance); quote -0.585. SAM vs relative central pressure
stays null (-0.031 to +0.015), so the independence is robust rather than
an artefact of construction. Seasonal: DJF -0.691, MAM -0.395,
JJA -0.603, SON -0.625.

ZW3 phase vs ASL longitude: r = +0.240 (p 1e-8) but regression slope only
+0.328 +/- 0.056, twelve standard errors below unity. The wave modulates
the position; it does not carry the low. ZW3 amplitude vs relative central
pressure r = -0.177 (3% of variance), correct sign but weak.

Seasonal: wave control on longitude strongest in MAM (+0.335) and DJF
(+0.314), weakest in JJA (+0.131), which is the inverse of the amplitude
seasonality. Autumn is simultaneously the least SAM-controlled and most
wave-controlled season.

CAUTION: the climatological match (ZW3 trough 230.7E, ASL mean 230.8E) is
weaker evidence than it appears, since the sector midpoint is 234E and any
centrally distributed detection would land near there.

## ZW3-ASL: diagnostics resolve the pooled slope into two regimes

Range restriction refuted: phase-anomaly variance is flat across seasons
(sd 19-25 deg), and slopes separate with correlations, so the seasonal
contrast is physical.

Regression dilution confirmed independently of season: terciles on
deseasonalised amplitude give slopes 0.207, 0.398, 0.544 while being
seasonally balanced (strong tercile is if anything winter-enriched, which
works against the gradient).

Within-season terciles are decisive. DJF weak wave 0.311+/-0.179, DJF
strong wave 0.904+/-0.234 (consistent with unity, ~4 sigma from zero).
JJA weak 0.115+/-0.160, strong 0.139+/-0.275 (both null). The pooled 0.33
averages a near-complete summer coupling with an absent winter one.

Boundary contamination negligible: excluding the 10 sector-edge months
moves the slope 0.328 -> 0.351.

HYPOTHESIS REFUTED. Predicted that coupling weakens with meridional
separation between the wave (49S) and the ASL centre (67.8S DJF, 71.6S
JJA). The latitude profile shows the opposite: slope peaks at 49S (0.328)
and decays poleward to 0.057 at 70S, i.e. strongest where they are
furthest apart. Steering is a mid-latitude waveguide phenomenon, not a
local one. Note also that ZW3 phase at 65S does not track ASL longitude
(DJF slope -0.009), so the high-latitude field is not well described by
wave 3.

Alternative interpretation, untested: the ASL as a time-average of
transient cyclones (cf. Fogt et al. 2012) would have its winter position
set by a strong storm track that swamps the stationary wave, with the wave
dominating only in the quieter summer. Testing needs daily fields and an
eddy-activity metric.

## Out-of-sample replication: the strong-wave summer result FAILS

Split 1979-2001 / 2002-2025, anomalies centred within each half.

REPLICATES:
  SAM vs ASL absolute central pressure: slope -1.773+/-0.149 (r -0.584)
  and -1.881+/-0.157 (r -0.579). Exceptional stability. Headline result.
  ZW3 phase vs ASL longitude, all months: +0.308+/-0.084 and
  +0.207+/-0.075. Both significant, overlapping. Modest but real.
  DJF phase vs longitude: +0.412+/-0.205 and +0.377+/-0.177.

DOES NOT REPLICATE:
  DJF strong-wave slope. Full record gave +0.904+/-0.234; halves give
  +0.198+/-0.520 and +0.469+/-0.331, both null. JJA strong-wave sign flips
  (-0.289, +0.624). At n~23 nothing is estimable. The amplitude-dilution
  gradient is now also in doubt, sharing the same construction.

Likely cause: full-record calendar-month centring retains shared
low-frequency covariance between phase and longitude that within-half
centring removes. Tested in scripts/test_lowfreq_inflation.py.

METHODOLOGICAL LESSON: >40 tests were run on 564 months across four
scripts. The one finding selected as most striking failed replication.
Split-sample validation should be run before, not after, interpretation.
