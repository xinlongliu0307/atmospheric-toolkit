# What a Withdrawn Result Taught Me

A worked case study, with a public commit history, in which a
statistically significant, out-of-sample-replicated finding turned out to
be an artefact of the instrument that measured it.

None of the statistical lessons below are novel. They are standard, and
that is the point: they were all known to me in the abstract and I
violated each one in practice. The value here is the audit trail.

## The finding, and its withdrawal

Pre-registered hypothesis: eastward displacement of the Amundsen Sea Low
shifts the node of a sea ice concentration dipole eastward. Measured
slope +0.321 +/- 0.063 on 564 months of ERA5, replicated at +0.309 and
+0.354 in two independent halves, all below a Bonferroni threshold.

Withdrawn. Two alternative node estimators gave approximately zero, the
three estimators intercorrelated at only +0.13, -0.09 and +0.15, and the
original argmax estimator placed 43.4% of nodes exactly on a search
boundary. A purpose-built replacement, tested against a null matched to
the field's own longitudinal autocorrelation, passed real data at 55.0%
and matched red noise at 53.1%. The construct does not exist in the data.

## Lesson 1: replication does not protect against a biased instrument

Split-sample replication tests for sampling noise. A systematic artefact
replicates perfectly, because it is present in every sample. The argmax
estimator's boundary saturation was as stable across halves as any real
signal would have been. Robustness to methodological choices must be
tested separately from, and in addition to, out-of-sample replication.

## Lesson 2: white noise is the wrong null for a smooth field

The estimator rejected 3000 white-noise draws absolutely, with a maximum
r-squared of 0.120 against a 0.5 threshold. Against AR(1) profiles
matched to the data's own lag-1 autocorrelation of 0.985, it passed 53.1%
of draws. The window spanned one fitted period while the field
decorrelated over half of it, leaving roughly two independent points, so
almost any smooth profile was a wavenumber-1 pattern. The null must be
matched to the autocorrelation of the field being tested.

## Lesson 3: test a statistic on degenerate cases before using it

An antisymmetry statistic was introduced to separate dipoles from
monopoles. It anchored on the fitted node rather than the structure's own
centre, so monopoles scored +0.63 to +1.00, higher than ideal dipoles.
It was invalid on exactly the case it was built to exclude, and this was
discovered only because the reference cases were computed alongside the
observations rather than assumed.

## Lesson 4: state the degenerate case before running the test

Several proposed follow-ups were abandoned after arithmetic rather than
after data. A chain test was found to have an expected correlation of
0.046 against a detectable threshold of 0.089, and was not run. A
sensitivity test using a geostrophically related variable was framed as
partially rather than fully independent, because the two predictors are
related by a constant phase offset. Computing the expected effect size
before collecting evidence prevents uninterpretable nulls.

## What survived

Two findings, both robust to the checks above.

SAM controls Amundsen Sea Low absolute central pressure at -1.8 hPa per
standard deviation, r = -0.58, replicating at -1.773 +/- 0.149 and
-1.881 +/- 0.157 across independent halves. It depends on no estimator
choice.

ZW3 phase modestly steers ASL longitude, slope near 0.30, stable across
three filtering treatments and replicating at +0.308 and +0.207.

Also two validations rather than claims: Fourier-derived ZW3 ridge
longitudes recover reference points chosen independently two decades
earlier, and the ASL seasonal cycle reproduces the published contrast
between a semiannual oscillation in absolute central pressure and a
winter minimum in relative central pressure.

## Scoreboard

Findings withdrawn: 1. Mechanistic hypotheses refuted: 2.
Pre-registered predictions unsupported: 2. Predictions by the analyst
that proved incomplete or wrong: 6. Findings surviving all checks: 2.

Every item above is recorded in OBSERVATIONS.md and
PREREGISTRATION_seaice.md, timestamped in the commit history, including
the pre-registrations written before the analyses that tested them.

## Lesson 5: the method working prospectively

The four lessons above are post-mortems. This one is not.

A sensitivity test was planned using a meridional-wind formulation of the
same circulation index, on the reasoning that a different variable and a
different physical quantity would provide independent corroboration of a
surviving finding. Before running it, the registration required two things
to be stated in advance: an analytic prediction for the phase relationship
between the two indices, and a degeneracy diagnostic to be reported before
the primary test.

The analytic prediction was confirmed to a tenth of a degree, which
validated both implementations. The degeneracy diagnostic came out at
0.998. The two predictors were the same quantity plus a constant offset,
and the primary test was therefore uninformative by construction.

Its result: slope +0.320 +/- 0.055, p = 1.25e-08, replicating in both
independent halves, from a different variable. Reported without the
degeneracy check, that is a compelling independent confirmation. It is
nothing of the kind. The registration caught it before the data was
examined rather than after a reviewer asked.

The general form: when a test uses a predictor derived from or physically
related to the original one, compute their correlation and state the
threshold at which the test becomes uninformative BEFORE running it.
