# Pre-registered hypotheses: ASL, ZW3, and Antarctic sea ice

Written before examining the sea ice data. Amendments must be dated and
justified below, not silently edited.

## Predictions

H1. A deeper ASL relative central pressure produces a meridional dipole
    in sea ice concentration: negative anomalies in the Bellingshausen
    sector (east flank, warm northerly advection) and positive anomalies
    in the Ross sector (west flank, cold southerly advection).

H2. Eastward displacement of the ASL longitude shifts the node of that
    dipole eastward by a comparable amount.

H3. SAM's sea ice signature is circumpolar rather than sector-confined,
    and is therefore distinguishable from H1 in spatial pattern.

## Analysis specified in advance

- Monthly ERA5 sea ice cover, 1979-2025, deseasonalised by calendar month.
- Sector means: Ross 160E-230E, Amundsen 230E-260E,
  Bellingshausen 260E-290E, all 60-75S.
- Single primary test per hypothesis. Significance with effective sample
  size adjusted for lag-1 autocorrelation.
- Bonferroni threshold p < 0.017 for the three primary tests.
- Split-sample replication (1979-2001, 2002-2025) run at the same time as
  the full-record analysis, not afterwards.
- Any test beyond the three above is exploratory and labelled as such.

## Rationale

The preceding ZW3-ASL analysis ran >40 tests before validating, and its
most striking finding failed replication. This registration exists to
prevent a repeat.

## Outcomes (recorded 2026-07-26, after analysis)

H1  NOT SUPPORTED. r = +0.043 (full), +0.039, +0.045 (halves); predicted
    negative. Stable across samples, so a genuine null rather than an
    underpowered test. Limitation: relative central pressure may be the
    wrong predictor for an advective mechanism, which depends on the
    pressure gradient rather than depth relative to a sector mean.

H2  SUPPORTED, with replication. Slope +0.321+/-0.063 (full),
    +0.294+/-0.095 and +0.360+/-0.087 (halves), all p < 0.017. A ten
    degree eastward ASL displacement moves the ice dipole node about
    three degrees east. Slope well below unity indicates the ice pattern
    is anchored by factors beyond the low's position.

H3  NOT SUPPORTED, and mis-specified. Full record r = +0.094 (p 0.067);
    halves +0.296 (significant) and -0.019 (null), opposite signs
    cancelling. The descriptive profiles show the prediction was wrong in
    kind: SAM's ice signature is zonally structured (mean +0.035,
    sd 0.113, sign reversals by sector) while the ASL relative-pressure
    signature is weak and uniform (mean +0.042, sd 0.037), the inverse of
    what was predicted. Pattern correlation between the two profiles is
    -0.048, confirming the indices act differently on the ice.
    The between-half non-stationarity is recorded without explanation; a
    mechanism would require its own pre-registration.

Note: one dead-code block in node_longitude was removed after the run;
H2 values verified unchanged.

## Post-hoc robustness check on H2 (recorded 2026-07-26)

H2 was supported under the pre-registered node estimator (argmax of a
west-minus-east contrast): slope +0.32, replicated in both halves. An
alternative sign-weighted-centroid estimator gives approximately zero
(-0.003 full record; -0.035 and +0.029 in the halves, opposite signs).
The finding is therefore NOT robust to the estimator choice, and the
pre-registered slope should not be quoted without this caveat.

Suspected cause: the argmax estimator has a node-anomaly sd of 44.2 deg
within a 110 deg candidate window, exceeding the 31.8 deg a uniform
distribution would give, which indicates boundary saturation and
near-binary behaviour rather than position estimation. Diagnosed in
scripts/diagnose_node_estimator.py.

METHODOLOGICAL LESSON, more important than the physics: the argmax result
replicated cleanly across two independent halves (+0.309, +0.354). Split-
sample replication protects against sampling noise but not against a
systematically biased instrument, because a stable artefact replicates
perfectly. Robustness to methodological choices must be tested separately
from, and in addition to, out-of-sample replication.

The ZW3 chain test is UNDERPOWERED under both estimators (expected r
approximately 0.046 against a detectable threshold of 0.089 at n_eff
~490) and will not be pre-registered. This dataset cannot resolve the
chain at monthly resolution.

## H2 WITHDRAWN (recorded 2026-07-26)

Diagnostic verdict: the pre-registered argmax node estimator is
boundary-saturated. 43.4% of months fall exactly on a candidate boundary
(23.8% at 170E, 19.7% at 280E); the decile distribution is U-shaped; the
node sd of 44.7 deg exceeds the 31.8 deg of a uniform distribution over
the window. A median-split binary flag reproduces r = +0.182 against the
estimator's own +0.211, so the apparent signal is a west/east dichotomy,
not a position relationship.

Two estimators without those pathologies give null results:
  B sign-weighted centroid: -0.003 +/- 0.013 (halves -0.035, +0.029)
  C fitted-sinusoid phase:  -0.007 +/- 0.055 (halves +0.043, -0.035)

The three estimators intercorrelate at only +0.131, -0.093, and +0.152.
Three methods targeting one physical quantity would agree; these do not,
indicating the target is not well defined. The sector sea ice anomaly
profile frequently lacks a single locatable dipole node, so the
pre-registered question presupposed a structure that is often absent.
That presupposition was never itself tested.

CONCLUSION: no reliable evidence that ASL longitude displaces a sea ice
dipole node. The +0.32 slope is withdrawn and must not be cited.

A well-posed reformulation, for a future pre-registration only: regress
the SIC anomaly at each longitude on the ASL longitude anomaly and assess
the resulting spatial pattern. This requires no node to exist. It was not
run here, to avoid post-hoc testing on the same data.

## H2 construct measured non-existent (recorded 2026-07-26)

A validity-aware estimator (windtools/dipole.py, 11 analytic tests) was
built to replace the boundary-saturated argmax. It fixed the pathology
(9.4% of nodes near a window edge, against 43.4%) and rejects white noise
absolutely (3000 draws, max r2 0.120, none above the 0.5 threshold).

It does not, however, discriminate structure. Against the correct null,
AR(1) profiles matched to the data's own longitudinal autocorrelation of
0.985, red noise passes at 53.1% with median r2 0.521; the observations
pass at 55.0% with median r2 0.547. Indistinguishable.

Cause: the 130 deg window spans one fitted period while the field
decorrelates over ~66 deg, leaving ~2 independent points across the
window. Any smooth profile is a wavenumber-1 pattern at that resolution.

An antisymmetry discriminator was attempted and is INVALID: it anchors on
the fitted node rather than the structure's own centre, so monopoles
score +0.63 to +1.00, higher than ideal dipoles. It measures nothing and
its observed distribution should be disregarded.

CONCLUSION: the sector sea ice anomaly does not contain a locatable
dipole node at monthly resolution. The withdrawn H2 asked about a
structure that is not present in this field. This closes the question by
measurement rather than by a failed correlation, and explains why three
estimators disagreed: there was nothing for them to agree about.

The one real signal in the descriptive pass is seasonal: 28% valid in
DJF against 70-75% in JJA-SON, consistent with the ice edge lying inside
the 60-75S band in winter. This reflects anomaly magnitude and smoothness,
not dipole structure.
