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
