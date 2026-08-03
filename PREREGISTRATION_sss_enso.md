# Pre-registered: does the freshwater pool edge lead Nino 3.4 SST?

Written before any salinity data was examined.

## Mechanism

Weakening easterly trades release the western Pacific freshwater pool,
whose eastern edge migrates east. ESA report that this displacement
accompanies the eastward shift of deep convection during El Nino, and
that researchers are developing an ENSO-specific salinity index to test
whether it provides an earlier signal than SST alone.

## Primary hypothesis

H1. The fresh pool edge longitude anomaly correlates positively with the
    Nino 3.4 SST anomaly at a lag of 3 months, salinity leading.
    Criterion: supported if the full-record correlation is positive with
    p < 0.05 using an autocorrelation-adjusted effective sample size, AND
    both split halves agree in sign.

The 3-month lag is specified in advance from the mechanism, not chosen
from a correlogram. The full lag correlogram will be plotted and reported
as EXPLORATORY; its peak is not a test result.

## Descriptive quantities reported before the test

D1. Fraction of months where fresh_pool_edge returns valid.
    If below 0.8, the index is not well defined on this product and the
    test is reported as uninformative.
D2. Fraction of valid months flagged at_boundary. If above 0.1 the
    longitude window is too narrow and must be widened before testing.
D3. Correlation between the edge longitude and Nino 3.4 at zero lag,
    to establish how much of any lagged signal is simple contemporaneity.

## Specification

- Isohaline threshold 34.8 psu, standard in the fresh pool literature.
- Meridional band 5S-5N, averaged to a longitude profile before the edge
  is located.
- Longitude window 130E-270E.
- Monthly means; anomalies by calendar-month climatology.
- Significance via windtools.stats.corr_neff.
- Split-sample replication computed in the same run as the full record.
- No other test will be run. Anything further is exploratory and labelled.

## Rationale

The preceding sea ice analysis ran more than forty tests before
validating, and its most striking finding failed replication. A lag
analysis is unusually prone to the same error because every lag is a
free parameter.
