# Sampling designs

A design describes \*how\* to sample, independently of \*what\* you are
sampling. Build one with a \`design\_\*()\` constructor, then hand it to
\[draw()\] with the data. The same design can be reused across data
sets, stored, printed and passed around.

## The designs

- \[design_simple()\]:

  Rows uniformly at random.

- \[design_stratified()\]:

  A share of each stratum.

- \[design_systematic()\]:

  Every \*k\*-th row.

- \[design_cluster()\]:

  Whole clusters.

- \[design_multistage()\]:

  Clusters, then rows within them.

- \[design_weighted()\]:

  Rows with probability governed by a weight.

- \[design_reservoir()\]:

  A fixed-size sample from a stream.

- \[design_bootstrap()\]:

  Resampled replicates.

- \[design_temporal()\]:

  A share of each time interval.

- \[design_spatial()\]:

  Rows inside a region.

## What every design guarantees

Arguments mean the same thing everywhere they appear. \`n\` is always
the total number of rows drawn, never a per-group figure; \`allocation\`
always says how a total is split across groups; \`na_rm\` always decides
whether missing keys are dropped or raise an error; and \[draw()\]
always restores the caller's random number stream before returning.

## See also

\[draw()\]
