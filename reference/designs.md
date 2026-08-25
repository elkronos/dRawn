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

Arguments mean the same thing everywhere they appear:

\* \`n\` is always the \*\*total\*\* number of rows drawn, never a
per-group figure. \[design_temporal()\] uses \`per_interval\` precisely
because that one \*is\* per group, and \[design_cluster()\] has no \`n\`
at all because the row count follows from which clusters were selected.
\* \`allocation\` always says how a total is split across groups, either
\`"proportional"\` or \`"equal"\`. Splitting is exact: the
largest-remainder method is used, so a request for 60 rows returns 60
rather than drifting with per-group rounding. \* \`na_rm\` always
decides whether rows with a missing key are dropped or raise an error.
It is never silently assumed either way. \* \`replace\` always means
sampling with replacement within whatever group the design works on. \*
\[draw()\] always restores the caller's random number stream before
returning, and always gives back a data frame with the input's class and
column order.

## Estimating from a sample

Because a design is a value, it can report its own inclusion
probabilities before any sampling happens. \[inclusion_prob()\] gives
first-order probabilities, \[joint_prob()\] second-order ones, and
\[ht_total()\] combines them into a population total with a standard
error.

Not every design has a closed form for these — see \[inclusion_prob()\]
for which, why, and what to do instead. Those that do not say so rather
than returning an approximation.

## Choosing one

- Every unit equally likely:

  \[design_simple()\], or \[design_systematic()\] when the frame has a
  useful order, or \[design_reservoir()\] when it does not fit in
  memory.

- Guaranteed coverage of subgroups:

  \[design_stratified()\], with \`min_per_stratum\` if rare groups must
  appear.

- Fieldwork cost matters more than efficiency:

  \[design_cluster()\] or \[design_multistage()\] — visiting five sites
  is cheaper than visiting fifty, at the cost of precision.

- Large units matter more:

  \[design_weighted()\], with \`method = "systematic"\` or \`"poisson"\`
  if you intend to estimate.

- Coverage across time or space:

  \[design_temporal()\], \[design_spatial()\].

- Uncertainty of a statistic, not a population total:

  \[design_bootstrap()\].

## See also

\[draw()\], \[inclusion_prob()\], \[ht_total()\]
