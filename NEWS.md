# drawn (development version)

## Wrong answers, fixed

* **Reordering, subsetting or binding a weighted sample gave wrong standard
  errors.** The sample matched its rows to the frame by position, so a sort
  moved each row onto another row's stratum, cluster or systematic position.
  Sorting one stratified sample by `y` cut its standard error from 3,015 to
  1,332, and `head()`, `s[cond, ]` and `rbind()` were estimated as if drawn.
  Each row now finds its frame row through the new `.row` column, so
  reordering by `[`, `dplyr::arrange()` or `data.table::setorder()` changes
  nothing, and `ht_total()`, `ht_mean()`, `as_svydesign()` and
  `sample_summary()` refuse a subset, a sample bound to another, or a join
  that repeats rows.
* **Editing the frame in place after drawing changed the sample's
  estimates.** The sample kept a reference to the frame, so a data.table
  `setkey()`, `setorder()` or `:=` on it afterwards changed the standard
  error, interval, degrees of freedom and design effect of a sample already
  drawn. In one example `setkey()` halved a cluster sample's standard error;
  in another `setorder()` took a certainty sample's design effect from 1.04
  to 171. `draw()` now stores a copy of the frame, column names included, so
  `setnames()`, `setcolorder()` or deleting a column with `:=` does not
  reach the sample either.
* **A certainty design over a with-replacement `rest` was weighted anyway.**
  `draw(weights = TRUE)` checked `replace` only on the outer design, so
  `design_certainty()` over `design_simple(replace = TRUE)` or
  `design_stratified(replace = TRUE)` returned duplicates weighted as distinct
  units, and in one example a total a third too high. It is now refused like
  any other with-replacement design. Drawing it without weights is still
  allowed.
* **Renaming, recoding or dropping a design column in a sample changed its
  standard errors.** The degrees of freedom, the jackknife, `as_svydesign()`
  and `sample_summary()` read strata and cluster labels from the sample's
  columns, while the analytic variance read them from the frame. Dropping a
  cluster sample's cluster column set the degrees of freedom to 0, which
  removed the interval, and made the jackknife delete rows instead of
  clusters, 61% low in one example. Recoding stratum labels broke the
  jackknife and `as_svydesign()`, or gave them the wrong strata. Every path
  now reads design columns from the copy of the frame stored with the sample,
  so only `y` and `by` come from the sample's own columns, and
  `as_svydesign()` passes the clusters to `survey` in a `.cluster` column.
* **A spatial sample analysed under another `sf::sf_use_s2()` setting got
  the wrong standard error.** Which rows lie inside the region was worked out
  again in the current session. A sample drawn under `sf::sf_use_s2(FALSE)`,
  the documented workaround for a region that crosses the antimeridian, and
  analysed after switching back or after `readRDS()` in a new session got
  standard errors for a different region: 5.7 times too large in one example,
  and 0 with a zero-width interval in another. The point estimate was not
  affected. `ht_total()`, `ht_mean()`, `as_svydesign()` and
  `sample_summary()` now check that the design still gives the frame the
  inclusion probabilities it gave at draw time, and refuse the sample when it
  does not, naming `sf::sf_use_s2()` for a spatial design. `?design_spatial`
  says the setting must be the same at draw time and after.
* **The jackknife and the local-mean variance treated rows taken with
  certainty as random.** A row that is in every possible sample was deleted,
  or compared with its neighbours, as if it could have been left out. With
  the units `design_weighted(method = "systematic")` takes with certainty,
  the jackknife averaged 181 times the true variance on one frame, and
  `variance = "auto"` fell back to it when one sampled row was below
  certainty. `design_spread(size =)` with dominant sizes averaged 1,270 times
  the true variance. Both estimators now hold certain units fixed, and the
  jackknife's finite population correction counts only the other units. With
  one unit below certainty the variance is reported as missing instead.
* **A census got a nonzero jackknife or local-mean standard error.**
  `variance = "jackknife"` on a sample of every row or every cluster applied
  no finite population correction, giving 114 for a six-cluster census in one
  example, and `design_spread()` with `n` equal to the frame reported a
  standard error of about 5% of the exact total. Both now give 0, as the
  analytic estimators did.
* **`design_weighted()` could give a certain unit a probability a rounding
  error below 1.** Depending on the scale of the sizes, a unit whose share came to
  exactly 1 got `1 - 1e-16`, so `.prob == 1` missed it and a census reported
  a small standard error under a "Deville approximation" note. A share within
  `1e-9` of 1 is now taken with certainty at exactly 1, in
  `design_weighted()` and `design_spread(size =)`.
* **A certainty design whose `rest` drew no rows reported an exact total.**
  When a `rest` of random size, such as Poisson sampling with a small `n`,
  drew none of the rows below the threshold, `ht_total()` and `ht_mean()`
  gave a standard error of 0 and a zero-width interval although those rows
  went unobserved. The total is exact only when no row lies below the
  threshold. Otherwise the variance is now missing, with a note naming the
  rows `rest` did not reach.
* **Units taken with certainty counted as degrees of freedom.** The t
  intervals of `ht_total()` and `ht_mean()` counted every sampled unit,
  including the ones `design_weighted()` and `design_spread(size =)` take
  with certainty, strata and `design_temporal()` intervals taken whole, and
  such units inside a certainty design's `rest`. With eight certain units
  among ten, systematic PPS reported 9 degrees of freedom where the same
  selection built with `design_certainty()` had 1, and its interval covered
  79%; a temporal design whose sparse days were taken whole covered 83%.
  Only units below certainty count now, and a stratum or interval holding
  none of them drops out. A sample with no unit below certainty has
  `df = Inf`, and a certainty design whose `rest` drew nothing has `df` 0.
* **`as_svydesign()` left units taken with certainty in `survey`'s Brewer
  stratum.** For `design_weighted(method = "systematic")`,
  `design_spread(size =)` and a certainty design over a systematic PPS
  `rest`, `survey` still centred and scaled Brewer's approximation over the
  certain units, which put its standard error at 23 times `ht_total()`'s on
  a frame with four dominant units. They now go to `survey` in a stratum of
  their own, taken whole, and the two standard errors agreed within 0.02%
  there.
* **A blank stratum or cluster value broke stratified and multistage
  designs.** An empty string, such as a blank field read from a file, was
  looked up by name, which never matches it. `inclusion_prob()`,
  `joint_prob()` and `draw(weights = TRUE)` stopped with "subscript out of
  bounds" for a stratified design or a certainty design over one. Neyman
  allocation read the blank stratum's spread as 0 and, with
  `min_per_stratum` of 1 or more, gave its rows to the other strata without
  a word. A multistage sample got `NA` probabilities, weights and totals, and
  `sample_summary()` gave a blank group an `NA` frame count. A blank value is
  now a key like any other, shown as `(blank)` in messages and summaries.
* **Two strata could be counted as one.** Several `strata` columns were
  joined into one key with a carriage return between the values, so
  `("x\ry", "z")` and `("x", "y\rz")` made one stratum of two, and the
  allocation, the probabilities and the variance were all worked out over
  the merged stratum. Two numbers that print alike, such as `0.3` and
  `0.1 + 0.2`, were one stratum too. Keys are now built from each column's
  integer codes, which cannot collide.
* **A sample estimated under another collation could get the wrong standard
  error.** Text strata and a text `order_by` were sorted by the session's
  collation, which differs between the C locale of a server and the en_US
  locale of a desktop. Where two strata tied in the largest-remainder split,
  the extra row went to a different stratum in each, and a stratified
  sample drawn in one and estimated in the other got a standard error 5
  times too large in one example. A systematic sample walked on text got
  the wrong neighbours for the successive-difference variance, and
  `joint_prob()` gave 0 for pairs drawn together. Text is now sorted byte by
  byte, as in the C locale, in every session.
* **Two cluster labels that print alike were one cluster.** `0.3` and
  `0.1 + 0.2` print the same, and the cluster variance, `joint_prob()` for
  cluster and multistage designs, and the jackknife compared cluster labels
  as text. A `design_cluster()` sample holding both got no variance, and two
  rows of the two clusters got the joint probability of rows in one cluster.
  A multistage sample got a Sen-Yates-Grundy variance built on those joint
  probabilities. Clusters are now compared by integer codes everywhere.
* **`design_multistage()` refused a factor cluster column with unused
  levels.** An unused level counted as a cluster of no rows, so
  `inclusion_prob()`, `sampling_weight()`, `joint_prob()`,
  `draw(weights = TRUE)` and `plot(type = "probability")` stopped with a
  false "no closed-form inclusion probability". Unused levels are now
  ignored, as `droplevels()` would drop them.
* **`joint_prob()` gave `NA` for a multistage row with no cluster.** With
  `na_rm = TRUE`, a row whose cluster label is missing is never drawn, so its
  joint probabilities are 0, as `design_cluster()` gives. `design_multistage()`
  gave `NA`, which filled the default matrix over every row with `NA`.
* **`inclusion_prob()` overstated the probabilities of
  `design_multistage(replace = TRUE)`.** It gave `(a / A) * m / N_h`, the
  expected number of times a row is drawn, where the probability of drawing
  it at least once is `(a / A) * (1 - (1 - 1 / N_h)^m)`. In one example, 3
  draws from a cluster of 4 rows gave 0.375 against 0.289. A cluster smaller
  than the per-cluster take was refused as capped at its own size, which it
  is not with replacement. `draw(weights = TRUE)` refuses with-replacement designs,
  so no estimate used these figures.
* **A Hajek mean over one primary sampling unit had a standard error of 0.**
  A `by` domain with one sampled row, or with all its sampled rows in one
  cluster, and a Poisson sample of one row got a standard error of 0 and a
  zero-width interval labelled "analytic", because the linearised residuals
  of a mean sum to zero inside the one unit. `survey` gives the same 0.
  `ht_mean()` now gives such a mean no standard error, with a note, unless
  every unit of the domain was taken with certainty, which makes it exact.
  `?ht_total` adds that intervals for domains with few sampled rows cover
  less often than `level` says.
* **A systematic sample understated the variance when `interval` did not
  divide the frame.** The sample size is then random, 30 or 31 rows for 244
  rows at an interval of 8, and successive differences cannot see what that
  adds to a total. A constant `y` got a standard error of 0, and a variable
  with a coefficient of variation of 0.02 got 0.04 of its true variance and
  an interval that never covered the truth. `ht_total()`, domain totals and
  `ht_mean(estimator = "ht")` now add a term for the random size, and the
  note says so. The default Hajek mean was not affected.
* **`deff()` measured a design against a simple random sample of the whole
  frame, rows it could never reach included.** A total estimates over the
  rows with a probability above 0, so an equal-probability sample of 20 from
  the first 100 of 400 rows, `design_reservoir(max_items = 100)`, reported a
  design effect of 0.05, "better than simple random sampling", and a
  temporal window or a spatial region over part of the frame read low the
  same way. The simple random sample is now drawn from the rows the design
  can reach, as `survey` measures it, and that sample gets 1. A
  Horvitz-Thompson mean gets its total's design effect. `?deff` says so.
* **A constant `y` got a design effect and a Hajek standard error made of
  rounding error.** With unequal weights, a total's design effect came out
  as `1.36e+32`, "worse than simple random sampling", and the Hajek mean's
  standard error as about `1e-16`, with a design effect that changed from
  seed to seed. A Hajek mean of a `y` that takes one value in every sampled
  row, overall or in a `by` domain, now has a variance of exactly 0, with a
  note, and the design effect of a constant `y` is `NA`.
* **Rows a design drops were counted in its finite population correction.**
  `as_svydesign()` gave `survey` the frame's row count as the `fpc` of
  `design_systematic(na_rm = TRUE)`, and the jackknife counted every row of
  the frame for `design_spread()` and `design_weighted()` with
  `na_rm = TRUE`. Both now count the rows the design can reach, and agree
  with the same draw from the frame without the dropped rows. The standard
  errors moved by under 1% in the examples checked.
* **Size measures whose total overflowed gave every row a probability of
  0.** Positive weights near `1e306` sum to `Inf`, so
  `design_weighted()` with `method = "systematic"` or `"poisson"` and
  `design_spread(size =)` gave every row a probability of 0 and drew empty
  samples without a word, and a certainty design over them took only the
  rows above its threshold. Such weights are now divided by the largest of
  them first. Weights with a finite total give the same probabilities, and
  the same rows for a seed, as before.
* **`joint_prob()` treated a row given twice in `rows` as two rows.** The
  pair a repeated row made with itself got the joint probability of two
  different rows, 0.015 for a simple random sample of 30 from 240 where the
  row's own probability is 0.125, and `simulate = TRUE` gave it 0. It now
  gets the row's own probability on both paths. `rows` that are not whole
  numbers are refused, where 1.7 was read as row 1.

## Refusals and clearer errors

* `as_svydesign()` and `sample_summary()` give the same reason as
  `ht_total()` when a sample has lost its design, and all three name the
  steps that drop it (`merge()`, column selection, `subset()`, `transform()`,
  `cbind()`) and how to add fieldwork columns instead.
* The refusal of `weights = TRUE` for a with-replacement design, and
  `?draw`, gave the bias as "about 15%". It is a factor of
  `(n / N) / (1 - (1 - 1 / N)^n)`, about 5% at `n = N / 10`, 58% at `n = N`
  and more beyond, applied within each stratum, and both now say so and give
  the per-stratum Hansen-Hurwitz form.
* `ht_total()`, `ht_mean()`, `as_svydesign()` and `sample_summary()` refuse a
  sample whose `.prob` column was removed, is not numeric, holds a value
  outside (0, 1], or differs from the probabilities recorded at draw time. A
  removed `.prob` gave a total of 0, and a missing value gave `NA` with a
  note saying the jackknife had run. The message says that adjusted weights
  belong in `.weight`.
* `ht_total()` and `ht_mean()` refuse an infinite or `NaN` `y`, with the
  count and the rows, instead of returning an infinite total or standard
  error. A variance that overflows is reported as missing, and the note no
  longer says the jackknife ran when it gave no figure.
* `ht_total()` and `ht_mean()` warn, with class `drawn_weight_edited`, when
  `.weight` is no longer `1 / .prob`. They use the design probabilities in
  `.prob`, while `as_svydesign()` hands `.weight` to `survey`, so after a
  nonresponse adjustment the two gave different totals without a word.
  `?ht_total` and `?as_svydesign` say which column each uses.
* `ht_total()` and `ht_mean()` warn, with class `drawn_grouped_df`, when the
  sample is grouped with `dplyr::group_by()`, which they ignore, and give the
  `by` that estimates per group.
* A `by` column that the sample does not have is named in the error.
* `draw(weights = TRUE)` refuses a `design_certainty()` whose `rest` gives
  every row below the threshold a probability of 0, such as one with
  `n = 0`, because the weighted total would leave those rows out. Drawing
  without weights still gives the cut-off sample of the rows above the
  threshold.
* An error raised by a certainty design's `rest` in `draw()` or
  `inclusion_prob()` now says which rows `rest` was given, as in "In `rest`,
  applied to the 230 rows with `w` below 100". `rest` counts only those rows,
  so its own message spoke of 230 rows on a frame of 240.
* `design_multistage()` refuses an `n` smaller than `n_clusters`, `n = 0`
  included, which selected clusters and then took no rows from some of them,
  and a `min_per_cluster` that `n` cannot cover across `n_clusters`, which
  `draw()` refused only later. `draw()` warns, with class
  `drawn_empty_cluster`, when proportional allocation gives a selected
  cluster no rows, and names it.
* A multistage sample from a single cluster gets no variance, with the
  reason "A variance needs at least two clusters", as a `design_cluster()`
  sample does. The note blamed an unlucky sample for a Sen-Yates-Grundy
  estimate that is negative on every such sample, and a constant `y` got an
  exact variance of 0.
* An `allocation`, `method`, `unit`, `variance`, `estimator`, `target` or
  `type` outside its choices is refused with a message that names the
  argument, the value passed and the choices. The message from `match.arg()`
  named neither the argument nor the value. A missing cluster label is called
  a cluster label by `inclusion_prob()` as well as by `draw()`.
* `ht_total()` and `ht_mean()` stop with the error when something fails
  inside a variance estimator, such as running out of memory. Every error was
  read as the estimator declining, so `variance = "auto"` fell back to the
  jackknife and gave R's message, "cannot allocate vector of size 244.1 Mb"
  in one example, as the reason. Only a decline, a condition of class
  `drawn_declined`, leads to the jackknife now. `inclusion_prob()` and
  `joint_prob()` raise that class where a design has no closed form.
* A variance too large to represent, from values of `y / .prob` above about
  `1e154`, is reported as missing with a note to rescale `y`, whichever
  estimator ran, and the note says it once.
* `variance = "jackknife"` reports `method = "none"` when the jackknife gives
  no figure, as for a systematic or Poisson sample, which is what `?ht_total`
  says. It reported `"jackknife"`, in `by` tables too. Its note for strata
  with one unit says how many there are, where it read "needs two units per
  stratum as well" with nothing before it.
* A `df` below 1 is refused. It was accepted and the interval dropped without
  a word. When the design leaves 0 degrees of freedom, as a Poisson sample of
  one row does, the note says that, and `print()` shows it. It blamed one
  unit per stratum, and was not printed beside an analytic standard error.
* `inclusion_prob()`, `sampling_weight()`, `joint_prob()` and
  `plot(type = "probability")` refuse a design that `draw()` refuses, with
  the same message. They answered for it: an `n` of 300 from 240 rows gave
  every row 1.25 under `design_simple()`, and a census of 1s under
  `design_stratified()` and the systematic and Poisson methods of
  `design_weighted()`; an `n_clusters` above the clusters available gave
  probabilities above 1; and a window shorter than one interval, or a key
  missing in every row under `na_rm = TRUE`, gave every row 0. The same
  holds for an `n` above the rows inside a spatial region, for a region
  with none, and for the `rest` of a certainty design, checked against the
  rows below its threshold. `design_multistage()` refuses an `n` that even
  its largest `n_clusters` clusters cannot hold, on both paths.
* The weights of `design_weighted()`, and of a certainty design over it,
  must be numeric on every path. `inclusion_prob()`, `joint_prob()`,
  `draw(weights = TRUE)` and `plot(type = "probability")` took a column of
  `TRUE` as weights of 1 and stopped with a misleading message on a factor
  or text one, where `draw()` refused all three. Every message about the
  weights names the column, and the one for a weight of 0 or below says to
  drop those rows from the frame, since every design that selects in
  proportion to size refuses them. `?inclusion_prob` and `?ht_mean` said
  such rows get a probability of 0.
* `inclusion_prob()`, `sampling_weight()`, `joint_prob()` and the
  successive-difference variance refuse a `design_systematic()` whose
  `order_by` column is not in the frame, as `draw()` did. A misspelt
  `order_by` gave every row a probability of 0.
* A missing coordinate under `design_spatial()` is reported in the same
  words by `draw()` and `inclusion_prob()`, naming the two `coords`
  columns.

## Behaviour changes

* `draw(weights = TRUE)` adds a third column, `.row`, after `.prob` and
  `.weight`. It holds the row of `data` each sampled row came from. A frame
  that already has a `.row` column is refused, as one with `.prob` or
  `.weight` is.
* A sample drawn from a data.table takes `:=` without the "Invalid
  .internal.selfref" warning. data.table joins `Suggests`.
* `ht_total()` of an empty sample from a design of random size, which Poisson
  sampling or a systematic start past the last row can give, returns a total
  of 0 with no variance and a note, where it stopped with an error. Over
  repeated draws that 0 is what keeps the total unbiased, so a simulation loop
  no longer has to catch it. `ht_mean()` refuses the Hajek mean of an empty
  sample, which is 0 / 0, and gives 0 for `estimator = "ht"`. An empty sample
  from a design that reaches no row, such as one with `n = 0`, is still
  refused, and the message says so. `?ht_total` describes this.
* Text keys are sorted byte by byte, as in the C locale, whatever the
  session's locale: the order in which strata are drawn and in which they
  break allocation ties, the rows of `by` estimates and of
  `sample_summary()`, and the walk of a text `order_by`. Factor columns keep
  their level order. In a locale whose collation differs from C, a seed can
  therefore select different rows than it did in drawn 0.2.0, and a tie can
  go to another stratum. `?design_stratified` gives the tie rule.
* `by` estimates come out in key order: numbers in sequence (1, 2, 10, not
  1, 10, 2) and factors in level order, where they were sorted as text.
  Domains with no sampled row are still left out, as `survey::svyby()`
  leaves them out.
* An sf data frame can be passed as `data` to `design_stratified()` and
  `design_spread()`, and an sf sample to `by` in `ht_total()` and `ht_mean()`
  and to `sample_summary()` of a stratified or cluster sample. All of these
  stopped inside base R on the geometry column. A `by` estimate of an sf
  sample is a plain data frame without the geometry.
* When a multistage sample with one row per cluster falls back to the
  jackknife and more than a fifth of the clusters were selected, the note
  gives the factor by which the jackknife inflates the variance between
  clusters. When every cluster was selected, it says that variance is zero by
  design and points to `design_stratified()`.
* A `by` domain gets a variance when a stratum or interval with one sampled
  row holds none of its units, read from the `by` columns of the frame
  stored with the sample. A domain of strata or intervals taken whole is
  exact, with a zero-width interval even at 0 degrees of freedom. Such
  domains were refused with a note saying the standard error would be
  understated, which was false for them. A domain with units in the lonely
  stratum is still refused, and when the frame has no `by` columns the note
  says the frame cannot tell.
* `ht_mean()` decides whether `y` is a proportion, and so whether its
  interval is formed on the logit scale, from what `y` is rather than from
  the values the sample drew: a logical `y`, or a column of the frame whose
  values are all 0 or 1. A count that came out 0 or 1 in every sampled row
  got a logit interval, which cannot pass 1 and so never covered a mean of 1
  or more; it now gets the symmetric interval. A 0/1 column added to the
  sample after the draw, or a numeric vector, gets the symmetric interval
  unless `proportion = TRUE` says it is a proportion. `by` estimates follow
  the same rule.
* A proportion estimated at exactly 0 or 1, such as no deviations in 60
  sampled invoices, gets Korn and Graubard's beta interval, with the sampled
  rows below certainty as its effective sample size, and a note saying so.
  It got a zero-width interval of 0 to 0. For 0 in 60 from a simple random
  sample the 95% upper bound is now 4.87%. A census, or a domain of rows
  taken with certainty, is still exact.
* Estimates from `ht_total()` and `ht_mean()` hold `unreachable`, the number
  of frame rows with inclusion probability 0, and print it on a line of its
  own when there are any, whether or not there is a standard error. `by`
  tables carry it as an attribute and print it too. `draw()` warns about
  such rows, but a sample saved and read back, or drawn with the warning
  muffled, gave a total over part of the frame with no sign of the rest.
* The printout of `ht_total()` and `ht_mean()` shows the note whenever the
  standard error is 0, such as the one saying every sampled row was taken
  with certainty.
* `sample_summary()` gains `N_reachable`, the rows the design gives a
  probability above 0, and prints the sampling fraction over them. `N` is
  still every row of the frame.

## New

* `ht_total()` and `ht_mean()` take `alternative = "less"` for a one-sided
  upper bound and `"greater"` for a lower bound, labelled as such when
  printed. The open end of `ci` is `-Inf` or `Inf`, or 0 or 1 for a
  proportion. A one-sided bound is the end of a two-sided interval at
  `1 - 2 * (1 - level)`, which is how it had to be asked for before.
* `ht_mean()` takes `proportion`, to say whether `y` is a proportion rather
  than leave it to the rule above, and `ci_method = "beta"` for Korn and
  Graubard's interval, the one `survey::svyciprop(method = "beta")` gives,
  with its effective sample size capped at the number of sampled rows below
  certainty.
* `ht_total()` and `ht_mean()` take `contrast`, weights over the `by`
  domains such as `c(a = 1, b = -1)`, and return the weighted sum of the
  domain estimates with its standard error, interval and degrees of freedom,
  as a `drawn_contrast`. The variance comes from one run of the design's
  estimator, so the covariance between domains is counted, and for simple
  random, stratified, cluster and multistage samples the standard error
  equals that of `survey::svycontrast()` on `survey::svyby(covmat = TRUE)`.
  A difference of domain means needed linearised residuals built by hand,
  and the `survey` route stopped with an error for systematic PPS samples.
  The interval is symmetric, for proportions too.
* `certainty_set()` gives the rows of the frame that `design_weighted()`
  with `method = "systematic"` or `"poisson"`, or `design_spread()` with
  `size`, takes with certainty, before anything is drawn. Its `threshold`
  attribute, the smallest size among them, can be passed to
  `design_certainty()` to take the same rows with the cut-off made
  explicit.

## Speed and memory

* Stratum and `by` keys are built from each column's integer codes, so
  their cost grows with the number of rows, however many combinations the
  columns' levels could form. On 200,000 rows keyed by two columns of 3,000
  values each, the key took 8 seconds where it now takes 0.2, and at 6,000
  values each it ran out of memory.
* `inclusion_prob()`, `draw(weights = TRUE)` and the estimators of a
  stratified design take time linear in the number of strata: 0.02 seconds
  at 40,000 strata, where they took 9.
* `design_cluster()` and `design_multistage()` find the rows of the
  selected clusters in one pass over the frame instead of one pass per
  cluster. Selecting 1,000 clusters from a million rows took 16 seconds and
  now takes 0.2. The clusters chosen for a seed are unchanged.
* The jackknife sums each unit with `rowsum()`, so its time grows with the
  number of rows rather than with rows times units. 20,000 rows took 2.6
  seconds and now take 0.03.
* The analytic variance of simple random, reservoir, spatial, stratified,
  temporal and equal-take multistage samples, alone or as the `rest` of
  `design_certainty()`, is computed from the variance within each stratum,
  interval or cluster, in time and memory that grow with the number of
  sampled rows. The Sen-Yates-Grundy sum over every pair of sampled rows,
  which gives the same figures, held several matrices with a row and a
  column per sampled row: a stratified sample of 4,000 rows took 4.8 seconds
  and 1.4 GB, growing with the square of the sample, and every `by` domain
  built them again. The same sample now takes 0.01 seconds, and 10,000 rows
  take 0.04. A census of one row gets its variance of 0 this way rather than
  from the jackknife.
* `joint_prob()` refuses a matrix that would need more than about 1 GB,
  from about 11,200 rows, and gives the size and `rows = s$.row` as the way
  to ask for the rows of a sample `s`. It used to go on and could exhaust
  memory. The full matrix is still the default. Stratified, temporal,
  cluster, multistage and systematic joint probabilities are filled one
  group at a time and need little memory beyond the matrix itself, where they
  needed 4 to 12 times its size: a stratified matrix over 5,000 rows took
  2.4 GB and 6.7 seconds, and now takes 0.2 GB and 0.4 seconds.
* `by` estimates work out the strata, clusters, systematic order, certainty
  split, spread neighbours and joint probabilities once per call rather than
  once per domain, and the variances of simple random, reservoir, spatial,
  stratified and temporal samples, alone or as the `rest` of
  `design_certainty()`, come from one pass over all the domains. 200 domains
  of a stratified sample of 1,000 from 100,000 rows took 3.4 seconds and now
  take 0.1; 200 domains of a spread sample of 1,000 from 20,000 rows took 16
  seconds and now take 0.3. The standard errors agree with the old ones to
  rounding error.

## Handing off to survey

* `as_svydesign()` hands `survey` a copy of the sample, so an in-place sort
  of the sample afterwards no longer changes the survey design's standard
  errors.
* Units taken with certainty go to `survey` as one primary sampling unit
  per stratum with an `fpc` of 1, through a new `.psu` id column, so
  `survey::degf()` gives the same degrees of freedom as `ht_total()` for
  certainty designs, strata and intervals taken whole, and PPS samples with
  certain units, and `survey::svyciprop()` gives the same interval as
  `ht_mean()`. `survey` had counted each certain row as a unit.
* A census, a sample whose every `.prob` is 1, goes to `survey` as
  `ids = ~1` with `fpc` the number of rows. Under systematic PPS, alone or
  as a certainty design's `rest`, Poisson sampling and
  `design_spread(size =)` it stopped `svydesign()` with "subscript out of
  bounds". An empty sample is refused with a message
  instead of failing inside `survey`, and so is a certainty sample whose
  `rest` drew none of the rows below the threshold, which `survey` would
  have reported as exact.
* `?as_svydesign`, the estimation article and `cran-comments.md` no longer
  say that `survey` agrees with systematic PPS within 0.2%. Standard errors
  of totals usually agree within 0.3% and within about 1%; those of means
  differ by up to 7% at a sampling fraction of one half.
* `as_svydesign()` labels strata built from several columns as `a . b`.
  The labels held a carriage return, which broke a CSV export of the design
  variables. Where two labels would read alike, the strata are numbered
  instead, so `survey` cannot merge them.
* `as_svydesign()` numbers the clusters in `.cluster`. `survey` made a
  factor of the labels, which merged two that print alike, such as `0.3` and
  `0.1 + 0.2`, into one primary sampling unit.

## Documentation and site

* `?ht_total` said the exact multistage estimator understates by about a
  third with one row per cluster. It understates by `f1` times the
  within-cluster component, where `f1` is the share of clusters selected, and
  the jackknife used in its place overstates by `f1` times the
  between-cluster component, which came to 3.3 times the true variance in one
  example with 30 of 40 clusters selected. `?design_multistage` and
  `?inclusion_prob` give the probability under `replace = TRUE`.

# drawn 0.2.0

An adversarial review of 0.1.0 found estimates that were wrong without saying
so. Every fix below was reproduced by simulation first, and each now has a
test that fails on the old behaviour.

## Wrong answers, fixed

* **One sampled row per stratum reported a variance of 0.** Sen-Yates-Grundy
  has no within-stratum pair to sum over, so a stratified sample with one row
  in every stratum printed a standard error of exactly zero and a zero-width
  interval. `ht_total()` and `ht_mean()` now decline with the reason, as
  `survey` does by default, for strata and for `design_temporal()` intervals.
* **A stratum allocated no rows silently biased the total.** Proportional
  allocation gives nothing to a stratum smaller than about `N / n`, and Neyman
  allocation gives nothing to one whose auxiliary variable is constant. Its
  rows then have probability 0 and drop out of every estimate; one example was
  84% low. `draw()` now warns, naming the strata and the fix. The
  documentation previously claimed `min_per_stratum = 0` "leaves allocation
  unbiased" and that a floor biases a total; both were wrong, and are
  corrected.
* **The jackknife ignored strata.** It deleted rows across the whole sample,
  counting the differences *between* strata as sampling variance: sixty times
  too large on a well-stratified frame. It is now the stratified (JKn)
  jackknife, and equals the analytic figure exactly for `design_stratified()`
  and `design_temporal()`.
* **One row per cluster understated a multistage variance by ~30%.** The
  exact two-stage estimator cannot see within-cluster variation when only one
  row is taken per cluster. It now declines, and `variance = "auto"` falls back
  to the delete-a-cluster (ultimate-cluster) jackknife.
* **Neyman allocation abandoned its own rule when a stratum was capped.** The
  surplus from a stratum too small for its share was dealt out round-robin;
  it is now re-split in proportion to `size * sd`, as in Cochran (1977, 5.9).
* `ht_mean(estimator = "ht")` printed itself as a Hajek mean.
* `design_systematic(order_by = )` returned rows in sorted order, breaking the
  documented promise that rows come back in frame order.
* The documentation for `design_weighted(method = "systematic")` said some
  pairs of rows can never be drawn together and pointed to
  `sampling::UPsystematicpi2()` for joint probabilities. The rows are shuffled
  before the walk, so every pair can co-occur and that function answers a
  different design.

## Variance estimators from the literature

Every estimate now reports which estimator produced its variance in `method`,
and names any approximation when it prints.

* **Confidence intervals use t on the design's degrees of freedom** — primary
  sampling units minus strata (Korn and Graubard 1999), matching
  `survey::degf()`. At four clusters the normal interval covered about 87%;
  the t interval restores nominal coverage. `df = Inf` gives the normal
  interval.
* **Systematic PPS uses Deville's (1999) approximation**, which needs only
  first-order probabilities and scales each row by its own `1 - pi`. It
  equals `sampling::varest()` exactly and tracks the empirical variance at
  0.98, where the previous jackknife ran at 0.89.
* **Systematic samples get the successive-difference approximation** (Wolter
  2007) rather than no variance at all. On a frame sorted by something related
  to `y` it stays near the truth where the simple random formula overstates
  two- to thirtyfold. Its limits — a cycle matching the interval, or a trend
  so smooth the random start is the only variation — are documented.

## New

* `design_spread()` draws a spatially balanced sample by the local pivotal
  method (Grafström, Lundström and Schelin 2012): exact first-order
  inclusion probabilities, equal or proportional to `size`, with neighbours
  competing so the sample covers the space evenly. On a smooth surface its
  variance was 14–25% of a simple random sample's. It spreads across any
  numeric columns, so it also balances a sample on covariates without
  choosing strata. `ht_total()` uses the local mean variance estimator
  (Stevens and Olsen 2003) for it.
* `ht_total()` and `ht_mean()` take `by =` for domain estimates computed from
  the whole sample rather than a subset of it. Standard errors match
  `survey::svyby()` exactly.
* `ht_mean()` accepts a logical `y`, so a proportion is one call.
* `plot(design, data, type = "map")` places every row by two columns and
  fills in the selected ones — the view that shows whether a sample covers the
  ground.
* Key functions carry `@references` to the methods they implement.

## Handing off to survey

* Multistage designs are expressed as two stages, so `survey`'s standard error
  now matches `ht_total()` exactly; it previously differed by up to 5–10%.
* Systematic PPS, alone or as a certainty design's `rest`, maps to
  `pps = "brewer"` and agrees with `ht_total()` within 0.2%.
* A reservoir capped by `max_items` uses the rows the stream reached as its
  population, not the whole frame.

## Documentation and site

* The README is generated from `README.Rmd`, so every output shown is real.
* New articles: a getting-started walkthrough, how the estimates work (every
  estimator, its reference and its simulation check), an audit-sampling
  walkthrough, and a spatially balanced field survey.
* `graphics` and `utils` are declared in `Imports`.

# drawn 0.1.0

First release.

## Designs are objects

Sampling is expressed as a design you build and then apply, rather than a
function call with the design buried in its arguments:

```r
plan <- design_stratified(strata = "site", n = 500)
draw(data, plan, seed = 1)
```

A design can be printed, stored, reused across data sets, and — the part that
matters for estimation — asked what it does before you draw anything.

Eleven designs ship: `design_simple()`, `design_stratified()`,
`design_systematic()`, `design_cluster()`, `design_multistage()`,
`design_weighted()`, `design_certainty()`, `design_reservoir()`,
`design_bootstrap()`, `design_temporal()` and `design_spatial()`.

`design_certainty()` composes: rows at or above a threshold are all taken, and
any other design is applied to what remains. That is the standard shape of an
audit sample, and the cheapest way to cut variance when a few units dominate a
total. Certainty rows have inclusion probability exactly 1, so they contribute
their own value with a weight of 1 and add nothing to the variance; inclusion
and joint probabilities compose exactly across the two parts.

## One contract across all of them

* `n` is always the total number of rows drawn, never a per-group figure.
  `design_temporal()` uses `per_interval` precisely because that one *is* per
  group.
* `allocation` always says how a total is split across groups.
* `na_rm` always decides whether missing keys are dropped or raise an error,
  and means the same thing whether you are drawing or asking for probabilities.
* `replace` always rules out `weights = TRUE`: an inclusion probability
  describes distinct units, and a sample holding duplicates cannot be weighted
  by one.
* `draw()` restores the caller's random number stream before returning, so
  sampling inside a simulation does not shift the simulation's own draws.
* What comes back has the same class and the same columns, in the same order,
  as what went in — tibble in, tibble out, with or without weights. Rows come
  back in frame order except where draw order is meaningful: `design_simple()`
  and `design_weighted()` return draw order, and `design_bootstrap()` returns
  replicates behind a leading `.replicate` column.

## Inclusion probabilities and estimation

`inclusion_prob()` gives the first-order inclusion probability of every
population row under a design, without drawing anything — a quick way to check
a plan before committing to it. `draw(weights = TRUE)` attaches those
probabilities and the corresponding design weights to a sample as `.prob` and
`.weight`.

`joint_prob()` gives second-order probabilities, and `ht_total()` and
`ht_mean()` use them to return an estimate with a standard error and confidence
interval:

```r
s <- draw(pop, design_stratified("site", n = 40), seed = 1, weights = TRUE)
ht_total(s, "spend")
ht_mean(s, "spend")
```

`ht_mean()` defaults to the Hajek estimator, `sum(y/pi) / sum(1/pi)`, which
divides by the estimated population size rather than the known one;
`estimator = "ht"` divides by the true `N`. The two coincide exactly whenever
the weights of the drawn rows sum to `N` — every fixed-size equal-probability
design, and every stratified design without replacement — and part company when
the sample size is random or the weights vary within a fixed-size sample. Hajek
is usually the steadier of the two, which is why it is the default; `?ht_mean`
sets out when it is not, and why `"ht"` is unbiased only where every row is
reachable.

`deff()` reports the design effect on either result — the design's variance
against simple random sampling of the same size, which reads as an exchange rate
on sample size. Stratifying on a variable that matters can put it far below 1;
clustering on the same variable puts it well above.

The variance estimator matches how each design randomises: Sen-Yates-Grundy for
fixed-size designs, the independent-units form for Poisson sampling, the
cluster-level form for cluster designs — whose row count is random whenever
clusters differ in size — and, for a certainty design, whatever `rest` uses,
since certainty rows are in every possible sample and contribute nothing.

Where a design has inclusion probabilities but no closed-form joint ones,
`ht_total()` falls back to a delete-a-group jackknife and reports which method
it used, or `"none"` where neither could produce a figure. The jackknife is
declined for systematic sampling, which has a single primary sampling unit, and
for Poisson sampling, which has an exact variance already.

Every analytic estimator is checked in the test suite against the empirical
sampling variance of its own estimator, every inclusion and joint probability
against the observed frequency over thousands of draws, and every design that
`survey` can express against `survey::svytotal()`.

## Planning a sample, and checking one

`plan_size()` solves for `n` given the precision wanted, rather than asking for
a guess. It takes a margin of error and a measure of spread, and applies three
corrections: a finite population correction from `N`, an inflation for the
design from `deff`, and an inflation for non-response from `response`. It works
for a mean, a proportion or a total, and says so plainly when the margin asked
for is unreachable by sampling.

```r
plan_size(margin = 5, sd = 40, N = 20000, deff = 2.5, response = 0.7)
```

`sample_summary()` reports what was drawn against what was there: the sampling
fraction, the range and coefficient of variation of the design weights, the
per-stratum or per-cluster counts against the frame, and the number of rows the
design could never have reached.

## Handing off to the survey package

`as_svydesign()` builds a `survey::svydesign()` object from a drawn sample, so
the analysis this package does not do — subpopulation estimates, regression,
calibration, quantiles with proper standard errors — can be done by the package
that does.

The two packages compute variance from different starting points: this one from
the design's joint inclusion probabilities, `survey` from the design's shape.
Each design is therefore expressed in `survey`'s own terms — strata for a
stratified or temporal design, cluster ids for a cluster design,
`survey::poisson_sampling()` for Poisson, and a taken-whole stratum for the
certainty rows. Totals and standard errors then agree to floating point for
every design but two: `survey` uses the ultimate-cluster approximation for
multistage designs, and returns the conservative simple-random figure for
systematic sampling where `ht_total()` declines to return one. `?as_svydesign`
has the table. Compositions with no single `survey` equivalent — a certainty
design over a cluster, multistage or Poisson `rest` — are refused rather than
approximated.

## Seeing a design

`plot(design, data)` shows which rows a design selects — every frame row as a
dot, selected ones filled in — and `type = "probability"` plots inclusion
probability against frame position. Base graphics, no plotting dependency.

## Allocation

`design_stratified()` accepts `allocation = "neyman"` with `allocation_by`,
giving strata shares proportional to `size * sd`. That minimises the variance of
a total for a fixed `n` by putting more rows where the values vary most.

## Saying when it cannot

Not every design has a closed-form inclusion probability, and the ones that
don't report that rather than returning an approximation:

* `design_cluster(balanced = TRUE)` — the per-cluster take is the smallest
  *selected* cluster's size, which is random.
* `design_multistage(allocation = "proportional")` — the second-stage
  allocation depends on which clusters were selected.
* `design_multistage()` where the per-cluster take is not constant — `n` not
  divisible by `n_clusters`, or a cluster smaller than `n / n_clusters`. The
  allocation runs over the *selected* clusters, so a row's probability depends
  on which others were drawn.
* `design_weighted(method = "successive")` — successive sampling has no closed
  form. This is the default, and the reason the other two methods exist.
* `design_bootstrap()` — resampling a sample is not a probability sample of a
  finite population.

`simulate = TRUE` estimates any of them by Monte Carlo, for first-order
probabilities and, in `joint_prob()`, for second-order ones — which is the
general answer where no formula exists, including
`design_weighted(method = "systematic")`, whose joint probabilities depend on
the order units are visited. `design_bootstrap()` is the one exception, and
simulating it is refused: every row appears in some replicate, so the count
converges to 1 for all of them. Systematic sampling has first-order
probabilities but no design-unbiased variance at all, because most pairs of rows
can never co-occur; `ht_total()` says so instead of returning a number.

## Probability proportional to size

`design_weighted()` offers three methods, and they are not interchangeable:

| `method` | Inclusion probabilities | Sample size |
|---|---|---|
| `"successive"` (default) | Not proportional to weight, no closed form | Fixed |
| `"systematic"` | Exactly `n * p_i` | Fixed |
| `"poisson"` | Exactly `n * p_i` | Random, mean `n` |

The default is what `base::sample(prob = )` does: a good way to bias selection
toward heavy units, but the weights govern each sequential draw rather than the
probability of ending up in the sample. Use `"systematic"` or `"poisson"` when
the sample will be estimated from. Units heavy enough that `n * p_i > 1` are
taken with certainty and the remainder rescaled, iterating until every
probability is valid.

## Notes

* `sf` is a suggested dependency, needed only by `design_spatial()`, and
  `survey` only by `as_svydesign()`. Nothing else requires either.
* Under spherical geometry a longitude/latitude polygon whose edge spans more
  than 180 degrees is drawn the short way across the antimeridian, so a "whole
  world" rectangle collapses to a narrow strip. `design_spatial()` warns when a
  region has such an edge.
* Reservoir sampling takes a vectorised path for a data frame, whose length is
  already known, and Algorithm L for a genuine stream — a connection or a
  generator function.
* Calendar units in `design_temporal()` step by calendar months and years, not
  by fixed 30.44-day durations.
