# drawn: Design-Based Sampling with Known Inclusion Probabilities

Draws probability samples from data frames under an explicit, reusable
sampling design: simple random, stratified, systematic, cluster,
multi-stage, weighted (including probability-proportional-to-size),
certainty-plus-sample, reservoir, bootstrap, temporal and spatial. A
design is a value that can be stored, printed and reused, and that knows
its own first- and second-order inclusion probabilities, so a drawn
sample carries the weights needed for Horvitz-Thompson estimation:
'ht_total()' and 'ht_mean()' return a population total or mean with a
standard error, confidence interval and design effect. 'plan_size()'
solves for the sample size a target margin of error requires, and
'as_svydesign()' hands a sample to the 'survey' package for analysis.
Designs whose inclusion probabilities have no closed form report that
rather than supplying an approximation, and can be estimated by
simulation instead.

## See also

Useful links:

- <https://github.com/elkronos/dRawn>

- <https://elkronos.github.io/dRawn/>

- Report bugs at <https://github.com/elkronos/dRawn/issues>

## Author

**Maintainer**: Justin Chase <jchase.msu@gmail.com>
