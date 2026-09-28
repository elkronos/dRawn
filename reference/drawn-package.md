# drawn: Design-Based Sampling with Known Inclusion Probabilities

Draws probability samples from data frames under an explicit, reusable
sampling design: simple random, stratified (proportional, equal or
Neyman), systematic, cluster, multi-stage, probability-proportional-to-
size, certainty-plus-sample, spatially balanced (the local pivotal
method of Grafstrom, Lundstrom and Schelin, 2012), reservoir, bootstrap,
temporal and region-restricted. A design is a value that can be stored,
printed and reused, and that knows its own first- and second-order
inclusion probabilities, so a drawn sample carries the weights needed
for Horvitz-Thompson estimation: 'ht_total()' and 'ht_mean()' return a
population total, mean or proportion, overall or by domain, with a
standard error, a confidence interval on the design's degrees of freedom
and a design effect. Each design gets the variance estimator the
literature recommends for it, and says which: Sen-Yates-Grundy where
joint probabilities are known, Deville's approximation for
probability-proportional-to-size, successive differences for systematic
samples, local means for spatially balanced ones, and a stratified
jackknife otherwise. 'plan_size()' solves for the sample size a target
margin of error requires, and 'as_svydesign()' hands a sample to the
'survey' package for further analysis. Designs whose inclusion
probabilities have no closed form report that rather than supplying an
approximation, and can be estimated by simulation instead.

## See also

Useful links:

- <https://github.com/elkronos/dRawn>

- <https://elkronos.github.io/dRawn/>

- Report bugs at <https://github.com/elkronos/dRawn/issues>

## Author

**Maintainer**: Justin Chase <jchase.msu@gmail.com>
