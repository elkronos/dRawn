# drawn: Design-Based Sampling with Known Inclusion Probabilities

Draws probability samples from data frames under an explicit, reusable
sampling design: simple random, stratified, systematic, cluster,
multi-stage, weighted (including probability-proportional-to-size),
reservoir, bootstrap, temporal and spatial. A design is a value that can
be stored, printed and reused, and that knows its own first- and
second-order inclusion probabilities, so a drawn sample carries the
weights needed for Horvitz-Thompson estimation and 'ht_total()' returns
a population total with a standard error and confidence interval.
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
