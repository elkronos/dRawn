# Fixture for `Rscript tools/check-duplicates.R --self-test`. Not package
# code. Every definition below except `chained`, `once` and `not_a_function`
# is repeated in b.R, mostly in forms a line-based match misses.
plain <- function(x) x
`%or%` <- function(a, b) if (is.null(a)) b else a
equals = function(x) x
"quoted" <- function(x) x
split_def <-
  function(x) x
chained <- also_chained <- function(x) x
once <- function(x) x
not_a_function <- 1
