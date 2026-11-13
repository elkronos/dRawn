# Fixture for `Rscript tools/check-duplicates.R --self-test`. Not package
# code: every definition below is repeated in b.R, in a form a line-based
# match misses, except `chained`, `once` and `not_a_function`.
plain <- function(x) x
`%or%` <- function(a, b) if (is.null(a)) b else a
equals = function(x) x
"quoted" <- function(x) x
split_def <-
  function(x) x
chained <- also_chained <- function(x) x
once <- function(x) x
not_a_function <- 1
