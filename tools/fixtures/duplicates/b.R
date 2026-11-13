# Fixture for `Rscript tools/check-duplicates.R --self-test`. Not package
# code.
plain <- function(x) x + 1
`%or%` = function(a, b) a
equals <- function(x) x
quoted = function(x) x
split_def = function(
  x
) x
also_chained <- function(x) x
not_a_function <- function(x) x
not_a_function_either <- 2
