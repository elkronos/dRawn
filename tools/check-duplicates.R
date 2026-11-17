#!/usr/bin/env Rscript
# Guard against a function being defined in more than one file under R/.
# R sources R/ alphabetically, so a stale copy left behind by a refactor
# silently wins over the intended one -- which is how a fixed allocate() went
# on using its unfixed twin.
#
# Each file is read with parse(), and every top-level assignment of a
# function counts: `<-` or `=`, a backtick or quoted name such as `%||%`, a
# definition split over several lines, and each name in a chained
# assignment.
#
#   Rscript tools/check-duplicates.R              checks R/
#   Rscript tools/check-duplicates.R --self-test  checks the checker on
#                                                 tools/fixtures/duplicates

is_assignment <- function(e) {
  is.call(e) && length(e) == 3L &&
    (identical(e[[1]], as.name("<-")) || identical(e[[1]], as.name("=")))
}

assigned_name <- function(lhs) {
  if (is.name(lhs)) return(as.character(lhs))
  if (is.character(lhs) && length(lhs) == 1L) return(lhs)
  NA_character_
}

# The names a top-level expression binds to a function, if any.
function_names <- function(e) {
  out <- character(0)
  while (is_assignment(e)) {
    out <- c(out, assigned_name(e[[2]]))
    e <- e[[3]]
  }
  if (is.call(e) && identical(e[[1]], as.name("function"))) {
    out[!is.na(out)]
  } else {
    character(0)
  }
}

function_defs <- function(files) {
  rows <- lapply(files, function(f) {
    exprs <- parse(f, keep.source = TRUE)
    refs <- attr(exprs, "srcref")
    found <- lapply(seq_along(exprs), function(i) {
      nm <- function_names(exprs[[i]])
      if (!length(nm)) return(NULL)
      data.frame(fn = nm, file = basename(f), line = refs[[i]][1])
    })
    do.call(rbind, found)
  })
  defs <- do.call(rbind, rows)
  if (is.null(defs)) defs <- data.frame(fn = character(0), file = character(0),
                                        line = integer(0))
  defs
}

duplicated_defs <- function(defs) {
  dupes <- defs[defs$fn %in% defs$fn[duplicated(defs$fn)], ]
  dupes[order(dupes$fn, dupes$file, dupes$line), ]
}

if ("--self-test" %in% commandArgs(trailingOnly = TRUE)) {
  files <- list.files("tools/fixtures/duplicates", pattern = "\\.R$",
                      full.names = TRUE)
  found <- sort(unique(duplicated_defs(function_defs(files))$fn))
  expected <- sort(c("%or%", "also_chained", "equals", "plain", "quoted",
                     "split_def"))
  if (!identical(found, expected)) {
    cat("Self-test failed.\n  expected: ", paste(expected, collapse = ", "),
        "\n  found:    ", paste(found, collapse = ", "), "\n", sep = "")
    quit(status = 1)
  }
  cat("Self-test passed: ", length(expected), " duplicates found in the ",
      "fixtures.\n", sep = "")
  quit(status = 0)
}

defs <- function_defs(list.files("R", pattern = "\\.[Rr]$", full.names = TRUE))
dupes <- duplicated_defs(defs)
if (nrow(dupes)) {
  cat("Duplicate definitions:\n"); print(dupes, row.names = FALSE)
  quit(status = 1)
}
cat("No duplicate definitions among ", nrow(defs), " functions across ",
    length(unique(defs$file)), " files.\n", sep = "")
