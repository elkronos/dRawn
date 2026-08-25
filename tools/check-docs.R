#!/usr/bin/env Rscript
# Documentation guards, run in CI.
#
#  1. Every argument of an exported function is documented, and nothing is
#     documented that is not an argument.
#  2. Every [foo()] cross-reference in R/, README and NEWS resolves.
#  3. Nothing references an identifier the package no longer has.
#  4. Every design constructor appears in the overview topic, README and NEWS.
#
# R CMD check catches undocumented objects but not documentation that has
# quietly gone stale, which is the failure mode this covers.

suppressMessages(library(drawn))
bad <- 0L
note <- function(...) { cat("  ", ..., "\n", sep = ""); bad <<- bad + 1L }

rd_files <- list.files("man", pattern = "\\.Rd$", full.names = TRUE)
lines_of <- function(f) readLines(f, warn = FALSE)

# \item{} entries inside \arguments{} only -- \describe{} blocks in a
# @section also use \item{} and are not parameters.
rd_args <- function(f) {
  txt <- paste(lines_of(f), collapse = "\n")
  m <- regexpr("\\\\arguments\\{", txt)
  if (m < 0) return(character(0))
  rest <- substring(txt, m + attr(m, "match.length"))
  depth <- 1L; end <- 0L
  ch <- strsplit(rest, "")[[1]]
  for (i in seq_along(ch)) {
    if (ch[i] == "{") depth <- depth + 1L
    if (ch[i] == "}") depth <- depth - 1L
    if (depth == 0L) { end <- i; break }
  }
  block <- substring(rest, 1, max(end - 1L, 0L))
  items <- regmatches(block, gregexpr("\\\\item\\{([^}]*)\\}", block))[[1]]
  unique(trimws(unlist(strsplit(gsub("\\\\item\\{|\\}", "", items), ","))))
}
rd_alias <- function(f) {
  gsub("^\\\\alias\\{|\\}$", "", grep("^\\\\alias\\{", lines_of(f), value = TRUE))
}

exports <- getNamespaceExports("drawn")

cat("1. arguments documented\n")
for (f in rd_files) {
  aliases <- intersect(rd_alias(f), exports)
  if (!length(aliases)) next
  documented <- rd_args(f)
  all_formals <- unique(unlist(lapply(aliases, function(a)
    setdiff(names(formals(get(a, asNamespace("drawn")))), "..."))))
  for (a in aliases) {
    fm <- setdiff(names(formals(get(a, asNamespace("drawn")))), "...")
    miss <- setdiff(fm, documented)
    if (length(miss)) note(a, "(): undocumented -> ", paste(miss, collapse = ", "))
  }
  extra <- setdiff(documented, all_formals)
  if (length(extra)) note(basename(f), ": documents non-arguments -> ",
                          paste(extra, collapse = ", "))
}

cat("2. cross-references resolve\n")
all_alias <- unlist(lapply(rd_files, rd_alias))
for (f in c(list.files("R", full.names = TRUE), "README.md", "NEWS.md")) {
  if (!file.exists(f)) next
  txt <- lines_of(f)
  links <- unique(gsub("\\[|\\(\\)\\]", "",
    unlist(regmatches(txt, gregexpr("\\[[a-zA-Z_.][a-zA-Z0-9_.]*\\(\\)\\]", txt)))))
  for (l in links) {
    if (!(l %in% all_alias) && !exists(l, where = asNamespace("drawn")) &&
        !exists(l, where = baseenv())) {
      note(basename(f), ": [", l, "()] does not resolve")
    }
  }
}

cat("3. no references to removed API\n")
removed <- c("max_rows", "normalization", "equal_samples", "design_weight",
             "proportional_stage_two", "stage_two_sample_size", "data_stream",
             "simple_random_sampling", "stratified_sampling", "cluster_sampling",
             "multi_stage_sampling", "weighted_sampling", "temporal_sampling",
             "spatial_sampling", "reservoir_sampling", "bootstrap_sampling",
             "systematic_sampling")
for (f in c(list.files("R", full.names = TRUE), rd_files, "README.md",
            "NEWS.md", "DESCRIPTION", "_pkgdown.yml")) {
  if (!file.exists(f)) next
  txt <- lines_of(f)
  for (r in removed) {
    hit <- grep(paste0("\\b", r, "\\b(?!ed)"), txt, perl = TRUE)
    if (length(hit)) note(basename(f), ":", hit[1], " mentions removed `", r, "`")
  }
}

cat("4. every design is documented everywhere\n")
ctors <- sort(grep("^design_", exports, value = TRUE))
for (f in c("man/designs.Rd", "README.md", "NEWS.md")) {
  if (!file.exists(f)) next
  txt <- paste(lines_of(f), collapse = " ")
  miss <- ctors[!vapply(ctors, grepl, logical(1), x = txt, fixed = TRUE)]
  if (length(miss)) note(basename(f), " omits ", paste(miss, collapse = ", "))
}

if (bad) {
  cat("\n", bad, " documentation problem(s)\n", sep = "")
  quit(status = 1)
}
cat("\nDocumentation consistent: ", length(ctors), " designs, ",
    length(exports), " exported objects, ", length(rd_files), " topics.\n", sep = "")
