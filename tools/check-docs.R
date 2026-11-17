#!/usr/bin/env Rscript
# Documentation guards, run in CI.
#
#  1. Every argument of an exported function is documented, and nothing is
#     documented that is not an argument.
#  2. Every [foo()] cross-reference in R/, README and NEWS resolves.
#  3. Nothing references an identifier the package no longer has.
#  4. Every design constructor appears in the overview topic, README and NEWS.
#  5. Every exported function is linked from the overview topic, and the
#     generated Rd is real Rd -- not markdown that was never processed.
#
# R CMD check catches undocumented objects but not documentation that has
# quietly gone stale, which is the failure mode this covers.
#
#   Rscript tools/check-docs.R              checks the installed package
#   Rscript tools/check-docs.R --self-test  checks rules 4 and 5 on
#                                           tools/fixtures/docs

# A name as a whole word: `draw` is not found in `drawn` or `draw_design`.
mentions <- function(text, name) {
  pattern <- paste0("(?<![\\w.])", gsub("(\\W)", "\\\\\\1", name, perl = TRUE),
                    "(?!\\w)")
  grepl(pattern, paste(text, collapse = " "), perl = TRUE)
}

# Plain text of a parsed Rd file, outside \usage and \examples. Code leaves
# are left out too, because `x[f()]` is valid R there.
rd_prose <- function(rd) {
  out <- character(0)
  walk <- function(x) {
    tag <- attr(x, "Rd_tag")
    if (!is.null(tag) && tag %in% c("\\usage", "\\examples")) return()
    if (is.list(x)) {
      for (el in x) walk(el)
    } else if (identical(tag, "TEXT")) {
      out <<- c(out, as.character(x))
    }
  }
  walk(rd)
  paste(out, collapse = "")
}

# A markdown link roxygen never processed: [foo()] left as text.
unprocessed_link <- function(rd) {
  grepl("\\[[a-zA-Z_.][a-zA-Z0-9_.]*\\(\\)\\]", rd_prose(rd))
}

# Topics a parsed Rd file links to: \link{f}, \link[=f]{...},
# \link[pkg:f]{...} and \link[pkg]{f}.
rd_links <- function(rd) {
  out <- character(0)
  walk <- function(x) {
    if (identical(attr(x, "Rd_tag"), "\\link")) {
      opt <- paste(unlist(attr(x, "Rd_option")), collapse = "")
      out <<- c(out, if (grepl("^=|:", opt)) {
        sub("^=|^[^:]*:", "", opt)
      } else {
        paste(unlist(x), collapse = "")
      })
    }
    if (is.list(x)) for (el in x) walk(el)
  }
  walk(rd)
  unique(out)
}

if ("--self-test" %in% commandArgs(trailingOnly = TRUE)) {
  failed <- character(0)
  expect <- function(ok, what) if (!ok) failed <<- c(failed, what)
  fx <- function(f) tools::parse_Rd(file.path("tools/fixtures/docs", f))
  expect(unprocessed_link(fx("unprocessed.Rd")),
         "a [foo()] link in \\details is found")
  expect(!unprocessed_link(fx("processed.Rd")),
         "[foo()] in \\usage, \\examples and \\code is not reported")
  hub <- fx("hub.Rd")
  expect(identical(sort(rd_links(hub)),
                   c("ht_mean", "ht_total", "median", "ratio")),
         "the targets of every form of \\link are read")
  expect(!"draw" %in% rd_links(hub),
         "the words drawn and draw_design do not link draw")
  prose <- rd_prose(hub)
  expect(!mentions(prose, "draw") && mentions(prose, "drawn"),
         "draw is not mentioned through drawn or draw_design")
  expect(mentions("see design_simple().", "design_simple") &&
           !mentions("design_simple_x", "design_simple") &&
           !mentions("x.design_simple", "design_simple"),
         "a name followed by ( counts, a longer or dotted name does not")
  if (length(failed)) {
    cat("Self-test failed:\n", paste0("  ", failed, "\n"), sep = "")
    quit(status = 1)
  }
  cat("Self-test passed.\n")
  quit(status = 0)
}

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
  documented <- setdiff(rd_args(f), "...")
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
  txt <- lines_of(f)
  miss <- ctors[!vapply(ctors, mentions, logical(1), text = txt)]
  if (length(miss)) note(basename(f), " omits ", paste(miss, collapse = ", "))
}

cat("5. Rd is processed, and every export is linked from the overview\n")
# Without `Roxygen: list(markdown = TRUE)` in DESCRIPTION every [foo()] and
# `code` span is written to Rd verbatim, so the help pages render as raw
# markdown with no working links. Cheap to check, silent to miss.
if (!any(grepl("markdown\\s*=\\s*TRUE", lines_of("DESCRIPTION")))) {
  note("DESCRIPTION: needs Roxygen: list(markdown = TRUE)")
}
for (f in rd_files) {
  if (unprocessed_link(tools::parse_Rd(f))) {
    note(basename(f), ": unprocessed markdown link -- see DESCRIPTION")
  }
}
hub <- rd_links(tools::parse_Rd("man/designs.Rd"))
for (e in setdiff(exports, c(grep("^design_", exports, value = TRUE),
                             "draw_design", "is_design"))) {
  if (!e %in% hub) {
    note("designs.Rd: no link to ", e, "() from the overview topic")
  }
}

if (bad) {
  cat("\n", bad, " documentation problem(s)\n", sep = "")
  quit(status = 1)
}
cat("\nDocumentation consistent: ", length(ctors), " designs, ",
    length(exports), " exported objects, ", length(rd_files), " topics.\n", sep = "")
