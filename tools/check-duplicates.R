#!/usr/bin/env Rscript
# Guard against a function being defined in more than one file under R/.
# R sources R/ alphabetically, so a stale copy left behind by a refactor
# silently wins over the intended one -- which is how a fixed allocate() went
# on using its unfixed twin.
defs <- do.call(rbind, lapply(list.files("R", full.names = TRUE), function(f) {
  nm <- sub(" *<- *function.*", "", grep("^[a-zA-Z_.][a-zA-Z0-9_.]* *<- *function",
                                         readLines(f, warn = FALSE), value = TRUE))
  if (length(nm) == 0) NULL else data.frame(fn = trimws(nm), file = basename(f))
}))
dupes <- defs[defs$fn %in% defs$fn[duplicated(defs$fn)], ]
if (nrow(dupes)) {
  cat("Duplicate definitions:\n"); print(dupes[order(dupes$fn), ], row.names = FALSE)
  quit(status = 1)
}
cat("No duplicate definitions across", length(unique(defs$file)), "files.\n")
