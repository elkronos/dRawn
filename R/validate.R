# Argument validation. Every design constructor and draw method routes through
# these, so a given argument fails the same way and with the same wording
# wherever it appears.

#' Validate a count-like argument
#'
#' Passing these straight to [base::sample()] means `0`, `2.7`, `NA` and
#' length-2 vectors either fail cryptically or succeed silently.
#'
#' @noRd
check_count <- function(x, arg, allow_zero = TRUE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x != trunc(x)) {
    stop("`", arg, "` must be a single whole number, not ",
         format_bad(x, value = TRUE), ".", call. = FALSE)
  }
  if (x < 0 || (!allow_zero && x < 1)) {
    stop("`", arg, "` must be ",
         if (allow_zero) "non-negative" else "at least 1", ", not ", x, ".",
         call. = FALSE)
  }
  # as.integer() turns anything past 2^31 - 1 into NA with a warning, and the
  # design then prints `n NA` and fails somewhere far away.
  if (x > .Machine$integer.max) {
    stop("`", arg, "` must be at most ", .Machine$integer.max, ", not ", x, ".",
         call. = FALSE)
  }
  as.integer(x)
}

#' `n` given as a plan from `plan_size()`
#'
#' A design's `n` can be the plan itself, whose `n` is the number of rows to
#' draw. Only `n` arguments take a plan, since `check_count()` also checks
#' counts such as `n_clusters` and `reserve`, where a plan means nothing.
#'
#' @noRd
plan_n <- function(n) {
  if (inherits(n, "drawn_plan")) n[["n"]] else n
}

#' Validate a character vector of column names given at design time
#' @noRd
check_columns <- function(x, arg, max_len = Inf) {
  if (!is.character(x) || length(x) == 0L || anyNA(x) || any(!nzchar(x))) {
    stop("`", arg, "` must be one or more non-empty column names.",
         call. = FALSE)
  }
  if (length(x) > max_len) {
    stop("`", arg, "` must name ",
         if (max_len == 1L) "a single column" else paste(max_len, "columns"),
         ", not ", length(x), ".", call. = FALSE)
  }
  x
}

#' @noRd
check_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", arg, "` must be TRUE or FALSE.", call. = FALSE)
  }
  x
}

#' Match an argument against its allowed values
#'
#' The user-facing stand-in for [base::match.arg()], whose message reads
#' "'arg' should be one of" without naming the argument. `x` left at its
#' default, the whole of `choices`, or `NULL` gives the first choice, and a
#' unique abbreviation is completed, as with `match.arg()`. Code that passes
#' `NULL` on for an argument it was not given gets the default.
#'
#' @noRd
check_choice <- function(x, choices, arg) {
  if (is.null(x) || identical(x, choices)) return(choices[[1L]])
  i <- if (is.character(x) && length(x) == 1L && !is.na(x)) {
    pmatch(x, choices)
  } else {
    NA_integer_
  }
  if (is.na(i)) {
    shown <- paste0("\"", choices, "\"")
    stop("`", arg, "` must be one of ",
         paste(utils::head(shown, -1L), collapse = ", "), " or ",
         utils::tail(shown, 1L), ", not ",
         if (is.character(x) && length(x) == 1L && !is.na(x)) {
           paste0("\"", x, "\"")
         } else {
           format_bad(x)
         }, ".", call. = FALSE)
  }
  choices[[i]]
}

#' A plot title: one string, or NULL for the default
#'
#' `main` also takes an expression, a call such as `quote(alpha^2)` or a
#' name such as `quote(alpha)`, which `mtext()` draws as plotmath. `sub` is
#' added to the line under the title as text, so it takes a string only.
#'
#' @noRd
check_title <- function(x, arg) {
  string <- is.character(x) && length(x) == 1L && !is.na(x)
  math <- arg == "main" && is.language(x)
  if (is.null(x) || string || math) return(x)
  stop("`", arg, "` must be a single string",
       if (arg == "main") ", an expression", " or NULL, not ",
       format_bad(x), ".", call. = FALSE)
}

#' The plot palette, with the user's entries in place of the defaults
#'
#' Checked before anything is drawn: an unknown name, a NULL or a value that
#' is not one colour would otherwise stop the plot halfway through, or be
#' ignored. A named character vector is read as a list.
#'
#' @noRd
check_palette <- function(palette) {
  base <- drawn_pal()
  if (is.null(palette)) return(base)
  if (is.character(palette) && !is.null(names(palette))) {
    palette <- as.list(palette)
  }
  if (!is.list(palette) || is.data.frame(palette)) {
    stop("`palette` must be a named list of colours, such as ",
         "list(accent = \"darkred\"), not ",
         if (is.character(palette)) "an unnamed character vector"
         else a_class(palette), ".", call. = FALSE)
  }
  nm <- names(palette)
  unnamed <- if (is.null(nm)) {
    seq_along(palette)
  } else {
    which(is.na(nm) | !nzchar(nm))
  }
  if (length(unnamed)) {
    stop("`palette` ", if (length(unnamed) == 1L) "entry " else "entries ",
         paste(unnamed, collapse = ", "), " ",
         if (length(unnamed) == 1L) "has" else "have", " no name. Name ",
         "each entry after the colour it replaces, such as ",
         "list(accent = \"darkred\").", call. = FALSE)
  }
  unknown <- setdiff(nm, names(base))
  if (length(unknown)) {
    stop("`palette` has ", if (length(unknown) == 1L) "an entry" else
           "entries", " ", quote_names(unknown), ", which the plot does not ",
         "use. Name entries from ", quote_names(names(base)), ".",
         call. = FALSE)
  }
  for (k in nm) {
    v <- palette[[k]]
    ok <- (is.character(v) || is.numeric(v)) && length(v) == 1L &&
      !is.na(v) &&
      !inherits(tryCatch(grDevices::col2rgb(v), error = identity), "error")
    if (!ok) {
      shown <- if (is.null(v)) {
        "NULL"
      } else if (is.character(v) && length(v) == 1L && !is.na(v)) {
        paste0("\"", v, "\"")
      } else {
        format_bad(v)
      }
      stop("`palette$", k, "` must be one colour, a name from colors() or ",
           "a code such as \"#1B6B52\", not ", shown, ".", call. = FALSE)
    }
  }
  utils::modifyList(base, palette)
}

#' Refuse `...` arguments to plot() that are not graphical parameters
#'
#' `...` reaches the call that sets up the panel. A name it does not take,
#' such as `simulate`, was dropped with a warning from deep inside base
#' graphics. Read by name only, so nothing in `...` is evaluated before the
#' refusal.
#'
#' @noRd
check_panel_args <- function(...) {
  if (...length() == 0L) return(invisible(TRUE))
  nm <- ...names()
  if (is.null(nm)) nm <- rep("", ...length())
  nm[is.na(nm)] <- ""
  if (!all(nzchar(nm))) {
    stop("Arguments in `...` must be named graphical parameters, such as ",
         "`xlim = c(0, 10)`; ", sum(!nzchar(nm)), " ",
         if (sum(!nzchar(nm)) == 1L) "has" else "have", " no name.",
         call. = FALSE)
  }
  bad <- setdiff(nm, panel_parameters())
  if (length(bad)) {
    bad <- paste0("`", bad, "`")
    if (length(bad) > 1L) {
      bad <- paste(paste(bad[-length(bad)], collapse = ", "), "or",
                   bad[length(bad)])
    }
    stop("plot() for a design does not take ", bad, ". Its ",
         "`...` takes graphical parameters for the panel, such as `xlim`, ",
         "`ylim`, `asp` or `log`.", call. = FALSE)
  }
  invisible(TRUE)
}

#' The names plot() of a design passes on to the panel
#'
#' The arguments of `plot.default()` that set up a panel, and the graphical
#' parameters of `par()` that can be set. `panel.first` and `panel.last` are
#' left out, since `...` is evaluated before the panel exists.
#'
#' @noRd
panel_parameters <- function() {
  c("xlim", "ylim", "log", "xlab", "ylab", "ann", "axes", "frame.plot",
    "asp", "xgap.axis", "ygap.axis",
    "adj", "ask", "bg", "bty", "cex", "cex.axis", "cex.lab", "cex.main",
    "cex.sub", "col", "col.axis", "col.lab", "col.main", "col.sub", "crt",
    "err", "family", "fg", "fig", "fin", "font", "font.axis", "font.lab",
    "font.main", "font.sub", "lab", "las", "lend", "lheight", "ljoin",
    "lmitre", "lty", "lwd", "mai", "mar", "mex", "mfcol", "mfg", "mfrow",
    "mgp", "mkh", "new", "oma", "omd", "omi", "pch", "pin", "plt", "ps",
    "pty", "smo", "srt", "tck", "tcl", "usr", "xaxp", "xaxs", "xaxt",
    "xlog", "xpd", "yaxp", "yaxs", "yaxt", "ylbias", "ylog")
}

#' A value that failed a check, for the message
#'
#' Text is shown in backticks, as a name, by default. Where a number or a
#' flag was wanted, `value = TRUE` shows a single string as one, so a seed
#' given as "1" reads `the string "1"` rather than the `1` of a number.
#'
#' @noRd
format_bad <- function(x, value = FALSE) {
  if (value && is.character(x) && length(x) == 1L && !is.na(x)) {
    return(paste0("the string \"", x, "\""))
  }
  if (is.character(x)) {
    return(paste0("`", x, "`", collapse = ", "))
  }
  # A function, a name or an environment has no value to show, and is.na()
  # warns on it. is.na() of a data frame gives one answer per cell.
  if (is.data.frame(x) || !(is.null(x) || is.atomic(x) || is.list(x))) {
    return(a_class(x))
  }
  if (length(x) != 1L) {
    return(paste0("a length-", length(x), " ", class(x)[1]))
  }
  if (is.na(x)) return("NA")
  if (!is.numeric(x)) return(a_class(x))
  paste0(x)
}

#' A word with its indefinite article
#'
#' "an" before a vowel, so a class reads "an integer" and "a list".
#'
#' @noRd
with_article <- function(word) {
  paste(if (grepl("^[aeiouAEIOU]", word)) "an" else "a", word)
}

#' The class of `x`, with its article, for messages
#' @noRd
a_class <- function(x) with_article(class(x)[1])

#' Require a design, and say so when it and the data are swapped
#'
#' `data |> draw(design)` written the other way round, `design |> draw(data)`,
#' puts the design where the data go.
#'
#' @noRd
check_design <- function(design, data, fun) {
  if (is_design(data) && !is_design(design)) {
    stop("`data` is a sampling design and `design` is ", a_class(design),
         ", so the two look swapped. Call ", fun, "(data, design), with the ",
         "data first.", call. = FALSE)
  }
  if (!is_design(design)) {
    stop("`design` must come from one of the design_*() constructors, not ",
         a_class(design), ". See ?designs.", call. = FALSE)
  }
  invisible(design)
}

#' Validate a data frame and, optionally, that it holds the given columns
#' @noRd
validate_data <- function(data, required_columns = NULL, arg = "data") {
  if (!is.data.frame(data)) {
    stop("`", arg, "` must be a data frame, not ", a_class(data), ".",
         call. = FALSE)
  }
  if (nrow(data) == 0L) {
    stop("`", arg, "` has no rows.", call. = FALSE)
  }
  if (!is.null(required_columns)) {
    missing_cols <- setdiff(required_columns, names(data))
    if (length(missing_cols) > 0L) {
      stop("`", arg, "` is missing ",
           if (length(missing_cols) == 1L) "column: " else "columns: ",
           paste0("`", missing_cols, "`", collapse = ", "), ".", call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' Refuse to draw more rows than exist, without replacement
#'
#' Shared so the message is identical across every design that needs it.
#'
#' @noRd
check_draw_size <- function(n, n_rows, replace, arg = "n") {
  if (!isTRUE(replace) && n > n_rows) {
    stop("`", arg, "` (", n, ") cannot exceed the number of rows (", n_rows,
         ") when sampling without replacement.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Drop rows flagged by `bad`, or stop if `na_rm` is FALSE
#'
#' Every design that reads a key column routes through here, so `na_rm` means
#' the same thing everywhere and the error names the offending column.
#'
#' @noRd
drop_na_rows <- function(data, bad, na_rm, what) {
  bad[is.na(bad)] <- TRUE
  if (!any(bad)) {
    return(data)
  }
  if (!isTRUE(na_rm)) {
    stop(sum(bad), " row(s) have ", what, ". Set na_rm = TRUE to drop them.",
         call. = FALSE)
  }
  out <- reindex(data, which(!bad))
  if (nrow(out) == 0L) {
    stop("Every row has ", what, ".", call. = FALSE)
  }
  out
}

#' Require that a column can serve as a grouping or ordering key
#'
#' A data frame column can itself be a list or a matrix. Neither groups or
#' orders the way a key must: `unique()` on a matrix column works element-wise,
#' so cluster sampling silently returned one row from a twenty-row frame, and
#' `order()` on a list column fails inside base R with "unimplemented type
#' 'list' in 'orderVector1'".
#'
#' @noRd
check_key_columns <- function(data, columns, arg) {
  for (nm in columns) {
    col <- data[[nm]]
    what <- if (!is.null(dim(col))) {
      paste0("a ", length(dim(col)), "-dimensional ", class(col)[1], " column")
    } else if (is.list(col) && !is.factor(col)) {
      "a list column"
    } else {
      next
    }
    stop("`", arg, "` names `", nm, "`, which is ", what,
         ". Grouping and ordering need a plain vector or factor; ",
         "flatten it first.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Apply the design's missing-data policy outside of `draw()`
#'
#' `drop_na_rows()` both raises the error and does the dropping, which suits
#' `draw_design()`. The inclusion-probability methods need only the error: they
#' report 0 for a row the design cannot reach, so they have nothing to drop.
#' Without this they answer for a design that `draw()` would refuse, and `na_rm`
#' stops meaning the same thing everywhere.
#'
#' @noRd
check_na_policy <- function(bad, na_rm, what) {
  bad[is.na(bad)] <- TRUE
  if (any(bad) && !isTRUE(na_rm)) {
    stop(sum(bad), " row(s) have ", what, ". Set na_rm = TRUE to drop them.",
         call. = FALSE)
  }
  invisible(bad)
}

#' A grouping key built from integer codes
#'
#' Each column is coded on its own: a factor by its levels, anything else by
#' its distinct values in byte order (`method = "radix"`), which is the same
#' in every locale. The tuples of codes that occur are numbered in
#' lexicographic order, first column slowest. The work grows with the number
#' of rows, however many combinations the columns' levels could form, and no
#' value can make two tuples share a key. A row with a missing value in any
#' column gets `NA`.
#'
#' The levels are the values joined by a carriage return, which
#' `group_label()` shows as ` . `. When that would give two tuples the same
#' level, as a value holding a carriage return can, the levels are the code
#' tuples instead. The codes depend on the data keyed, so keys that are
#' compared across two frames come from `group_key2()`.
#'
#' @noRd
group_key <- function(data, cols) {
  key_from_columns(lapply(cols, function(col) data[[col]]))
}

#' Keys for two frames on one set of codes
#'
#' A group has the same code in `x` and in `table`, whatever else either
#' holds, so the two can be compared by position.
#'
#' @noRd
group_key2 <- function(x, table, cols) {
  n_x <- nrow(x)
  key <- key_from_columns(lapply(cols, function(col) {
    stack_values(x[[col]], table[[col]])
  }))
  list(x = key[seq_len(n_x)], table = key[n_x + seq_len(nrow(table))])
}

#' @noRd
key_from_columns <- function(columns) {
  coded <- lapply(columns, column_codes)
  codes <- lapply(coded, `[[`, "code")
  n <- length(codes[[1L]])
  ok <- which(!Reduce(`|`, lapply(codes, is.na)))
  id <- rep(NA_integer_, n)
  sub <- lapply(codes, `[`, ok)
  o <- do.call(order, c(unname(sub), list(method = "radix")))
  sorted <- lapply(sub, `[`, o)
  m <- length(o)
  starts <- if (m) {
    c(TRUE, Reduce(`|`, lapply(sorted, function(s) s[-1L] != s[-m])))
  } else {
    logical(0)
  }
  id[ok[o]] <- cumsum(starts)
  tuples <- lapply(sorted, `[`, which(starts))
  lev <- do.call(paste, c(Map(function(cc, tt) cc$labels[tt], coded, tuples),
                          sep = "\r"))
  if (anyDuplicated(lev)) {
    lev <- do.call(paste, c(unname(tuples), sep = "\r"))
  }
  structure(id, levels = as.character(lev), class = "factor")
}

#' Integer codes and labels for one key column
#'
#' Classed vectors such as dates are coded on their underlying values, since
#' `match()` would compare them as text. A double whose text repeats another's
#' (`0.3` and `0.1 + 0.2`) is labelled with all its digits.
#'
#' @noRd
column_codes <- function(x) {
  if (is.factor(x)) {
    return(list(code = as.integer(x), labels = levels(x)))
  }
  v <- key_text(if (is.object(x) && is.atomic(x)) unclass(x) else x)
  u <- sort(unique(v), method = sort_method(v))
  code <- match(v, u)
  labels <- as.character(x[match(u, v)])
  if (is.double(v) && anyDuplicated(labels)) labels <- sprintf("%.17g", u)
  list(code = code, labels = labels)
}

#' Order a key column the same way in every locale
#'
#' `method = "radix"` sorts text by byte, as the C locale does, where the
#' default method follows the session's collation. It does not take complex or
#' raw vectors, which have no collation to follow.
#'
#' @noRd
key_order <- function(x) order(key_text(x), method = sort_method(x))

#' The distinct values of a key column in key order
#'
#' The values themselves, with their encodings, so they can be matched back
#' to the column.
#'
#' @noRd
key_values <- function(x) {
  u <- unique(x[!is.na(x)])
  u[key_order(u)]
}

#' Text that a radix sort takes, with the same bytes
#'
#' `method = "radix"` refuses a non-ASCII string of unknown encoding, which
#' is what `read.csv()`, `readLines()` and `data.table::fread()` give unless
#' told the file's encoding. In the copy returned, such a string is marked
#' UTF-8 when its bytes are valid UTF-8 and "bytes" otherwise. Its bytes are
#' unchanged, so text sorts byte by byte, in the same order in a UTF-8
#' session as in the C locale. `enc2utf8()` is not used because in the C
#' locale it turns each non-ASCII byte into a `<xx>` escape, which sorts
#' differently. In the copy, an unmarked string and a UTF-8 string with the
#' same bytes are one value in any session. Callers take labels from the
#' original values.
#'
#' @noRd
key_text <- function(x) {
  if (!is.character(x)) return(x)
  odd <- which(Encoding(x) == "unknown" &
                 grepl("[^\001-\177]", x, useBytes = TRUE))
  if (length(odd)) {
    v <- x[odd]
    Encoding(v) <- ifelse(validUTF8(v), "UTF-8", "bytes")
    x[odd] <- v
  }
  x
}

#' @noRd
sort_method <- function(x) {
  if (is.factor(x) ||
      typeof(x) %in% c("logical", "integer", "double", "character")) {
    "radix"
  } else {
    "auto"
  }
}

#' One vector holding the values of two key columns
#'
#' Factors keep their level order. Plain vectors combine as `c()` combines
#' them, so an integer column matches a double one. Other columns of different
#' kinds are compared as text.
#'
#' @noRd
stack_values <- function(a, b) {
  if (is.factor(a) && is.factor(b)) {
    return(factor(c(as.character(a), as.character(b)),
                  levels = union(levels(a), levels(b))))
  }
  plain <- !is.object(a) && !is.object(b)
  if (!is.factor(a) && !is.factor(b) &&
      (plain || identical(class(a), class(b)))) {
    return(c(a, b))
  }
  c(as.character(a), as.character(b))
}

#' A key's levels for display
#'
#' Columns are separated by ` . ` and an empty value is shown as `(blank)`.
#' Two groups can share a display label; `unique_labels()` gives labels that
#' cannot.
#'
#' @noRd
group_label <- function(x) {
  x <- gsub("(?<![^\r])(?![^\r])", "(blank)", as.character(x), perl = TRUE)
  gsub("\r", " . ", x, fixed = TRUE)
}

#' Display labels for a key's levels, or their codes where two would match
#'
#' For output read as identifiers, such as the strata handed to `survey`,
#' which would merge two strata given one label.
#'
#' @noRd
unique_labels <- function(key) {
  lab <- group_label(levels(key))
  if (anyDuplicated(lab)) lab <- as.character(seq_along(lab))
  lab
}
