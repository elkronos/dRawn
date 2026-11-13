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
    stop("`", arg, "` must be a single whole number, not ", format_bad(x), ".",
         call. = FALSE)
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
#' default, the whole of `choices`, gives the first choice, and a unique
#' abbreviation is completed, as with `match.arg()`.
#'
#' @noRd
check_choice <- function(x, choices, arg) {
  if (identical(x, choices)) return(choices[[1L]])
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

#' @noRd
format_bad <- function(x) {
  if (is.character(x)) {
    return(paste0("`", x, "`", collapse = ", "))
  }
  if (length(x) != 1L) {
    return(paste0("a length-", length(x), " ", class(x)[1]))
  }
  if (is.na(x)) return("NA")
  if (!is.numeric(x)) return(paste0("a ", class(x)[1]))
  paste0(x)
}

#' Validate a data frame and, optionally, that it holds the given columns
#' @noRd
validate_data <- function(data, required_columns = NULL, arg = "data") {
  if (!is.data.frame(data)) {
    stop("`", arg, "` must be a data frame, not a ", class(data)[1], ".",
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
  out <- data[!bad, , drop = FALSE]
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
  v <- if (is.object(x) && is.atomic(x)) unclass(x) else x
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
key_order <- function(x) order(x, method = sort_method(x))

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
