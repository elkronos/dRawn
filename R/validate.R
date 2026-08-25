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
  as.integer(x)
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

#' @noRd
format_bad <- function(x) {
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
