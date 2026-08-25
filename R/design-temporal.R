#' Temporal sampling
#'
#' Divides a time window into equal intervals and samples within each one.
#' `per_interval` is named for what it is: unlike `n` elsewhere in the package,
#' it is a per-group figure, not a total.
#'
#' Buckets are assigned arithmetically in a single pass rather than by
#' re-filtering the data per interval, so cost scales with the data rather than
#' with the number of intervals.
#'
#' @param time Column of timestamps. `POSIXct` and `Date` columns are used as
#'   they are; character columns are parsed with [lubridate::ymd_hms()].
#' @param from,to Window bounds, as timestamps or strings. The window is
#'   half-open: `[from, to)`. Trailing partial intervals are dropped.
#' @param interval Interval width, in units of `unit`.
#' @param per_interval Rows to draw from each interval. Intervals holding fewer
#'   rows contribute all of them.
#' @param unit One of `"seconds"`, `"minutes"`, `"hours"`, `"days"`, `"weeks"`,
#'   `"months"` or `"years"`. `"months"` and `"years"` step by calendar units,
#'   not by fixed 30.44-day durations.
#' @param tz Time zone the timestamps are expressed in. Applied when the column
#'   *and* both bounds are parsed, so it genuinely shifts the window rather than
#'   relabelling the output.
#' @param na_rm Drop rows whose timestamp is missing or unparseable instead of
#'   raising an error.
#'
#' @return A design object, for use with [draw()].
#'
#' @examples
#' df <- data.frame(
#'   id = 1:48,
#'   ts = seq(as.POSIXct("2020-01-01", tz = "UTC"), by = "hour", length.out = 48)
#' )
#' d <- design_temporal("ts", from = "2020-01-01", to = "2020-01-02",
#'                      interval = 6, per_interval = 2, unit = "hours")
#' draw(df, d, seed = 1)$ts
#'
#' @family designs
#' @seealso [draw()]
#' @export
design_temporal <- function(time, from, to, interval, per_interval,
                            unit = c("hours", "seconds", "minutes", "days",
                                     "weeks", "months", "years"),
                            tz = NULL, na_rm = FALSE) {
  if (!is.null(tz) && (!is.character(tz) || length(tz) != 1L)) {
    stop("`tz` must be a single time zone name, or NULL.", call. = FALSE)
  }
  new_design("temporal", list(
    time = check_columns(time, "time", 1L),
    from = from,
    to = to,
    interval = check_count(interval, "interval", allow_zero = FALSE),
    per_interval = check_count(per_interval, "per_interval"),
    unit = match.arg(unit),
    tz = tz,
    na_rm = check_flag(na_rm, "na_rm")
  ))
}

#' @export
draw_design.drawn_design_temporal <- function(design, data) {
  validate_data(data, required_columns = design$time)
  col <- design$time
  tz <- design$tz %||% "UTC"

  tv <- parse_time(data[[col]], tz, col)
  if (anyNA(tv)) {
    data <- drop_na_rows(
      data, is.na(tv), design$na_rm,
      paste0("a missing or unparseable timestamp in `", col, "`")
    )
    tv <- parse_time(data[[col]], tz, col)
  }

  from_dt <- parse_time(design$from, tz, "from")
  to_dt <- parse_time(design$to, tz, "to")
  if (length(from_dt) != 1L || length(to_dt) != 1L ||
      is.na(from_dt) || is.na(to_dt)) {
    stop("`from` and `to` must each be a single parseable date-time.",
         call. = FALSE)
  }
  if (to_dt <= from_dt) {
    stop("`to` must be after `from`.", call. = FALSE)
  }

  breaks <- make_breaks(from_dt, to_dt, design$interval, design$unit)
  if (length(breaks) < 2L) {
    stop("The window (", format(difftime(to_dt, from_dt)),
         ") is shorter than one interval of ", design$interval, " ",
         design$unit, "; no complete interval exists.", call. = FALSE)
  }

  in_window <- tv >= breaks[1] & tv < breaks[length(breaks)]
  if (!any(in_window)) {
    return(empty_like(data))
  }

  keep <- which(in_window)
  bucket <- findInterval(tv[keep], breaks)
  idx_by_bucket <- split(keep, bucket)
  n_alloc <- pmin(lengths(idx_by_bucket), design$per_interval)

  reindex(data, take_within(idx_by_bucket, n_alloc, FALSE), sort = TRUE)
}

#' Parse a value to POSIXct, leaving date-time classes alone
#'
#' Running ymd_hms() over an already-parsed column round-trips it through
#' as.character(), and midnight values silently become NA.
#'
#' @noRd
parse_time <- function(x, tz, arg) {
  if (inherits(x, "POSIXct")) {
    return(lubridate::with_tz(x, tzone = tz))
  }
  if (inherits(x, "POSIXlt")) {
    return(lubridate::with_tz(as.POSIXct(x), tzone = tz))
  }
  if (inherits(x, "Date")) {
    return(as.POSIXct(format(x, "%Y-%m-%d 00:00:00"), tz = tz))
  }
  if (is.character(x) || is.factor(x)) {
    x <- as.character(x)
    parsed <- suppressWarnings(lubridate::ymd_hms(x, tz = tz, quiet = TRUE))
    if (all(is.na(parsed) | is.na(x))) {
      parsed <- suppressWarnings(lubridate::ymd(x, tz = tz, quiet = TRUE))
    }
    return(parsed)
  }
  stop("`", arg, "` must hold date-times, dates or parseable strings, not ",
       class(x)[1], ".", call. = FALSE)
}

#' Interval boundaries, calendar-aware for months and years
#'
#' lubridate::dmonths() is a fixed 30.4375-day duration, so "monthly" buckets
#' drift off month boundaries immediately.
#'
#' @noRd
make_breaks <- function(from_dt, to_dt, interval, unit) {
  breaks <- if (unit %in% c("months", "years")) {
    seq(from_dt, to_dt, by = paste(interval, sub("s$", "", unit)))
  } else {
    step <- interval * switch(unit,
      seconds = 1, minutes = 60, hours = 3600,
      days = 86400, weeks = 604800
    )
    seq(from_dt, to_dt, by = step)
  }
  # Keep complete intervals only.
  breaks[breaks <= to_dt]
}

# ---- inclusion probability ------------------------------------------------

#' @noRd
exact_inclusion.drawn_design_temporal <- function(design, data) {
  validate_data(data, required_columns = design$time)
  tz <- design$tz %||% "UTC"
  tv <- parse_time(data[[design$time]], tz, design$time)
  from_dt <- parse_time(design$from, tz, "from")
  to_dt <- parse_time(design$to, tz, "to")
  breaks <- make_breaks(from_dt, to_dt, design$interval, design$unit)

  out <- numeric(nrow(data))
  if (length(breaks) < 2L) return(out)

  in_window <- !is.na(tv) & tv >= breaks[1] & tv < breaks[length(breaks)]
  if (!any(in_window)) return(out)

  keep <- which(in_window)
  bucket <- findInterval(tv[keep], breaks)
  for (b in unique(bucket)) {
    idx <- keep[bucket == b]
    out[idx] <- min(1, design$per_interval / length(idx))
  }
  out
}
