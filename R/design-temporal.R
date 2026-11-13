#' Temporal sampling
#'
#' Divides a time window into intervals of `interval` units and samples within
#' each one. `per_interval` is named for what it is: unlike `n` elsewhere in
#' the package, it is a per-group figure, not a total.
#'
#' Each row's interval is found in a single pass over the data, so the draw
#' and the inclusion probabilities take time that grows with the number of
#' rows rather than with rows times intervals.
#'
#' @section Time zones and calendar steps:
#' The window and the intervals are laid out in one time zone: `tz` when it
#' is given, otherwise the zone of a `POSIXct` `time` column, otherwise UTC.
#' A `POSIXct` column with an empty `tzone`, which R shows in the session's
#' zone, also gets UTC, so a design gives the same intervals on every
#' machine.
#'
#' A `Date`, in the column or as a bound, is its midnight in that zone, and a
#' string without an offset is a clock time there. A clock time that a
#' daylight-saving change skips moves to the moment the clocks jump, and one
#' that happens twice takes its first occurrence, so a `Date` stands for the
#' first instant of its day: 01:00 where the change skips midnight, as in
#' Santiago or Cairo, and the first of two midnights where it repeats one, as
#' in the Azores. A string ending in `Z` or an offset such as `+01:00` keeps
#' the instant it names.
#'
#' `"seconds"`, `"minutes"` and `"hours"` are fixed durations counted from
#' `from`, so hourly intervals ignore the clock change. `"days"` and
#' `"weeks"` step by calendar days, so each interval starts at the same
#' clock time as `from` and lasts 23 or 25 hours across a change.
#' `"months"` and `"years"` step by calendar units counted from `from`, and a
#' day the month lacks moves to its last day: a window starting on 31
#' January 2024 has intervals starting on 29 February, 31 March and 30 April.
#'
#' @param time Column of timestamps: `POSIXct`, `POSIXlt`, `Date`, character
#'   or factor. Strings are read as year-month-day with an optional time, in
#'   the forms [lubridate::ymd_hms()] accepts, such as `"2024-03-01 14:30:00"`
#'   or `"2024-03-01T14:30:00Z"`, and also `"2024-03-01 14:30"`,
#'   `"2024-03-01 14"` and `"2024-03-01"`.
#' @param from,to Window bounds, each a single `POSIXct`, `Date` or string
#'   read as `time` is. The window is half-open: `[from, to)`. Trailing
#'   partial intervals are dropped. Both are checked when the design is
#'   built, and `to` must come after `from`.
#' @param interval Interval width, in units of `unit`.
#' @param per_interval Rows to draw from each interval. Intervals holding fewer
#'   rows contribute all of them.
#' @param unit One of `"seconds"`, `"minutes"`, `"hours"`, `"days"`, `"weeks"`,
#'   `"months"` or `"years"`. See "Time zones and calendar steps".
#' @param tz Time zone, a name from [OlsonNames()], in which the bounds, any
#'   `Date` or string timestamps and the calendar steps are read. `NULL`, the
#'   default, takes the zone of a `POSIXct` `time` column and UTC otherwise.
#'   It shifts the window rather than relabelling the output.
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
  if (!is.null(tz)) check_zone(tz)
  time <- check_columns(time, "time", 1L)
  interval <- check_count(interval, "interval", allow_zero = FALSE)
  per_interval <- check_count(per_interval, "per_interval")
  unit <- check_choice(unit, c("hours", "seconds", "minutes", "days",
                               "weeks", "months", "years"), "unit")
  na_rm <- check_flag(na_rm, "na_rm")
  # The bounds are read in `tz`, or in UTC until a column supplies the zone.
  # Whether they parse does not depend on the zone; their order can, for a
  # pair that mixes an offset with a clock time, so the window is checked
  # again when it is laid out.
  zone <- tz %||% "UTC"
  check_window(time_bound(from, zone, "from")$at,
               time_bound(to, zone, "to")$at)
  # `from` and `to` are stored as given, so the design prints what was
  # written, and are read again in the zone in force when it is used.
  new_design("temporal", list(
    time = time,
    from = from,
    to = to,
    interval = interval,
    per_interval = per_interval,
    unit = unit,
    tz = tz,
    na_rm = na_rm
  ))
}

#' @export
draw_design.drawn_design_temporal <- function(design, data) {
  tm <- temporal_setup(design, data)
  keep <- which(!is.na(tm$bucket))
  if (!length(keep)) {
    return(empty_like(data))
  }
  idx_by_bucket <- split(keep, tm$bucket[keep])
  n_alloc <- pmin(lengths(idx_by_bucket), design$per_interval)

  reindex(data, take_within(idx_by_bucket, n_alloc, FALSE), sort = TRUE)
}

#' The window, its intervals and each row's interval
#'
#' The checks and the parsing the draw, the probabilities and the variance
#' share, so [inclusion_prob()] refuses what [draw()] refuses, with the same
#' message: a timestamp missing under `na_rm = FALSE` or in every row, a
#' bound that is not one parseable date-time, a window that runs backwards,
#' and a window shorter than one interval. Returns `bucket`, each row's
#' interval as an integer, `NA` for a row outside the window or with no
#' timestamp, `breaks` and `tz`.
#'
#' @noRd
temporal_setup <- function(design, data) {
  validate_data(data, required_columns = design$time)
  col <- design$time
  win <- temporal_window(design, data)

  tv <- parse_time(data[[col]], win$tz, col)
  what <- paste0("a missing or unparseable timestamp in `", col, "`")
  if (is.character(data[[col]]) || is.factor(data[[col]])) {
    what <- paste0(what, " (strings must be year-month-day with an optional ",
                   "time, such as \"2024-03-01\", \"2024-03-01 14:30\" or ",
                   "\"2024-03-01T14:30:00Z\")")
  }
  bad <- check_na_policy(is.na(tv), design$na_rm, what)
  if (all(bad)) stop("Every row has ", what, ".", call. = FALSE)

  breaks <- win$breaks
  bucket <- rep(NA_integer_, nrow(data))
  inw <- !bad & tv >= breaks[1] & tv < breaks[length(breaks)]
  bucket[inw] <- findInterval(tv[inw], breaks)
  list(bucket = bucket, breaks = breaks, tz = win$tz)
}

#' The time zone, bounds and interval boundaries of a temporal design
#'
#' The one place the zone is resolved: `tz` when the design has one, else the
#' zone of a `POSIXct` `time` column, else UTC. The bounds are read in that
#' zone and the breaks laid out from them.
#'
#' @noRd
temporal_window <- function(design, data) {
  tz <- temporal_tz(design, data[[design$time]])
  from <- time_bound(design$from, tz, "from")
  to <- time_bound(design$to, tz, "to")
  check_window(from$at, to$at)

  breaks <- make_breaks(from$at, to$at, design$interval, design$unit,
                        clock = from$clock)
  if (length(breaks) < 2L) {
    stop("The window (", format(difftime(to$at, from$at)),
         ") is shorter than one interval of ", design$interval, " ",
         design$unit, "; no complete interval exists.", call. = FALSE)
  }
  list(tz = tz, from = from$at, to = to$at, breaks = breaks)
}

#' The zone a design is laid out in
#'
#' `x` is the `time` column. A `POSIXct` column carries its zone in the
#' `tzone` attribute; an empty one means the session's zone, which differs
#' between machines, so it gets UTC like a column with no zone at all.
#'
#' @noRd
temporal_tz <- function(design, x) {
  if (!is.null(design$tz)) return(design$tz)
  zone <- if (inherits(x, "POSIXt")) attr(x, "tzone")[1L]
  if (is.null(zone) || is.na(zone) || !nzchar(zone)) return("UTC")
  if (!known_zone(zone)) {
    stop("`", design$time, "` is in time zone \"", zone, "\", which R does ",
         "not know. Set `tz` in design_temporal() to a name from ",
         "OlsonNames().", call. = FALSE)
  }
  zone
}

#' @noRd
known_zone <- function(tz) tz %in% c("UTC", "GMT") || tz %in% OlsonNames()

#' @noRd
check_zone <- function(tz) {
  if (!is.character(tz) || length(tz) != 1L || is.na(tz)) {
    stop("`tz` must be a single time zone name, or NULL.", call. = FALSE)
  }
  if (!known_zone(tz)) {
    stop("`tz` is \"", tz, "\", which is not a time zone R knows. Use a name ",
         "from OlsonNames(), such as \"Europe/London\" or \"UTC\".",
         call. = FALSE)
  }
  invisible(tz)
}

#' One window bound, as an instant and as a clock reading
#'
#' `at` is the instant in `tz`. `clock` is the local clock time the bound
#' was written as, held in UTC so that calendar steps from it are free of
#' daylight-saving changes. A `Date` or a string without an offset keeps the
#' clock time as written, which can differ from `at` when that time does not
#' exist in `tz`; anything else reads its clock from `at`.
#'
#' @noRd
time_bound <- function(x, tz, arg) {
  hint <- paste0("Give a POSIXct, a Date, or a string written year-month-day ",
                 "with an optional time, such as \"2024-03-01\" or ",
                 "\"2024-03-01 14:30\".")
  if (length(x) != 1L) {
    stop("`", arg, "` must be a single date-time, not ", length(x),
         " values. ", hint, call. = FALSE)
  }
  read <- if (!is.na(x)) read_time(x, tz, arg)
  if (is.null(read) || is.na(read$at)) {
    shown <- if (is.na(x)) {
      "NA"
    } else if (is.character(x) || is.factor(x)) {
      paste0("\"", as.character(x), "\"")
    } else {
      format(x)
    }
    stop("`", arg, "` is ", shown, ", which does not read as a date-time. ",
         hint, call. = FALSE)
  }
  clock <- read$clock
  if (is.null(clock) || is.na(clock)) {
    clock <- lubridate::force_tz(read$at, "UTC")
  }
  list(at = read$at, clock = clock)
}

#' Parse timestamps to instants in a time zone
#'
#' @noRd
parse_time <- function(x, tz, arg) read_time(x, tz, arg)$at

#' Read timestamps as instants, keeping the clock times they were written as
#'
#' Returns `at`, `POSIXct` instants in `tz`, and `clock`: for a `Date` or a
#' string without an offset, the clock time as written, in UTC; `NA` for a
#' string with its own offset; `NULL` for a date-time, which is an instant
#' already. A clock time is placed in `tz` by `zone_time()`.
#'
#' Base `as.POSIXct()` reads an ISO 8601 string such as
#' "2024-03-01T14:30:00Z" by its date alone, and drops the PM of
#' "2024-03-01 02:30:00 PM", so strings go through lubridate.
#'
#' @noRd
read_time <- function(x, tz, arg) {
  if (inherits(x, "POSIXt")) {
    # A date-time is an instant; only its display zone changes. Running
    # ymd_hms() over it would round-trip it through as.character(), and
    # midnight values would become NA.
    at <- as.POSIXct(x)
    attr(at, "tzone") <- tz
    return(list(at = at, clock = NULL))
  }
  if (inherits(x, "Date")) {
    clock <- .POSIXct(floor(unclass(x)) * 86400, tz = "UTC")
    return(list(at = zone_time(clock, tz), clock = clock))
  }
  if (is.character(x) || is.factor(x)) {
    read <- read_strings(as.character(x))
    at <- zone_time(read$clock, tz)
    if (any(read$own)) {
      at[read$own] <- read$clock[read$own]
      read$clock[read$own] <- NA
    }
    return(list(at = at, clock = read$clock))
  }
  stop("`", arg, "` must hold date-times, dates or parseable strings, not ",
       class(x)[1], ".", call. = FALSE)
}

#' Parse strings as clock times in UTC
#'
#' Returns `clock`, each string read in UTC, and `own`, `TRUE` for a string
#' that names its own instant with `Z` or an offset. Such a string reads as
#' the same instant in any zone, which is how it is told apart: it is read a
#' second time an hour east of UTC.
#'
#' Full date-times go through ymd_hms() and dates through ymd(). Whatever is
#' left is tried against "YYYY-mm-dd HH:MM" and "YYYY-mm-dd HH" exactly, with
#' a space or a `T` before the time. ymd_hms(truncated = 3) would also read
#' those, but it reads "01/02/2024" as 21 February 2001.
#'
#' @noRd
read_strings <- function(x) {
  clock <- suppressWarnings(lubridate::ymd_hms(x, tz = "UTC", quiet = TRUE))
  east <- suppressWarnings(
    lubridate::ymd_hms(x, tz = "Etc/GMT-1", quiet = TRUE))
  own <- !is.na(clock) & clock == east

  # Retry the date-only parser on whatever is still missing, rather than only
  # when nothing parsed at all. A column mixing "2024-01-01" with
  # "2024-01-01 05:00:00" otherwise loses every date-only row.
  todo <- is.na(clock) & !is.na(x)
  if (any(todo)) {
    clock[todo] <- suppressWarnings(
      lubridate::ymd(x[todo], tz = "UTC", quiet = TRUE))
  }
  todo <- is.na(clock) & !is.na(x)
  if (any(todo)) {
    clock[todo] <- suppressWarnings(lubridate::parse_date_time(
      x[todo], c("%Y-%m-%d %H:%M", "%Y-%m-%dT%H:%M", "%Y-%m-%d %H",
                 "%Y-%m-%dT%H"),
      tz = "UTC", exact = TRUE, quiet = TRUE))
  }
  list(clock = clock, own = own)
}

#' Place clock times, held in UTC, in a time zone
#'
#' A clock time that a daylight-saving change skips moves to the moment the
#' clocks jump, and one that happens twice takes its first occurrence. A
#' `Date` thus becomes the first instant of its day, whatever the other rows
#' hold, and a day that a zone skipped altogether, as Samoa skipped
#' 2011-12-30, becomes the first instant of the next one.
#'
#' @noRd
zone_time <- function(clock, tz) {
  if (identical(tz, "UTC")) {
    attr(clock, "tzone") <- tz
    return(clock)
  }
  lubridate::force_tz(clock, tzone = tz, roll_dst = c("boundary", "pre"))
}

#' Interval boundaries
#'
#' Seconds, minutes and hours are fixed durations from `from_dt`. Days and
#' weeks step by calendar days and months and years by calendar units, each
#' counted from `clock`, the clock time of `from_dt` held in UTC, and then
#' placed in the zone of `from_dt`. Stepping the clock keeps each break at
#' the same local time across a daylight-saving change, where a fixed 86,400
#' seconds would drift off midnight. Months and years are added to the
#' anchor rather than to the previous break, and a day the month lacks moves
#' to its last day, so a 31 January anchor gives 29 February, 31 March and
#' 30 April. lubridate::dmonths() would be a fixed 30.4375 days.
#'
#' Only complete intervals are used: a row counts only while
#' `tv < breaks[length(breaks)]`.
#'
#' @noRd
make_breaks <- function(from_dt, to_dt, interval, unit, clock = NULL) {
  if (unit %in% c("seconds", "minutes", "hours")) {
    step <- interval * c(seconds = 1, minutes = 60, hours = 3600)[[unit]]
    return(seq(from_dt, to_dt, by = step))
  }
  tz <- attr(from_dt, "tzone")[1L] %||% ""
  if (!nzchar(tz)) tz <- "UTC"
  clock <- clock %||% lubridate::force_tz(from_dt, "UTC")
  end <- lubridate::force_tz(to_dt, "UTC")

  steps <- if (unit %in% c("days", "weeks")) {
    days <- interval * c(days = 1, weeks = 7)[[unit]]
    span <- as.numeric(difftime(end, clock, units = "days"))
    clock + 86400 * days * seq_len(max(0, ceiling(span / days) + 1))
  } else {
    a <- as.POSIXlt(clock)
    b <- as.POSIXlt(end)
    gap <- 12 * (b$year - a$year) + (b$mon - a$mon)
    per <- interval * c(months = 1, years = 12)[[unit]]
    k <- seq_len(max(0, gap %/% per + 1))
    lubridate::`%m+%`(clock, lubridate::period(month = per * k))
  }
  at <- zone_time(steps, tz)
  # A day the zone skipped places two clock times at one instant; the empty
  # interval between them is dropped.
  unique(c(from_dt, at[at > from_dt & at <= to_dt]))
}

# ---- inclusion probability ------------------------------------------------

#' @noRd
exact_inclusion.drawn_design_temporal <- function(design, data) {
  temporal_pi(design, temporal_setup(design, data)$bucket)
}

#' Inclusion probabilities from each row's interval
#'
#' One count per interval, so the work grows with the number of rows rather
#' than with rows times intervals.
#'
#' @noRd
temporal_pi <- function(design, bucket) {
  out <- numeric(length(bucket))
  inw <- !is.na(bucket)
  if (!any(inw)) return(out)
  sizes <- tabulate(bucket[inw])
  out[inw] <- pmin(1, design$per_interval / sizes[bucket[inw]])
  out
}

#' Which sampling interval each row falls in
#'
#' `NA` for a row outside the window, which the design can never select. Shared
#' by the joint-probability and `survey` mappings so the three views of a
#' temporal design cannot drift apart.
#'
#' @noRd
temporal_bucket <- function(design, data) {
  as.character(temporal_setup(design, data)$bucket)
}

#' The window has to run forwards
#'
#' `make_breaks()` hands a backwards window to `seq()`, which fails with "wrong
#' sign in 'by' argument" -- a base-R message that says nothing about the design.
#'
#' @noRd
check_window <- function(from_dt, to_dt) {
  if (!is.na(from_dt) && !is.na(to_dt) && to_dt <= from_dt) {
    stop("`to` must be after `from`.", call. = FALSE)
  }
  invisible(NULL)
}
