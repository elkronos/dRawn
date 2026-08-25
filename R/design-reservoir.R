#' Reservoir sampling
#'
#' Draws a uniform sample of fixed size from a stream of unknown length in a
#' single pass, using Algorithm L.
#'
#' A data frame is not a stream: its length is already known and it already fits
#' in memory, so [draw()] takes a direct vectorised path for one. Reach for the
#' streaming path when the data genuinely does not fit — pass a connection, or a
#' function that returns the next item and `NULL` at end of stream. For those,
#' [draw()] returns a list rather than a data frame.
#'
#' @param n Reservoir size. Fewer items come back if the stream is shorter; the
#'   result is never padded.
#' @param max_items Stop after reading this many items. Rows past it are never
#'   read, so [inclusion_prob()] gives them `0` and gives the rows before it
#'   `n / max_items` rather than `n / N`. `Inf` reads the whole
#'   stream. Warns when the cap actually truncates.
#'
#' @return A design object, for use with [draw()].
#'
#' @examples
#' nrow(draw(data.frame(id = 1:1000), design_reservoir(n = 10), seed = 1))
#'
#' # A real stream: a generator that yields 1..50 then stops
#' i <- 0
#' gen <- function() {
#'   i <<- i + 1
#'   if (i > 50) NULL else i
#' }
#' unlist(draw(gen, design_reservoir(n = 5), seed = 1))
#'
#' @family designs
#' @seealso [draw()]
#' @export
design_reservoir <- function(n, max_items = Inf) {
  if (!(identical(max_items, Inf) ||
        (is.numeric(max_items) && length(max_items) == 1L &&
         !is.na(max_items) && max_items >= 1 &&
         max_items == trunc(max_items)))) {
    stop("`max_items` must be Inf or a single whole number of 1 or more.",
         call. = FALSE)
  }
  new_design("reservoir", list(
    n = check_count(n, "n"),
    max_items = max_items
  ))
}

#' @export
draw_design.drawn_design_reservoir <- function(design, data) {
  n <- design$n
  max_items <- design$max_items

  if (is.data.frame(data)) {
    validate_data(data)
    n_rows <- nrow(data)
    if (is.finite(max_items) && max_items < n_rows) {
      warning("`max_items` (", max_items, ") is below the ", n_rows,
              " available rows; sampling only the first ", max_items, ".",
              call. = FALSE)
      n_rows <- as.integer(max_items)
    }
    idx <- if (n >= n_rows) seq_len(n_rows) else sort(sample.int(n_rows, n))
    return(reindex(data, idx))
  }

  next_item <- as_stream(data)
  con <- attr(next_item, "close_when_done")
  if (!is.null(con)) on.exit(close(con), add = TRUE)
  reservoir <- vector("list", n)
  count <- 0L
  truncated <- FALSE

  if (n > 0L) {
    # Algorithm L: once the reservoir fills, skip ahead geometrically rather
    # than drawing once per item.
    w <- exp(log(stats::runif(1)) / n)
    next_swap <- n + floor(log(stats::runif(1)) / log(1 - w)) + 1

    repeat {
      item <- next_item()
      if (is.null(item)) break
      # Only now is it known that the stream held more than `max_items`. Testing
      # the count before reading warned on a stream of exactly `max_items`,
      # where nothing was left behind.
      if (count >= max_items) {
        truncated <- TRUE
        break
      }
      count <- count + 1L

      if (count <= n) {
        reservoir[[count]] <- item
      } else if (count == next_swap) {
        reservoir[[sample.int(n, 1L)]] <- item
        w <- w * exp(log(stats::runif(1)) / n)
        next_swap <- next_swap + floor(log(stats::runif(1)) / log(1 - w)) + 1
      }
    }
  }

  if (truncated) {
    warning("Stopped after `max_items` = ", max_items, " item(s); the stream ",
            "may have held more.", call. = FALSE)
  }

  # Trim rather than pad with NULLs when the stream was short.
  reservoir[seq_len(min(count, n))]
}

#' Turn the supported inputs into a next-item function
#' @noRd
as_stream <- function(x) {
  if (is.function(x)) {
    if (length(formals(x)) > 0L) {
      stop("A stream function must take no arguments.", call. = FALSE)
    }
    return(x)
  }
  if (inherits(x, "connection")) {
    # An on.exit() here would close the connection when *this* function
    # returns, long before the reader below has read a byte from it. The
    # caller's own open connections are left alone; one opened here is closed
    # by draw_design(), which knows when the reading is done.
    opened_here <- !isOpen(x)
    if (opened_here) open(x, "r")
    reader <- function() {
      line <- readLines(x, n = 1L, warn = FALSE)
      if (length(line) == 0L) NULL else line
    }
    attr(reader, "close_when_done") <- if (opened_here) x else NULL
    return(reader)
  }
  if (is.list(x)) {
    i <- 0L
    n <- length(x)
    return(function() {
      i <<- i + 1L
      if (i > n) NULL else x[[i]]
    })
  }
  stop("`data` must be a data frame, list, connection, or zero-argument ",
       "function for a reservoir design, not a ", class(x)[1], ".",
       call. = FALSE)
}

# ---- inclusion probability ------------------------------------------------

#' @noRd
exact_inclusion.drawn_design_reservoir <- function(design, data) {
  # `max_items` stops the read, so rows past it are never seen -- and the rows
  # before it are competing against a shorter frame, not the whole of `data`.
  reach <- reservoir_reach(design, nrow(data))
  out <- numeric(nrow(data))
  if (reach > 0L) out[seq_len(reach)] <- min(1, design$n / reach)
  out
}

#' How many rows of a data frame a reservoir design can actually see
#' @noRd
reservoir_reach <- function(design, n_rows) {
  as.integer(min(n_rows, design$max_items))
}
