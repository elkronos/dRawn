#' Reservoir sampling
#'
#' Draws a uniform sample of fixed size from a stream of unknown length in a
#' single pass, using Algorithm L (Li 1994), which skips ahead geometrically
#' instead of drawing a random number for every item (Vitter 1985).
#'
#' A data frame or a list is not a stream: its length is already known and it
#' already fits in memory, so [draw()] picks `n` positions in it at random and
#' returns them in their original order. Every element of a list is an item,
#' `NULL` included, and the elements keep their names. The streaming path is
#' for data that does not fit in memory. Pass it a connection, or a function
#' that returns the next item and `NULL` at end of stream. It holds at most `n`
#' items at a time, plus the block it is reading, and returns them in the order
#' they arrived. For a list, a connection or a function, [draw()] returns a
#' list rather than a data frame, except with `batches = TRUE`.
#'
#' The two paths use different algorithms, so for one seed a stream gives a
#' different sample from a data frame or list holding the same items. Each
#' path on its own gives the same sample for the same seed.
#'
#' @section Items in a stream:
#' Each line of a connection is one item, a character string as [readLines()]
#' returns it, and a header line is an item like any other. To leave a CSV
#' header out of the sample, open the connection, read the header with
#' `readLines(con, n = 1)` and then draw from the connection, which is read
#' from where the header left off. `read.csv(text = c(header, unlist(lines)))`
#' turns the sampled lines into a data frame (see the examples). [draw()]
#' reads a connection in blocks of 10,000 lines, which gives the same sample
#' for a seed as reading it line by line.
#'
#' A function is called once for each item, and whatever it returns is that
#' item, a data frame included, so a function that returns blocks of rows
#' gives a sample of whole blocks. With `batches = TRUE` the rows are the
#' items instead: each call returns a data frame of rows, every batch with
#' the same columns in the same order, and [draw()] returns `n` of those rows
#' as one data frame, bound with [rbind()] in the order they arrived. Cutting
#' the same rows into batches of other sizes gives the same sample for a
#' seed. Between batches it keeps at most `max(2 * n, 1024)` rows, some of
#' them already dropped from the sample, besides the batch being read.
#'
#' @section What is read:
#' A connection that is not open is opened, read to its end and closed, as
#' [read.csv()] does, so it cannot be drawn from twice. A connection that is
#' already open is read from where it stands and left open.
#'
#' With a finite `max_items`, [draw()] reads one item past the cap, because
#' only that read tells a stream of exactly `max_items` items from a longer
#' one. On an open connection in text mode the line goes back with
#' [pushBack()], so the next read starts with it. A connection opened in
#' binary mode (`"rb"`) cannot take a line back, and loses it. A function is
#' called once past the cap and what it returns is dropped; with
#' `batches = TRUE`, so are the rows of a batch that run past the cap.
#'
#' @param n Reservoir size, or a plan from [plan_size()]. Fewer items come
#'   back if the stream is shorter; the result is never padded.
#' @param max_items Stop after this many items, or rows with
#'   `batches = TRUE`. Items past the cap are never sampled, so
#'   [inclusion_prob()] gives rows past it `0` and gives the rows before it
#'   `n / max_items` rather than `n / N`. A stream is read one item past the
#'   cap; see "What is read". `Inf` reads the whole stream. Warns when the cap
#'   truncates.
#' @param batches If `TRUE`, `data` is a function that returns a data frame of
#'   rows on each call and `NULL` at the end, and the reservoir samples rows
#'   rather than whole data frames. [draw()] then returns a data frame. A data
#'   frame given as `data` is sampled by rows either way; a list or a
#'   connection is refused.
#'
#' @return A design object, for use with [draw()].
#'
#' @references
#' Vitter, J. S. (1985). Random sampling with a reservoir. *ACM Transactions
#' on Mathematical Software*, 11, 37–57.
#'
#' Li, K.-H. (1994). Reservoir-sampling algorithms of time complexity
#' O(n(1 + log(N/n))). *ACM Transactions on Mathematical Software*, 20,
#' 481–493.
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
#' # A CSV file read line by line, with its header left out of the sample
#' path <- tempfile(fileext = ".csv")
#' write.csv(data.frame(id = 1:200, value = (1:200) %% 7), path,
#'           row.names = FALSE)
#' con <- file(path, "r")
#' header <- readLines(con, n = 1)
#' lines <- draw(con, design_reservoir(n = 5), seed = 1)
#' close(con)
#' read.csv(text = c(header, unlist(lines)))
#'
#' # The same file parsed 50 rows at a time, sampling rows: the same five
#' # rows as above for the same seed
#' con <- file(path, "r")
#' header <- readLines(con, n = 1)
#' next_rows <- function() {
#'   lines <- readLines(con, n = 50)
#'   if (length(lines) == 0) NULL else read.csv(text = c(header, lines))
#' }
#' draw(next_rows, design_reservoir(n = 5, batches = TRUE), seed = 1)
#' close(con)
#' unlink(path)
#'
#' @family designs
#' @seealso [draw()]
#' @export
design_reservoir <- function(n, max_items = Inf, batches = FALSE) {
  if (!(identical(max_items, Inf) ||
        (is.numeric(max_items) && length(max_items) == 1L &&
         !is.na(max_items) && max_items >= 1 &&
         max_items == trunc(max_items)))) {
    stop("`max_items` must be Inf or a single whole number of 1 or more.",
         call. = FALSE)
  }
  new_design("reservoir", list(
    n = check_count(plan_n(n), "n"),
    max_items = max_items,
    batches = check_flag(batches, "batches")
  ))
}

#' @export
draw_design.drawn_design_reservoir <- function(design, data) {
  n <- design$n
  max_items <- design$max_items
  # A design saved before `batches` existed has none, and reads items.
  batches <- isTRUE(design$batches)

  if (is.data.frame(data)) {
    validate_data(data)
    n_rows <- nrow(data)
    if (is.finite(max_items) && max_items < n_rows) {
      warn_drawn(paste0("`max_items` (", max_items, ") is below the ",
                        n_rows, " available rows; sampling only the first ",
                        max_items, "."), "drawn_stream_truncated")
      n_rows <- as.integer(max_items)
    }
    idx <- if (n >= n_rows) seq_len(n_rows) else sort(sample.int(n_rows, n))
    return(reindex(data, idx))
  }

  # A list's length is known, so its positions are sampled as a data frame's
  # rows are. Reading it as a stream would end it at its first NULL element.
  if (is.list(data)) {
    if (batches) {
      stop("`batches = TRUE` needs a function that returns data frames, and ",
           "`data` is a list. Bind its data frames with do.call(rbind, data) ",
           "and draw from the result, or leave `batches` FALSE to sample ",
           "whole elements.", call. = FALSE)
    }
    n_items <- length(data)
    if (is.finite(max_items) && max_items < n_items) {
      truncation_warning(max_items)
      n_items <- as.integer(max_items)
    }
    idx <- if (n >= n_items) seq_len(n_items) else sort(sample.int(n_items, n))
    return(data[idx])
  }

  stream <- as_stream(data, batches)
  if (!is.null(stream$close)) on.exit(close(stream$close), add = TRUE)
  out <- if (is.null(stream$read)) {
    reservoir_stream(data, n, max_items)
  } else {
    reservoir_blocks(stream$read, n, max_items, rows = batches,
                     give_back = stream$give_back)
  }
  if (out$truncated) truncation_warning(max_items)
  out$items
}

#' Algorithm L over a next-item function
#'
#' Returns `items`, the sampled items in the order they arrived, `count`, the
#' number of items read, and `truncated`, whether the stream held more than
#' `max_items`. `count` is a double, which counts exactly past
#' `.Machine$integer.max`. The reservoir starts with at most 1,024 slots and
#' doubles until it holds `n`, so a large `n` costs memory only for the items
#' a stream supplies.
#'
#' A function is called once per item, so this loop handles one item at a
#' time; [reservoir_blocks()] makes the same draws for a stream read in
#' blocks.
#'
#' @noRd
reservoir_stream <- function(next_item, n, max_items) {
  slots <- min(n, 1024)
  reservoir <- vector("list", slots)
  arrived <- numeric(slots)
  count <- 0
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
      count <- count + 1

      if (count <= n) {
        if (count > slots) {
          slots <- min(n, 2 * slots)
          length(reservoir) <- slots
          length(arrived) <- slots
        }
        reservoir[[count]] <- item
        arrived[count] <- count
      } else if (count == next_swap) {
        slot <- sample.int(n, 1L)
        reservoir[[slot]] <- item
        arrived[slot] <- count
        w <- w * exp(log(stats::runif(1)) / n)
        next_swap <- next_swap + floor(log(stats::runif(1)) / log(1 - w)) + 1
      }
    }
  }

  # Trim rather than pad with NULLs when the stream was short.
  held <- seq_len(min(count, n))
  list(items = reservoir[held][order(arrived[held])], count = count,
       truncated = truncated)
}

#' Algorithm L over a stream read in blocks
#'
#' `read(want)` returns the next block of items, or `NULL` at the end of the
#' stream. A block is a list or a character vector of items, or with
#' `rows = TRUE` a data frame whose rows are the items. `want` is the most
#' items worth reading, one past `max_items`; a reader may return fewer or,
#' for a function's batches, more. `give_back(block, k)` receives a block
#' that ran past `max_items` and the number of its items that were within
#' it.
#'
#' The swaps Algorithm L makes are worked out block by block, as positions
#' within the block, and draw the same random numbers in the same order as
#' reading one item at a time would. A seed therefore gives the same sample
#' however the stream is cut into blocks.
#'
#' Returns `items`, the sampled items in the order they arrived (a list, or a
#' data frame with `rows = TRUE`), `count`, the number of items read, and
#' `truncated`, whether the stream held more than `max_items`. `count` is a
#' double, which counts exactly past `.Machine$integer.max`.
#'
#' @noRd
reservoir_blocks <- function(read, n, max_items, rows = FALSE,
                             give_back = NULL) {
  held <- if (rows) row_reservoir(n) else item_reservoir(n)
  count <- 0
  truncated <- FALSE

  if (n > 0L) {
    # Algorithm L: once the reservoir fills, skip ahead geometrically rather
    # than drawing once per item.
    w <- exp(log(stats::runif(1)) / n)
    next_swap <- n + floor(log(stats::runif(1)) / log(1 - w)) + 1

    repeat {
      block <- read(max_items - count + 1)
      if (is.null(block)) break
      held$check(block)
      k <- NROW(block)
      # Only an item past `max_items` shows that the stream held more.
      # Testing the count before reading warned on a stream of exactly
      # `max_items`, where nothing was left behind.
      if (count + k > max_items) {
        truncated <- TRUE
        k <- max_items - count
        if (!is.null(give_back)) give_back(block, k)
      }

      pos <- numeric(0)
      slot <- integer(0)
      if (count < n) {
        pos <- seq_len(min(k, n - count))
        slot <- count + pos
      }
      while (next_swap <= count + k) {
        m <- length(pos) + 1L
        pos[m] <- next_swap - count
        slot[m] <- sample.int(n, 1L)
        w <- w * exp(log(stats::runif(1)) / n)
        next_swap <- next_swap + floor(log(stats::runif(1)) / log(1 - w)) + 1
      }
      if (length(pos) > 0L) held$put(block, pos, slot, count + pos)

      count <- count + k
      if (truncated) break
    }
  }

  list(items = held$items(), count = count, truncated = truncated)
}

#' The reservoir of a stream of items
#'
#' `put(block, pos, slot, at)` stores the items at positions `pos` of a block
#' in `slot`, in order, so a slot written twice keeps the later item. The
#' reservoir starts with at most 1,024 slots and doubles until it holds `n`,
#' so a large `n` costs memory only for the items a stream supplies.
#'
#' @noRd
item_reservoir <- function(n) {
  slots <- min(n, 1024)
  kept <- vector("list", slots)
  arrived <- numeric(slots)
  list(
    check = function(block) NULL,
    put = function(block, pos, slot, at) {
      top <- max(slot)
      if (top > slots) {
        slots <<- min(n, max(2 * slots, top))
        length(kept) <<- slots
        length(arrived) <<- slots
      }
      kept[slot] <<- block[pos]
      arrived[slot] <<- at
    },
    items = function() {
      # Trim rather than pad with NULLs when the stream was short.
      filled <- which(arrived > 0)
      kept[filled[order(arrived[filled])]]
    }
  )
}

#' The reservoir of a stream of data frames, sampled by row
#'
#' Each block that adds rows keeps a piece holding only those rows, and every
#' slot points at a row of a piece. A piece is dropped once none of its rows
#' is held, and the pieces are bound into one when they hold more than twice
#' the rows in the reservoir, or 1,024, so between blocks they hold at most
#' `max(2 * n, 1024)` rows. `check(block)` refuses a block that is not a data
#' frame or whose columns differ from the first block's.
#'
#' @noRd
row_reservoir <- function(n) {
  slots <- min(n, 1024)
  piece_of <- integer(slots)
  row_of <- integer(slots)
  arrived <- numeric(slots)
  pieces <- list()
  live <- integer(0)
  stored <- 0
  filled <- 0
  first <- NULL
  seen <- 0

  # The held rows, as one data frame in the order they arrived. Pieces are
  # numbered in arrival order and rows within a piece too, so taking each
  # piece's rows in turn keeps that order.
  bind_held <- function() {
    held <- seq_len(filled)
    held <- held[order(arrived[held])]
    from <- piece_of[held]
    parts <- lapply(unique(from), function(p) {
      reindex(pieces[[p]], row_of[held][from == p])
    })
    out <- if (length(parts) == 1L) parts[[1L]] else do.call(rbind, parts)
    # The batches' row names number rows within each batch, not rows of a
    # frame, so the sample numbers its rows afresh.
    rownames(out) <- NULL
    list(data = out, slots = held)
  }

  list(
    check = function(block) {
      seen <<- seen + 1
      if (!is.data.frame(block)) {
        stop("`data` returned ", a_class(block), " on call ", seen,
             ", where `batches = TRUE` needs a data frame of rows. Return ",
             "each batch as a data frame, and NULL at the end of the stream.",
             call. = FALSE)
      }
      if (is.null(first)) {
        first <<- empty_like(block)
      } else if (!identical(names(block), names(first))) {
        stop("Batch ", seen, " from `data` ",
             column_change(names(block), names(first)),
             " where batch 1 has ", quote_labels(names(first)), ". Return ",
             "the same columns, in the same order, from every call.",
             call. = FALSE)
      }
      NULL
    },
    put = function(block, pos, slot, at) {
      last <- !duplicated(slot, fromLast = TRUE)
      pos <- pos[last]
      slot <- slot[last]
      top <- max(slot)
      if (top > slots) {
        slots <<- min(n, max(2 * slots, top))
        length(piece_of) <<- slots
        length(row_of) <<- slots
        length(arrived) <<- slots
      }
      # Slots fill in order, so the highest slot written is the count held.
      filled <<- max(filled, top)
      old <- piece_of[slot]
      old <- old[!is.na(old) & old > 0L]
      id <- length(pieces) + 1L
      pieces[[id]] <<- reindex(block, pos)
      live[id] <<- length(pos)
      stored <<- stored + length(pos)
      piece_of[slot] <<- id
      row_of[slot] <<- seq_along(pos)
      arrived[slot] <<- at[last]
      if (length(old) > 0L) {
        live <<- live - tabulate(old, nbins = id)
        for (p in unique(old)[live[unique(old)] == 0L]) {
          stored <<- stored - nrow(pieces[[p]])
          pieces[p] <<- list(NULL)
        }
      }
      if (stored > 2 * max(filled, 512)) {
        one <- bind_held()
        pieces <<- list(one$data)
        live <<- length(one$slots)
        stored <<- length(one$slots)
        piece_of[one$slots] <<- 1L
        row_of[one$slots] <<- seq_along(one$slots)
      }
    },
    items = function() {
      if (is.null(first)) return(data.frame())
      if (filled == 0) return(first)
      bind_held()$data
    }
  )
}

#' How a batch's columns differ from the first batch's
#' @noRd
column_change <- function(now, first) {
  added <- setdiff(now, first)
  lost <- setdiff(first, now)
  parts <- c(
    if (length(added) > 0L) paste("adds", quote_labels(added)),
    if (length(lost) > 0L) paste("lacks", quote_labels(lost))
  )
  if (length(parts) == 0L) {
    return(paste0("has the columns ", quote_labels(now), ","))
  }
  paste0(paste(parts, collapse = " and "), ",")
}

#' @noRd
truncation_warning <- function(max_items) {
  warn_drawn(paste0("Stopped after `max_items` = ", max_items, " item(s); ",
                    "the stream may have held more."), "drawn_stream_truncated")
}

#' Turn a connection or function into a block reader
#'
#' Returns `read(want)`, which gives the next block of items or `NULL`,
#' `close`, a connection opened here for the caller to close when the reading
#' is done, and `give_back(block, k)`, which pushes the items of `block` past
#' the first `k` back onto a connection the caller opened. A function of
#' items, read one call at a time, gets none of them.
#'
#' @noRd
as_stream <- function(x, batches = FALSE, block = 10000L) {
  if (is.function(x)) {
    if (length(formals(x)) > 0L) {
      stop("A stream function must take no arguments.", call. = FALSE)
    }
    # A function of items is read one call at a time by reservoir_stream().
    return(if (batches) list(read = function(want) x()) else list())
  }
  if (inherits(x, "connection")) {
    if (batches) {
      stop("`batches = TRUE` needs a function that returns data frames, and ",
           "`data` is a connection, which gives lines. Read the connection ",
           "in a function that returns each block of lines as a data frame; ",
           "see the examples in ?design_reservoir.", call. = FALSE)
    }
    # An on.exit() here would close the connection when *this* function
    # returns, long before the reader below has read a byte from it. The
    # caller's own open connections are left alone; one opened here is closed
    # by draw_design(), which knows when the reading is done.
    opened_here <- !isOpen(x)
    if (opened_here) open(x, "r")
    read <- function(want) {
      lines <- readLines(x, n = min(block, want), warn = FALSE)
      if (length(lines) == 0L) NULL else lines
    }
    # The line read past `max_items` goes back onto a connection the caller
    # opened, so their next read starts with it. Only a text-mode connection
    # takes lines back, and one opened here is closed anyway.
    give_back <- if (!opened_here && identical(summary(x)$text, "text")) {
      function(lines, k) pushBack(lines[-seq_len(k)], x)
    }
    return(list(read = read, close = if (opened_here) x,
                give_back = give_back))
  }
  stop("`data` must be a data frame, list, connection, or zero-argument ",
       "function for a reservoir design, not ", a_class(x), ".",
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
