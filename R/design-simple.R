#' Simple random sampling
#'
#' Draws `n` rows uniformly at random.
#'
#' @param n Number of rows to draw. A single non-negative whole number.
#' @param replace Sample with replacement? Note that `draw(weights = TRUE)`
#'   declines a with-replacement design: `.prob` would be the chance of
#'   appearing at least once, which does not weight a sample holding
#'   duplicates. See [draw()].
#'
#' @return A design object, for use with [draw()].
#'
#' @examples
#' df <- data.frame(id = 1:100)
#' draw(df, design_simple(n = 10), seed = 1)
#' draw(df, design_simple(n = 150, replace = TRUE), seed = 1)
#'
#' @family designs
#' @seealso [draw()]
#' @export
design_simple <- function(n, replace = FALSE) {
  new_design("simple", list(
    n = check_count(n, "n"),
    replace = check_flag(replace, "replace")
  ))
}

#' @export
draw_design.drawn_design_simple <- function(design, data) {
  validate_data(data)

  check_draw_size(design$n, nrow(data), design$replace)

  idx <- sample.int(nrow(data), size = design$n, replace = design$replace)
  reindex(data, idx)
}

# ---- inclusion probability ------------------------------------------------

#' @noRd
exact_inclusion.drawn_design_simple <- function(design, data) {
  if (design$replace) {
    # P(selected at least once) over n independent draws.
    return(rep(1 - (1 - 1 / nrow(data))^design$n, nrow(data)))
  }
  rep(design$n / nrow(data), nrow(data))
}
