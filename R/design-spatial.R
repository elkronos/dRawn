#' Spatial sampling
#'
#' Samples rows whose coordinates fall inside a region.
#'
#' @section Coordinate order:
#' `coords` is `c(x, y)` — longitude first, then latitude — matching
#' [sf::st_as_sf()]. Getting this backwards puts your points somewhere else
#' entirely, so it is worth checking against a known landmark once.
#'
#' @section Regions that cross the antimeridian:
#' With spherical geometry enabled ([sf::sf_use_s2()], the default), consecutive
#' polygon vertices are joined by the *shortest* great-circle path. An edge from
#' longitude -179 to +179 therefore spans the 2 degrees across the antimeridian,
#' not the 358 degrees the coordinates suggest. A "whole world" rectangle
#' collapses into a 2-degree-wide pole-to-pole strip of 2.8 million km2 --
#' against the roughly 510 million km2 of the globe -- and contains almost
#' nothing.
#'
#' Ring orientation is not the cause and reversing it does not help: `sf`
#' normalises winding, so both directions give the same strip. Adding
#' intermediate vertices does not reliably help either, because the intermediate
#' edges still bow along geodesics.
#'
#' [draw()] warns whenever any edge of `region` spans more than 180 degrees of
#' longitude. If you hit it, either split the region at the antimeridian into two
#' polygons, or switch to planar interpretation with `sf::sf_use_s2(FALSE)`.
#'
#' @param coords Two column names, `c(x, y)` — longitude then latitude.
#' @param region An `sf` or `sfc` geometry, or a list of them, which is unioned.
#'   Reprojected to `crs` when it differs.
#' @param n Number of rows to draw from inside the region.
#' @param crs Coordinate reference system of the coordinate columns. Defaults to
#'   EPSG:4326.
#' @param na_rm Drop rows with missing coordinates instead of raising an error.
#'
#' @return A design object, for use with [draw()].
#'
#' @examplesIf requireNamespace("sf", quietly = TRUE)
#' box <- sf::st_sfc(
#'   sf::st_polygon(list(cbind(c(0, 10, 10, 0, 0), c(0, 0, 10, 10, 0)))),
#'   crs = 4326
#' )
#' df <- data.frame(id = 1:20,
#'                  lon = seq(1, 9, length.out = 20),
#'                  lat = seq(1, 9, length.out = 20))
#' draw(df, design_spatial(c("lon", "lat"), region = box, n = 5), seed = 1)
#'
#' @family designs
#' @seealso [draw()]
#' @export
design_spatial <- function(coords, region, n, crs = 4326, na_rm = FALSE) {
  coords <- check_columns(coords, "coords")
  if (length(coords) != 2L) {
    stop("`coords` must name exactly two columns, c(x, y) -- longitude then ",
         "latitude. Got ", length(coords), ".", call. = FALSE)
  }
  check_region(region, crs)
  new_design("spatial", list(
    coords = coords,
    region = region,
    n = check_count(n, "n"),
    crs = crs,
    na_rm = check_flag(na_rm, "na_rm")
  ))
}

#' @export
draw_design.drawn_design_spatial <- function(design, data) {
  require_suggested("sf", "spatial sampling")
  validate_data(data, required_columns = design$coords)

  bad <- is.na(data[[design$coords[1]]]) | is.na(data[[design$coords[2]]])
  if (any(bad)) {
    data <- drop_na_rows(data, bad, design$na_rm, "missing coordinates")
  }

  inside <- which(spatial_inside(design, data, warn = TRUE))

  if (length(inside) == 0L) {
    stop("No rows fall inside `region`. If you expected a much larger area, ",
         "see the \"Regions that cross the antimeridian\" section of ",
         "?design_spatial.", call. = FALSE)
  }
  if (design$n > length(inside)) {
    stop("`n` (", design$n, ") exceeds the ", length(inside),
         " row(s) inside `region`.", call. = FALSE)
  }

  reindex(data, sort(inside[sample.int(length(inside), design$n)]))
}

#' Which rows fall inside the design's region
#'
#' Shared by the draw method and [inclusion_prob()], so the two can never
#' disagree about what "inside" means.
#'
#' @noRd
spatial_inside <- function(design, data, warn = FALSE) {
  x_col <- design$coords[1]
  y_col <- design$coords[2]
  x <- data[[x_col]]
  y <- data[[y_col]]
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("`", x_col, "` and `", y_col, "` must both be numeric.", call. = FALSE)
  }

  region <- design$region
  if (is.list(region) && !inherits(region, c("sf", "sfc"))) {
    # do.call(st_union, region) maps the third and later elements onto
    # st_union()'s own arguments.
    region <- sf::st_union(do.call(c, lapply(region, sf::st_geometry)))
  }
  region <- sf::st_geometry(region)

  if (is.na(sf::st_crs(region))) {
    stop("`region` has no CRS. Set one with sf::st_set_crs().", call. = FALSE)
  }

  ok <- !is.na(x) & !is.na(y)
  check_na_policy(!ok, design$na_rm,
                  paste0("a missing `", x_col, "` or `", y_col, "`"))
  data_sf <- sf::st_as_sf(data[ok, , drop = FALSE], coords = c(x_col, y_col),
                          crs = design$crs, remove = FALSE)
  if (sf::st_crs(region) != sf::st_crs(data_sf)) {
    region <- sf::st_transform(region, sf::st_crs(data_sf))
  }

  if (isTRUE(warn)) warn_if_antimeridian(region)

  out <- rep(FALSE, nrow(data))
  out[ok] <- lengths(sf::st_intersects(data_sf, region)) > 0L
  out
}

#' Warn when a ring edge spans more than 180 degrees of longitude
#'
#' Under s2 such an edge is drawn the short way, across the antimeridian, so the
#' polygon is not the region the coordinates suggest. Checking the raw
#' coordinates is exact and cheap; comparing areas is not, because a collapsed
#' ring and its bounding box collapse together.
#'
#' @noRd
warn_if_antimeridian <- function(region) {
  if (!isTRUE(sf::sf_use_s2())) {
    return(invisible(NULL))
  }
  if (is.na(sf::st_crs(region)) || !isTRUE(sf::st_is_longlat(region))) {
    return(invisible(NULL))
  }

  span <- tryCatch({
    cs <- sf::st_coordinates(region)
    ring_cols <- intersect(c("L1", "L2", "L3"), colnames(cs))
    rings <- split(
      cs[, "X"],
      interaction(as.data.frame(cs[, ring_cols, drop = FALSE]), drop = TRUE)
    )
    max(vapply(rings, function(v) if (length(v) < 2L) 0 else max(abs(diff(v))),
               numeric(1)))
  }, error = function(e) NA_real_)

  if (is.finite(span) && span > 180) {
    warning("An edge of `region` spans ", round(span),
            " degrees of longitude. Under spherical geometry that edge is ",
            "drawn the short way, across the antimeridian, so `region` is much ",
            "smaller than its coordinates suggest. See the \"Regions that cross ",
            "the antimeridian\" section of ?design_spatial.", call. = FALSE)
  }
  invisible(NULL)
}

# ---- inclusion probability ------------------------------------------------

#' @noRd
exact_inclusion.drawn_design_spatial <- function(design, data) {
  require_suggested("sf", "spatial sampling")
  validate_data(data, required_columns = design$coords)
  inside <- spatial_inside(design, data)
  out <- numeric(nrow(data))
  if (!any(inside)) return(out)
  out[inside] <- min(1, design$n / sum(inside))
  out
}

#' A region has to be geometry `sf` can read, with a CRS
#'
#' Both were stored unchecked, so a NULL region or a nonsense CRS surfaced much
#' later as a raw `sf` message that named neither argument.
#'
#' @noRd
check_region <- function(region, crs) {
  if (!requireNamespace("sf", quietly = TRUE)) return(invisible(NULL))
  geom <- tryCatch({
    if (is.list(region) && !inherits(region, c("sf", "sfc"))) {
      if (!length(region)) {
        stop("`region` is an empty list.", call. = FALSE)
      }
      sf::st_union(do.call(c, lapply(region, sf::st_geometry)))
    } else {
      sf::st_geometry(region)
    }
  }, error = function(e) {
    stop("`region` must be an sf or sfc geometry, or a list of them, not ",
         if (is.null(region)) "NULL" else class(region)[1], ". (", 
         conditionMessage(e), ")", call. = FALSE)
  })
  if (is.na(sf::st_crs(geom))) {
    stop("`region` has no CRS. Set one with sf::st_set_crs().", call. = FALSE)
  }
  tryCatch(sf::st_crs(crs), error = function(e) {
    stop("`crs` is not a coordinate reference system sf recognises: ",
         conditionMessage(e), call. = FALSE)
  })
  if (is.na(sf::st_crs(crs))) {
    stop("`crs` is not a coordinate reference system sf recognises.",
         call. = FALSE)
  }
  invisible(NULL)
}
