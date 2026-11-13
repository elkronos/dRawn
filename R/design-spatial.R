#' Sampling within a region
#'
#' Draws a simple random sample of the rows whose coordinates fall inside a
#' region. It decides *where* sampling happens, not how evenly it covers the
#' ground: for a sample spread evenly across space — usually far more precise
#' when what you measure varies smoothly — use [design_spread()], on the frame
#' restricted to the region if need be.
#'
#' @section Coordinate order:
#' `coords` is `c(x, y)` — longitude first, then latitude — matching
#' [sf::st_as_sf()]. Getting this backwards puts your points somewhere else
#' entirely, so it is worth checking against a known landmark once. When the
#' coordinates are longitude and latitude, a latitude outside -90 to 90 or a
#' longitude outside -180 to 360 is refused, since it means the columns are
#' the wrong way round or hold projected coordinates.
#'
#' A POINT geometry column, such as the geometry of an `sf` frame, can stand
#' in for the two columns: name it alone, as in `coords = "geometry"`, and
#' each row's position is read with [sf::st_coordinates()]. The column's own
#' CRS is used when it has one, and `region` is reprojected to it.
#'
#' @section Edges of a longitude/latitude region:
#' With spherical geometry ([sf::sf_use_s2()], the default), each edge of a
#' longitude/latitude `region` follows the great circle between its two
#' corners, the shortest path on the globe. On a flat longitude/latitude map
#' that path bows towards the nearer pole wherever the edge runs east-west
#' away from the equator; edges along a meridian or the equator do not bow.
#' A wide box is therefore a different shape from the rectangle it looks like
#' on such a map. The top edge of a box from -50 to 50 degrees of longitude
#' and -30 to 30 of latitude reaches latitude 41.9 at longitude 0, so a row
#' at (0, 41) is inside it, and its bottom edge reaches -41.9. A box from -5
#' to 5 degrees of longitude and 40 to 50 of latitude takes in a row at
#' (0, 50.05) and leaves out one at (0, 40.05).
#'
#' [draw()], [inclusion_prob()], [joint_prob()] and
#' [plot()][plot.drawn_design] warn, with class `drawn_geodesic_edges`, when
#' an edge of `region` passes more than 0.1 degrees (about 11 km) from the
#' straight line between its corners on a longitude/latitude map. To make the
#' edges follow those straight lines, add vertices along them in degree space
#' before you build the design:
#'
#' ```
#' flat <- sf::st_set_crs(region, NA)    # coordinates as plain numbers
#' flat <- sf::st_segmentize(flat, 0.1)  # a vertex at least every 0.1 degrees
#' region <- sf::st_set_crs(flat, 4326)
#' ```
#'
#' An edge 0.1 degrees long bows by about a metre at most. The CRS has to come
#' off for the middle step: on a region that keeps it,
#' [sf::st_segmentize()] places the new vertices along the great circles,
#' under either s2 setting, and the bow stays. A region larger than a
#' hemisphere needs splitting first; see the next section.
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
#' normalises winding, so both directions give the same strip. Adding vertices
#' along the edges in degree space, as in the previous section, does not
#' rescue a region this large either. s2 takes the smaller of the two areas a
#' ring divides the globe into as its inside, so a ring around more than half
#' the globe encloses the rest of the globe instead.
#'
#' The short way is the right reading for a region meant to cross the
#' antimeridian, such as a box from longitude 179 to -179. [draw()],
#' [inclusion_prob()], [joint_prob()] and [plot()][plot.drawn_design] warn,
#' with class `drawn_antimeridian`, whenever an edge of `region` joins two
#' longitudes more than 180 degrees apart, and the warning gives both
#' readings. If you meant the short way, `region` is what you drew. If you
#' meant the long way, split the region into parts smaller than a hemisphere,
#' add vertices to each as in the previous section, and pass the parts as a
#' list; or switch to planar interpretation with `sf::sf_use_s2(FALSE)`.
#' Splitting works under the default setting. Switching has to happen before
#' [draw()] and stay in force for every estimate from the sample; see the next
#' section.
#'
#' @section The s2 setting at draw time and after:
#' Which rows fall inside `region` depends on [sf::sf_use_s2()]. With spherical
#' geometry the edges of a longitude/latitude region are great-circle arcs, and
#' without it they are straight lines in longitude and latitude, so rows near an
#' edge can be inside under one setting and outside under the other. Use the
#' same setting when you draw a sample and whenever you analyse it.
#'
#' [ht_total()], [ht_mean()], [as_svydesign()] and [sample_summary()] find the
#' region's rows again each time they run, because the variance, the design
#' effect and `survey`'s population count all depend on them. They compare the
#' result with the inclusion probabilities stored at draw time and refuse the
#' sample when the two differ, naming `sf::sf_use_s2()`. The point estimate
#' would not have moved, since it uses the stored `.prob`, but the standard
#' error would.
#'
#' `sf_use_s2()` is a session setting, and a new R session starts with it
#' `TRUE` unless the `_SF_USE_S2` environment variable says otherwise. A sample
#' drawn under `sf::sf_use_s2(FALSE)` and saved with [saveRDS()] needs
#' `sf::sf_use_s2(FALSE)` again after [readRDS()].
#' [inclusion_prob()] and [joint_prob()] on the frame follow the current
#' setting too.
#'
#' @param coords Two column names, `c(x, y)` with longitude first, or the
#'   name of one POINT geometry column. See "Coordinate order".
#' @param region An `sf` or `sfc` geometry of POLYGON or MULTIPOLYGON
#'   features, or a list of them. A row inside any of them is inside the
#'   region. Each is reprojected to the coordinates' CRS when its own differs.
#' @param n Number of rows to draw from inside the region.
#' @param crs Coordinate reference system of the coordinate columns. Defaults to
#'   EPSG:4326. A geometry column named in `coords` keeps its own CRS when it
#'   has one.
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
  if (length(coords) > 2L) {
    stop("`coords` must name exactly two columns, c(x, y) -- longitude then ",
         "latitude -- or one POINT geometry column. Got ", length(coords), ".",
         call. = FALSE)
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
  inside <- which(spatial_setup(design, data))
  reindex(data, sort(inside[sample.int(length(inside), design$n)]))
}

#' Rows inside the region, refusing a region that cannot supply `n`
#'
#' The checks the draw and the probabilities share, so [inclusion_prob()]
#' refuses what [draw()] refuses, with the same message: no row inside
#' `region`, or fewer than `n`. Returns `spatial_inside()`.
#'
#' @noRd
spatial_setup <- function(design, data) {
  require_suggested("sf", "spatial sampling")
  inside <- spatial_inside(design, data)
  k <- sum(inside)
  if (k == 0L) {
    stop("No rows fall inside `region`. If you expected a much larger area, ",
         "see the \"Regions that cross the antimeridian\" section of ",
         "?design_spatial.", call. = FALSE)
  }
  if (design$n > k) {
    stop("`n` (", design$n, ") exceeds the ", k,
         " row(s) inside `region`.", call. = FALSE)
  }
  inside
}

#' Which rows fall inside the design's region
#'
#' Shared by the draw method and [inclusion_prob()], so the two can never
#' disagree about what "inside" means. A row with a missing coordinate is
#' outside, once the design's `na_rm` allows it.
#'
#' @noRd
spatial_inside <- function(design, data) {
  pts <- spatial_points(design, data)
  x <- pts$x
  y <- pts$y

  ok <- !is.na(x) & !is.na(y)
  what <- paste0("missing coordinates in ", pts$label)
  check_na_policy(!ok, design$na_rm, what)
  if (length(ok) && !any(ok)) {
    stop("Every row has ", what, ".", call. = FALSE)
  }
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("`", design$coords[1], "` and `", design$coords[2],
         "` must both be numeric.", call. = FALSE)
  }
  check_lonlat(x[ok], y[ok], pts)

  region <- region_geometry(design$region, pts$crs)
  data_sf <- sf::st_as_sf(data.frame(x = x[ok], y = y[ok]),
                          coords = c("x", "y"), crs = pts$crs)
  out <- rep(FALSE, nrow(data))
  out[ok] <- lengths(sf::st_intersects(data_sf, region)) > 0L
  out
}

#' Each row's two coordinates, and the CRS they are in
#'
#' `coords` names two columns, or one POINT geometry column whose X and Y
#' are read with `sf::st_coordinates()`, and `geometry` says which. `label`
#' names the columns for messages, and `names` the two coordinates.
#'
#' @noRd
spatial_points <- function(design, data) {
  coords <- design$coords
  validate_data(data, required_columns = coords)
  cols <- coordinate_columns(data, coords, "coords")
  if (length(cols) != 2L) {
    stop("`coords` must name two coordinate columns, c(x, y) -- longitude ",
         "then latitude -- or one POINT geometry column, not ",
         format_bad(coords), ".", call. = FALSE)
  }
  list(
    x = cols[[1L]], y = cols[[2L]], names = names(cols),
    crs = spatial_crs(design, data), geometry = length(coords) == 1L,
    label = paste0("`", coords, "`", collapse = " or ")
  )
}

#' The CRS of the coordinates: a geometry column's own, else `crs`
#' @noRd
spatial_crs <- function(design, data) {
  geom <- if (length(design$coords) == 1L && is.data.frame(data)) {
    data[[design$coords]]
  }
  if (inherits(geom, "sfc") && !is.na(sf::st_crs(geom))) {
    return(sf::st_crs(geom))
  }
  sf::st_crs(design$crs)
}

#' Numeric coordinates behind column names
#'
#' A POINT geometry column, such as an `sf` frame's geometry, stands for its
#' X and Y coordinates, read with `sf::st_coordinates()`; an empty point
#' gives `NA`. Any other geometry is refused. With `time = TRUE` a Date column
#' becomes days and a date-time column seconds since 1970. Other columns come
#' back as they are, for the caller to check. The result is a named list of
#' vectors, two for a geometry column and one for any other.
#'
#' @noRd
coordinate_columns <- function(data, cols, arg, time = FALSE) {
  out <- lapply(cols, function(nm) {
    col <- data[[nm]]
    if (inherits(col, "sfc")) {
      return(point_coordinates(col, nm, arg))
    }
    if (time && inherits(col, "Date")) {
      col <- as.numeric(col)
    } else if (time && inherits(col, "POSIXt")) {
      col <- as.numeric(as.POSIXct(col))
    }
    stats::setNames(list(col), nm)
  })
  do.call(c, out)
}

#' The X and Y coordinates of a POINT geometry column
#' @noRd
point_coordinates <- function(col, nm, arg) {
  require_suggested("sf", "coordinates from a geometry column")
  type <- unique(as.character(sf::st_geometry_type(col, by_geometry = TRUE)))
  other <- setdiff(type, "POINT")
  if (length(other)) {
    stop("`", arg, "` names `", nm, "`, a geometry column holding ", other[1],
         " geometry. Only POINT geometry gives each row one position; ",
         "reduce each feature to a point first, for example with ",
         "sf::st_point_on_surface().", call. = FALSE)
  }
  xy <- sf::st_coordinates(col)
  stats::setNames(list(unname(xy[, "X"]), unname(xy[, "Y"])),
                  paste(nm, c("X", "Y")))
}

#' Refuse coordinates that cannot be longitude and latitude
#'
#' Under a longitude/latitude CRS, sf reads longitudes from -180 to 360 and
#' latitudes from -90 to 90, and only warns about values outside them. A
#' latitude past 90 is the sign of columns in (lat, lon) order, and a value
#' in the hundreds of thousands the sign of projected coordinates under the
#' default `crs`.
#'
#' @noRd
check_lonlat <- function(x, y, pts) {
  if (!isTRUE(sf::st_is_longlat(pts$crs))) return(invisible(NULL))
  bad_y <- y[!(y >= -90 & y <= 90)]
  bad_x <- x[!(x >= -180 & x <= 360)]
  if (!length(bad_y) && !length(bad_x)) return(invisible(NULL))
  k <- if (length(bad_y)) 2L else 1L
  bad <- if (k == 2L) bad_y else bad_x
  v <- bad[which.max(abs(bad))]
  what <- if (k == 2L) {
    "-90 to 90, the range of a latitude"
  } else {
    "-180 to 360, the widest range sf reads as longitude"
  }
  col <- if (pts$geometry) {
    paste0("The ", c("X", "Y")[k], " coordinate of ", pts$label)
  } else {
    paste0("`", pts$names[k], "`")
  }
  fix <- if (pts$geometry) {
    "If they are projected, set the column's CRS with sf::st_set_crs()."
  } else {
    paste0("`coords` is c(x, y), longitude first: if it names the columns ",
           "in (lat, lon) order, swap them. If the coordinates are ",
           "projected, set `crs` to their system.")
  }
  stop(col, " holds ", format(v, big.mark = ",", scientific = FALSE,
                              trim = TRUE),
       ", outside ", what, ", and the coordinates are in a ",
       "longitude/latitude CRS (", crs_label(pts$crs), "). ", fix,
       call. = FALSE)
}

#' @noRd
crs_label <- function(crs) {
  if (!is.na(crs$epsg)) paste0("EPSG:", crs$epsg) else crs$input
}

#' `region` as one geometry in the coordinates' CRS
#'
#' Each member of a list is checked and reprojected on its own before they
#' are combined, since `c()` refuses geometries in different CRSs. The
#' members stay separate features rather than being unioned: a row inside
#' any of them is inside the region either way, and under s2 the union of
#' two parts that together cover more than a hemisphere comes back as the
#' rest of the globe.
#'
#' @noRd
region_geometry <- function(region, crs) {
  crs <- sf::st_crs(crs)
  listed <- is.list(region) && !inherits(region, c("sf", "sfc"))
  parts <- if (listed) region else list(region)
  if (!length(parts)) {
    stop("`region` is an empty list.", call. = FALSE)
  }
  geoms <- lapply(seq_along(parts), function(i) {
    part <- parts[[i]]
    arg <- if (listed) paste0("`region[[", i, "]]`") else "`region`"
    g <- tryCatch(sf::st_geometry(part), error = function(e) {
      stop(arg, " must be an sf or sfc geometry",
           if (!listed) ", or a list of them", ", not ",
           if (is.null(part)) "NULL" else class(part)[1], ".", call. = FALSE)
    })
    if (is.na(sf::st_crs(g))) {
      stop(arg, " has no CRS. Set one with sf::st_set_crs().", call. = FALSE)
    }
    check_region_type(g, arg)
    if (sf::st_crs(g) != crs) g <- sf::st_transform(g, crs)
    g
  })
  if (length(geoms) == 1L) geoms[[1L]] else do.call(c, geoms)
}

#' A region has to enclose an area
#'
#' A point or a line has no inside, so no row could fall within it, and a
#' sample drawn from it would be refused only at [draw()] with "No rows fall
#' inside `region`".
#'
#' @noRd
check_region_type <- function(g, arg) {
  type <- unique(as.character(sf::st_geometry_type(g, by_geometry = TRUE)))
  other <- setdiff(type, c("POLYGON", "MULTIPOLYGON"))
  if (!length(other)) return(invisible(NULL))
  fix <- if (any(other %in% c("GEOMETRYCOLLECTION", "GEOMETRY"))) {
    "Take its polygons with sf::st_collection_extract(region, \"POLYGON\")."
  } else {
    paste0("A point or a line encloses no area; for the area around it, ",
           "use sf::st_buffer().")
  }
  stop(arg, " must be POLYGON or MULTIPOLYGON geometry, not ", other[1],
       ". ", fix, call. = FALSE)
}

#' The outline of `region`, as `plot(type = "map")` draws it
#'
#' The rings of `region` in the coordinates' CRS, as `x` and `y` with an `NA`
#' between rings. Under spherical geometry each edge of a longitude/latitude
#' region is the great-circle arc between its corners, which bows away from
#' the straight line a flat map would draw, so vertices are added along each
#' arc at most 50 km apart before the outline is drawn. Under planar
#' geometry, or with projected coordinates, the edges are the straight lines
#' between the corners.
#'
#' @noRd
region_outline <- function(design, data) {
  crs <- spatial_crs(design, data)
  region <- region_geometry(design$region, crs)
  if (isTRUE(sf::st_is_longlat(crs)) && isTRUE(sf::sf_use_s2())) {
    region <- sf::st_segmentize(region, 50000)
  }
  cs <- sf::st_coordinates(sf::st_cast(region, "MULTIPOLYGON"))
  ring <- cs[, intersect(c("L1", "L2", "L3"), colnames(cs)), drop = FALSE]
  id <- do.call(paste, c(unname(as.data.frame(ring)), sep = "."))
  pieces <- split(seq_len(nrow(cs)), factor(id, levels = unique(id)))
  at <- unlist(lapply(pieces, function(i) c(i, NA)), use.names = FALSE)
  at <- at[-length(at)]
  list(x = unname(cs[at, "X"]), y = unname(cs[at, "Y"]))
}

# ---- warnings about the region's edges --------------------------------------

#' Warn about edges that s2 draws differently from a flat map
#'
#' Called once by [draw()], [inclusion_prob()] and [joint_prob()], and through
#' them by `plot()`, and not by the draw method, so a simulation of `R` draws
#' warns once. Only spherical geometry on longitude/latitude coordinates is
#' affected. Anything that fails here fails again, with its own message, in
#' the call that follows, so failures are passed over in silence.
#'
#' @noRd
warn_region <- function(design, data) {
  while (inherits(design, "drawn_design_certainty")) design <- design$rest
  if (!inherits(design, "drawn_design_spatial") ||
      !requireNamespace("sf", quietly = TRUE) || !isTRUE(sf::sf_use_s2())) {
    return(invisible(NULL))
  }
  edges <- tryCatch({
    crs <- spatial_crs(design, data)
    if (isTRUE(sf::st_is_longlat(crs))) {
      region_edges(region_geometry(design$region, crs))
    }
  }, error = function(e) NULL)
  if (is.null(edges) || !length(edges$x0)) return(invisible(NULL))
  warn_antimeridian(edges)
  warn_geodesic_edges(edges)
  invisible(NULL)
}

#' Every edge of every ring, as the coordinates of its two ends
#' @noRd
region_edges <- function(region) {
  cs <- sf::st_coordinates(sf::st_cast(region, "MULTIPOLYGON"))
  ring <- cs[, intersect(c("L1", "L2", "L3"), colnames(cs)), drop = FALSE]
  n <- nrow(cs)
  from <- if (n > 1L) {
    which(rowSums(ring[-1L, , drop = FALSE] != ring[-n, , drop = FALSE]) == 0)
  } else {
    integer(0)
  }
  list(x0 = unname(cs[from, "X"]), y0 = unname(cs[from, "Y"]),
       x1 = unname(cs[from + 1L, "X"]), y1 = unname(cs[from + 1L, "Y"]))
}

#' Warn when an edge joins longitudes more than 180 degrees apart
#'
#' s2 draws such an edge the short way, across the antimeridian. That is
#' what a region meant to cross it needs, and the opposite of what a region
#' meant to go the long way round needs, so the warning gives both readings.
#' Checking the raw coordinates is exact and cheap; comparing areas is not,
#' because a collapsed ring and its bounding box collapse together. An edge
#' whose ends are the same point, 360 degrees apart, is not counted.
#'
#' @noRd
warn_antimeridian <- function(e) {
  dx <- e$x1 - e$x0
  span <- abs(dx)
  span[abs(dx - 360 * round(dx / 360)) < 1e-9] <- 0
  i <- which.max(span)
  if (!length(i) || !(span[i] > 180)) return(invisible(NULL))
  warn_drawn(paste0(
    "An edge of `region` joins longitudes ", fmt_coord(e$x0[i]), " and ",
    fmt_coord(e$x1[i]), ". Under spherical geometry it takes the short way, ",
    "across the antimeridian, a span of ", fmt_coord(360 - span[i]),
    " degrees. If you meant it to span the ", fmt_coord(span[i]),
    " degrees of longitude between the two values, split `region`; see the ",
    "\"Regions that cross the antimeridian\" section of ?design_spatial."),
    "drawn_antimeridian")
}

#' Warn when an edge strays far from the straight line in degrees
#'
#' The distance is measured from points a quarter, a half and three quarters
#' of the way along the straight line between the edge's ends, in longitude
#' and latitude, to the great circle through those ends, as an angle at the
#' centre of the globe. An edge longer than 180 degrees of longitude is taken
#' the short way, as s2 takes it.
#'
#' @noRd
warn_geodesic_edges <- function(e, limit = 0.1) {
  dx <- e$x1 - e$x0
  dx <- dx - 360 * round(dx / 360)
  a <- unit_vector(e$x0, e$y0)
  b <- unit_vector(e$x0 + dx, e$y1)
  normal <- cbind(a[, 2] * b[, 3] - a[, 3] * b[, 2],
                  a[, 3] * b[, 1] - a[, 1] * b[, 3],
                  a[, 1] * b[, 2] - a[, 2] * b[, 1])
  len <- sqrt(rowSums(normal^2))
  # Ends that coincide, or lie opposite each other, fix no great circle.
  flat <- len < 1e-12
  len[flat] <- 1
  bow <- numeric(length(len))
  for (t in c(0.25, 0.5, 0.75)) {
    p <- unit_vector(e$x0 + t * dx, e$y0 + t * (e$y1 - e$y0))
    off <- abs(rowSums(p * normal)) / len
    bow <- pmax(bow, asin(pmin(off, 1)) * 180 / pi)
  }
  bow[flat] <- 0
  i <- which.max(bow)
  if (!length(i) || !(bow[i] > limit)) return(invisible(NULL))
  km <- bow[i] * pi / 180 * 6371
  warn_drawn(paste0(
    "An edge of `region`, from (", fmt_coord(e$x0[i]), ", ",
    fmt_coord(e$y0[i]), ") to (", fmt_coord(e$x1[i]), ", ",
    fmt_coord(e$y1[i]), "), follows a great circle under spherical geometry ",
    "and passes ", format(signif(bow[i], 3)), " degrees (about ",
    format(signif(km, 2), big.mark = ","), " km) from the straight line ",
    "between its ends on a longitude/latitude map. Rows near it can be ",
    "inside `region` though outside the polygon as drawn on such a map, or ",
    "the reverse. See the \"Edges of a longitude/latitude region\" section ",
    "of ?design_spatial to make the edges follow those lines."),
    "drawn_geodesic_edges")
}

#' Points on the unit sphere, one row per longitude and latitude in degrees
#' @noRd
unit_vector <- function(lon, lat) {
  lon <- lon * pi / 180
  lat <- lat * pi / 180
  cbind(cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat))
}

#' @noRd
fmt_coord <- function(x) format(signif(x, 6))

# ---- inclusion probability ------------------------------------------------

#' @noRd
exact_inclusion.drawn_design_spatial <- function(design, data) {
  inside <- spatial_setup(design, data)
  out <- numeric(nrow(data))
  out[inside] <- design$n / sum(inside)
  out
}

#' A region has to be polygons `sf` can read, with a CRS
#'
#' Both were stored unchecked, so a NULL region or a nonsense CRS surfaced much
#' later as a raw `sf` message that named neither argument. The region is
#' built once here, as `region_geometry()` builds it at draw time, so each
#' member of a list is checked and reprojected to `crs` before it is stored.
#'
#' @noRd
check_region <- function(region, crs) {
  if (!requireNamespace("sf", quietly = TRUE)) return(invisible(NULL))
  tryCatch(sf::st_crs(crs), error = function(e) {
    stop("`crs` is not a coordinate reference system sf recognises: ",
         conditionMessage(e), call. = FALSE)
  })
  if (is.na(sf::st_crs(crs))) {
    stop("`crs` is not a coordinate reference system sf recognises.",
         call. = FALSE)
  }
  region_geometry(region, crs)
  invisible(NULL)
}
