# Every design in design_rows() down every path a sample takes: draw() with
# and without weights, inclusion_prob() exact and simulated,
# sampling_weight(), joint_prob(), ht_total() under each variance, ht_mean(),
# domains, degrees of freedom, as_svydesign(), sample_summary(), plot() and
# deff(). A path an entry refuses on purpose is checked for its message, as
# any other result is; every other path has to run and agree with the
# others. The bugs a path once had are pinned in the tests for its area.
# This file checks that the paths agree for every design, including the ones
# those tests leave out, and a new design added to design_rows() is checked
# the same way.

rows <- design_rows()

# Runs `code` for `row`, muffling the warning the design gives on every draw.
quietly <- function(row, code) {
  withCallingHandlers(code, warning = function(w) {
    if (!is.null(row$warns) && inherits(w, row$warns)) {
      invokeRestart("muffleWarning")
    }
  })
}

# The value of `code`, or, for a path `row` refuses, NULL after checking that
# `code` stops with the message the entry gives.
cell <- function(row, path, code, label) {
  why <- row$refuse[[path]]
  if (is.null(why)) return(quietly(row, code))
  expect_error(quietly(row, code), why, fixed = TRUE, label = label)
  NULL
}

# The weighted sample of `row`, or NULL for a design that refuses one.
grid_sample <- function(row) {
  if (!is.null(row$refuse$weights)) return(NULL)
  quietly(row, draw(row$frame, row$design, seed = 1, weights = TRUE))
}

# Sampling units below certainty less the strata that hold them, worked out
# from the entry's own description of its units and strata.
reference_df <- function(row, s) {
  free <- s$.prob < 1 - 1e-9
  if (!any(free)) return(Inf)
  units <- row$frame[s$.row[free], , drop = FALSE]
  psu <- if (is.null(row$psu)) seq_len(nrow(units)) else units[[row$psu]]
  stratum <- if (is.null(row$strata)) {
    rep(1, nrow(units))
  } else if (is.function(row$strata)) {
    row$strata(units)
  } else {
    do.call(paste, c(unname(as.list(units[row$strata])), sep = "\r"))
  }
  length(unique(psu)) - length(unique(stratum))
}

test_that("every design draws the same rows for a seed, with the frame's columns", {
  for (nm in names(rows)) {
    row <- rows[[nm]]
    if (is.null(row$warns)) {
      a <- draw(row$frame, row$design, seed = 1)
    } else {
      expect_warning(a <- draw(row$frame, row$design, seed = 1),
                     class = row$warns)
    }
    expect_identical(quietly(row, draw(row$frame, row$design, seed = 1)), a,
                     label = nm)
    if (is.data.frame(a)) {
      expect_identical(setdiff(names(a), ".replicate"), names(row$frame),
                       label = nm)
      size <- nrow(a)
    } else {
      size <- length(a)
    }
    if (row$fixed) expect_equal(size, row$n, label = nm)
  }
})

test_that("a weighted sample holds the rows drawn and their probabilities", {
  for (nm in names(rows)) {
    row <- rows[[nm]]
    s <- cell(row, "weights",
              draw(row$frame, row$design, seed = 1, weights = TRUE), nm)
    if (is.null(s)) next
    plain <- quietly(row, draw(row$frame, row$design, seed = 1))
    expect_identical(names(s), c(".prob", ".weight", ".row",
                                 names(row$frame)), label = nm)
    expect_identical(s$id, plain$id, label = nm)
    expect_identical(s$id, row$frame$id[s$.row], label = nm)
    expect_equal(s$.weight, 1 / s$.prob, label = nm)
    if (is.null(row$refuse$inclusion)) {
      p <- quietly(row, inclusion_prob(row$frame, row$design))
      expect_equal(s$.prob, p[s$.row], label = nm)
    } else {
      # Weighted on the clusters selected, which no frame-wide probability
      # gives.
      expect_true(all(s$.prob > 0 & s$.prob <= 1), label = nm)
    }
  }
})

test_that("inclusion probabilities sum to the expected size and give the weights", {
  for (nm in names(rows)) {
    row <- rows[[nm]]
    p <- cell(row, "inclusion", inclusion_prob(row$frame, row$design), nm)
    w <- cell(row, "inclusion", sampling_weight(row$frame, row$design), nm)
    if (is.null(p)) next
    expect_length(p, nrow(row$frame))
    expect_true(all(p >= 0 & p <= 1), label = nm)
    expect_equal(sum(p), row$expected, label = nm)
    expect_equal(w, ifelse(p > 0, 1 / p, NA_real_), label = nm)
  }
})

test_that("simulated inclusion agrees with the exact probabilities", {
  skip_on_cran()
  R <- 100
  for (nm in names(rows)) {
    row <- rows[[nm]]
    if (!is.null(row$refuse$inclusion)) next
    p <- quietly(row, inclusion_prob(row$frame, row$design))
    sim <- quietly(row, inclusion_prob(row$frame, row$design, simulate = TRUE,
                                       R = R, seed = 1))
    # Each row's count of draws is binomial, exactly, so its smaller tail
    # probability is a test of that row whatever its probability. A row at 0
    # or 1 fails on a single draw against it.
    hits <- round(sim * R)
    tail <- pmin(stats::pbinom(hits, R, p),
                 stats::pbinom(hits - 1, R, p, lower.tail = FALSE))
    expect_gt(min(tail), 1e-7, label = nm)
  }
})

test_that("joint probabilities agree with the inclusion probabilities", {
  for (nm in names(rows)) {
    row <- rows[[nm]]
    m <- cell(row, "joint", joint_prob(row$frame, row$design), nm)
    if (is.null(m)) next
    m <- unname(m)
    p <- quietly(row, inclusion_prob(row$frame, row$design))
    expect_equal(diag(m), p, label = nm)
    expect_equal(m, t(m), label = nm)
    expect_true(all(m >= 0 & m <= outer(p, p, pmin) + 1e-12), label = nm)
    if (row$fixed) {
      # A row's joint probabilities with every row, itself included, add up
      # to n times its own when every draw has n rows.
      expect_equal(rowSums(m), row$n * p, label = nm)
    } else {
      # Rows of a design of random size are taken independently.
      want <- outer(p, p)
      diag(want) <- p
      expect_equal(m, want, label = nm)
    }
  }
})

test_that("ht_total() gives the Horvitz-Thompson total under every variance", {
  for (nm in names(rows)) {
    row <- rows[[nm]]
    s <- grid_sample(row)
    if (is.null(s)) next
    auto <- ht_total(s, "y")
    expect_equal(auto$total, sum(s$y * s$.weight), label = nm)
    expect_equal(auto$method, row$method, label = nm)
    expect_true(is.finite(auto$se) && auto$se > 0, label = nm)
    expect_equal(ht_total(s, "y", variance = "analytic")$se, auto$se,
                 label = nm)
    expect_true(is.na(ht_total(s, "y", variance = "none")$se), label = nm)

    jk <- ht_total(s, "y", variance = "jackknife")
    if (identical(row$jackknife, "same")) {
      expect_equal(jk$se, auto$se, label = nm)
    } else if (identical(row$jackknife, "own")) {
      expect_equal(jk$method, "jackknife", label = nm)
      expect_gt(jk$se / auto$se, 2 / 3, label = nm)
      expect_lt(jk$se / auto$se, 3 / 2, label = nm)
    } else {
      expect_true(is.na(jk$se), label = nm)
      expect_equal(jk$method, "none", label = nm)
      expect_match(jk$note, row$jackknife, fixed = TRUE, label = nm)
    }
  }
})

test_that("ht_mean() is the Hajek mean, with a logit interval for a 0/1 variable", {
  for (nm in names(rows)) {
    row <- rows[[nm]]
    s <- grid_sample(row)
    if (is.null(s)) next
    m <- ht_mean(s, "y")
    expect_equal(m$mean, sum(s$y * s$.weight) / sum(s$.weight), label = nm)
    expect_true(is.finite(m$se) && m$se > 0, label = nm)

    b <- ht_mean(s, "flag")
    expect_equal(b$ci_scale, "logit", label = nm)
    half <- stats::qt(0.975, b$df) * b$se / (b$mean * (1 - b$mean))
    expect_equal(b$ci, stats::plogis(stats::qlogis(b$mean) + c(-1, 1) * half),
                 label = nm)
  }
})

test_that("domain totals add up to the whole sample's", {
  for (nm in names(rows)) {
    row <- rows[[nm]]
    s <- grid_sample(row)
    if (is.null(s)) next
    by <- ht_total(s, "y", by = "site")
    expect_equal(sum(by$total), ht_total(s, "y")$total, label = nm)
    expect_equal(sum(by$n), nrow(s), label = nm)
    expect_true(all(is.finite(by$se) & by$se > 0), label = nm)
  }
})

test_that("degrees of freedom are the units below certainty less their strata", {
  for (nm in names(rows)) {
    row <- rows[[nm]]
    s <- grid_sample(row)
    if (is.null(s)) next
    want <- reference_df(row, s)
    expect_equal(ht_total(s, "y")$df, want, label = nm)
    expect_equal(ht_mean(s, "y")$df, want, label = nm)
    expect_equal(attr(ht_total(s, "y", by = "site"), "df"), want, label = nm)
  }
})

test_that("survey gets the total, the degrees of freedom and the same estimator", {
  skip_if_not_installed("survey")
  for (nm in names(rows)) {
    row <- rows[[nm]]
    s <- grid_sample(row)
    if (is.null(s)) next
    ours <- ht_total(s, "y")
    sv <- as_svydesign(s)
    theirs <- survey::svytotal(~y, sv)
    expect_equal(as.numeric(theirs), ours$total, label = nm)
    expect_equal(survey::degf(sv), ours$df, label = nm)
    se <- as.numeric(survey::SE(theirs))
    if (row$survey == "same") {
      expect_equal(se, ours$se, tolerance = 1e-6, label = nm)
    } else if (row$survey == "close") {
      # Brewer's approximation in survey, Deville's here.
      expect_equal(se, ours$se, tolerance = 0.01, label = nm)
    } else {
      # Successive differences and the local mean have no survey estimator;
      # survey treats these samples as simple random ones.
      expect_true(is.finite(se) && se > 0, label = nm)
    }
  }
})

test_that("sample_summary() counts the sample and the rows the design reaches", {
  for (nm in names(rows)) {
    row <- rows[[nm]]
    s <- grid_sample(row)
    if (is.null(s)) next
    x <- sample_summary(s)
    expect_equal(x$n, nrow(s), label = nm)
    expect_equal(x$N, nrow(row$frame), label = nm)
    expect_equal(x$n_certain, sum(s$.prob >= 1 - 1e-9), label = nm)
    if (row$conditional) {
      # Every row of a cluster left out could have been drawn, at a
      # probability the clusters selected do not give.
      expect_equal(x$N_reachable, nrow(row$frame), label = nm)
      expect_true(is.na(x$expected_n), label = nm)
      expect_equal(x$conditional, length(unique(s$cl)), label = nm)
    } else {
      p <- quietly(row, inclusion_prob(row$frame, row$design))
      expect_equal(x$N_reachable, sum(p > 0), label = nm)
      expect_equal(x$expected_n, row$expected, label = nm)
      expect_true(is.na(x$conditional), label = nm)
    }
  }
})

test_that("every view of every design plots, or says why it cannot", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (nm in names(rows)) {
    row <- rows[[nm]]
    for (type in c("selection", "probability", "map")) {
      path <- if (!is.null(row$refuse$plot)) "plot" else type
      out <- cell(row, path,
                  plot(row$design, row$frame, type = type, seed = 1,
                       coords = c("x", "yy")),
                  paste(nm, type))
      if (!is.null(out)) expect_identical(out, row$design, label = nm)
    }
  }
})

test_that("deff() reads the design effect under drawn's generic and survey's", {
  survey <- requireNamespace("survey", quietly = TRUE)
  for (nm in names(rows)) {
    row <- rows[[nm]]
    s <- grid_sample(row)
    if (is.null(s)) next
    r <- ht_total(s, "y")
    m <- ht_mean(s, "y")
    expect_gt(r$deff, 0, label = nm)
    expect_equal(deff(r), r$deff, label = nm)
    expect_equal(deff(m), m$deff, label = nm)
    if (survey) expect_equal(survey::deff(r), r$deff, label = nm)
  }
})
