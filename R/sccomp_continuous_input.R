#' Priors for continuous-outcome sccomp
#'
#' @param location_intercept,location_coefficient Normal prior mean and SD for
#'   location intercepts and other coefficients, respectively.
#' @param log_variance_sd Normal prior SD for log variance coefficients.
#' @param log_s0,log_prior_df,log_observation_df Normal prior mean and SD for
#'   log baseline SD and log(df - 2) of variance and observation distributions.
#' @param trend_sd Normal prior SD for the optional variance trend.
#' @param random_log_sd Normal prior mean and SD for random-SD hypermeans.
#' @param random_log_sd_sd Half-normal SD for random-SD hyperdispersion.
#' @param lkj LKJ correlation shape, at least one.
#' @return A validated continuous prior object. Defaults assume global scaling.
#' @export
sccomp_continuous_prior <- function(
    location_intercept = c(0, 5), location_coefficient = c(0, 2.5),
    log_variance_sd = 1, log_s0 = c(0, 1),
    log_prior_df = c(log(8), 1), log_observation_df = c(log(8), 1),
    trend_sd = 0.5, random_log_sd = c(log(0.5), 1),
    random_log_sd_sd = 0.5, lkj = 2) {
  prior <- list(
    location_intercept = location_intercept,
    location_coefficient = location_coefficient,
    log_variance_sd = log_variance_sd, log_s0 = log_s0,
    log_prior_df = log_prior_df, log_observation_df = log_observation_df,
    trend_sd = trend_sd, random_log_sd = random_log_sd,
    random_log_sd_sd = random_log_sd_sd, lkj = lkj
  )
  validate_continuous_prior(prior)
  structure(prior, class = "sccomp_continuous_prior")
}

validate_continuous_prior <- function(prior) {
  pairs <- c("location_intercept", "location_coefficient", "log_s0",
             "log_prior_df", "log_observation_df", "random_log_sd")
  sds <- c("log_variance_sd", "trend_sd", "random_log_sd_sd")
  if (!is.list(prior) || !setequal(names(prior), c(pairs, sds, "lkj"))) {
    stop("Supply a complete prior from sccomp_continuous_prior().", call. = FALSE)
  }
  for (name in pairs) {
    x <- prior[[name]]
    if (!is.numeric(x) || length(x) != 2L || any(!is.finite(x)) || x[2] <= 0) {
      stop("Prior `", name, "` must be a finite mean and positive SD.",
           call. = FALSE)
    }
  }
  for (name in sds) {
    x <- prior[[name]]
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
      stop("Prior `", name, "` must be one finite positive SD.", call. = FALSE)
    }
  }
  if (!is.numeric(prior$lkj) || length(prior$lkj) != 1L ||
      !is.finite(prior$lkj) || prior$lkj < 1) {
    stop("Prior `lkj` must be finite and at least one.", call. = FALSE)
  }
  invisible(prior)
}

continuous_design <- function(formula, data, variance = FALSE) {
  mf <- stats::model.frame(formula, data, na.action = stats::na.fail)
  factors <- names(mf)[vapply(mf, function(x) is.factor(x) || is.character(x),
                             logical(1))]
  for (name in factors) {
    if (is.character(mf[[name]])) mf[[name]] <- factor(mf[[name]])
    if (nlevels(mf[[name]]) < 2L) {
      stop("Covariate `", name, "` has fewer than two levels.", call. = FALSE)
    }
  }
  contrasts <- lapply(mf[factors], function(x) {
    if (variance) stats::contr.sum(levels(x)) else {
      stats::contr.treatment(levels(x))
    }
  })
  x <- stats::model.matrix(formula, mf, contrasts.arg = contrasts)
  center <- rep(0, ncol(x))
  if (variance) {
    x <- x[, colnames(x) != "(Intercept)", drop = FALSE]
    center <- colMeans(x)
    x <- sweep(x, 2L, center, "-")
  }
  if (any(!is.finite(x)) || (ncol(x) && qr(x)$rank != ncol(x))) {
    stop("The ", if (variance) "variance" else "location",
         " design must be finite and full rank.", call. = FALSE)
  }
  list(matrix = x, terms = stats::terms(mf),
       xlevels = lapply(mf[factors], levels), contrasts = contrasts,
       center = center)
}

continuous_model_input <- function(.data, formula_location, formula_variance,
                                   sample, feature, value, family,
                                   variance_trend, scaling, student_t_df, prior,
                                   use_data = TRUE, grainsize = 1L) {
  validate_continuous_prior(prior)
  if (length(unique(c(sample, feature, value))) != 3L) {
    stop("Sample, feature and value must name different columns.", call. = FALSE)
  }
  if (!is.logical(use_data) || length(use_data) != 1L || is.na(use_data) ||
      !is.numeric(grainsize) || length(grainsize) != 1L ||
      !is.finite(grainsize) || grainsize < 1 || grainsize != as.integer(grainsize)) {
    stop("use_data must be logical and grainsize a positive integer.",
         call. = FALSE)
  }
  if (!is.data.frame(.data) || !nrow(.data)) {
    stop("Continuous input must be a nonempty data frame.", call. = FALSE)
  }
  if (!inherits(formula_location, "formula") || length(formula_location) != 2L ||
      !inherits(formula_variance, "formula") || length(formula_variance) != 2L) {
    stop("Location and variance formulas must be one-sided.", call. = FALSE)
  }
  if (any(grepl("|", attr(stats::terms(formula_variance), "term.labels"),
                fixed = TRUE))) {
    stop("Random effects are unavailable in formula_variance.", call. = FALSE)
  }
  required <- unique(c(sample, feature, value, all.vars(formula_location),
                       all.vars(formula_variance)))
  if (!all(required %in% names(.data))) {
    stop("Missing input columns: ", paste(setdiff(required, names(.data)),
                                          collapse = ", "), call. = FALSE)
  }
  if (anyNA(.data[required]) || !is.numeric(.data[[value]]) ||
      any(!is.finite(.data[[value]]))) {
    stop("Required columns must be observed and values finite numeric.",
         call. = FALSE)
  }
  if (anyDuplicated(.data[c(sample, feature)])) {
    stop("Duplicate sample-feature pairs are unsupported.", call. = FALSE)
  }
  samples <- unique(as.character(.data[[sample]]))
  features <- unique(as.character(.data[[feature]]))
  if (any(!nzchar(samples)) || any(!nzchar(features))) {
    stop("Sample and feature identifiers must be nonempty.", call. = FALSE)
  }
  n <- length(samples)
  g <- length(features)
  if (g < 2L) stop("At least two features are required.", call. = FALSE)
  if (g < 10L) warning("Fewer than ten features: variance learning may be weak.",
                       call. = FALSE)
  if (nrow(.data) != n * g) {
    stop("Require one observed value for every sample-feature pair.",
         call. = FALSE)
  }
  sample_index <- match(as.character(.data[[sample]]), samples)
  for (name in setdiff(required, c(sample, feature, value))) {
    groups <- split(.data[[name]], sample_index)
    if (any(vapply(groups, function(x) length(unique(x)) != 1L, logical(1)))) {
      stop("Covariate `", name, "` must be constant within sample.", call. = FALSE)
    }
  }
  sample_data <- .data[match(samples, as.character(.data[[sample]])),
                       unique(c(sample, setdiff(required, c(feature, value)))),
                       drop = FALSE]
  rownames(sample_data) <- samples
  y <- matrix(NA_real_, n, g, dimnames = list(samples, features))
  y[cbind(sample_index, match(as.character(.data[[feature]]), features))] <-
    .data[[value]]
  center <- 0
  scale <- 1
  if (scaling == "global") {
    center <- stats::median(y)
    mads <- apply(y, 2L, stats::mad, constant = 1.4826)
    positive <- mads[is.finite(mads) & mads > 0]
    if (!length(positive)) {
      stop("Global scaling requires a positive within-feature MAD.", call. = FALSE)
    }
    scale <- stats::median(positive)
  }
  y <- (y - center) / scale
  if (any(!is.finite(y))) {
    stop("The global transformation produced nonfinite values.", call. = FALSE)
  }
  location_terms <- stats::terms(formula_location)
  labels <- attr(location_terms, "term.labels")
  fixed_labels <- labels[!grepl("|", labels, fixed = TRUE)]
  fixed_formula <- if (length(fixed_labels)) {
    stats::reformulate(fixed_labels,
                      intercept = attr(location_terms, "intercept"),
                      env = environment(formula_location))
  } else stats::as.formula(
    if (attr(location_terms, "intercept")) "~ 1" else "~ 0",
    env = environment(formula_location)
  )
  variance_labels <- attr(stats::terms(formula_variance), "term.labels")
  variance_formula <- if (length(variance_labels)) {
    stats::reformulate(variance_labels, env = environment(formula_variance))
  } else stats::as.formula("~ 1", env = environment(formula_variance))
  location <- continuous_design(fixed_formula, sample_data)
  if (!ncol(location$matrix)) {
    stop("The location design needs at least one fixed coefficient.",
         call. = FALSE)
  }
  variance <- continuous_design(variance_formula, sample_data, variance = TRUE)
  trend <- rep(0, g)
  if (variance_trend == "linear") {
    means <- colMeans(y)
    if (!is.finite(stats::sd(means)) || stats::sd(means) == 0) {
      stop("Linear variance trend requires differing feature means.", call. = FALSE)
    }
    trend <- as.numeric(base::scale(means))
  }
  estimate_df <- identical(student_t_df, "estimate") && family == "student_t"
  if (!identical(student_t_df, "estimate") &&
      (!is.numeric(student_t_df) || length(student_t_df) != 1L ||
       !is.finite(student_t_df) || student_t_df <= 2)) {
    stop("student_t_df must be 'estimate' or one finite number greater than 2.",
         call. = FALSE)
  }
  input <- list(
    N = n, G = g, C = ncol(location$matrix), A = ncol(variance$matrix),
    y = y, X = location$matrix, W = variance$matrix,
    location_intercept = as.integer(colnames(location$matrix) == "(Intercept)"),
    student_t = as.integer(family == "student_t"),
    estimate_df = as.integer(estimate_df),
    fixed_df = if (is.numeric(student_t_df)) student_t_df else 4,
    has_trend = as.integer(variance_trend == "linear"), trend_covariate = trend,
    use_data = as.integer(use_data), grainsize = as.integer(grainsize),
    prior_intercept = prior$location_intercept,
    prior_coefficient = prior$location_coefficient,
    prior_variance_sd = prior$log_variance_sd, prior_log_s0 = prior$log_s0,
    prior_eta0 = prior$log_prior_df, prior_eta_obs = prior$log_observation_df,
    prior_trend_sd = prior$trend_sd, prior_mu_omega = prior$random_log_sd,
    prior_tau_omega_sd = prior$random_log_sd_sd, prior_lkj = prior$lkj
  )
  clauses <- formula_to_random_effect_formulae(formula_location)
  if (nrow(clauses) > 4L) stop("At most four random clauses are supported.",
                               call. = FALSE)
  random <- vector("list", nrow(clauses))
  for (b in seq_len(4L)) {
    z <- matrix(numeric(), n, 0L)
    levels <- character()
    group <- rep(1L, n)
    if (b <= nrow(clauses)) {
      grouping <- clauses$grouping[b]
      if (!grouping %in% names(sample_data)) {
        stop("A random grouping term must name one input column.", call. = FALSE)
      }
      levels <- unique(as.character(sample_data[[grouping]]))
      if (length(levels) < 2L) stop("Random effects need two grouping levels.",
                                   call. = FALSE)
      if (length(levels) < 5L) warning("Fewer than five random grouping levels.",
                                      call. = FALSE)
      group <- match(as.character(sample_data[[grouping]]), levels)
      random_formula <- clauses$formula[[b]]
      environment(random_formula) <- environment(formula_location)
      design <- continuous_design(random_formula, sample_data)
      z <- design$matrix
      if (!ncol(z)) stop("A random clause must contain a coefficient.",
                         call. = FALSE)
      slopes <- which(colnames(z) != "(Intercept)")
      for (q in slopes) {
        varying <- vapply(split(z[, q], group), function(x) {
          length(unique(x)) > 1L
        }, logical(1))
        if (!any(varying)) {
          stop("Random slope `", colnames(z)[q],
               "` has no within-group variation.", call. = FALSE)
        }
      }
      random[[b]] <- list(grouping = grouping, levels = levels,
                           design = design, formula = random_formula)
    }
    input[[paste0("K", b)]] <- ncol(z)
    input[[paste0("J", b)]] <- length(levels)
    input[[paste0("Z", b)]] <- z
    input[[paste0("group", b)]] <- group
  }
  continuous_check_covariance_design(input)
  list(stan_data = input, sample_data = sample_data, original_data = .data,
       samples = samples, features = features, sample = sample, feature = feature,
       value = value, location = location, variance = variance, random = random,
       fixed_formula = fixed_formula, variance_formula = variance_formula,
       global_center = center, global_scale = scale, prior = prior)
}

# Necessary identifiability check: distinguish random covariance from residual variance.
continuous_check_covariance_design <- function(input) {
  blocks <- which(vapply(1:4, function(b) input[[paste0("K", b)]] > 0L, logical(1)))
  if (!length(blocks)) return(invisible(NULL))
  pairs <- cbind(seq_len(input$N), seq_len(input$N))
  for (b in blocks) {
    groups <- split(seq_len(input$N), input[[paste0("group", b)]])
    for (indices in groups) if (length(indices) >= 2L) {
      pairs <- rbind(pairs, t(utils::combn(indices, 2L)))
    }
  }
  pairs <- unique(pairs)
  diagonal <- pairs[, 1L] == pairs[, 2L]
  residual <- cbind(1, input$W)
  basis <- residual[pairs[, 1L], , drop = FALSE] * diagonal
  labels <- c("residual variance", colnames(input$W))
  for (b in blocks) {
    z <- input[[paste0("Z", b)]]
    group <- input[[paste0("group", b)]]
    same_group <- group[pairs[, 1L]] == group[pairs[, 2L]]
    for (a in seq_len(ncol(z))) for (q in seq.int(a, ncol(z))) {
      column <- z[pairs[, 1L], a] * z[pairs[, 2L], q]
      if (a != q) column <- column + z[pairs[, 1L], q] * z[pairs[, 2L], a]
      basis <- cbind(basis, column * same_group)
      labels <- c(labels, sprintf("random clause %d covariance %s/%s", b, colnames(z)[a], colnames(z)[q]))
    }
  }
  norms <- sqrt(colSums(basis^2))
  if (any(!is.finite(basis)) || any(!is.finite(norms)) || any(norms == 0)) {
    stop("Random/residual covariance design is not identifiable or finite.", call. = FALSE)
  }
  decomposition <- qr(sweep(basis, 2L, norms, "/"), tol = 1e-8)
  if (decomposition$rank < ncol(basis)) {
    dependent <- labels[decomposition$pivot[seq.int(decomposition$rank + 1L, ncol(basis))]]
    stop("Random/residual covariance design is not identifiable: ",
         paste(dependent, collapse = ", "),
         ". Add repeated observations or simplify random slopes/variance terms.", call. = FALSE)
  }
  invisible(NULL)
}
