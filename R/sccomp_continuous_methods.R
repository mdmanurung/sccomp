continuous_fit_parts <- function(fit) {
  metadata <- attr(fit, "continuous")
  draws <- attr(fit, "fit")
  if (is.null(metadata$input) || !inherits(draws, "draws_array")) {
    stop("This continuous table has no retained model metadata/posterior draws.",
         call. = FALSE)
  }
  list(metadata = metadata, input = metadata$input, draws = draws)
}

#' Test continuous effects in original input units
#'
#' @param .data A continuous fit.
#' @param contrasts Location expressions, or a list with `location` and/or
#'   `log_variance` expression vectors. Names label contrasts. Only linear
#'   arithmetic of coefficient names is accepted; backtick nonstandard names.
#'   NULL tests location coefficients except the intercept and, when requested,
#'   log-variance coefficients. Variance coefficients use centered sum coding;
#'   inspect the fit's parameter names before specifying contrasts.
#' @param percent_false_positive Not applicable to continuous fits; retained
#'   only for S3 consistency with the `sccomp_test` generic, and supplying it
#'   is an error. Use `effect_threshold` instead.
#' @param test_composition_above_logit_fold_change Not applicable to continuous
#'   fits; retained only for S3 consistency with the `sccomp_test` generic, and
#'   supplying it is an error. Use `effect_threshold` instead.
#' @param pass_fit Retain portable posterior draws in the tested table.
#' @param ... Unused; passing arguments through `...` is an error.
#' @param effect_threshold Required positive minimum absolute location effect
#'   in original input units.
#' @param variance_ratio_threshold Optional ratio greater than one; the
#'   corresponding log-variance threshold is its natural logarithm.
#' @param posterior_fdr Posterior expected false discovery probability cutoff.
#' @param allow_exploratory Explicit opt-in for variational selection.
#' @return Component/contrast/feature statistics. `p_null` is posterior mass
#'   within the threshold; `p_error` includes that mass and the opposite
#'   direction. `posterior_fdr` is cumulative mean `p_error` across features within each
#'   component and contrast; `selected` applies the requested cutoff. This is a
#'   ties-conservative rule: equal `p_error` values share the cumulative mean at
#'   the end of their tied group, so feature order does not determine selection.
#'   posterior expected error rule, not a universal frequentist FDR guarantee.
#' @export
sccomp_test.sccomp_continuous_tbl <- function(
    .data, contrasts = NULL, percent_false_positive = 5,
    test_composition_above_logit_fold_change = 0.1, pass_fit = TRUE,
    ..., effect_threshold, variance_ratio_threshold = NULL,
    posterior_fdr = .05, allow_exploratory = FALSE) {
  if (!missing(percent_false_positive) ||
      !missing(test_composition_above_logit_fold_change)) {
    stop("sccomp says: continuous fits use `effect_threshold` and `variance_ratio_threshold`, not `percent_false_positive` and `test_composition_above_logit_fold_change`.",
         call. = FALSE)
  }
  if (length(list(...))) stop("Unused continuous test arguments.", call. = FALSE)
  if (!is.logical(pass_fit) || length(pass_fit) != 1L || is.na(pass_fit)) {
    stop("pass_fit must be logical.", call. = FALSE)
  }
  parts <- continuous_fit_parts(.data)
  meta <- parts$metadata
  if (missing(effect_threshold) || !is.numeric(effect_threshold) ||
      length(effect_threshold) != 1L || !is.finite(effect_threshold) ||
      effect_threshold <= 0) {
    stop("Supply a finite positive effect_threshold in original input units.",
         call. = FALSE)
  }
  if (!is.null(variance_ratio_threshold) &&
      (!is.numeric(variance_ratio_threshold) || length(variance_ratio_threshold) != 1L ||
       !is.finite(variance_ratio_threshold) || variance_ratio_threshold <= 1)) {
    stop("variance_ratio_threshold must be finite and greater than one.", call. = FALSE)
  }
  if (!is.numeric(posterior_fdr) || length(posterior_fdr) != 1L ||
      !is.finite(posterior_fdr) || posterior_fdr <= 0 || posterior_fdr >= 1) {
    stop("posterior_fdr must be between zero and one.", call. = FALSE)
  }
  if (!is.logical(allow_exploratory) || length(allow_exploratory) != 1L ||
      is.na(allow_exploratory)) stop("allow_exploratory must be logical.", call. = FALSE)
  if (identical(meta$inference_status, "prior_only")) {
    stop("Prior-only fits cannot be used for discovery.", call. = FALSE)
  }
  if (meta$inference_method == "variational") {
    if (!allow_exploratory) stop("VI selection requires allow_exploratory = TRUE.",
                                call. = FALSE)
  } else if (!isTRUE(meta$diagnostics$hmc_qualified)) {
    stop("HMC diagnostics have not passed the continuous inference gates.",
         call. = FALSE)
  }
  if (is.character(contrasts)) contrasts <- list(location = contrasts)
  if (!is.null(contrasts) && (!is.list(contrasts) ||
      is.null(names(contrasts)) || anyDuplicated(names(contrasts)) ||
      !all(names(contrasts) %in% c("location", "log_variance")))) {
    stop("contrasts must be expressions or a component-named list.", call. = FALSE)
  }
  if (!is.null(contrasts$log_variance) && is.null(variance_ratio_threshold)) {
    stop("Variance contrasts require variance_ratio_threshold.", call. = FALSE)
  }
  rows <- list()
  for (component in c("location", "log_variance")) {
    if (component == "log_variance" && is.null(variance_ratio_threshold)) next
    parameters <- if (component == "location") colnames(parts$input$stan_data$X) else {
      colnames(parts$input$stan_data$W)
    }
    expressions <- contrasts[[component]]
    if (is.null(expressions)) {
      expressions <- parameters[parameters != "(Intercept)"]
      names(expressions) <- expressions
      expressions <- paste0("`", expressions, "`") |>
        stats::setNames(names(expressions))
    }
    if (!is.character(expressions) || anyNA(expressions) || any(!nzchar(expressions))) {
      stop("Each contrast must be a nonempty linear expression.", call. = FALSE)
    }
    labels <- names(expressions)
    if (is.null(labels)) labels <- expressions
    if (any(!nzchar(labels)) || anyDuplicated(labels)) {
      stop("Contrast names must be nonempty and unique per component.", call. = FALSE)
    }
    threshold <- if (component == "location") effect_threshold else log(variance_ratio_threshold)
    for (j in seq_along(expressions)) {
      weights <- continuous_contrast_weights(expressions[j], parameters)
      for (g in seq_along(parts$input$features)) {
        prefix <- if (component == "location") "beta" else "gamma"
        values <- matrix(0, dim(parts$draws)[1L], dim(parts$draws)[2L])
        for (q in which(weights != 0)) {
          x <- unclass(parts$draws)[, , sprintf("%s[%d,%d]", prefix, q, g)]
          if (component == "location") {
            x <- x * parts$input$global_scale + parts$input$global_center *
              parts$input$coefficient_offset[q]
          }
          values <- values + weights[q] * x
        }
        if (any(!is.finite(values))) {
          stop("Contrast values must be finite; check coefficient weights.", call. = FALSE)
        }
        diagnostics <- c(rhat = posterior::rhat(values),
                         ess_bulk = posterior::ess_bulk(values),
                         ess_tail = posterior::ess_tail(values))
        if (meta$inference_method == "hmc" && sum(weights != 0) > 1L &&
            (any(!is.finite(diagnostics)) || diagnostics["rhat"] >= 1.01 ||
             diagnostics["ess_bulk"] < 400 || diagnostics["ess_tail"] < 400)) {
          stop("HMC diagnostics failed for requested contrast `", labels[j],
               "` in feature `", parts$input$features[g], "`.", call. = FALSE)
        }
        positive <- mean(values > threshold)
        negative <- mean(values < -threshold)
        direction <- if (positive >= negative) "positive" else "negative"
        rows[[length(rows) + 1L]] <- data.frame(
          feature = parts$input$features[g], component = component,
          parameter = labels[j], estimate = mean(values),
          lower = unname(stats::quantile(values, .025)),
          upper = unname(stats::quantile(values, .975)),
          rhat = unname(diagnostics["rhat"]), ess_bulk = unname(diagnostics["ess_bulk"]),
          ess_tail = unname(diagnostics["ess_tail"]),
          threshold = threshold, direction = direction,
          p_null = mean(abs(values) <= threshold),
          p_false_sign = if (direction == "positive") mean(values <= 0) else mean(values >= 0),
          p_error = 1 - max(positive, negative),
          probability_mcse = if (meta$inference_method == "hmc") {
            posterior::mcse_mean(if (direction == "positive") values > threshold else values < -threshold)
          } else NA_real_
        )
      }
    }
  }
  if (!length(rows)) stop("No location or variance effects are available to test.", call. = FALSE)
  result <- tibble::as_tibble(do.call(rbind, rows))
  result$posterior_fdr <- NA_real_
  families <- interaction(result$component, result$parameter, drop = TRUE)
  for (indices in split(seq_len(nrow(result)), families)) {
    ordered <- indices[order(result$p_error[indices])]
    probability <- result$p_error[ordered]
    cumulative <- cumsum(probability) / seq_along(ordered)
    # Equal posterior error probabilities must not select features by input order.
    ends <- which(!duplicated(probability, fromLast = TRUE))
    result$posterior_fdr[ordered] <- cumulative[ends][match(probability, probability[ends])]
  }
  result$selected <- result$posterior_fdr <= posterior_fdr
  result$inference_status <- meta$inference_status
  class(result) <- c("sccomp_continuous_tbl", "sccomp_tbl", class(result))
  meta$testing <- list(contrasts = contrasts, effect_threshold = effect_threshold,
                        variance_ratio_threshold = variance_ratio_threshold,
                        posterior_fdr = posterior_fdr, allow_exploratory = allow_exploratory)
  attr(result, "continuous") <- meta
  attr(result, "inference_method") <- meta$inference_method
  attr(result, "fit") <- if (pass_fit) parts$draws else NULL
  result
}

continuous_contrast_weights <- function(expression, parameters) {
  parsed <- tryCatch(parse(text = expression, keep.source = FALSE),
                     error = function(e) stop("Invalid contrast expression.", call. = FALSE))
  if (length(parsed) != 1L) stop("Supply one expression per contrast.", call. = FALSE)
  zero <- stats::setNames(rep(0, length(parameters) + 1L), c(parameters, ".constant"))
  walk <- function(node) {
    if (is.numeric(node) && length(node) == 1L && is.finite(node)) {
      x <- zero
      x[length(x)] <- node
      return(x)
    }
    if (is.symbol(node)) {
      name <- as.character(node)
      if (!name %in% parameters) stop("Unknown contrast coefficient: ", name, call. = FALSE)
      x <- zero
      x[name] <- 1
      return(x)
    }
    if (!is.call(node)) stop("Only linear arithmetic is accepted.", call. = FALSE)
    operator <- as.character(node[[1L]])
    if (!operator %in% c("(", "+", "-", "*", "/") || length(node) > 3L) {
      stop("Only linear arithmetic is accepted.", call. = FALSE)
    }
    a <- walk(node[[2L]])
    if (length(node) == 2L) {
      if (operator %in% c("(", "+")) return(a)
      if (operator == "-") return(-a)
      stop("Invalid contrast operator.", call. = FALSE)
    }
    b <- walk(node[[3L]])
    if (operator == "+") return(a + b)
    if (operator == "-") return(a - b)
    if (operator == "*" && all(a[-length(a)] == 0)) return(a[length(a)] * b)
    if (operator == "*" && all(b[-length(b)] == 0)) return(b[length(b)] * a)
    if (operator == "/" && all(b[-length(b)] == 0) && b[length(b)] != 0) return(a / b[length(b)])
    stop("Contrasts must be linear; division requires a nonzero constant.", call. = FALSE)
  }
  weights <- walk(parsed[[1L]])
  if (weights[length(weights)] != 0 || any(!is.finite(weights)) || all(weights == 0)) {
    stop("Contrasts require finite nonzero coefficient weights and no constant offset.", call. = FALSE)
  }
  weights[-length(weights)]
}

continuous_new_design <- function(design, data, variance = FALSE) {
  mf <- stats::model.frame(design$terms, data, xlev = design$xlevels,
                           na.action = stats::na.fail)
  x <- stats::model.matrix(design$terms, mf, contrasts.arg = design$contrasts)
  if (variance) x <- x[, colnames(x) != "(Intercept)", drop = FALSE]
  if (!identical(colnames(x), colnames(design$matrix)) || any(!is.finite(x))) {
    stop("New data do not reproduce the fitted design columns.", call. = FALSE)
  }
  if (variance) x <- sweep(x, 2L, design$center, "-")
  x
}

continuous_term_mask <- function(formula, fitted_formula, design, random = FALSE) {
  if (is.null(formula)) return(rep(TRUE, ncol(design)))
  if (!inherits(formula, "formula") || length(formula) != 2L) {
    stop("Prediction formulas must be one-sided subsets of fitted terms.", call. = FALSE)
  }
  terms <- stats::terms(formula)
  requested <- attr(terms, "term.labels")
  requested <- requested[!grepl("|", requested, fixed = TRUE)]
  fitted <- attr(stats::terms(fitted_formula), "term.labels")
  if (!all(requested %in% fitted) ||
      (!random && any(grepl("|", attr(terms, "term.labels"), fixed = TRUE)))) {
    stop("Prediction formulas may contain only fitted terms.", call. = FALSE)
  }
  assignments <- attr(design, "assign")
  assignments %in% c(if (attr(terms, "intercept") == 1L) 0L,
                     match(requested, fitted))
}

continuous_predict_draws <- function(fit, new_data, number_of_draws, seed,
                                     summary = FALSE, generate = TRUE,
                                     formula_location = NULL, formula_variance = NULL) {
  parts <- continuous_fit_parts(fit)
  input <- parts$input
  if (!is.numeric(number_of_draws) || length(number_of_draws) != 1L ||
      !is.finite(number_of_draws) || number_of_draws < 1 || number_of_draws != floor(number_of_draws) ||
      !is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || seed < 1 || seed != floor(seed)) {
    stop("number_of_draws and mcmc_seed must be positive integers.", call. = FALSE)
  }
  data <- if (is.null(new_data)) input$sample_data else new_data
  if (!is.data.frame(data) || !nrow(data) || !input$sample %in% names(data) ||
      anyNA(data[[input$sample]]) || anyDuplicated(data[[input$sample]]) ||
      any(!nzchar(as.character(data[[input$sample]])))) {
    stop("New data require one row per nonempty sample identifier.", call. = FALSE)
  }
  x <- continuous_new_design(input$location, data)
  w <- continuous_new_design(input$variance, data, variance = TRUE)
  location_mask <- continuous_term_mask(formula_location, input$fixed_formula,
                                         input$location$matrix, random = TRUE)
  x[, !location_mask] <- 0
  # Centering/sweeping strips assign; reconstruct it from the fitted model terms.
  variance_design <- stats::model.matrix(input$variance$terms,
    stats::model.frame(input$variance$terms, input$sample_data,
                       xlev = input$variance$xlevels), contrasts.arg = input$variance$contrasts)
  variance_mask <- continuous_term_mask(formula_variance, input$variance_formula,
                                         variance_design)
  variance_mask <- variance_mask[colnames(variance_design) != "(Intercept)"]
  w[, !variance_mask] <- 0
  random_blocks <- seq_along(input$random)
  if (!is.null(formula_location)) {
    requested <- formula_to_random_effect_formulae(formula_location)
    keys <- vapply(input$random, function(block) {
      paste(block$grouping, paste(deparse(block$formula), collapse = ""))
    }, character(1))
    requested_keys <- if (nrow(requested)) vapply(seq_len(nrow(requested)), function(b) {
      paste(requested$grouping[b], paste(deparse(requested$formula[[b]]), collapse = ""))
    }, character(1)) else character()
    random_blocks <- match(requested_keys, keys)
    if (anyNA(random_blocks) || anyDuplicated(random_blocks)) {
      stop("Prediction random clauses must match fitted clauses.", call. = FALSE)
    }
  }
  random <- lapply(random_blocks, function(b) {
    block <- input$random[[b]]
    group <- match(as.character(data[[block$grouping]]), block$levels)
    if (length(group) != nrow(data) || anyNA(group)) {
      stop("Prediction requires observed random-effect grouping levels.", call. = FALSE)
    }
    list(z = continuous_new_design(block$design, data), group = group, block = b)
  })
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  draws <- posterior::as_draws_matrix(parts$draws)
  indices <- sort(sample.int(nrow(draws), min(number_of_draws, nrow(draws))))
  draws <- draws[indices, , drop = FALSE]
  rows <- vector("list", length(input$features))
  for (g in seq_along(input$features)) {
    mu <- x %*% t(draws[, sprintf("beta[%d,%d]", seq_len(ncol(x)), g), drop = FALSE])
    for (b in seq_along(random)) {
      for (q in seq_len(ncol(random[[b]]$z))) {
        coefficients <- draws[, sprintf("random_effect_%d[%d,%d,%d]",
                                         random[[b]]$block, g, q, random[[b]]$group), drop = FALSE]
        mu <- mu + random[[b]]$z[, q] * t(coefficients)
      }
    }
    log_v <- matrix(log(draws[, sprintf("sigma_sq[%d]", g)]),
                     nrow(data), nrow(draws), byrow = TRUE)
    if (ncol(w)) {
      log_v <- log_v + w %*% t(draws[, sprintf("gamma[%d,%d]", seq_len(ncol(w)), g), drop = FALSE])
    }
    variance <- exp(log_v) * input$global_scale^2
    location <- as.numeric(x %*% input$coefficient_offset) * input$global_center +
      input$global_scale * mu
    generated <- NULL
    if (generate) {
      noise <- if (parts$metadata$family == "gaussian") stats::rnorm(length(mu)) else {
        nu <- rep(draws[, "nu_obs"], each = nrow(data))
        stats::rt(length(mu), nu) * sqrt((nu - 2) / nu)
      }
      generated <- location + sqrt(variance) * noise
    }
    if (summary) {
      rows[[g]] <- data.frame(
        sample = as.character(data[[input$sample]]), feature = input$features[g],
        location = rowMeans(location), residual_variance = rowMeans(variance)
      )
      if (generate) {
        rows[[g]]$generated_value <- rowMeans(generated)
        rows[[g]]$lower <- apply(generated, 1L, stats::quantile, .025)
        rows[[g]]$upper <- apply(generated, 1L, stats::quantile, .975)
      }
    } else {
      rows[[g]] <- data.frame(
        sample = rep(as.character(data[[input$sample]]), nrow(draws)),
        feature = input$features[g], .draw = rep(indices, each = nrow(data)),
        generated_value = as.vector(generated), location = as.vector(location),
        residual_variance = as.vector(variance)
      )
    }
  }
  result <- tibble::as_tibble(do.call(rbind, rows))
  result$inference_status <- parts$metadata$inference_status
  result
}

#' @export
sccomp_predict.sccomp_continuous_tbl <- function(
    fit, formula_composition = NULL, new_data = NULL, number_of_draws = 500,
    mcmc_seed = sample_seed(), summary_instead_of_draws = TRUE, robust = FALSE,
    ..., formula_location = NULL) {
  if (!is.null(formula_composition)) {
    stop("sccomp says: continuous fits use `formula_location` instead of `formula_composition`.",
         call. = FALSE)
  }
  if (!identical(robust, FALSE)) {
    stop("sccomp says: robust summaries are not available for continuous fits.",
         call. = FALSE)
  }
  if (length(list(...))) {
    stop("Unused continuous prediction arguments.", call. = FALSE)
  }
  if (!is.logical(summary_instead_of_draws) || length(summary_instead_of_draws) != 1L ||
      is.na(summary_instead_of_draws)) stop("summary_instead_of_draws must be logical.", call. = FALSE)
  continuous_predict_draws(fit, new_data, number_of_draws, mcmc_seed,
                            summary = summary_instead_of_draws,
                            formula_location = formula_location)
}

#' @export
sccomp_replicate.sccomp_continuous_tbl <- function(
    fit, formula_composition = NULL, formula_variability = NULL,
    number_of_draws = 1, mcmc_seed = sample_seed(),
    cache_stan_model = sccomp_stan_models_cache_dir, ...,
    formula_location = NULL, formula_variance = NULL, new_data = NULL) {
  if (!is.null(formula_composition) || !is.null(formula_variability)) {
    stop("sccomp says: continuous fits use `formula_location` and `formula_variance`, not `formula_composition` and `formula_variability`.",
         call. = FALSE)
  }
  if (length(list(...))) {
    stop("Unused continuous replication arguments.", call. = FALSE)
  }
  continuous_predict_draws(fit, new_data, number_of_draws, mcmc_seed,
                            formula_location = formula_location, formula_variance = formula_variance)
}

#' @export
sccomp_calculate_residuals.sccomp_continuous_tbl <- function(.data) {
  parts <- continuous_fit_parts(.data)
  prediction <- continuous_predict_draws(
    .data, NULL, prod(dim(parts$draws)[1:2]), parts$metadata$seed,
    summary = TRUE, generate = FALSE
  )
  input <- parts$input
  observed <- input$original_data[c(input$sample, input$feature, input$value)]
  names(observed) <- c("sample", "feature", "observed_value")
  observed$sample <- as.character(observed$sample)
  observed$feature <- as.character(observed$feature)
  dplyr::left_join(prediction, observed, by = c("sample", "feature")) |>
    dplyr::mutate(residual = .data$observed_value - .data$location,
                   standardized_residual = .data$residual / sqrt(.data$residual_variance))
}

#' @export
sccomp_proportional_fold_change.sccomp_continuous_tbl <- function(.data, ...) {
  stop("Proportional fold change is unavailable for continuous outcomes; use original-unit location contrasts.",
       call. = FALSE)
}

#' @export
print.sccomp_continuous_tbl <- function(x, ...) {
  meta <- attr(x, "continuous")
  cat("sccomp continuous:", meta$family, "|", meta$inference_status,
      "| HMC diagnostic gate:", isTRUE(meta$diagnostics$hmc_qualified), "\n")
  cat("Location:", deparse(meta$formula_location), "| log variance:",
      deparse(meta$formula_variance), "\n")
  print(tibble::as_tibble(x), ...)
  invisible(x)
}

#' @export
plot.sccomp_continuous_tbl <- function(x, ...) {
  if (length(list(...))) stop("Unused continuous plot arguments.", call. = FALSE)
  meta <- attr(x, "continuous")
  table <- x[x$component %in% c("location", "log_variance"), ]
  if (!nrow(table)) stop("No location/log-variance effects to plot.", call. = FALSE)
  ggplot2::ggplot(table, ggplot2::aes(x = .data$estimate, y = .data$feature)) +
    ggplot2::geom_segment(ggplot2::aes(x = .data$lower, xend = .data$upper,
                                       yend = .data$feature)) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(ggplot2::vars(.data$component, .data$parameter), scales = "free") +
    ggplot2::labs(x = "Location effect (input units) / log variance effect", y = "Feature",
                   subtitle = paste(meta$family, meta$inference_status,
                     if (isTRUE(meta$diagnostics$hmc_qualified)) "diagnostics passed" else "diagnostics not qualified")) +
    ggplot2::theme_bw()
}
