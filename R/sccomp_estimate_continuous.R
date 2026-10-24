#' Fit continuous feature values with Bayesian variance moderation
#'
#' @description Fit transformed continuous assays using Gaussian or robust
#'   Student-t observations. Features share residual-variance hyperparameters;
#'   biological coefficients have fixed priors. This development API requires
#'   independent qualification for the intended assay and design.
#' @param .data Complete long data frame, one value per sample-feature pair.
#' @param formula_location,formula_variance One-sided formulas for location and
#'   log variance. Random intercepts/slopes are supported only for location.
#'   Random and residual covariance designs must be distinguishable. Correlated
#'   binary intercept/slopes with just two observations per group fail this
#'   check; repeat each condition within groups or simplify the random clause.
#' @param sample,feature,value Column names, as strings or bare names.
#' @param family Observation family; Student-t is the default.
#' @param variance_trend Optional linear observed-mean variance trend.
#' @param scaling Global robust affine scaling, or none with explicit priors.
#' @param student_t_df Shared learned df, or a fixed numeric df greater than two.
#' @param prior Priors from [sccomp_continuous_prior()].
#' @param inference_method HMC, or explicitly exploratory variational inference.
#' @param cores Total CPU threads available to sampling.
#' @param chains,iter_warmup,iter_sampling HMC chain count and iterations per chain.
#' @param mcmc_seed Reproducible integer seed.
#' @param output_directory Parent directory retaining a unique subdirectory of
#'   sampler CSVs and metadata for each fit, preventing result overwrites.
#' @param cache_stan_model Writable compiled-model cache.
#' @param use_data Whether to condition on observations; FALSE samples priors.
#' @param grainsize Feature chunk size for Stan reduce_sum.
#' @param ... Additional arguments to CmdStan sampling or variational inference.
#'   HMC defaults to `init = 0.5` on the unconstrained scale to avoid extreme
#'   initial log-normal random scales; supply `init` to override this choice.
#' @return A long component table inheriting from `sccomp_continuous_tbl`, with
#'   original-unit estimates, uncertainty, diagnostics, portable posterior draws,
#'   and complete model metadata. No publication-readiness claim is implied.
#' @export
sccomp_estimate_continuous <- function(
    .data, formula_location = ~ 1, formula_variance = ~ 1,
    sample, feature, value, family = c("student_t", "gaussian"),
    variance_trend = c("none", "linear"), scaling = c("global", "none"),
    student_t_df = "estimate", prior = sccomp_continuous_prior(),
    inference_method = c("hmc", "variational"), cores = 4L, chains = 4L,
    iter_warmup = 1000L, iter_sampling = 1000L, mcmc_seed = sample_seed(),
    output_directory = "sccomp_continuous_draws",
    cache_stan_model = sccomp_stan_models_cache_dir,
    use_data = TRUE, grainsize = 1L, ...) {
  family <- match.arg(family)
  variance_trend <- match.arg(variance_trend)
  scaling <- match.arg(scaling)
  inference_method <- match.arg(inference_method)
  if (scaling == "none" && missing(prior)) {
    stop("scaling = 'none' requires explicitly supplied assay-unit priors.",
         call. = FALSE)
  }
  sample <- rlang::as_name(rlang::ensym(sample))
  feature <- rlang::as_name(rlang::ensym(feature))
  value <- rlang::as_name(rlang::ensym(value))
  for (name in c("cores", "chains", "iter_warmup", "iter_sampling", "mcmc_seed")) {
    x <- get(name)
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x < 1 || x != floor(x)) {
      stop("`", name, "` must be a positive integer.", call. = FALSE)
    }
  }
  input <- continuous_model_input(
    .data, formula_location, formula_variance, sample, feature, value, family,
    variance_trend, scaling, student_t_df, prior, use_data, grainsize
  )
  offset <- base::qr.coef(qr(input$stan_data$X), rep(1, input$stan_data$N))
  constant_supported <- max(abs(input$stan_data$X %*% offset - 1)) < 1e-8
  if (scaling == "global" && !constant_supported) {
    stop("Global centering requires a location design spanning a constant; ",
         "use an intercept or scaling = 'none' with explicit priors.",
         call. = FALSE)
  }
  input$coefficient_offset <- if (constant_supported) offset else rep(0, length(offset))
  check_and_install_cmdstanr()
  compilation_time <- system.time(compiled <- continuous_load_model(
    "glm_multi_continuous", family, cache_stan_model))
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  fit_directory <- tempfile("fit-", tmpdir = normalizePath(output_directory))
  if (!dir.create(fit_directory)) stop("Cannot create a new fit output directory.", call. = FALSE)
  parallel_chains <- min(chains, cores)
  threads <- max(1L, floor(cores / parallel_chains))
  sampler_arguments <- list(...)
  if (inference_method == "hmc") {
    # Wide unconstrained starts compound through both log-normal scale levels.
    if (!"init" %in% names(sampler_arguments)) sampler_arguments$init <- 0.5
    fit <- do.call(compiled$model$sample, c(list(
      data = input$stan_data, chains = chains, parallel_chains = parallel_chains,
      threads_per_chain = threads, iter_warmup = iter_warmup,
      iter_sampling = iter_sampling, seed = mcmc_seed,
      output_dir = fit_directory), sampler_arguments))
  } else {
    warning("Variational inference is exploratory; selection requires opt-in.",
            call. = FALSE)
    fit <- compiled$model$variational(
      data = input$stan_data, seed = mcmc_seed,
      threads = cores,
      output_samples = chains * iter_sampling,
      output_dir = fit_directory, ...
    )
  }
  draws <- fit$draws(format = "draws_array")
  status <- if (inference_method == "variational") "exploratory" else "hmc"
  if (!use_data) status <- "prior_only"
  result <- continuous_estimate_table(draws, input)
  hyper_variables <- c("log_s0", "nu0", "trend_slope")
  for (b in seq_along(input$random)) {
    hyper_variables <- c(hyper_variables, paste0("mu_omega_", b),
                         paste0("tau_omega_", b))
  }
  hyper_variables <- posterior::variables(draws)[
    sub("\\[.*", "", posterior::variables(draws)) %in% hyper_variables
  ]
  hyper <- continuous_hyperparameter_table(draws, hyper_variables, input$global_scale)
  chain_diagnostics <- if (inference_method == "hmc") fit$diagnostic_summary() else NULL
  variables <- posterior::variables(draws)
  base_names <- sub("\\[.*", "", variables)
  free_names <- c("lp__", "beta", "gamma", "log_variance_raw", "log_s0",
                  "eta0", "eta_obs", "trend_slope",
                  unlist(lapply(1:4, function(b) paste0(
                    c("mu_omega_", "tau_omega_", "z_omega_", "random_raw_"), b))))
  free <- base_names %in% free_names
  for (j in which(grepl("^L_[1-4]\\[", variables))) {
    indices <- as.integer(strsplit(sub("\\]$", "", sub("^L_[1-4]\\[", "", variables[j])), ",")[[1]])
    free[j] <- indices[2L] > 1L && indices[3L] <= indices[2L]
  }
  internal <- posterior::summarise_draws(draws[, , variables[free], drop = FALSE])
  decision_diagnostics <- rbind(
    as.data.frame(result[!(result$component == "student_t_df" &
                            input$stan_data$estimate_df == 0L),
                         c("rhat", "ess_bulk", "ess_tail")]),
    as.data.frame(hyper[c("rhat", "ess_bulk", "ess_tail")]),
    as.data.frame(internal[c("rhat", "ess_bulk", "ess_tail")])
  )
  decision_diagnostics$variable <- c(
    with(result[!(result$component == "student_t_df" & input$stan_data$estimate_df == 0L), ],
         paste(feature, component, parameter, sep = ":")),
    hyper$variable, internal$variable
  )
  qualified <- inference_method == "hmc" && use_data && chains >= 4L &&
    dim(draws)[2L] == chains && length(chain_diagnostics$num_divergent) == chains &&
    all(is.finite(as.matrix(decision_diagnostics[c("rhat", "ess_bulk", "ess_tail")]))) &&
    all(decision_diagnostics$rhat < 1.01) &&
    all(decision_diagnostics$ess_bulk >= 400) &&
    all(decision_diagnostics$ess_tail >= 400) &&
    all(chain_diagnostics$num_divergent == 0) &&
    all(chain_diagnostics$num_max_treedepth == 0) &&
    all(chain_diagnostics$ebfmi > 0.3)
  metadata <- list(
    family = family, variance_trend = variance_trend, scaling = scaling,
    student_t_df = student_t_df, formula_location = formula_location,
    formula_variance = formula_variance, input = input,
    inference_method = inference_method, inference_status = status,
    seed = mcmc_seed, sampler_arguments = list(cores = cores, chains = chains,
      iter_warmup = iter_warmup, iter_sampling = iter_sampling, extra = sampler_arguments),
    source_identity = compiled$identity, runtime = utils::sessionInfo(),
    model_load_elapsed_seconds = unname(compilation_time["elapsed"]),
    sampler_files = fit$output_files(), hyperparameters = hyper,
    hyperparameter_units = "log_s0 is log prior SD scale and mu_omega is mean log random SD, in original input units; other hyperparameters are unit invariant",
    diagnostics = list(chains = chain_diagnostics,
                       parameters = decision_diagnostics, hmc_qualified = qualified)
  )
  result$inference_status <- status
  class(result) <- c("sccomp_continuous_tbl", "sccomp_tbl", class(result))
  attr(result, "continuous") <- metadata
  attr(result, "fit") <- draws
  attr(result, "inference_method") <- inference_method
  attr(result, "model_input") <- input$stan_data
  saveRDS(metadata, file.path(fit_directory,
                               paste0(basename(fit$output_files()[1L]),
                                      "-metadata.rds")))
  result
}

continuous_hyperparameter_table <- function(draws, variables, global_scale) {
  values <- unclass(draws[, , variables, drop = FALSE])
  log_sd <- sub("\\[.*", "", variables) %in% c("log_s0", paste0("mu_omega_", 1:4))
  values[, , log_sd] <- values[, , log_sd] + log(global_scale)
  posterior::summarise_draws(posterior::as_draws_array(values))
}

continuous_load_model <- function(name, family, cache_dir) {
  stan_root <- system.file("stan", package = "sccomp")
  source <- file.path(stan_root, paste0(name, ".stan"))
  includes <- list.files(file.path(stan_root, "continuous"), full.names = TRUE)
  if (!file.exists(source) || !length(includes)) {
    stop("Continuous Stan sources are unavailable in this package.", call. = FALSE)
  }
  identity <- list(package_version = as.character(utils::packageVersion("sccomp")),
                   source = basename(source),
                   source_md5 = tools::md5sum(c(source, includes)),
                   family = family, cmdstan = as.character(cmdstanr::cmdstan_version()),
                   cpp_options = list(stan_threads = TRUE),
                   makeflags = Sys.getenv("MAKEFLAGS"),
                   build_environment = Sys.getenv(c("CXX", "CXXFLAGS", "LDFLAGS", "STANCFLAGS")),
                   make_local_md5 = tools::md5sum(file.path(cmdstanr::cmdstan_path(), "make/local")))
  names(identity$source_md5) <- c(basename(source),
                                paste0("continuous/", basename(includes)))
  fingerprint <- tempfile("sccomp_source_identity_")
  on.exit(unlink(fingerprint), add = TRUE)
  saveRDS(identity, fingerprint, version = 2)
  identity$hash <- unname(tools::md5sum(fingerprint))
  directory <- file.path(cache_dir, identity$package_version, identity$hash)
  dir.create(file.path(directory, "continuous"), showWarnings = FALSE,
             recursive = TRUE)
  destination <- file.path(directory, basename(source))
  cached <- file.path(directory, names(identity$source_md5))
  valid <- all(file.exists(cached)) &&
    identical(unname(tools::md5sum(cached)), unname(identity$source_md5))
  if (!valid) {
    stopifnot(file.copy(source, destination, overwrite = TRUE),
              all(file.copy(includes, file.path(directory, "continuous"),
                            overwrite = TRUE)))
    stopifnot(identical(unname(tools::md5sum(cached)), unname(identity$source_md5)))
    saveRDS(identity, file.path(directory, "identity.rds"))
  }
  model <- cmdstanr::cmdstan_model(destination, include_paths = directory,
                                   cpp_options = identity$cpp_options)
  list(model = model, identity = identity)
}

continuous_summary <- function(draws, mapping) {
  summary <- posterior::summarise_draws(
    draws, estimate = mean,
    lower = function(x) unname(stats::quantile(x, 0.025)),
    upper = function(x) unname(stats::quantile(x, 0.975)),
    rhat = posterior::rhat, ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail
  )
  index <- match(mapping$variable, summary$variable)
  stopifnot(!anyNA(index))
  tibble::as_tibble(cbind(mapping[c("feature", "component", "parameter")],
                          as.data.frame(summary[index, -1L, drop = FALSE])))
}

continuous_estimate_table <- function(draws, input) {
  data <- input$stan_data
  values <- unclass(draws)
  mapping <- data.frame(feature = character(), component = character(),
                         parameter = character(), variable = character())
  transformed <- list()
  add <- function(variable, feature, component, parameter, scale = 1, offset = 0) {
    transformed[[variable]] <<- values[, , variable] * scale + offset
    mapping <<- rbind(mapping, data.frame(feature = feature, component = component,
                                           parameter = parameter, variable = variable))
  }
  for (g in seq_len(data$G)) {
    for (c in seq_len(data$C)) {
      add(sprintf("beta[%d,%d]", c, g), input$features[g], "location",
          colnames(data$X)[c], input$global_scale,
          input$global_center * input$coefficient_offset[c])
    }
    for (a in seq_len(data$A)) {
      add(sprintf("gamma[%d,%d]", a, g), input$features[g], "log_variance",
          colnames(data$W)[a])
    }
    add(sprintf("sigma_sq[%d]", g), input$features[g], "residual_variance",
        "baseline", input$global_scale^2)
    for (b in seq_along(input$random)) {
      block <- input$random[[b]]
      labels <- colnames(block$design$matrix)
      for (q in seq_along(labels)) {
        add(sprintf("omega_%d[%d,%d]", b, g, q), input$features[g], "random_sd",
            paste(b, block$grouping, labels[q], sep = ":"), input$global_scale)
        if (q > 1L) {
          for (r in seq_len(q - 1L)) {
            variable <- sprintf("correlation_%d[%d,%d,%d]", b, g, q, r)
            correlation <- values[, , sprintf("L_%d[%d,%d,1]", b, g, q)] *
              values[, , sprintf("L_%d[%d,%d,1]", b, g, r)]
            if (r > 1L) {
              for (k in 2:r) {
                correlation <- correlation +
                  values[, , sprintf("L_%d[%d,%d,%d]", b, g, q, k)] *
                  values[, , sprintf("L_%d[%d,%d,%d]", b, g, r, k)]
              }
            }
            transformed[[variable]] <- correlation
            mapping <- rbind(mapping, data.frame(
              feature = input$features[g], component = "random_correlation",
              parameter = paste(b, block$grouping, labels[q], labels[r], sep = ":"),
              variable = variable
            ))
          }
        }
      }
    }
  }
  if (data$student_t == 1L) add("nu_obs", NA_character_, "student_t_df", "nu_obs")
  array <- array(unlist(transformed),
                  dim = c(dim(values)[1:2], length(transformed)),
                  dimnames = list(dimnames(values)[[1]], dimnames(values)[[2]],
                                  names(transformed)))
  continuous_summary(posterior::as_draws_array(array), mapping)
}
