continuous_input_fixture <- function() {
  set.seed(1001)
  input <- expand.grid(sample = sprintf("sample_%02d", seq_len(12)),
                       feature = sprintf("feature_%02d", seq_len(12)),
                       stringsAsFactors = FALSE)
  input$value <- stats::rnorm(nrow(input)) - 5
  input$condition <- factor(rep(c("control", "treated"), 72L))
  input$donor <- factor(rep(rep(seq_len(6), each = 2L), 12L))
  input$time <- rep(seq_len(12), 12L)
  input
}

make_continuous_input <- function(data = continuous_input_fixture(),
                                 formula_location = ~ condition,
                                 formula_variance = ~ 1,
                                 family = "student_t", scaling = "global",
                                 variance_trend = "none",
                                 student_t_df = "estimate") {
  continuous_model_input(
    data, formula_location, formula_variance, "sample", "feature", "value",
    family, variance_trend, scaling, student_t_df, sccomp_continuous_prior()
  )
}

test_that("continuous prior objects reject invalid widths and shapes", {
  expect_s3_class(sccomp_continuous_prior(), "sccomp_continuous_prior")
  expect_error(sccomp_continuous_prior(log_prior_df = c(0, 0)), "positive SD")
  expect_error(sccomp_continuous_prior(log_observation_df = c(NA, 1)), "finite")
  expect_error(sccomp_continuous_prior(log_variance_sd = Inf), "finite")
  expect_error(sccomp_continuous_prior(lkj = 0.5), "at least one")
  prior <- sccomp_continuous_prior()
  prior$random_log_sd_sd <- -1
  expect_error(validate_continuous_prior(prior), "positive SD")
})

test_that("one global transform preserves relative variances and original data", {
  data <- continuous_input_fixture()
  input <- make_continuous_input(data)
  original <- matrix(data$value, 12, 12)
  expected_scale <- stats::median(apply(original, 2, stats::mad))
  expect_equal(input$global_center, stats::median(original))
  expect_equal(input$global_scale, expected_scale)
  expect_equal(unname(input$stan_data$y * input$global_scale +
                       input$global_center), original)
  expect_equal(apply(input$stan_data$y, 2, stats::var) * input$global_scale^2,
               apply(original, 2, stats::var), ignore_attr = TRUE)
  expect_identical(input$original_data, data)
  expect_true(any(input$original_data$value < 0))
  expect_equal(make_continuous_input(data, scaling = "none")$global_scale, 1)
  constant <- data
  constant$value <- rep(seq_len(12), each = 12L)
  expect_error(make_continuous_input(constant), "positive within-feature MAD")
})

test_that("continuous data fail closed on incomplete or inconsistent input", {
  data <- continuous_input_fixture()
  expect_error(make_continuous_input(data[-1, ]), "every sample-feature pair")
  expect_error(make_continuous_input(rbind(data, data[1, ])), "Duplicate")
  data$value[1] <- Inf
  expect_error(make_continuous_input(data), "finite numeric")
  data <- continuous_input_fixture()
  data$condition[1] <- "treated"
  expect_error(make_continuous_input(data), "constant within sample")
  expect_error(make_continuous_input(formula_location = ~ missing),
               "Missing input columns")
  expect_error(make_continuous_input(formula_location = ~ 0), "one fixed coefficient")
  expect_error(sccomp_estimate_continuous(data, sample = sample, feature = feature,
                                        value = value, scaling = "none"), "explicitly supplied")
  expect_error(sccomp_estimate_continuous(data, sample = sample, feature = sample,
                                        value = value), "different columns")
  expect_error(sccomp_estimate_continuous(data, sample = sample, feature = feature,
                                        value = value, use_data = NA), "use_data must")
  expect_error(sccomp_estimate_continuous(data, sample = sample, feature = feature,
                                        value = value, mcmc_seed = Inf), "positive integer")
})

test_that("continuous designs have centered variance effects and rank checks", {
  input <- make_continuous_input(formula_variance = ~ condition + time)
  expect_equal(colMeans(input$stan_data$W), c(condition1 = 0, time = 0))
  expect_equal(input$stan_data$W[, 1], rep(c(1, -1), 6L), ignore_attr = TRUE)
  expect_equal(make_continuous_input()$stan_data$A, 0L)
  expect_equal(make_continuous_input(formula_location = ~ 1)$stan_data$C, 1L)
  no_intercept <- make_continuous_input(formula_variance = ~ 0 + condition)
  expect_equal(ncol(no_intercept$stan_data$W), 1L)
  expect_error(make_continuous_input(formula_location = ~ time + I(2 * time)),
               "full rank")
  expect_error(make_continuous_input(formula_variance = ~ condition + (1 | donor)),
               "unavailable in formula_variance")
  trend <- make_continuous_input(variance_trend = "linear")$stan_data$trend_covariate
  expect_equal(mean(trend), 0, tolerance = 1e-12)
  expect_equal(stats::sd(trend), 1)
})

test_that("continuous random designs preserve pairing and identifiable slopes", {
  input <- make_continuous_input(formula_location = ~ condition + (1 | donor))
  expect_equal(input$stan_data$K1, 1L)
  expect_equal(input$stan_data$J1, 6L)
  expect_equal(input$stan_data$group1, rep(seq_len(6), each = 2L))
  expect_error(make_continuous_input(formula_location = ~ condition +
                                     (1 + condition | donor)), "covariance design is not identifiable")
  repeated <- continuous_input_fixture()
  repeated$donor <- factor(rep(rep(seq_len(3), each = 4L), 12L))
  slope <- suppressWarnings(make_continuous_input(repeated, formula_location = ~ condition +
                                 (1 + condition | donor)))
  expect_equal(slope$stan_data$K1, 2L)
  aliased <- continuous_input_fixture()
  aliased$donor_copy <- aliased$donor
  expect_error(make_continuous_input(aliased, formula_location = ~ condition +
                                     (1 | donor) + (1 | donor_copy)), "covariance design is not identifiable")
  singleton <- continuous_input_fixture()
  singleton$donor <- singleton$sample
  expect_error(make_continuous_input(singleton, formula_location = ~ condition +
                                     (1 | donor)), "covariance design is not identifiable")
  between <- continuous_input_fixture()
  between$between <- rep(rep(seq_len(6), each = 2L), 12L)
  expect_error(make_continuous_input(between,
                                    formula_location = ~ condition +
                                      (1 + between | donor)),
               "no within-group variation")
  one <- continuous_input_fixture()
  one$donor <- 1L
  expect_error(make_continuous_input(one,
                                    formula_location = ~ condition + (1 | donor)),
               "two grouping levels")
})

test_that("Gaussian and fixed-df branches have no learned observation df", {
  expect_equal(make_continuous_input()$stan_data$estimate_df, 1L)
  gaussian <- make_continuous_input(family = "gaussian")$stan_data
  expect_equal(gaussian$estimate_df, 0L)
  expect_equal(gaussian$student_t, 0L)
  fixed <- make_continuous_input(student_t_df = 100)$stan_data
  expect_equal(fixed$estimate_df, 0L)
  expect_equal(fixed$fixed_df, 100)
  expect_error(make_continuous_input(student_t_df = 2), "greater than 2")
})
