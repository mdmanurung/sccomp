library(dplyr)
library(tidyr)

student_t_fixture <- function() {
  tidyr::expand_grid(
    sample = paste0("sample_", seq_len(8)),
    feature = paste0("feature_", seq_len(3))
  ) |>
    dplyr::mutate(
      condition = factor(
        rep(rep(c("control", "treated"), each = 4), each = 3),
        levels = c("control", "treated")
      ),
      donor = factor(rep(rep(paste0("donor_", seq_len(4)), 2), each = 3)),
      value = c(
        1.0, 2.0, 3.0, 1.1, 2.1, 3.2, 0.9, 1.9, 2.8, 1.2, 2.2, 3.1,
        2.0, 2.4, 3.0, 2.1, 2.5, 3.1, 1.9, 2.6, 2.9, 2.2, 2.5, 3.2
      )
    )
}

test_that("Student-t data preparation preserves continuous values", {
  input <- student_t_fixture() |>
    dplyr::mutate(random_effect = "1")

  spread <- sccomp:::data_to_spread(
    input,
    ~ condition,
    sample,
    feature,
    value,
    "random_effect",
    noise_model = "student_t"
  )

  expect_true(all(spread$exposure == 1L))
  expect_equal(spread$feature_1, input$value[input$feature == "feature_1"])

  model_input <- sccomp:::data_spread_to_model_input(
    spread,
    ~ condition,
    sample,
    feature,
    value,
    approximate_posterior_inference = TRUE,
    formula_variability = ~ 1,
    random_effect_elements = sccomp:::parse_formula_random_effect(~ condition),
    noise_model = "student_t"
  )

  expect_false(model_input$is_proportion)
  expect_equal(dim(model_input$y), c(8L, 3L))
  expect_identical(rownames(model_input$y), spread$sample)
  expect_equal(unname(model_input$y[, "feature_1"]), spread$feature_1)
})

test_that("Student-t source uses an unconstrained robust likelihood", {
  stan_file <- system.file("stan", "glm_multi_student_t.stan", package = "sccomp")
  source <- paste(readLines(stan_file, warn = FALSE), collapse = "\n")

  expect_match(source, "matrix\\[C, M\\] beta", fixed = FALSE)
  expect_match(source, "student_t_lpdf", fixed = TRUE)
  expect_match(source, "sqrt((student_t_df - 2) / student_t_df)", fixed = TRUE)
  expect_false(grepl("softmax", source, fixed = TRUE))
  expect_false(grepl("sum_to_zero_vector", source, fixed = TRUE))
})

test_that("Student-t rejects missing sample-feature pairs", {
  input <- student_t_fixture()[-1, ]

  expect_error(
    sccomp:::sccomp_glm_data_frame_counts(
      input,
      formula_composition = ~ condition,
      formula_variability = ~ 1,
      sample = "sample",
      cell_group = "feature",
      abundance = "value",
      noise_model = "student_t",
      student_t_df = 4
    ),
    "requires one observed value for every sample/feature pair",
    fixed = TRUE
  )
})

test_that("Student-t fits, tests, predicts, and calculates residuals", {
  skip_cmdstan()
  input <- student_t_fixture()
  output_dir <- tempfile("sccomp_student_t_draws_")
  cache_dir <- tempfile("sccomp_student_t_cache_")
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  fit <- sccomp_estimate(
    input,
    formula_composition = ~ condition + (1 | donor),
    formula_variability = ~ 1,
    sample = "sample",
    cell_group = "feature",
    abundance = "value",
    noise_model = "student_t",
    student_t_df = 4,
    inference_method = "hmc",
    cores = 1,
    max_sampling_iterations = 50,
    mcmc_seed = 42,
    verbose = FALSE,
    output_directory = output_dir,
    cache_stan_model = cache_dir
  )

  expect_identical(attr(fit, "noise_model"), "student_t")
  expect_identical(attr(fit, "student_t_df"), 4)
  expect_true(all(is.finite(fit$c_effect)))

  tested <- sccomp_test(fit, test_composition_above_logit_fold_change = 0.1)
  expect_true(all(c("c_pH0", "c_FDR") %in% colnames(tested)))

  predicted <- sccomp_predict(fit, number_of_draws = 20, mcmc_seed = 42)
  expect_true(all(c(
    "location_mean", "location_lower", "location_upper",
    "value_mean", "value_lower", "value_upper"
  ) %in% colnames(predicted)))

  residuals <- sccomp_calculate_residuals(fit)
  expect_equal(nrow(residuals), nrow(input))
  expect_true(all(is.finite(residuals$residuals)))
})
