# One row per independent dataset/component/contrast; bootstrap these rows later.
continuous_discovery_metrics <- function(table, truth, effect_threshold) {
  stopifnot(!anyDuplicated(table$feature), !anyDuplicated(names(truth)),
            all(table$feature %in% names(truth)),
            is.finite(effect_threshold), effect_threshold > 0)
  actual <- unname(truth[table$feature])
  stopifnot(all(is.finite(actual)), !anyNA(table$selected),
            all(table$direction %in% c("positive", "negative")))
  selected <- table$selected
  null <- abs(actual) <= effect_threshold
  wrong_sign <- actual * ifelse(table$direction == "positive", 1, -1) <= 0
  discoveries <- sum(selected)
  false <- sum(selected & null)
  sign_errors <- sum(selected & wrong_sign)
  directional_errors <- sum(selected & (null | wrong_sign))
  divisor <- max(1L, discoveries)
  data.frame(features = nrow(table), discoveries = discoveries,
    true_discoveries = sum(selected & !null & !wrong_sign),
    threshold_false_discoveries = false, false_signs = sign_errors,
    directional_false_discoveries = directional_errors,
    no_discoveries = discoveries == 0L, fdp = false / divisor,
    false_sign_proportion = sign_errors / divisor,
    directional_fdp = directional_errors / divisor,
    any_discovery = discoveries > 0L,
    bias = mean(table$estimate - actual),
    rmse = sqrt(mean((table$estimate - actual)^2)),
    coverage_95 = mean(actual >= table$lower & actual <= table$upper),
    mean_width_95 = mean(table$upper - table$lower))
}

if (identical(commandArgs(trailingOnly = TRUE), "--self-check")) {
  truth <- c(a = 0, b = .5, c = 2, d = -2)
  table <- data.frame(feature = names(truth), estimate = c(.1, .4, 1.9, 2),
    lower = c(-.1, .1, 1.8, 1.5), upper = c(.2, .8, 2.2, 2.5),
    selected = TRUE, direction = "positive")
  x <- continuous_discovery_metrics(table, truth, 1)
  stopifnot(x$discoveries == 4L, x$threshold_false_discoveries == 2L,
    x$false_signs == 2L, x$directional_false_discoveries == 3L,
    x$true_discoveries == 1L, x$fdp == .5, x$directional_fdp == .75,
    x$coverage_95 == .75)
  table$selected <- FALSE
  x <- continuous_discovery_metrics(table, truth, 1)
  stopifnot(x$no_discoveries, x$fdp == 0, x$false_sign_proportion == 0)
  reordered <- table[c(4, 2, 3, 1), ]
  stopifnot(identical(x, continuous_discovery_metrics(reordered, truth, 1)))
  cat("Composite null/sign errors, coverage, no-discovery and truth alignment checks passed.\n")
}
