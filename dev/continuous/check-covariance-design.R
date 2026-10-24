# Exact binary-pair covariance alias, independent of the R input validator.
group <- rep(1:20, each = 2L)
condition <- rep(0:1, 20L)
n <- length(condition)
same_group <- outer(group, group, "==")
intercept <- same_group * 1
cross <- (outer(condition, rep(1, n)) + outer(rep(1, n), condition)) * same_group
slope <- outer(condition, condition) * same_group
stopifnot(identical(unname(intercept - cross + 2 * slope), unname(diag(n))))
cat("Two-observation binary pairs: I = intercept covariance - cross covariance + 2 * slope covariance.\n")
