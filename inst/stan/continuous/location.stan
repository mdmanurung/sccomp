if (K1 > 0) {
  for (g in 1:G)
    mu[, g] += rows_dot_product(Z1, random_effect_1[g][, group1]');
}
if (K2 > 0) {
  for (g in 1:G)
    mu[, g] += rows_dot_product(Z2, random_effect_2[g][, group2]');
}
if (K3 > 0) {
  for (g in 1:G)
    mu[, g] += rows_dot_product(Z3, random_effect_3[g][, group3]');
}
if (K4 > 0) {
  for (g in 1:G)
    mu[, g] += rows_dot_product(Z4, random_effect_4[g][, group4]');
}
