# ============================================================
# 05_simulation_results.R
# Main finite-sample simulation results
# ============================================================

# ------------------------------------------------------------
# 1. Helper
# ------------------------------------------------------------

compute_ess <- function(weights) {
  sum(weights)^2 / sum(weights^2)
}

R <- length(updated_samp_list)

# ------------------------------------------------------------
# 2. Population target
# ------------------------------------------------------------

theta0 <- mean(df$y_bin)

# ------------------------------------------------------------
# 3. Raw design-based estimator and ESS
# ------------------------------------------------------------

est_design_raw <- numeric(R)
ess_design_raw <- numeric(R)

for (k in seq_len(R)) {

  samp <- updated_samp_list[[k]]

  est_design_raw[k] <-
    sum(samp$adj_w * samp$y_bin) /
    sum(samp$adj_w)

  ess_design_raw[k] <-
    compute_ess(samp$adj_w)
}

# ------------------------------------------------------------
# 4. Raw linear equivalent-weight estimator and ESS
# ------------------------------------------------------------

est_linear_raw <- numeric(R)
ess_linear_raw <- numeric(R)

for (k in seq_len(R)) {

  samp <- updated_samp_list[[k]]

  w_linear_raw <- samp$linear_raw_w

  est_linear_raw[k] <-
    sum(w_linear_raw * samp$y_bin) /
    sum(w_linear_raw)

  ess_linear_raw[k] <-
    compute_ess(w_linear_raw)
}

# ------------------------------------------------------------
# 5. Raw logistic equivalent-weight estimator and ESS
# ------------------------------------------------------------

est_logit_raw <- numeric(R)
ess_logit_raw <- numeric(R)

for (k in seq_len(R)) {

  samp <- updated_samp_list[[k]]

  w_logit_raw <- samp$logit_raw_w

  est_logit_raw[k] <-
    sum(w_logit_raw * samp$y_bin) /
    sum(w_logit_raw)

  ess_logit_raw[k] <-
    compute_ess(w_logit_raw)
}

# ------------------------------------------------------------
# 6. Main raw-weight comparison
# ------------------------------------------------------------

results_raw <- data.frame(
  Method = c(
    "Design",
    "Linear",
    "Logistic"
  ),

  Population_Mean = rep(theta0, 3),

  Estimate = c(
    mean(est_design_raw),
    mean(est_linear_raw),
    mean(est_logit_raw)
  ),

  Bias = c(
    mean(est_design_raw) - theta0,
    mean(est_linear_raw) - theta0,
    mean(est_logit_raw) - theta0
  ),

  MSE = c(
    mean((est_design_raw - theta0)^2),
    mean((est_linear_raw - theta0)^2),
    mean((est_logit_raw - theta0)^2)
  ),

  Empirical_SD = c(
    sd(est_design_raw),
    sd(est_linear_raw),
    sd(est_logit_raw)
  ),

  Mean_ESS = c(
    mean(ess_design_raw),
    mean(ess_linear_raw),
    mean(ess_logit_raw)
  )
)

print(results_raw)

# ------------------------------------------------------------
# 7. Optional raked-weight comparison
# ------------------------------------------------------------

est_design_raked <- numeric(R)
est_linear_raked <- numeric(R)
est_logit_raked <- numeric(R)

ess_design_raked <- numeric(R)
ess_linear_raked <- numeric(R)
ess_logit_raked <- numeric(R)

for (k in seq_len(R)) {

  samp <- updated_samp_list[[k]]

  est_design_raked[k] <-
    sum(samp$raked_w * samp$y_bin) /
    sum(samp$raked_w)

  est_linear_raked[k] <-
    sum(samp$raked_linear_w * samp$y_bin) /
    sum(samp$raked_linear_w)

  est_logit_raked[k] <-
    sum(samp$raked_logit_w * samp$y_bin) /
    sum(samp$raked_logit_w)

  ess_design_raked[k] <-
    compute_ess(samp$raked_w)

  ess_linear_raked[k] <-
    compute_ess(samp$raked_linear_w)

  ess_logit_raked[k] <-
    compute_ess(samp$raked_logit_w)
}

results_raked <- data.frame(
  Method = c(
    "Design",
    "Linear",
    "Logistic"
  ),

  Population_Mean = rep(theta0, 3),

  Estimate = c(
    mean(est_design_raked),
    mean(est_linear_raked),
    mean(est_logit_raked)
  ),

  Bias = c(
    mean(est_design_raked) - theta0,
    mean(est_linear_raked) - theta0,
    mean(est_logit_raked) - theta0
  ),

  MSE = c(
    mean((est_design_raked - theta0)^2),
    mean((est_linear_raked - theta0)^2),
    mean((est_logit_raked - theta0)^2)
  ),

  Empirical_SD = c(
    sd(est_design_raked),
    sd(est_linear_raked),
    sd(est_logit_raked)
  ),

  Mean_ESS = c(
    mean(ess_design_raked),
    mean(ess_linear_raked),
    mean(ess_logit_raked)
  )
)

print(results_raked)

# ------------------------------------------------------------
# 8. Linear asymptotic variance
# ------------------------------------------------------------

V_hat_linear <- numeric(R)

for (k in seq_len(R)) {

  samp <- updated_samp_list[[k]]

  n_k <- nrow(samp)
  N <- sum(post_strat_df$n)
  N_pop <- post_strat_df$n

  X <- model.matrix(
    ~ mar + ethn + edu + age + city,
    data = samp
  )

  X_pop <- model.matrix(
    ~ mar + ethn + edu + age + city,
    data = post_strat_df
  )

  X_pop <- X_pop[, colnames(X), drop = FALSE]

  q_hat <- as.numeric(
    crossprod(X_pop, N_pop) / N
  )

  y <- samp$y_bin

  beta_hat <- solve(
    crossprod(X),
    crossprod(X, y)
  )

  residuals <- as.numeric(
    y - X %*% beta_hat
  )

  Q_hat <- crossprod(X) / n_k

  Q_inv_q <- solve(
    Q_hat,
    q_hat
  )

  a_hat <- as.numeric(
    X %*% Q_inv_q
  )

  psi_hat <- a_hat * residuals

  V_hat_linear[k] <-
    mean(psi_hat^2)
}

V_emp_linear <-
  nrow(updated_samp_list[[1]]) *
  var(est_linear_raw)

V_mean_linear <-
  mean(V_hat_linear)

linear_variance_comparison <- data.frame(
  Empirical_Asymptotic_Variance = V_emp_linear,
  Mean_Plugin_Variance = V_mean_linear,
  Ratio = V_mean_linear / V_emp_linear
)

print(linear_variance_comparison)

# ------------------------------------------------------------
# 9. 95% Wald CI coverage for the raw linear estimator
# ------------------------------------------------------------

n <- nrow(updated_samp_list[[1]])

linear_se_raw <-
  sqrt(V_hat_linear / n)

z_975 <- qnorm(0.975)

linear_lower_raw <-
  est_linear_raw - z_975 * linear_se_raw

linear_upper_raw <-
  est_linear_raw + z_975 * linear_se_raw

covered_linear <-
  linear_lower_raw <= theta0 &
  theta0 <= linear_upper_raw

linear_CI_results <- data.frame(
  Method = "Linear",
  Population_Mean = theta0,
  Mean_Estimate = mean(est_linear_raw),
  Mean_Variance = mean(V_hat_linear),
  Mean_SE = mean(linear_se_raw),
  Coverage_95 = mean(covered_linear),
  Mean_CI_Length =
    mean(linear_upper_raw - linear_lower_raw)
)

print(linear_CI_results)
