# ============================================================
# 07_theorem2_variance.R
# Theorem 2: Plug-in asymptotic variance
# Logistic equivalent-weight estimator
# ============================================================

expit <- function(x) {
  1 / (1 + exp(-x))
}

# ------------------------------------------------------------
# 1. Population-cell raw equivalent weights
# ------------------------------------------------------------

get_raw_weights <- function(
    beta,
    pi,
    X_pop,
    q_pop
) {

  p_pop <- as.numeric(
    expit(X_pop %*% beta)
  )

  v_pop <- p_pop * (1 - p_pop)

  Q <- crossprod(
    X_pop,
    X_pop * (pi * v_pop)
  )

  C <- crossprod(
    X_pop,
    q_pop * v_pop
  )

  as.numeric(
    X_pop %*%
      solve(Q, C)
  )
}

# ------------------------------------------------------------
# 2. Population functional evaluated on the sample
# ------------------------------------------------------------

T_function <- function(
    beta,
    pi,
    X_pop,
    sample_cell,
    y,
    q_pop
) {

  w_cell <- get_raw_weights(
    beta = beta,
    pi = pi,
    X_pop = X_pop,
    q_pop = q_pop
  )

  w_i <- w_cell[sample_cell]

  sum(w_i * y) / sum(w_i)
}

# ------------------------------------------------------------
# 3. Storage
# ------------------------------------------------------------

R <- length(updated_samp_list)

V_hat_logit <- numeric(R)
theta_hat_logit <- numeric(R)

# ------------------------------------------------------------
# 4. Main simulation loop
# ------------------------------------------------------------

for (k in seq_len(R)) {

  samp <- updated_samp_list[[k]]

  n <- nrow(samp)
  N <- sum(post_strat_df$n)

  X <- model.matrix(
    ~ mar + ethn + edu + age + city,
    data = samp
  )

  X_pop <- model.matrix(
    ~ mar + ethn + edu + age + city,
    data = post_strat_df
  )

  X_pop <- X_pop[, colnames(X), drop = FALSE]

  y <- samp$y_bin

  q_pop <- post_strat_df$n / N

  # ----------------------------------------------------------
  # Fit logistic model
  # ----------------------------------------------------------

  logit_fit <- glm(
    y_bin ~ mar + ethn + edu + age + city,
    data = samp,
    family = binomial()
  )

  beta_hat <- coef(logit_fit)
  beta_hat <- beta_hat[colnames(X)]

  p_hat <- as.numeric(
    fitted(logit_fit)
  )

  # ----------------------------------------------------------
  # Sample cell proportions
  # ----------------------------------------------------------

  sample_cell <- match(
    interaction(
      samp$mar,
      samp$ethn,
      samp$edu,
      samp$age,
      samp$city,
      drop = TRUE
    ),
    interaction(
      post_strat_df$mar,
      post_strat_df$ethn,
      post_strat_df$edu,
      post_strat_df$age,
      post_strat_df$city,
      drop = TRUE
    )
  )

  if (anyNA(sample_cell)) {
    stop(
      "Some sample observations could not be matched ",
      "to poststratification cells."
    )
  }

  cell_counts <- tabulate(
    sample_cell,
    nbins = nrow(post_strat_df)
  )

  pi_hat <- cell_counts / n

  # ----------------------------------------------------------
  # Raw equivalent weights
  # ----------------------------------------------------------

  w_cell <- get_raw_weights(
    beta = beta_hat,
    pi = pi_hat,
    X_pop = X_pop,
    q_pop = q_pop
  )

  w_raw <- w_cell[sample_cell]

  # Normalize raw weights for the ratio estimator.
  w_tilde <- w_raw / mean(w_raw)

  theta_hat <- mean(
    w_tilde * y
  )

  theta_hat_logit[k] <- theta_hat

  # ----------------------------------------------------------
  # Logistic information matrix
  # ----------------------------------------------------------

  W_hat <- p_hat * (1 - p_hat)

  I_hat <- crossprod(
    X,
    X * W_hat
  ) / n

  # ----------------------------------------------------------
  # B_beta
  #
  # Numerical derivative of the complete ratio functional.
  # ----------------------------------------------------------

  p_beta <- length(beta_hat)
  B_beta_hat <- numeric(p_beta)

  eps <- 1e-5

  for (j in seq_len(p_beta)) {

    beta_plus <- beta_hat
    beta_minus <- beta_hat

    beta_plus[j] <- beta_plus[j] + eps
    beta_minus[j] <- beta_minus[j] - eps

    T_plus <- T_function(
      beta = beta_plus,
      pi = pi_hat,
      X_pop = X_pop,
      sample_cell = sample_cell,
      y = y,
      q_pop = q_pop
    )

    T_minus <- T_function(
      beta = beta_minus,
      pi = pi_hat,
      X_pop = X_pop,
      sample_cell = sample_cell,
      y = y,
      q_pop = q_pop
    )

    B_beta_hat[j] <-
      (T_plus - T_minus) / (2 * eps)
  }

  # ----------------------------------------------------------
  # B_pi
  # ----------------------------------------------------------

  S <- nrow(post_strat_df)

  p_pop <- as.numeric(
    expit(X_pop %*% beta_hat)
  )

  v_pop <- p_pop * (1 - p_pop)

  Q_hat <- crossprod(
    X_pop,
    X_pop * (pi_hat * v_pop)
  )

  Q_inv <- solve(Q_hat)

  C_hat <- crossprod(
    X_pop,
    q_pop * v_pop
  )

  w_cell <- as.numeric(
    X_pop %*%
      Q_inv %*%
      C_hat
  )

  w_sample <- w_cell[sample_cell]

  A_hat <- mean(
    w_sample * y
  )

  D_hat <- mean(
    w_sample
  )

  # M[s, j] = x_s' Q^{-1} x_j.
  M <- X_pop %*%
    Q_inv %*%
    t(X_pop)

  # Column j:
  # -M[, j] * v_j * w_j
  dw_dpi <- -sweep(
    M,
    2,
    v_pop * w_cell,
    "*"
  )

  dw_sample_dpi <- dw_dpi[
    sample_cell,
    ,
    drop = FALSE
  ]

  dA_dpi <- colMeans(
    dw_sample_dpi * y
  )

  dD_dpi <- colMeans(
    dw_sample_dpi
  )

  B_pi_hat <- (
    dA_dpi * D_hat -
      A_hat * dD_dpi
  ) / D_hat^2

  # ----------------------------------------------------------
  # Three influence-function components
  # ----------------------------------------------------------

  # Component 1: outcome variation conditional on X.
  IF_1 <- w_tilde * (
    y - theta_hat
  )

  # Component 2: uncertainty from estimating beta.
  beta_score <- X * (
    y - p_hat
  )

  IF_2 <- as.numeric(
    beta_score %*%
      solve(I_hat, B_beta_hat)
  )

  # Component 3: uncertainty from estimating pi.
  e_minus_pi <- matrix(
    0,
    nrow = n,
    ncol = S
  )

  e_minus_pi[
    cbind(seq_len(n), sample_cell)
  ] <- 1

  e_minus_pi <- sweep(
    e_minus_pi,
    2,
    pi_hat,
    "-"
  )

  IF_3 <- as.numeric(
    e_minus_pi %*%
      B_pi_hat
  )

  # ----------------------------------------------------------
  # Total influence function
  # ----------------------------------------------------------

  psi_hat <-
    IF_1 +
    IF_2 +
    IF_3

  # ----------------------------------------------------------
  # Plug-in asymptotic variance
  # ----------------------------------------------------------

  V_hat_logit[k] <-
    mean(psi_hat^2)

  if (k == 1) {

    cat("\n")
    cat("============================================\n")
    cat("Theorem 2 diagnostics -- replication 1\n")
    cat("============================================\n")
    cat("theta_hat =", theta_hat, "\n")
    cat("V_hat =", V_hat_logit[k], "\n")
    cat("Mean IF =", mean(psi_hat), "\n")
    cat("Mean IF1 =", mean(IF_1), "\n")
    cat("Mean IF2 =", mean(IF_2), "\n")
    cat("Mean IF3 =", mean(IF_3), "\n")
    cat("Mean raw weight =", mean(w_raw), "\n")
    cat("Sum normalized weights =", sum(w_tilde), "\n")
    cat("Mean B_beta =", mean(B_beta_hat), "\n")
    cat("Mean B_pi =", mean(B_pi_hat), "\n")
  }

  cat("\rFinished Theorem 2 replication", k, "of", R)
}

cat("\n")

# ------------------------------------------------------------
# 5. Plug-in variance summary
# ------------------------------------------------------------

logit_variance_summary <- data.frame(
  Mean = mean(V_hat_logit),
  Q1 = quantile(V_hat_logit, 0.25),
  Median = quantile(V_hat_logit, 0.50),
  Q3 = quantile(V_hat_logit, 0.75),
  SD = sd(V_hat_logit)
)

print(logit_variance_summary)

# ------------------------------------------------------------
# 6. Empirical Monte Carlo variance
# ------------------------------------------------------------

theta_true <- mean(df$prob)

n <- nrow(updated_samp_list[[1]])

empirical_variance <- var(
  sqrt(n) *
    (theta_hat_logit - theta_true)
)

mean_plugin_variance <-
  mean(V_hat_logit)

variance_comparison <- data.frame(
  Empirical_Variance = empirical_variance,
  Mean_Plugin_Variance = mean_plugin_variance,
  Variance_Ratio =
    mean_plugin_variance /
    empirical_variance
)

print(variance_comparison)

# ------------------------------------------------------------
# 7. Standardized estimator
# ------------------------------------------------------------

Z_theorem2 <-
  sqrt(n) *
  (theta_hat_logit - theta_true) /
  sqrt(V_hat_logit)

# ------------------------------------------------------------
# 8. 95% Wald CI coverage
# ------------------------------------------------------------

z_975 <- qnorm(0.975)

SE_hat_logit <-
  sqrt(V_hat_logit / n)

CI_lower_logit <-
  theta_hat_logit -
  z_975 * SE_hat_logit

CI_upper_logit <-
  theta_hat_logit +
  z_975 * SE_hat_logit

coverage_logit <-
  CI_lower_logit <= theta_true &
  theta_true <= CI_upper_logit

logistic_CI_results <- data.frame(
  Method = "Logistic",
  Empirical_Variance = empirical_variance,
  Mean_Plugin_Variance = mean(V_hat_logit),
  Mean_SE = mean(SE_hat_logit),
  Coverage_95 = mean(coverage_logit),
  Mean_CI_Length =
    mean(CI_upper_logit - CI_lower_logit)
)

print(logistic_CI_results)
