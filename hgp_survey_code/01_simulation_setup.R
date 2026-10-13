# ============================================================
# 01_simulation_setup.R
# Companion code for the HGP survey-weighting paper
# ============================================================

# Required packages:
# install.packages(c("dplyr", "tidyr", "ggplot2", "MASS"))
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

inv_logit <- function(x) 1 / (1 + exp(-x))

rbf_kernel <- function(Fmat, alpha = 1.5, rho = 0.8, jitter = 1e-6) {
  D2 <- as.matrix(dist(Fmat))^2
  K  <- alpha^2 * exp(-D2 / (2 * rho^2))
  K + diag(jitter, nrow(Fmat))
}

rbf_kernel <- function(Fmat, alpha = 1.5, rho = 0.8, jitter = 1e-6) {
  D2 <- as.matrix(dist(Fmat))^2
  K  <- alpha^2 * exp(-D2 / (2 * rho^2))
  K + diag(jitter, nrow(Fmat))
}

make_terms <- function(pop_df) {
  mf <- model.frame(~ city + age + ethn + edu + mar,
                    data = pop_df, drop.unused.levels = FALSE)
  terms(~ city + age + ethn + edu + mar, data = mf)
}

mm_from_terms <- function(tt, df) {
  xlev <- .getXlevels(tt, df)
  mf <- model.frame(tt, data = df, xlev = xlev,
                    na.action = na.pass, drop.unused.levels = FALSE)
  model.matrix(tt, mf)
}

build_poststrat_table <- function(C_city, A_age, E_ethn, D_edu, M_mar, N_pop) {
  cell_df <- expand_grid(
    city = factor(1:C_city),
    age  = factor(1:A_age),
    ethn = factor(1:E_ethn),
    edu  = factor(1:D_edu),
    mar  = factor(1:M_mar)
  ) %>%
    mutate(group_K_full = row_number())

  w_raw  <- rgamma(nrow(cell_df), shape = 1.0, rate = 1.0)
  p_cell <- w_raw / sum(w_raw)

  cell_df$N_k <- as.vector(rmultinom(1, size = N_pop, prob = p_cell))
  cell_df
}

make_design_dense <- function(df) {
  X <- model.matrix(~ city + age + ethn + edu + mar, data = df)
  apply(X, 2, as.numeric)
}

generate_truth <- function(X_pop_clean, X_pop,
                           scenario = c("logit_parametric", "logit_gp"),
                           beta_mean = -0.3, beta_sd = 0.5,
                           gp_alpha = 1.5, gp_rho = 0.8) {
  scenario <- match.arg(scenario)
  K <- nrow(X_pop_clean)

  beta <- rnorm(ncol(X_pop), mean = beta_mean, sd = beta_sd)
  eta_fixed <- as.vector(X_pop %*% beta)

  gp_eff <- rep(0, K)
  if (scenario == "logit_gp") {
    Fmat <- scale(cbind(
      as.numeric(X_pop_clean$city),
      as.numeric(X_pop_clean$age),
      as.numeric(X_pop_clean$ethn),
      as.numeric(X_pop_clean$edu),
      as.numeric(X_pop_clean$mar)
    ))
    Sigma <- rbf_kernel(Fmat, alpha = gp_alpha, rho = gp_rho)
    gp_eff <- as.vector(MASS::mvrnorm(1, mu = rep(0, K), Sigma = Sigma))
  }

  eta_cell <- eta_fixed + gp_eff
  p_cell   <- inv_logit(eta_cell)
  list(beta = beta, gp_eff = gp_eff, p_cell = p_cell)
}

pop_mean <- function(pred, N_k) sum(pred * N_k) / sum(N_k)

compute_metrics_vec <- function(pred, truth, N_k) {
  mse_cell <- mean((pred - truth)^2)
  err_pop  <- pop_mean(pred, N_k) - pop_mean(truth, N_k)
  c(MSE_cell = mse_cell, PopMean_SqErr = err_pop^2)
}

compute_metrics_vec <- function(pred, truth, N_k) {
  pred  <- as.numeric(pred)
  truth <- as.numeric(truth)
  N_k   <- as.numeric(N_k)
  err <- pred - truth
  c(
    MSE_unweighted = mean(err^2),
    MSE_weighted   = sum(N_k * err^2) / sum(N_k)
  )
}

plot_cell_preds <- function(truth, pred, title, max_points = 2500) {
  plot_df <- data.frame(truth = truth, pred = pred)
  if (nrow(plot_df) > max_points) {
    set.seed(1)
    plot_df <- plot_df[sample.int(nrow(plot_df), max_points), ]
  }
  ggplot(plot_df, aes(x = truth, y = pred)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_abline(intercept = 0, slope = 1) +
    labs(title = title, x = "True cell probability", y = "Predicted cell probability") +
    theme_minimal()
}
