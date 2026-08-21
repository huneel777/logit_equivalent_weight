# ============================================================
# 04_logistic_weights.R
# Logistic equivalent weights and raking
# ============================================================

library(dplyr)
library(survey)

expit <- function(x) {
  1 / (1 + exp(-x))
}

for (k in seq_along(updated_samp_list)) {

  samp <- updated_samp_list[[k]]

  n <- nrow(samp)
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

  # ----------------------------------------------------------
  # Fit logistic outcome model
  # ----------------------------------------------------------

  logit_fit <- glm(
    y_bin ~ mar + ethn + edu + age + city,
    data = samp,
    family = binomial()
  )

  beta_hat <- coef(logit_fit)
  beta_hat <- beta_hat[colnames(X)]

  p_hat <- as.numeric(fitted(logit_fit))

  # ----------------------------------------------------------
  # Logistic information matrix
  # ----------------------------------------------------------

  W_hat <- p_hat * (1 - p_hat)

  XtWX <- crossprod(
    X,
    X * W_hat
  )

  XtWX_inv <- solve(XtWX)

  # ----------------------------------------------------------
  # Population logistic variance
  # ----------------------------------------------------------

  p_pop <- as.numeric(
    expit(X_pop %*% beta_hat)
  )

  W_pop <- p_pop * (1 - p_pop)

  # ----------------------------------------------------------
  # Raw logistic equivalent weights
  #
  # w = (n/N) X (X'WX)^(-1)
  #     X_pop' (N_pop * W_pop)
  # ----------------------------------------------------------

  w_raw <- as.numeric(
    (n / N) *
      X %*%
      XtWX_inv %*%
      crossprod(
        X_pop,
        N_pop * W_pop
      )
  )

  samp$logit_raw_w <- w_raw

  # ----------------------------------------------------------
  # Practical normalized weights
  #
  # Match the original analysis: replace nonpositive
  # weights with the smallest positive weight, then
  # normalize to sum to n.
  # ----------------------------------------------------------

  w <- w_raw

  min_positive <- min(w[w > 0])

  w[w <= 0] <- min_positive

  samp$logit_w <-
    w * n / sum(w)

  # ----------------------------------------------------------
  # Rake logistic weights to population margins
  # ----------------------------------------------------------

  design <- svydesign(
    ids = ~1,
    data = samp,
    weights = ~logit_w
  )

  raked_design <- rake(
    design,
    sample.margins = list(
      ~mar, ~ethn, ~edu, ~age, ~city
    ),
    population.margins = pop_margins
  )

  samp$raked_logit_w <-
    weights(raked_design) * n / sum(weights(raked_design))

  updated_samp_list[[k]] <- samp

  cat("\rFinished logistic replication", k)
}

cat("\n")
