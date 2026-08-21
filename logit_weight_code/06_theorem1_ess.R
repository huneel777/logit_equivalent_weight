# ============================================================
# 06_theorem1_ess.R
# Theorem 1: Effective sample size
# ============================================================

# ------------------------------------------------------------
# 1. Effective sample size
# ------------------------------------------------------------

compute_ess <- function(weights) {
  sum(weights)^2 / sum(weights^2)
}


# ------------------------------------------------------------
# 2. Population-cell probabilities
# ------------------------------------------------------------

df_cells <- df |>
  dplyr::group_by(
    mar, ethn, edu, age, city
  ) |>
  dplyr::summarise(
    prob = mean(prob),
    N_pop = n(),
    .groups = "drop"
  )


# ------------------------------------------------------------
# 3. Logistic equivalent weights
# ------------------------------------------------------------

logit_eqw_from_probs <- function(
    prob,
    prob_pop,
    X,
    X_pop,
    N_pop,
    n,
    N
) {
  
  if (length(prob) != nrow(X)) {
    stop("length(prob) must equal nrow(X).")
  }
  
  if (length(prob_pop) != nrow(X_pop)) {
    stop("length(prob_pop) must equal nrow(X_pop).")
  }
  
  if (length(N_pop) != nrow(X_pop)) {
    stop("length(N_pop) must equal nrow(X_pop).")
  }
  
  if (n != nrow(X)) {
    stop("n must equal nrow(X).")
  }
  
  if (N != sum(N_pop)) {
    stop("N must equal sum(N_pop).")
  }
  
  # Logistic variance terms
  s <- as.numeric(
    prob * (1 - prob)
  )
  
  s_pop <- as.numeric(
    prob_pop * (1 - prob_pop)
  )
  
  # Sample information matrix
  XtWX <- crossprod(
    X,
    X * s
  )
  
  XtWXi <- solve(XtWX)
  
  # Population component
  r_vec <-
    (n / N) *
    crossprod(
      X_pop,
      N_pop * s_pop
    )
  
  # Equivalent weights for sample observations
  w_sample <- as.numeric(
    t(r_vec) %*%
      XtWXi %*%
      t(X)
  )
  
  # Equivalent weights evaluated at population cells
  w_pop <- as.numeric(
    t(r_vec) %*%
      XtWXi %*%
      t(X_pop)
  )
  
  list(
    w_sample = w_sample,
    w_pop = w_pop
  )
}


# ------------------------------------------------------------
# 4. Calculate Theorem 1 quantities
# ------------------------------------------------------------

R <- length(updated_samp_list)

ess_logit_reconstructed <- numeric(R)
ess_logit_cell_samp <- numeric(R)


for (i in seq_len(R)) {
  
  samp <- updated_samp_list[[i]]
  
  n <- nrow(samp)
  N <- sum(post_strat_df$n)
  N_pop <- post_strat_df$n
  
  # ----------------------------------------------------------
  # Design matrices
  # ----------------------------------------------------------
  
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
  # Population probabilities
  # ----------------------------------------------------------
  
  prob_pop <- post_strat_df |>
    dplyr::select(
      mar, ethn, edu, age, city
    ) |>
    dplyr::left_join(
      df_cells |>
        dplyr::select(
          mar, ethn, edu, age, city, prob
        ),
      by = c(
        "mar",
        "ethn",
        "edu",
        "age",
        "city"
      )
    ) |>
    dplyr::pull(prob)
  
  if (anyNA(prob_pop)) {
    stop("Some population cells have no matching probability.")
  }
  
  
  # ----------------------------------------------------------
  # Logistic equivalent weights
  # ----------------------------------------------------------
  
  eqw <- logit_eqw_from_probs(
    prob = samp$prob,
    prob_pop = prob_pop,
    X = X,
    X_pop = X_pop,
    N_pop = N_pop,
    n = n,
    N = N
  )
  
  w_pop <- eqw$w_pop
  
  
  # ----------------------------------------------------------
  # Reconstructed individual ESS / n
  # ----------------------------------------------------------
  
  ess_logit_reconstructed[i] <-
    compute_ess(eqw$w_sample) / n
  
  
  # ----------------------------------------------------------
  # Sample-cell proportions
  # ----------------------------------------------------------
  
  n_s <- post_strat_df |>
    dplyr::select(
      mar, ethn, edu, age, city
    ) |>
    dplyr::left_join(
      samp |>
        dplyr::count(
          mar, ethn, edu, age, city,
          name = "n_s"
        ),
      by = c(
        "mar",
        "ethn",
        "edu",
        "age",
        "city"
      )
    ) |>
    dplyr::mutate(
      n_s = dplyr::coalesce(n_s, 0L)
    ) |>
    dplyr::pull(n_s)
  
  pi_hat <- n_s / n
  
  
  # ----------------------------------------------------------
  # Sample-cell ESS / n
  # ----------------------------------------------------------
  
  ess_logit_cell_samp[i] <-
    (sum(pi_hat * w_pop)^2) /
    sum(pi_hat * w_pop^2)
}


# ------------------------------------------------------------
# 5. Theorem 1 summary
# ------------------------------------------------------------

theorem1_summary <- data.frame(
  
  Quantity = c(
    "Reconstructed individual ESS / n",
    "Sample-cell ESS / n"
  ),
  
  Mean = c(
    mean(ess_logit_reconstructed),
    mean(ess_logit_cell_samp)
  ),
  
  Q1 = c(
    quantile(
      ess_logit_reconstructed,
      0.25
    ),
    quantile(
      ess_logit_cell_samp,
      0.25
    )
  ),
  
  Median = c(
    median(ess_logit_reconstructed),
    median(ess_logit_cell_samp)
  ),
  
  Q3 = c(
    quantile(
      ess_logit_reconstructed,
      0.75
    ),
    quantile(
      ess_logit_cell_samp,
      0.75
    )
  )
)

print(theorem1_summary)


# ------------------------------------------------------------
# 6. Compare with stored raw logistic weights
# ------------------------------------------------------------

ess_logit_raw <- numeric(R)

for (i in seq_len(R)) {
  
  ess_logit_raw[i] <-
    compute_ess(
      updated_samp_list[[i]]$logit_raw_w
    )
}


theorem1_stored_weight_summary <- data.frame(
  
  Mean_ESS_over_n = mean(
    ess_logit_raw /
      sapply(
        updated_samp_list,
        nrow
      )
  ),
  
  Q1 = quantile(
    ess_logit_raw /
      sapply(
        updated_samp_list,
        nrow
      ),
    0.25
  ),
  
  Median = median(
    ess_logit_raw /
      sapply(
        updated_samp_list,
        nrow
      )
  ),
  
  Q3 = quantile(
    ess_logit_raw /
      sapply(
        updated_samp_list,
        nrow
      ),
    0.75
  )
)

print(theorem1_stored_weight_summary)