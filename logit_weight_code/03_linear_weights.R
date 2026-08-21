# ============================================================
# 03_linear_weights.R
# Linear equivalent weights and raking
# ============================================================

library(dplyr)
library(survey)

updated_samp_list <- vector(
  "list",
  length(updated_samp_list_0)
)

for (k in seq_along(updated_samp_list_0)) {

  samp <- updated_samp_list_0[[k]]

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
  # Raw linear equivalent weights
  # ----------------------------------------------------------

  unit_weight <-
    (n / N) *
    X %*%
    solve(crossprod(X)) %*%
    crossprod(X_pop, N_pop)

  unit_weight <- as.numeric(unit_weight)

  # Keep the unmodified equivalent weights for the
  # raw-weight analysis.
  samp$linear_raw_w <- unit_weight

  # Match the original analysis for the practical estimator:
  # replace nonpositive weights with the smallest
  # strictly positive weight.
  min_positive <- min(unit_weight[unit_weight > 0])

  unit_weight[unit_weight <= 0] <- min_positive

  samp$linear_w <-
    unit_weight * n / sum(unit_weight)

  # ----------------------------------------------------------
  # Rake linear weights to population margins
  # ----------------------------------------------------------

  design <- svydesign(
    ids = ~1,
    data = samp,
    weights = ~linear_w
  )

  raked_design <- rake(
    design,
    sample.margins = list(
      ~mar, ~ethn, ~edu, ~age, ~city
    ),
    population.margins = pop_margins
  )

  samp$raked_linear_w <-
    weights(raked_design) * n / N

  # Remove intermediate join variables.
  samp <- samp |>
    select(
      -any_of(c(
        "pop_count",
        "samp_count",
        "weight"
      ))
    )

  updated_samp_list[[k]] <- samp

  cat("\rFinished linear replication", k)
}

cat("\n")
