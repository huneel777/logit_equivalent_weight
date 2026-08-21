# ============================================================
# 01_simulation_setup.R
# Simulation population and repeated samples
# ============================================================

library(dplyr)

sigmoid <- function(x) {
  exp(x) / (1 + exp(x))
}

# ------------------------------------------------------------
# 1. Generate population
# ------------------------------------------------------------

set.seed(333)

N <- 1000000

x_mar <- rbinom(N, 1, 0.602) |> as.factor()

x_ethn <- sample(
  1:4,
  N,
  replace = TRUE,
  prob = c(0.382, 0.226, 0.313, 0.079)
) |> as.factor()

x_edu <- sample(
  1:5,
  N,
  replace = TRUE,
  prob = c(0.098, 0.188, 0.300, 0.190, 0.224)
) |> as.factor()

x_age <- sample(
  1:7,
  N,
  replace = TRUE,
  prob = c(0.045, 0.080, 0.252, 0.258, 0.237, 0.116, 0.012)
) |> as.factor()

x_city <- sample(
  1:16,
  N,
  replace = TRUE,
  prob = c(
    0.009, 0.007, 0.039, 0.005,
    0.061, 0.021, 0.167, 0.031,
    0.043, 0.067, 0.085, 0.057,
    0.013, 0.084, 0.125, 0.186
  )
) |> as.factor()

# ------------------------------------------------------------
# 2. Generate covariates and outcomes
# ------------------------------------------------------------

set.seed(77)

X1 <- rnorm(N, mean = 3, sd = 5)
X2 <- rbinom(N, size = 1, prob = 0.5)
X3 <- x_mar
X4 <- x_ethn
X5 <- x_edu
X6 <- x_age

beta_true <- rnorm(17, 0, 1)

df <- data.frame(
  mar = x_mar,
  ethn = x_ethn,
  edu = x_edu,
  age = x_age,
  city = x_city,
  X1, X2, X3, X4, X5, X6
)

X_new <- model.matrix(
  ~ X1 + X2 + X3 + X4 + X5 + X6,
  data = df
)

df$y_cont <- as.numeric(
  X_new %*% beta_true + rnorm(N, mean = 0, sd = 1)
)

prob <- abs(
  1 - 1 / (1 + exp(-X_new %*% beta_true))
)

df$prob <- as.numeric(prob)

df$y_bin <- rbinom(
  N,
  size = 1,
  prob = df$prob
)

# ------------------------------------------------------------
# 3. Sampling probabilities
# ------------------------------------------------------------

city_sample_probs <- c(
  0.095, 0.095, 0.098, 0.096,
  0.097, 0.095, 0.087, 0.094,
  0.028, 0.031, 0.038, 0.029,
  0.029, 0.029, 0.029, 0.030
)

mar_sample_probs <- c(0.73, 0.27)

ethn_sample_probs <- c(0.19, 0.43, 0.22, 0.16)

edu_sample_probs <- c(0.042, 0.27, 0.308, 0.25, 0.13)

age_sample_probs <- c(
  0.047, 0.14, 0.33, 0.20,
  0.14, 0.082, 0.061
)

# Population probabilities used to construct the base weight.
city_pop_probs <- c(
  0.009, 0.007, 0.039, 0.005,
  0.061, 0.021, 0.167, 0.031,
  0.043, 0.067, 0.085, 0.057,
  0.013, 0.084, 0.125, 0.186
)

mar_pop_probs <- c(0.398, 0.602)

ethn_pop_probs <- c(0.382, 0.226, 0.313, 0.079)

edu_pop_probs <- c(0.098, 0.188, 0.300, 0.190, 0.224)

age_pop_probs <- c(
  0.045, 0.080, 0.252, 0.258,
  0.237, 0.116, 0.012
)

# ------------------------------------------------------------
# 4. Generate repeated samples
# ------------------------------------------------------------

R <- 300
sample_size <- 3300

samp_list <- vector("list", R)

check_missing_levels <- function(sampled_df, full_df, variable) {
  setdiff(
    unique(full_df[[variable]]),
    unique(sampled_df[[variable]])
  )
}

for (k in seq_len(R)) {

  set.seed(k + 1)

  selection_prob <- with(
    df,
    city_sample_probs[city] *
      mar_sample_probs[mar] *
      ethn_sample_probs[ethn] *
      edu_sample_probs[edu] *
      age_sample_probs[age]
  )

  selection_prob <- selection_prob / sum(selection_prob)

  attempt <- 1

  repeat {

    df_sampled <- df[
      sample(
        seq_len(nrow(df)),
        sample_size,
        replace = FALSE,
        prob = selection_prob
      ),
    ]

    missing_levels <- lapply(
      c("city", "ethn", "mar", "edu", "age"),
      function(v) check_missing_levels(df_sampled, df, v)
    )

    names(missing_levels) <- c(
      "city", "ethn", "mar", "edu", "age"
    )

    if (all(lengths(missing_levels) == 0)) {
      break
    }

    for (variable in names(missing_levels)) {

      for (level in missing_levels[[variable]]) {

        extra_sample <- df[
          df[[variable]] == level,
        ]

        extra_sample <- extra_sample[
          sample(
            seq_len(nrow(extra_sample)),
            size = 1,
            replace = TRUE
          ),
        ]

        df_sampled <- bind_rows(
          df_sampled,
          extra_sample
        )
      }
    }

    attempt <- attempt + 1

    if (attempt > 3) {
      set.seed(k + 200 * k)
      attempt <- 1
    }
  }

  # ----------------------------------------------------------
  # Base weights
  # ----------------------------------------------------------

  df_sampled$city_prob <-
    city_pop_probs[as.numeric(df_sampled$city)]

  df_sampled$mar_prob <-
    mar_pop_probs[as.numeric(df_sampled$mar)]

  df_sampled$ethn_prob <-
    ethn_pop_probs[as.numeric(df_sampled$ethn)]

  df_sampled$edu_prob <-
    edu_pop_probs[as.numeric(df_sampled$edu)]

  df_sampled$age_prob <-
    age_pop_probs[as.numeric(df_sampled$age)]

  df_sampled$combined_prob <-
    df_sampled$city_prob *
    df_sampled$mar_prob *
    df_sampled$ethn_prob *
    df_sampled$edu_prob *
    df_sampled$age_prob

  df_sampled$w <-
    df_sampled$combined_prob /
    sum(df_sampled$combined_prob) * N

  df_sampled <- df_sampled |>
    select(-city_prob:-combined_prob)

  rownames(df_sampled) <- NULL

  # ----------------------------------------------------------
  # Response adjustment
  # ----------------------------------------------------------

  set.seed(777)

  resp <- rbinom(
    nrow(df_sampled),
    size = 1,
    prob = 0.1
  )

  df_sampled$id <- seq_len(nrow(df_sampled))
  df_sampled$resp <- resp

  df_sampled <- df_sampled |>
    group_by(city) |>
    mutate(
      total_response_weight =
        sum(w[resp == 0]),
      total_population_weight =
        sum(w),
      participation_rate =
        total_response_weight /
        total_population_weight,
      adj_w = if_else(
        resp == 0,
        w / participation_rate,
        w
      )
    ) |>
    ungroup() |>
    filter(resp == 0) |>
    select(-resp, -total_response_weight,
           -total_population_weight,
           -participation_rate) |>
    select(id, adj_w, everything())

  samp_list[[k]] <- df_sampled

  cat("\rFinished replication", k, "of", R)
}

cat("\n")
