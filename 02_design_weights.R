# ============================================================
# 02_design_weights.R
# Design-based weights and raking
# ============================================================

library(dplyr)
library(survey)

# ------------------------------------------------------------
# 1. Population poststratification cells
# ------------------------------------------------------------

post_strat_df <- df |>
  group_by(mar, ethn, edu, age, city) |>
  summarise(
    n = n(),
    .groups = "drop"
  )

N <- sum(post_strat_df$n)

# ------------------------------------------------------------
# 2. Population margins for raking
# ------------------------------------------------------------

pop_margins <- list(
  data.frame(
    mar = levels(df$mar),
    Freq = as.numeric(table(df$mar))
  ),
  data.frame(
    ethn = levels(df$ethn),
    Freq = as.numeric(table(df$ethn))
  ),
  data.frame(
    edu = levels(df$edu),
    Freq = as.numeric(table(df$edu))
  ),
  data.frame(
    age = levels(df$age),
    Freq = as.numeric(table(df$age))
  ),
  data.frame(
    city = levels(df$city),
    Freq = as.numeric(table(df$city))
  )
)

# ------------------------------------------------------------
# 3. Design-based poststratification and raking
# ------------------------------------------------------------

updated_samp_list_0 <- vector(
  "list",
  length(samp_list)
)

for (k in seq_along(samp_list)) {

  samp <- samp_list[[k]]

  pop_counts <- post_strat_df

  samp_counts <- samp |>
    count(
      mar, ethn, edu, age, city,
      name = "samp_count"
    )

  weights <- pop_counts |>
    rename(pop_count = n) |>
    left_join(
      samp_counts,
      by = c("mar", "ethn", "edu", "age", "city")
    ) |>
    mutate(
      weight = pop_count / samp_count
    )

  samp <- samp |>
    left_join(
      weights,
      by = c("mar", "ethn", "edu", "age", "city")
    )

  design <- svydesign(
    ids = ~1,
    data = samp,
    weights = ~adj_w
  )

  raked_design <- rake(
    design,
    sample.margins = list(
      ~mar, ~ethn, ~edu, ~age, ~city
    ),
    population.margins = pop_margins
  )

  samp$raked_w <- weights(raked_design)

  updated_samp_list_0[[k]] <- samp

  cat("\rFinished design-based replication", k)
}

cat("\n")
