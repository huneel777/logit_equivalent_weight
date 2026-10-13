# ============================================================
# 02_fit_models.R
# Companion code for the HGP survey-weighting paper
# ============================================================

# This script assumes 01_simulation_setup.R has been sourced.
# Requires CmdStanR and a compatible CmdStan installation.
# The Stan model file used in the paper must be available at stan_file.

# install.packages("cmdstanr")  # installation instructions: https://mc-stan.org/cmdstanr/

fit_logit_map_with_priors <- function(data, prior_mean, sigma_u = 1.0, s_beta = 2.5, maxit = 300) {
  X <- model.matrix(~ city + age + ethn + edu + mar, data = data)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  y <- as.integer(data$y)
  h <- as.integer(data$group_H)
  H <- length(prior_mean)
  p <- ncol(X)

  nll <- function(par) {
    beta <- par[1:p]
    u    <- par[(p + 1):(p + H)]
    eta <- as.vector(X %*% beta + u[h])
    p_i <- inv_logit(eta)
    ll <- sum(y * log(p_i + 1e-12) + (1 - y) * log(1 - p_i + 1e-12))

    lp_beta <- sum(dnorm(beta, 0, s_beta, log = TRUE))
    lp_u    <- sum(dnorm(u, prior_mean, sigma_u, log = TRUE))
    -(ll + lp_beta + lp_u)
  }

  par0 <- c(rep(0, ncol(X)), prior_mean)
  fit <- optim(par0, nll, method = "BFGS", control = list(maxit = maxit))

  list(
    beta = fit$par[1:ncol(X)],
    u = fit$par[(ncol(X) + 1):(ncol(X) + H)],
    converged = (fit$convergence == 0),
    X_cols = colnames(X),
    sigma_u = sigma_u
  )
}

predict_logit_map <- function(map_fit, newdata) {
  Xn <- model.matrix(~ city + age + ethn + edu + mar, data = newdata)
  Xn <- Xn[, colnames(Xn) != "(Intercept)", drop = FALSE]
  Xn <- Xn[, map_fit$X_cols, drop = FALSE]
  h <- as.integer(newdata$group_H)
  eta <- as.vector(Xn %*% map_fit$beta + map_fit$u[h])
  as.numeric(inv_logit(eta))
}

run_scenario <- function(
  scenario = c("logit_parametric","logit_gp"),
  stan_file = "Weight_GP_Binary_Big.stan",
  C_city = 10, A_age = 6, E_ethn = 4, D_edu = 5, M_mar = 2,
  N_pop = 100000, n_samp = 2000, K_obs_target = 1000,
  seed_pop = 777, seed_truth = 100, seed_sample = 2025,
  chains = 2, iter_warmup = 500, iter_sampling = 500,
  adapt_delta = 0.95, max_treedepth = 12
) {
  scenario <- match.arg(scenario)

  ## A) population cells
  set.seed(seed_pop)
  K_cells  <- C_city * A_age * E_ethn * D_edu * M_mar
  H_groups <- C_city * A_age * E_ethn

  X_pop_clean <- build_poststrat_table(C_city, A_age, E_ethn, D_edu, M_mar, N_pop) %>%
    mutate(group_H = as.integer(factor(paste(city, age, ethn, sep = "_")))) %>%
    arrange(group_K_full)

  stopifnot(nrow(X_pop_clean) == K_cells)
  stopifnot(sum(X_pop_clean$N_k) == N_pop)
  stopifnot(n_distinct(X_pop_clean$group_H) == H_groups)

  ## expand to unit-level population
  df_pop <- X_pop_clean %>%
    dplyr::filter(N_k > 0) %>%
    tidyr::uncount(weights = N_k, .id = "unit_in_cell") %>%
    dplyr::mutate(id = dplyr::row_number()) %>%
    dplyr::select(id, city, age, ethn, edu, mar, group_K_full, group_H)


  ## truth
  X_pop_main <- sparse.model.matrix(~ city + age + ethn + edu + mar, data = X_pop_clean)
  set.seed(seed_truth)
  truth <- generate_truth(X_pop_clean, X_pop_main, scenario = scenario)
  X_pop_clean$p_true <- truth$p_cell

  df_pop <- df_pop %>%
    left_join(X_pop_clean %>% select(group_K_full, p_true), by = "group_K_full") %>%
    mutate(y = rbinom(n(), 1, p_true))

  ## sample restricted to K_obs_target frame cells
  set.seed(seed_sample)
  cells_in_frame <- sample(X_pop_clean$group_K_full, size = K_obs_target, replace = FALSE)
  df_scope <- df_pop %>% filter(group_K_full %in% cells_in_frame)
  df_samp  <- df_scope %>% slice_sample(n = n_samp, replace = FALSE)

  cat("\nSCENARIO:", scenario, "\n")
  cat("Sample size:", nrow(df_samp), "\n")
  cat("Observed cells in sample:", n_distinct(df_samp$group_K_full), "out of", K_cells, "\n")

  ## re-index K to observed cells for Stan
  obs_cells  <- sort(unique(df_samp$group_K_full))
  K_obs_real <- length(obs_cells)
  cat("K_obs_real (used in Stan):", K_obs_real, "\n")

  K_map <- setNames(seq_len(K_obs_real), obs_cells)
  df_samp <- df_samp %>% mutate(group_K = as.integer(K_map[as.character(group_K_full)]))

  X_unique_df <- X_pop_clean %>%
    dplyr::filter(group_K_full %in% obs_cells) %>%
    arrange(match(group_K_full, obs_cells))

  X_pop_df <- X_pop_clean %>% arrange(group_K_full)

  group_H_sample <- as.integer(df_samp$group_H)
  group_H_pop    <- as.integer(X_pop_df$group_H)

  ## priors by H from UNIT pop (df_pop)
  prior_by_H <- df_pop %>%
    group_by(group_H) %>%
    summarize(rate = mean(y), .groups = "drop") %>%
    arrange(group_H)
  prior_mean <- qlogis(pmin(pmax(prior_by_H$rate, 1e-6), 1 - 1e-6))

  ## Stan data
  X        <- make_design_dense(df_samp)
  X_unique <- make_design_dense(X_unique_df)
  X_pop    <- make_design_dense(X_pop_df)

  stan_data <- list(
    n = nrow(X),
    P = ncol(X),
    group_H = group_H_sample,
    group_K = as.integer(df_samp$group_K),
    H = H_groups,
    K = K_obs_real,
    y = as.integer(df_samp$y),
    X = X,
    X_unique = X_unique,
    prior_mean = prior_mean,
    n_pop = nrow(X_pop),
    X_pop = X_pop,
    group_H_pop = group_H_pop
  )

  init_fun <- function() list(
    mu_raw = rep(0, H_groups),
    tau_mu = 1.0,
    alpha = 1.0,
    rho = 1.0,
    eta = rep(0, K_obs_real)
  )

  GP_model <- cmdstan_model(stan_file)

  stan_fit <- GP_model$sample(
    data = stan_data,
    chains = chains,
    parallel_chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = 2025,
    init = init_fun,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth
  )

  ## Return “raw” results only; everything else computed in postprocess()
  list(
    scenario = scenario,
    df_samp = df_samp,
    X_pop_clean = X_pop_clean,
    X_pop_df = X_pop_df,
    df_pop = df_pop,
    prior_mean = prior_mean,
    stan_fit = stan_fit
  )
}

postprocess_scenario_noW <- function(res, max_points = 2500) {
  stopifnot(!is.null(res$stan_fit))

  df_samp     <- res$df_samp
  X_pop_clean <- res$X_pop_clean
  X_pop_df    <- res$X_pop_df
  prior_mean  <- res$prior_mean
  stan_fit    <- res$stan_fit
  scenario    <- res$scenario

  ## ---- HGP predictions (from generated quantities E_y_star) ----
  y_star_draws <- stan_fit$draws("E_y_star", format = "df")
  y_star_cols  <- grep("^E_y_star\\[", names(y_star_draws), value = TRUE)
  if (length(y_star_cols) == 0) {
    stop("No columns found for E_y_star[...]. Check Stan generated quantities name.")
  }
  y_star_only <- y_star_draws[, y_star_cols, drop = FALSE]
  X_pop_df$pred_hgp <- as.numeric(colMeans(y_star_only))

  ## ---- GLM (FIXED: master terms from FULL population) ----
  tt <- make_terms(X_pop_clean)
  Xs <- mm_from_terms(tt, df_samp)   # includes intercept
  Xp <- mm_from_terms(tt, X_pop_df)  # includes intercept; same columns as Xs

  glm_fit <- glm.fit(x = Xs, y = df_samp$y, family = binomial())
  X_pop_df$pred_glm <- as.numeric(inv_logit(Xp %*% coef(glm_fit)))

  ## ---- MAP ----
  map_fit <- fit_logit_map_with_priors(df_samp, prior_mean,
                                       sigma_u = 1.0, s_beta = 2.5, maxit = 300)
  cat("MAP converged:", map_fit$converged, "\n")
  X_pop_df$pred_map <- predict_logit_map(map_fit, X_pop_df)

  ## ---- Metrics ----
  truth_vec <- X_pop_df$p_true
  N_vec     <- X_pop_df$N_k

  metrics <- rbind(
    HGP             = compute_metrics_vec(X_pop_df$pred_hgp, truth_vec, N_vec),
    GLM_no_prior    = compute_metrics_vec(X_pop_df$pred_glm, truth_vec, N_vec),
    Logit_withprior = compute_metrics_vec(X_pop_df$pred_map, truth_vec, N_vec)
  )

  ## ---- Plots ----
  p_hgp <- plot_cell_preds(truth_vec, X_pop_df$pred_hgp,
                           paste0("HGP (Stan) | ", scenario), max_points)
  p_glm <- plot_cell_preds(truth_vec, X_pop_df$pred_glm,
                           paste0("GLM (no priors) | ", scenario), max_points)
  p_map <- plot_cell_preds(truth_vec, X_pop_df$pred_map,
                           paste0("Logit MAP priors | ", scenario), max_points)

  ## Return everything needed downstream (including stan_fit for hyperparam means)
  list(
    scenario    = scenario,
    metrics     = metrics,
    plots       = list(p_hgp = p_hgp, p_glm = p_glm, p_map = p_map),
    df_samp     = df_samp,
    X_pop_df    = X_pop_df,
    X_pop_clean = X_pop_clean,
    prior_mean  = prior_mean,
    glm_fit     = glm_fit,
    map_fit     = map_fit,
    stan_fit    = stan_fit
  )
}
