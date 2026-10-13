# ============================================================
# 03_compute_weights.R
# Companion code for the HGP survey-weighting paper
# ============================================================

# This script assumes 01_simulation_setup.R and 02_fit_models.R have been sourced.
# It constructs the HGP, GLM, and MAP equivalent weights and reports ESS.

inv_logit <- function(x) 1 / (1 + exp(-x))

rbf_kernel <- function(Fmat, alpha = 1.5, rho = 0.8, jitter = 1e-6) {
  D2 <- as.matrix(dist(Fmat))^2
  K  <- alpha^2 * exp(-D2 / (2 * rho^2))
  K + diag(jitter, nrow(Fmat))
}

normalize_to_sum <- function(w, target_sum) {
  w * (target_sum / sum(w))
}

handle_negative_weights <- function(w,
                                    method = c("none", "min_positive"),
                                    eps = 1e-12) {
  method <- match.arg(method)
  if (method == "none") return(w)

  if (any(w <= 0)) {
    w_pos <- w[w > eps]
    if (length(w_pos) == 0) stop("All weights are non-positive; cannot apply min_positive.")
    w[w <= eps] <- min(w_pos)
  }
  w
}

ess_from_w <- function(w) (sum(w)^2) / sum(w^2)


## ============================================================
## (1) Logistic equivalent weights from GLM MLE
##     w = (n/N) * (N_pop^T W_pop X_pop) (X^T W X)^(-1) X^T
##     (then renormalize to sum n)
## ============================================================
weights_from_glm <- function(glm_fit,
                             X_samp,      # n x P (incl intercept if used in fit)
                             X_pop,       # n_pop x P (same columns as X_samp)
                             N_k,         # length n_pop (cell counts)
                             target_sum = nrow(X_samp),
                             neg_method = c("min_positive", "none"),
                             ridge = 0) {

  neg_method <- match.arg(neg_method)

  beta <- coef(glm_fit)
  if (length(beta) != ncol(X_samp)) {
    stop("coef(glm_fit) length != ncol(X_samp). Make sure X_samp matches the fitted design.")
  }

  N_tot <- sum(N_k)
  n     <- nrow(X_samp)

  # sample weights W = p(1-p)
  p_s   <- as.numeric(inv_logit(X_samp %*% beta))
  W_s   <- p_s * (1 - p_s)

  # pop weights W_pop = p_pop(1-p_pop)
  p_p   <- as.numeric(inv_logit(X_pop %*% beta))
  W_p   <- p_p * (1 - p_p)

  # H = X^T W X
  # compute X^T (W * X) without forming diag(W)
  WX_s  <- X_samp * W_s
  H     <- crossprod(X_samp, WX_s)

  if (ridge > 0) H <- H + diag(ridge, nrow(H))

  # solve H^{-1} X^T  (P x n)
  Hinv_Xt <- solve(H, t(X_samp))

  # left = (N_pop^T W_pop X_pop)  (1 x P)
  left <- as.numeric(crossprod(N_k, X_pop * W_p))  # vector length P

  # w_raw (length n): (n/N) * left %*% H^{-1} X^T
  w_raw <- as.numeric((n / N_tot) * (left %*% Hinv_Xt))

  # negative-weight handling + renormalize to sum target_sum
  w_raw <- handle_negative_weights(w_raw, method = neg_method)
  w     <- normalize_to_sum(w_raw, target_sum)

  list(w = w, ess = ess_from_w(w), w_raw = w_raw)
}

weights_from_glm <- function(glm_fit,
                             X_samp,      # n x P (incl intercept if used in fit)
                             X_pop,       # n_pop x P (same columns as X_samp)
                             N_k,         # length n_pop (cell counts)
                             target_sum = nrow(X_samp),
                             neg_method = c("min_positive", "none"),
                             ridge = 0) {

  neg_method <- match.arg(neg_method)

  beta <- coef(glm_fit)
  if (length(beta) != ncol(X_samp)) {
    stop("coef(glm_fit) length != ncol(X_samp). Make sure X_samp matches the fitted design.")
  }

  N_tot <- sum(N_k)
  n     <- nrow(X_samp)

  # sample weights W = p(1-p)
  p_s   <- as.numeric(inv_logit(X_samp %*% beta))
  W_s   <- p_s * (1 - p_s)

  # pop weights W_pop = p_pop(1-p_pop)
  p_p   <- as.numeric(inv_logit(X_pop %*% beta))
  W_p   <- p_p * (1 - p_p)

  # H = X^T W X
  # compute X^T (W * X) without forming diag(W)
  WX_s  <- X_samp * W_s
  H     <- crossprod(X_samp, WX_s)

  if (ridge > 0) H <- H + diag(ridge, nrow(H))

  # solve H^{-1} X^T  (P x n)
  Hinv_Xt <- solve(H, t(X_samp))

  # left = (N_pop^T W_pop X_pop)  (1 x P)
  left <- as.numeric(crossprod(N_k, X_pop * W_p))  # vector length P

  # w_raw (length n): (n/N) * left %*% H^{-1} X^T
  w_raw <- as.numeric((n / N_tot) * (left %*% Hinv_Xt))

  # negative-weight handling + renormalize to sum target_sum
  w_raw <- handle_negative_weights(w_raw, method = neg_method)
  w     <- normalize_to_sum(w_raw, target_sum)

  list(w = w, ess = ess_from_w(w), w_raw = w_raw)
}

weights_from_map_fit <- function(map_fit,
                                 df_samp,
                                 X_pop_df,
                                 N_k,
                                 target_sum = nrow(df_samp),
                                 neg_method = c("min_positive", "none"),
                                 ridge = 1e-6) {

  neg_method <- match.arg(neg_method)

  ## Sample design
  Xs <- model.matrix(~ city + age + ethn + edu + mar, data = df_samp)
  Xs <- Xs[, colnames(Xs) != "(Intercept)", drop = FALSE]
  Xs <- Xs[, map_fit$X_cols, drop = FALSE]

  ## Population design
  Xp <- model.matrix(~ city + age + ethn + edu + mar, data = X_pop_df)
  Xp <- Xp[, colnames(Xp) != "(Intercept)", drop = FALSE]
  Xp <- Xp[, map_fit$X_cols, drop = FALSE]

  h_s <- as.integer(df_samp$group_H)
  h_p <- as.integer(X_pop_df$group_H)

  ## Linear predictors (u_h treated as fixed offsets)
  eta_s <- as.numeric(Xs %*% map_fit$beta + map_fit$u[h_s])
  eta_p <- as.numeric(Xp %*% map_fit$beta + map_fit$u[h_p])

  p_s <- inv_logit(eta_s)
  p_p <- inv_logit(eta_p)

  W_s <- p_s * (1 - p_s)
  W_p <- p_p * (1 - p_p)

  ## Fisher information for β only
  WXs <- Xs * W_s
  H   <- crossprod(Xs, WXs)

  ## ridge is REQUIRED for numerical stability
  H <- H + diag(ridge, nrow(H))

  Hinv_Xt <- solve(H, t(Xs))

  ## population score
  left <- as.numeric(crossprod(N_k, Xp * W_p))

  n     <- nrow(Xs)
  N_tot <- sum(N_k)

  w_raw <- as.numeric((n / N_tot) * (left %*% Hinv_Xt))

  w_raw <- handle_negative_weights(w_raw, method = neg_method)
  w     <- normalize_to_sum(w_raw, target_sum)

  list(w = w, ess = ess_from_w(w), w_raw = w_raw)
}

gp_kernel_sqexp <- function(X1, X2 = NULL, alpha, rho, jitter = 0) {
  X1 <- as.matrix(X1)
  if (is.null(X2)) X2 <- X1
  X2 <- as.matrix(X2)

  # squared distances ||x - x'||^2 = ||x||^2 + ||x'||^2 - 2 x x'
  s1 <- rowSums(X1^2)
  s2 <- rowSums(X2^2)
  D2 <- outer(s1, s2, "+") - 2 * (X1 %*% t(X2))
  D2 <- pmax(D2, 0)  # numerical safety

  K <- (alpha^2) * exp(-0.5 * D2 / (rho^2))

  # only add jitter when it's a square kernel (X2 = X1)
  if (jitter > 0 && nrow(X1) == nrow(X2) && ncol(X1) == ncol(X2)) {
    diag(K) <- diag(K) + jitter
  }
  K
}

compute_hgp_w_cg <- function(res,
                             sigma2 = 1e-6,
                             tol = 1e-8,
                             maxit = 2000,
                             normalize = TRUE,
                             eps_deriv = 1e-12,
                             p_hat_samp = NULL) {

  df_samp  <- res$df_samp
  X_pop_df <- res$X_pop_df
  n        <- nrow(df_samp)

  if (is.null(X_pop_df$N_k)) stop("X_pop_df must contain N_k.")
  N_total <- sum(X_pop_df$N_k)

  ## ---- 1) posterior means of hyperparams ----
  draws <- res$stan_fit$draws(c("alpha","rho","tau_mu"), format="df")
  alpha_hat <- mean(draws$alpha)
  rho_hat   <- mean(draws$rho)
  tau_hat   <- mean(draws$tau_mu)

  ## ---- 2) observed K-cells in SAME order ----
  obs_cells_full <- sort(unique(df_samp$group_K_full))
  Kobs <- length(obs_cells_full)

  X_unique_df <- res$X_pop_clean %>%
    dplyr::filter(group_K_full %in% obs_cells_full) %>%
    dplyr::mutate(.ord = match(group_K_full, obs_cells_full)) %>%
    dplyr::arrange(.ord)

  X_unique <- make_design_dense(X_unique_df)   # Kobs x P
  X_pop    <- make_design_dense(X_pop_df)      # n_pop x P

  ## ---- 3) kernels ----
  K_KK   <- gp_kernel_sqexp(X_unique, NULL,  alpha_hat, rho_hat, jitter = 1e-9)
  K_Kpop <- gp_kernel_sqexp(X_unique, X_pop, alpha_hat, rho_hat, jitter = 0)

  ## ---- 4) D_* N^pop vector at population level (logit case) ----
  if (is.null(X_pop_df$pred_hgp)) {
    stop("X_pop_df must contain pred_hgp (population-cell fitted probabilities).")
  }
  p_star <- X_pop_df$pred_hgp
  d_vec  <- X_pop_df$N_k * p_star * (1 - p_star)  # = D_* N^pop

  ## ---- 5) b = Sigma_{y,f_*} (D_* N^pop) ----
  b_cell_gp <- as.numeric(K_Kpop %*% d_vec)

  k_index <- match(df_samp$group_K_full, obs_cells_full)
  if (anyNA(k_index)) stop("Some sampled units have group_K_full not found in obs_cells_full.")
  b_gp <- b_cell_gp[k_index]

  d_by_H_pop <- tapply(d_vec, X_pop_df$group_H, sum)
  d_by_H_pop[is.na(d_by_H_pop)] <- 0

  h_i <- df_samp$group_H
  b_mu <- (tau_hat^2) * as.numeric(d_by_H_pop[as.character(h_i)])
  b_mu[is.na(b_mu)] <- 0

  b <- b_gp + b_mu

  ## ---- 6) fast matvec: v -> Sigma_yy v ----
  Sigma_mv <- function(v) {
    v_by_K <- as.numeric(tapply(v, k_index, sum))
    if (length(v_by_K) != Kobs) {
      tmp <- numeric(Kobs); tmp[as.integer(names(v_by_K))] <- v_by_K
      v_by_K <- tmp
    }
    gp_by_K <- as.numeric(K_KK %*% v_by_K)
    gp_part <- gp_by_K[k_index]

    v_by_H <- tapply(v, h_i, sum)
    mu_part <- (tau_hat^2) * as.numeric(v_by_H[as.character(h_i)])
    mu_part[is.na(mu_part)] <- 0

    mu_part + gp_part + sigma2 * v
  }

  ## ---- 7) conjugate gradient solve u = Sigma_yy^{-1} b ----
  cg_solve <- function(A_mv, b, tol, maxit) {
    x <- rep(0, length(b))
    r <- b - A_mv(x)
    p <- r
    rsold <- sum(r*r)
    for (it in seq_len(maxit)) {
      Ap <- A_mv(p)
      denom <- sum(p*Ap)
      if (!is.finite(denom) || denom <= 0) break
      alpha <- rsold / denom
      x <- x + alpha*p
      r <- r - alpha*Ap
      rsnew <- sum(r*r)
      if (sqrt(rsnew) < tol) break
      p <- r + (rsnew/rsold)*p
      rsold <- rsnew
    }
    x
  }

  u <- cg_solve(Sigma_mv, b, tol = tol, maxit = maxit)
  ## u = Sigma_yy^{-1} Sigma_{y,f_*} D_* N^pop  (A-free piece)

  ## ---- 8) Get unit fitted probs p_hat (for A^{-1}) ----
  if (is.null(p_hat_samp)) {
    if (!is.null(df_samp$pred_hgp)) {
      p_hat_samp <- df_samp$pred_hgp
    } else if (!is.null(df_samp$pred_map)) {
      p_hat_samp <- df_samp$pred_map
    } else if (!is.null(res$glm_fit)) {
      # fallback: use GLM fit on the sample design matrix if available
      tt <- make_terms(res$X_pop_clean)
      X_samp <- mm_from_terms(tt, df_samp)
      eta_hat <- as.numeric(X_samp %*% coef(res$glm_fit))
      p_hat_samp <- inv_logit(eta_hat)
    } else {
      stop("Need unit fitted probabilities for A^{-1}. Provide p_hat_samp, or add df_samp$pred_map/pred_hgp, or provide res$glm_fit.")
    }
  }

  ## ---- 9) Apply Case-B correction: w = (n/N) A^{-1} u ----
  a_diag <- p_hat_samp * (1 - p_hat_samp)          # logit derivative
  a_diag <- pmax(a_diag, eps_deriv)

  w_raw <- (n / N_total) * (u / a_diag)

  w_raw <- handle_negative_weights(w_raw, method = "min_positive")
  w     <- if (normalize) normalize_to_sum(w_raw, n) else w_raw

  list(w = w, w_raw = w_raw, ess = ess_from_w(w),
       alpha_hat = alpha_hat, rho_hat = rho_hat, tau_hat = tau_hat,
       N_total = N_total)
}

compute_all_weights <- function(res) {
  N_k <- res$X_pop_df$N_k
  n   <- nrow(res$df_samp)

  ## HGP via CG
  hgp <- compute_hgp_w_cg(res)

  ## GLM matrices aligned
  tt <- make_terms(res$X_pop_clean)
  X_samp <- mm_from_terms(tt, res$df_samp)
  X_pop  <- mm_from_terms(tt, res$X_pop_df)

  glm <- weights_from_glm(res$glm_fit, X_samp, X_pop, N_k, target_sum = n)
  map <- weights_from_map_fit(res$map_fit, res$df_samp, res$X_pop_df, N_k, target_sum = n)

  list(
    w_hgp = hgp$w, ess_hgp = hgp$ess,
    w_glm = glm$w, ess_glm = glm$ess,
    w_map = map$w, ess_map = map$ess
  )
}
