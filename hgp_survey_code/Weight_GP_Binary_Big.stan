functions {
  // Squared exponential kernel
  matrix gp_group_kernel(matrix X_group, real alpha, real rho, real delta) {
    int K = rows(X_group);
    matrix[K, K] K_mat;

    for (i in 1:K) {
      for (j in i:K) {
        real sq_dist = squared_distance(X_group[i], X_group[j]);
        K_mat[i, j] = square(alpha) * exp(-0.5 * sq_dist / square(rho));
        if (i != j) K_mat[j, i] = K_mat[i, j];
      }
    }
    for (i in 1:K) K_mat[i, i] += delta;
    return K_mat;
  }

  matrix gp_group_kernel_cross(matrix X1, matrix X2, real alpha, real rho) {
    int n1 = rows(X1);
    int n2 = rows(X2);
    matrix[n1, n2] K_cross;

    for (i in 1:n1) {
      for (j in 1:n2) {
        real sq_dist = squared_distance(X1[i], X2[j]);
        K_cross[i, j] = square(alpha) * exp(-0.5 * sq_dist / square(rho));
      }
    }
    return K_cross;
  }
}

data {
  // Training data
  int<lower=1> n;
  int<lower=1> P;
  array[n] int<lower=1> group_H;
  array[n] int<lower=1> group_K;
  int<lower=1> H;
  int<lower=1> K;
  array[n] int<lower=0, upper=1> y;

  matrix[n, P] X;
  matrix[K, P] X_unique;
  vector[H] prior_mean;

  // Prediction data
  int<lower=1> n_pop;
  matrix[n_pop, P] X_pop;
  array[n_pop] int<lower=1> group_H_pop;
}

transformed data {
  real delta = 1e-9;
}

parameters {
  // H-level offsets around prior_mean
  vector[H] mu_raw;
  real<lower=0> tau_mu;

  // GP hyperparameters + latent cell GP
  real<lower=0> alpha;
  real<lower=0> rho;
  vector[K] eta;
}

transformed parameters {
  vector[H] mu = prior_mean + tau_mu * mu_raw;

  matrix[K, K] K_eta = gp_group_kernel(X_unique, alpha, rho, delta);
  matrix[K, K] L_K_eta = cholesky_decompose(K_eta);
}

model {
  // priors
  mu_raw ~ normal(0, 1);
  tau_mu ~ cauchy(0, 2.5);

  alpha ~ normal(0, 1);
  rho ~ inv_gamma(5, 5);

  target += multi_normal_cholesky_lpdf(eta | rep_vector(0, K), L_K_eta);

  // likelihood
  for (i in 1:n)
    y[i] ~ bernoulli_logit(mu[group_H[i]] + eta[group_K[i]]);
}
generated quantities {
  vector[n_pop] E_y_star;

  // KxK GP on observed cells
  matrix[K, K] K_obs = gp_group_kernel(X_unique, alpha, rho, delta);
  matrix[K, K] L_obs = cholesky_decompose(K_obs);

  // n_pop x K cross kernel (pop vs observed cells)
  matrix[n_pop, K] K_cross = gp_group_kernel_cross(X_pop, X_unique, alpha, rho);

  // K_obs^{-1} * eta via Cholesky
  vector[K] Kinv_eta;
  Kinv_eta = mdivide_left_tri_low(L_obs, eta);
  Kinv_eta = mdivide_right_tri_low(Kinv_eta', L_obs)';

  // latent GP mean at pop cells
  vector[n_pop] f_star_mean = K_cross * Kinv_eta;

  for (j in 1:n_pop)
    E_y_star[j] = inv_logit(mu[group_H_pop[j]] + f_star_mean[j]);
}

