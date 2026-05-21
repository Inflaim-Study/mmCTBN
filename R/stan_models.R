# =============================================================================
# stan_models.R  --  Stan model strings for the four prior families
# =============================================================================

.STAN_SPIKE_SLAB <- "
data {
  int<lower=1> N_rows;
  int<lower=1> Q;
  int<lower=1> P;
  array[N_rows] int<lower=0> N_obs;
  matrix[N_rows, Q] X;
  matrix[N_rows, P] Z;
  vector<lower=0>[N_rows] T;
  array[Q] int<lower=0> beta_order;
  real<lower=0, upper=1> pi0;
  real<lower=0>          theta;
  real<lower=0>          a0;
  real<lower=0>          b0;
  real<lower=0>          spike_var;
  int<lower=1>           max_order;
}
transformed data {
  vector[Q] pi_beta;
  for (j in 1:Q)
    pi_beta[j] = pi0 * exp(-theta * beta_order[j]);
  real pi_gamma = pi0;
}
parameters {
  vector[Q]                       beta;
  vector[P]                       gamma;
  vector<lower=0>[max_order + 1]  sigma2;
  real<lower=0>                   sigma2_gamma;
}
model {
  for (p in 0:max_order)
    sigma2[p + 1] ~ inv_gamma(a0, b0);
  sigma2_gamma ~ inv_gamma(a0, b0);

  for (j in 1:Q)
    target += log_mix(pi_beta[j],
                      normal_lpdf(beta[j]  | 0, sqrt(sigma2[beta_order[j] + 1])),
                      normal_lpdf(beta[j]  | 0, sqrt(spike_var)));

  for (k in 1:P)
    target += log_mix(pi_gamma,
                      normal_lpdf(gamma[k] | 0, sqrt(sigma2_gamma)),
                      normal_lpdf(gamma[k] | 0, sqrt(spike_var)));

  vector[N_rows] log_mu = X * beta + Z * gamma + log(T);
  target += poisson_log_lpmf(N_obs | log_mu);
}
generated quantities {
  vector[Q] pip_beta;
  vector[P] pip_gamma;
  for (j in 1:Q) {
    real ls = log(pi_beta[j])   + normal_lpdf(beta[j]  | 0, sqrt(sigma2[beta_order[j] + 1]));
    real lp = log1m(pi_beta[j]) + normal_lpdf(beta[j]  | 0, sqrt(spike_var));
    pip_beta[j] = exp(ls - log_sum_exp(ls, lp));
  }
  for (k in 1:P) {
    real ls = log(pi_gamma)   + normal_lpdf(gamma[k] | 0, sqrt(sigma2_gamma));
    real lp = log1m(pi_gamma) + normal_lpdf(gamma[k] | 0, sqrt(spike_var));
    pip_gamma[k] = exp(ls - log_sum_exp(ls, lp));
  }
}
"

.STAN_STRUCTURED <- "
data {
  int<lower=1> N_rows;
  int<lower=1> Q;
  int<lower=1> P;
  array[N_rows] int<lower=0> N_obs;
  matrix[N_rows, Q] X;
  matrix[N_rows, P] Z;
  vector<lower=0>[N_rows] T;
  array[Q] int<lower=0> beta_order;
  real<lower=0> theta;
  real<lower=0> a0;
  real<lower=0> b0;
  int<lower=1>  max_order;
}
parameters {
  vector[Q]                      beta;
  vector[P]                      gamma;
  vector<lower=0>[max_order + 1] sigma2;
  real<lower=0>                  sigma2_gamma;
}
model {
  for (p in 0:max_order)
    sigma2[p + 1] ~ inv_gamma(a0, b0);
  sigma2_gamma ~ inv_gamma(a0, b0);

  for (j in 1:Q) {
    real eff_var = sigma2[beta_order[j] + 1] * exp(-theta * beta_order[j]);
    beta[j] ~ normal(0, sqrt(eff_var));
  }
  gamma ~ normal(0, sqrt(sigma2_gamma));
  vector[N_rows] log_mu = X * beta + Z * gamma + log(T);
  target += poisson_log_lpmf(N_obs | log_mu);
}
generated quantities {
  vector[Q] post_sd_beta;
  for (j in 1:Q)
    post_sd_beta[j] = sqrt(sigma2[beta_order[j] + 1] * exp(-theta * beta_order[j]));
}
"

.STAN_LASSO <- "
data {
  int<lower=1> N_rows;
  int<lower=1> Q;
  int<lower=1> P;
  array[N_rows] int<lower=0> N_obs;
  matrix[N_rows, Q] X;
  matrix[N_rows, P] Z;
  vector<lower=0>[N_rows] T;
  array[Q] int<lower=0> beta_order;
  real<lower=0> theta;
  real<lower=0> a0;
  real<lower=0> b0;
  int<lower=1>  max_order;
}
parameters {
  vector[Q]                      beta;
  vector[P]                      gamma;
  vector<lower=0>[max_order + 1] lambda2;
  real<lower=0>                  lambda2_gamma;
}
model {
  for (p in 0:max_order)
    lambda2[p + 1] ~ inv_gamma(a0, b0);
  lambda2_gamma ~ inv_gamma(a0, b0);

  for (j in 1:Q) {
    real scale_j = sqrt(lambda2[beta_order[j] + 1]) * exp(-theta * beta_order[j] / 2.0);
    beta[j] ~ double_exponential(0, scale_j);
  }
  gamma ~ double_exponential(0, sqrt(lambda2_gamma));

  vector[N_rows] log_mu = X * beta + Z * gamma + log(T);
  target += poisson_log_lpmf(N_obs | log_mu);
}
generated quantities {
  vector[Q] kappa;
  for (j in 1:Q) {
    real lam2_eff = lambda2[beta_order[j] + 1] * exp(-theta * beta_order[j]);
    kappa[j] = 1.0 / (1.0 + 2.0 * lam2_eff);
  }
}
"

.STAN_HORSESHOE <- "
data {
  int<lower=1> N_rows;
  int<lower=1> Q;
  int<lower=1> P;
  array[N_rows] int<lower=0> N_obs;
  matrix[N_rows, Q] X;
  matrix[N_rows, P] Z;
  vector<lower=0>[N_rows] T;
  array[Q] int<lower=0> beta_order;
  real<lower=0> theta;
  real<lower=0> tau0;
  real<lower=0> slab_df;
  real<lower=0> slab_scale;
  int<lower=1>  max_order;
}
parameters {
  vector[Q]          z_beta;
  vector[P]          gamma;
  vector<lower=0>[Q] lambda;
  real<lower=0>      tau;
  real<lower=0>      c2;
}
transformed parameters {
  vector[Q] beta;
  {
    for (j in 1:Q) {
      real tau_j       = tau * exp(-theta * beta_order[j]);
      real lambda2_j   = square(lambda[j]);
      real lambda2_t   = c2 * lambda2_j / (c2 + square(tau_j) * lambda2_j);
      beta[j]          = z_beta[j] * sqrt(lambda2_t) * tau_j;
    }
  }
}
model {
  tau ~ cauchy(0, tau0);
  c2 ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df * square(slab_scale));
  lambda ~ cauchy(0, 1);
  z_beta ~ std_normal();
  gamma ~ normal(0, 1);
  vector[N_rows] log_mu = X * beta + Z * gamma + log(T);
  target += poisson_log_lpmf(N_obs | log_mu);
}
generated quantities {
  vector[Q] kappa;
  for (j in 1:Q) {
    real tau_j      = tau * exp(-theta * beta_order[j]);
    real lambda2_j  = square(lambda[j]);
    real lambda2_t  = c2 * lambda2_j / (c2 + square(tau_j) * lambda2_j);
    kappa[j] = 1.0 / (1.0 + N_rows * square(tau_j) * lambda2_t);
  }
}
"

# ── Stan model cache lookup ─────────────────────────────────────────────────

.get_stan_model <- function(prior) {
  if (!requireNamespace("rstan", quietly = TRUE))
    stop("ctbn_stan() requires the 'rstan' package.\n",
         "  install.packages('rstan')", call. = FALSE)

  if (!is.null(.mmctbn_env$stan_model_cache[[prior]]))
    return(.mmctbn_env$stan_model_cache[[prior]])

  code <- switch(prior,
    spike_slab = .STAN_SPIKE_SLAB,
    structured = .STAN_STRUCTURED,
    lasso      = .STAN_LASSO,
    horseshoe  = .STAN_HORSESHOE,
    stop("Unknown prior: '", prior, "'.")
  )
  message(sprintf("Compiling Stan model for prior='%s' (once per session)...",
                  prior))
  sm <- rstan::stan_model(model_code = code,
                            model_name = paste0("ctbn_", prior),
                            verbose    = FALSE)
  .mmctbn_env$stan_model_cache[[prior]] <- sm
  sm
}
