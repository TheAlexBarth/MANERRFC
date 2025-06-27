//script for chl-regression
// I try to orient similar to NIMBLE/JAGS with definition of const and data
data {
    // constant term:
    int<lower=0> N_obs; 
    int<lower=0> K_frac;
    int<lower=0> K_site;
    int<lower=0> K_x; // number of covariate
    // data terms
    vector[N_obs] y; // chl-a
    array[N_obs] int frac; // grouping factor
    array[N_obs] int site;
    matrix[N_obs, K_x] X; // covariates
}

parameters {
    vector<lower=0>[K_frac] sigma_f;
    array[K_frac, K_site, K_x] real beta;
    matrix[K_frac, K_x] mu;
    matrix<lower=0>[K_frac, K_x] sigma_beta;
}

model{
    // priors
    for(k in 1:K_x) {
        mu[,k] ~ normal(0, 2);
        sigma_beta[,k] ~ inv_gamma(2,1);
        for(s in 1:K_site) {
            beta[,s,k] ~ normal(mu[,k],sigma_beta[,k]);
        }
    }
    sigma_f[] ~ inv_gamma(2,1);
    //models
    for(i in 1:N_obs) {
         real mu_term = dot_product(X[i,], to_vector(beta[frac[i],site[i],]));
        y[i] ~ normal(mu_term, sigma_f[frac[i]]);
    }
}

generated quantities {
   real dev_obs = 0;
   real dev_sim = 0;

   for(i in 1:N_obs) {
    real mu_term = dot_product(X[i,], to_vector(beta[frac[i],site[i],]));
    real y_sim = normal_rng(mu_term, sigma_f[frac[i]]);
    real log_lik = normal_lpdf(y[i] | mu_term, sigma_f[frac[i]]);
    real log_lik_sim = normal_lpdf(y_sim | mu_term , sigma_f[frac[i]]);
    dev_obs += -2 * log_lik;
    dev_sim += -2 * log_lik_sim;
   }
}