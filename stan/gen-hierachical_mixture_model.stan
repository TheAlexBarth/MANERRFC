data {
    // format
    int<lower=0> N_obs; // sampling observations
    int<lower=0> N_mes; // individual measurements
    int<lower=0> N_groups; // number of groups

    // response inputs
    array[N_obs] int n; // n_counts
    vector[N_mes] log_b; // biomass of individuals (logged)

    // predictor inputs
    int<lower=0> K_count; // count preds
    vector[N_obs] img_vol; // imaging volume offset
    matrix[N_obs, K_count] X_count;

    // group factor
    array[N_obs] int group_counts; // factor id
    array[N_mes] int group_bio;

    // programmable priors
    real theta_b_mu; // prior on theta_b
    real theta_b_sig; // prior var for theta_b
    real tau_b_mu; // uniform lower
    real tau_b_sig; // uniform higher
}

parameters {
    // count model hyperparams
    vector[K_count] theta_count;
    vector<lower=0>[K_count] tau_count;
    // biomass mod hypers
    real theta_bio;
    real<lower=0> tau_bio;

    // group-level counts covars
    matrix[N_groups, K_count] beta_count;
    vector[N_groups] mu_bio;
    // biomass obs
    real<lower=1e-6> sigma_bio;
    // count vars
    real<lower=1e-6> phi;
}

model {
    // hyperpriors
    theta_count ~ normal(0, 1);
    tau_count ~ normal(1,2);
    theta_bio ~ normal(theta_b_mu, theta_b_sig);// need to set to scale of data
    tau_bio ~ normal(tau_b_mu, tau_b_sig); // set to scale of data
    phi ~ gamma(2, 0.1); 


    // group priors
    for (k in 1:K_count) {
        beta_count[, k] ~ normal(theta_count[k], tau_count[k]);
    }

    for(k in 1:N_groups) {
        mu_bio[k] ~ normal(theta_bio, tau_bio);
    }

    // count data model
    for (i in 1:N_obs) {
        real log_lambda = dot_product(beta_count[group_counts[i]], X_count[i,]) + log(img_vol[i]);
        n[i] ~ neg_binomial_2_log(log_lambda, phi);
    }

    // biomass data model
    for (i in 1:N_mes) {
        log_b[i] ~ normal(mu_bio[group_bio[i]], sigma_bio);
    }
}

generated quantities {
    real MSE_obs_n = 0;
    real MSE_mod_n = 0;
    real MSE_obs_b = 0;
    real MSE_mod_b = 0;
 
    for(i in 1:N_obs){

        real log_lambda = dot_product(beta_count[group_counts[i]], X_count[i]) + log(img_vol[i]);
        int n_tilde = neg_binomial_2_log_rng(log_lambda, phi);

        //note to self here - group_counts[i] needs to match the right index for
        // this all to work!!
        // noticed this error and should re-check the pvalue calcs with negbin vs pois

        MSE_obs_n += square(n[i] - exp(log_lambda));
        MSE_mod_n += square(n_tilde - exp(log_lambda));


    }

    for(i in 1:N_mes) {
        real b_tilde = normal_rng(mu_bio[group_bio[i]], sigma_bio);
        MSE_obs_b += square(log_b[i] - mu_bio[group_bio[i]]);
        MSE_mod_b += square(b_tilde - mu_bio[group_bio[i]]);
    }

    MSE_obs_n /= N_obs;
    MSE_mod_n /= N_obs;
    MSE_obs_b /= N_mes;
    MSE_mod_b /= N_mes;
}
