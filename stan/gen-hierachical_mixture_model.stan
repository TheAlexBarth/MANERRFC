
// generalized script to build mixture model of count/biomass models
// Future efforts should include some model checking and comparison to see if a NIB is more appropriate
// Structured as hierachcial with intent of species but may be useful for analysis with sites,etc
// hierarchy is one-level only where N_groups are compared 
// N_groups must be the same for both models but are fed in separately as numeric group_xx
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
    int<lower=0> K_bmass; // biomass preds
    vector[N_obs] img_vol; // imaging volume offset
    matrix[N_obs, K_count] X_count;
    matrix[N_mes, K_bmass] X_bio;

    // group factor
    array[N_obs] int group_counts; // factor id
    array[N_mes] int group_bio;

    // for prediction
    vector[N_obs] y_real;
    matrix[N_obs, K_bmass] Xpred_indv;
}

parameters {
    // count model hyperparams
    vector[K_count] mu_count;
    vector<lower=0>[K_count] sigma_count;
    // biomass mod hypers
    vector[K_bmass] mu_bio;
    vector<lower=0>[K_bmass] sigma_bio;

    // group-level counts covars
    matrix[N_groups, K_count] beta_count;
    // group-level biomass covars
    matrix[N_groups, K_bmass] beta_bio;
    real<lower=0> v; // observation-level sd
}

model {
    // hyperpriors
    mu_count ~ normal(0, 5);
    sigma_count ~ normal(0, 2);
    mu_bio ~ normal(0, 5);
    sigma_bio ~ normal(0, 2);

    // group priors
    for (k in 1:K_count) {
        beta_count[, k] ~ normal(mu_count[k], sigma_count[k]);
    }
    for (k in 1:K_bmass) {
        beta_bio[, k] ~ normal(mu_bio[k], sigma_bio[k]);
    }
    v ~ cauchy(0, 5);

    // count data model
    for (i in 1:N_obs) {
        real lambda = exp(dot_product(beta_count[group_counts[i]], X_count[i,]) + log(img_vol[i]));
        n[i] ~ poisson(lambda);
    }

    // biomass data model
    for (i in 1:N_mes) {
        log_b[i] ~ normal(dot_product(beta_bio[group_bio[i]], X_bio[i,]), v);
    }
}

generated quantities {
    real SSE=0;
    real RMSPE;

    for(i in 1:N_obs){
        real lambda = exp(
            dot_product(beta_count[group_counts[i]], X_count[i]) + log(img_vol[i]) 
        );
        int n_sim = poisson_rng(lambda);

        // this isn't direct but simulate using the posterior params 
        // from the group
        // this only works if all levels are the same!
        real mu = dot_product(
            beta_bio[group_counts[i]], Xpred_indv[i]
        );

        real logbmass_sim = normal_rng(mu, v);
        real y_sim = (n_sim * exp(logbmass_sim)) / img_vol[i];
        SSE += (y_sim - y_real[i])^2;
    }
    RMSPE = sqrt(SSE/N_obs);
}
