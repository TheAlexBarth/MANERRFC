// Simple Salinity Full Model
// need to clean up with matrices rather than all these mus

data {
    int<lower=0> N_obs;
    int<lower=0> N_mes; // number of measurements
    int<lower=0> N_groups;
    array[N_obs] int n; // counts
    vector[N_obs] img_vol; // imaging volume offset
    vector[N_obs] sal; // salinity
    vector[N_obs] temp;
    array[N_obs] int group_counts; // taxonomic group for count data
    array[N_mes] int group_bio; // group for biomass data
    vector[N_mes] log_b; // biomass concentration
    vector[N_mes] sal_b;
    vector[N_mes] temp_b;
    // for model evaluation
    vector[N_obs] obs_biomass_conc;
}

parameters {
    // poisson parms
    // hyperparms for count
    real mu0;
    real mu1;
    real mu2;
    real sigma0;
    real sigma1;
    real sigma2;

    // count model-group parms
    vector[N_groups] beta0;
    vector[N_groups] beta1;
    vector[N_groups] beta2;

    // biomass parms
    // hyperparms
    real theta0;
    real theta1;
    real theta2;
    real tau0;
    real tau1;
    real tau2;
    // group-level
    vector[N_groups] a0;
    vector[N_groups] a1;
    vector[N_groups] a2;
    real v; //individual level var common
}

model {

    // hyperpriors
    mu0 ~ normal(0,5);
    mu1 ~ normal(0,2);
    mu2 ~ normal(0,2);
    sigma0 ~ normal(0,2);
    sigma1 ~ normal(0,2);
    sigma2 ~ normal(0,2);
    theta0 ~ normal(0,5);
    theta1 ~ normal(0,2);
    theta2 ~ normal(0,2);
    tau0 ~ normal(0,5);
    tau1 ~ normal(0,5);
    tau2 ~ normal(0,5);
    

    // group priors
    beta0 ~ normal(mu0, sigma0);
    beta1 ~ normal(mu1, sigma1);
    beta2 ~ normal(mu2, sigma2);
    a0 ~ normal(theta0, tau0);
    a1 ~ normal(theta1, tau1);
    a2 ~ normal(theta2, tau2);
    v ~ cauchy(0,5);

    
    // count data model
    for(i in 1:N_obs) {
        real nu = exp(beta0[group_counts[i]] + beta1[group_counts[i]] * sal[i] + beta2[group_counts[i]] * temp[i] + log(img_vol[i]));
        n[i] ~ poisson(nu);
    }
    for(i in 1:N_mes){
        log_b[i] ~ normal(a0[group_bio[i]] + a1[group_bio[i]] * sal_b[i] + beta2[group_bio[i]] * temp_b[i], v);
    }
}

// generated quantities {
//    real rmse_pred;
//    real SSE = 0;
   
//    for(i in 1:N_obs) {
//         real n_pred = poisson_rng(exp(beta0 + beta1 * sal[i] + log(img_vol[i])));
//         real log_b_pred = exp(normal_rng(a0 + a_sal * sal[i], b_sig));
//         real pred_biomass_conc = n_pred * log_b_pred;
//         SSE += (obs_biomass_conc[i] - pred_biomass_conc) ^2;
//    }

//    rmse_pred = sqrt(SSE/N_obs);
// }