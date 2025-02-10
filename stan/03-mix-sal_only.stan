// Simple Salinity Full Model


data {
    int<lower=0> N_obs;
    int<lower=0> N_mes; // number of measurements
    array[N_obs] int n; // counts
    vector[N_obs] img_vol; // imaging volume offset
    vector[N_obs] sal; // salinity
    vector[N_mes] log_b; // biomass concentration
    vector[N_mes] sal_b;
    vector[N_mes] temp_b;
    // for model evaluation
    vector[N_obs] obs_biomass_conc;
}

parameters {
    // poisson parms
    real beta0;
    real beta1;
    // biomass params
    real a0;
    real a_sal;
    real b_sig;
}

model {
    beta0 ~ normal(0,5);
    beta1 ~ normal(0,2);
    a0 ~ normal(0,5);
    a_sal ~ normal(0,5);
    b_sig ~ cauchy(0,2);
    
    
    vector[N_obs] nu = exp(beta0 + beta1 * sal + log(img_vol));
    vector[N_mes] b_mu = a0 + a_sal * sal_b;

    for(i in 1:N_obs){
        n[i] ~ poisson(nu[i]);
    }
    for(i in 1:N_mes) {
        log_b[i] ~ normal(b_mu, b_sig);
    }
}

generated quantities {
   real rmse_pred;
   real SSE = 0;
   
   for(i in 1:N_obs) {
        real n_pred = poisson_rng(exp(beta0 + beta1 * sal[i] + log(img_vol[i])));
        real log_b_pred = exp(normal_rng(a0 + a_sal * sal[i], b_sig));
        real pred_biomass_conc = n_pred * log_b_pred;
        SSE += (obs_biomass_conc[i] - pred_biomass_conc) ^2;
   }

   rmse_pred = sqrt(SSE/N_obs);
}