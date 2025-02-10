// Simple Salinity Hurdle Model


data {
    int<lower=0> N_obs;
    array[N_obs] int z; // observations of 0/1
    vector[N_obs] y; // biomass concentration
    vector[N_obs] sal;
}

parameters {
    // hurdle parameters
    real alpha0;
    real alpha1;
    // regression params
    real beta0; 
    real beta1;
    real<lower=0> sigma;
}

model {
    alpha0 ~ normal(0,5);
    alpha1 ~ normal(0,5);
    beta0 ~ normal(0,5);
    beta1 ~ normal(0,2);
    sigma ~ cauchy(0,2);
 
    vector[N_obs] pi = alpha0 + alpha1 * sal;

    for(n in 1:N_obs){
        z[n] ~ bernoulli_logit(pi[n]);
        if(z[n] == 1) {
            y[n] ~ normal(beta0 + beta1 * sal[n], sigma);
        }
    }
}

generated quantities {
    real rmse_pred;
    real SSE = 0;
   
    for(i in 1:N_obs) {
        real pi_pred = 1/(1+exp(-(alpha0 + alpha1 * sal[i])));
        real y_pred = pi_pred * 10^(beta0 + beta1 * sal[i]);
        SSE += (10^y[i] - y_pred)^2;
    }

   rmse_pred = sqrt(SSE/N_obs);
}