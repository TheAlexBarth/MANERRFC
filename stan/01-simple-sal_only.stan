// log-normal modal for just Salinity

data {
    int<lower=0> N;
    vector[N] y;
    vector[N] sal;
}

parameters {
    real beta0;
    real beta1;
    real<lower=0> sigma;
}

model {
    beta0 ~ normal(0,5);
    beta1 ~ normal(0,2);
    sigma ~ cauchy(0,2);

    y ~ normal(beta0 + beta1 * sal, sigma);
}


generated quantities {
   real rmse_pred;
   real SSE = 0;
   
   for(i in 1:N) {
        real pred_biomass_conc = 10^normal_rng(beta0 + beta1 * sal[i], sigma);
        SSE += (10^y[i] - pred_biomass_conc)^2;
   }

   rmse_pred = sqrt(SSE/N);
}