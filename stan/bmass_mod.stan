data {
    int<lower=0> N_mes; // individual measurements
    int<lower=0> K_role; // number of functional roles
    int<lower=0> K_taxo; // number of taxonomic/morphologic groups
    vector[N_mes] b;
    array[K_taxo] int taxo_role; //mapping role to taxo for prior
    array[N_mes] int taxo;
}

parameters {
    // normal mod
    vector<lower=0>[K_taxo] sigma_g;
    vector<lower=0>[K_taxo] eta_g; // group mean
    vector<lower=0>[K_role] eta_r; // role mean
    vector<lower=0>[K_role] sigma_r;
}

model {
       // biomass model
    // this is where stan syntax sucks
    // -- prior ---
    sigma_g[] ~ exponential(1);;
    eta_r ~ normal(7, 2);
    sigma_r ~ exponential(0.5);;
    for(g in 1:K_taxo) { //yuck
        eta_g[g] ~ normal(eta_r[taxo_role[g]], sigma_r[taxo_role[g]]);
    }

    // -- data model ----
    for(i in 1:N_mes) {
        log(b[i]) ~ normal(eta_g[taxo[i]], sigma_g[taxo[i]]);
    } 
}

generated quantities {
    real dev_obs = 0;
    real dev_mod = 0;

    for(i in 1:N_mes) {
        real logb_sim = normal_rng(eta_g[taxo[i]], sigma_g[taxo[i]]);

        dev_obs += -2 * normal_lpdf(log(b[i]) | eta_g[taxo[i]], sigma_g[taxo[i]]);
        dev_mod += -2 * normal_lpdf(logb_sim | eta_g[taxo[i]], sigma_g[taxo[i]]);
    }

}