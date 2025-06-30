// Revised to only include count predictors
// assumes a possion regression for counts
// just a normal distribution to estimate biomass of trophic groups.

data {
    // constants
    int<lower=0> N_obs; // sampling observations
    int<lower=0> N_mes; // individual measurements
    int<lower=0> K_role; // number of functional roles
    int<lower=0> K_taxo; // number of taxonomic/morphologic groups
    int<lower=0> K_site;
    int<lower=0> K_x; // INCLUDE chlorophyll predictors
    // data terms
    array[N_obs] int n; // counts of individuals
    vector[N_obs] img_vol;
    matrix[N_obs, K_x] W; // chlorophyll MUST BE LAST
    vector[N_mes] b;
    // indexing
    array[N_obs] int role;
    array[N_obs] int site;
    array[N_obs] int incld_chl;
    array[K_taxo] int taxo_role; //mapping role to taxo for prior
    array[N_mes] int taxo;
}

parameters {
    // pois reg
    vector<lower=0>[K_role] gamma;
    array[K_role, K_site, K_x] real alpha;
    matrix[K_role, K_x] mu_alpha;
    matrix<lower=0>[K_role, K_x] sigma_alpha;

    // normal mod
    vector<lower=0>[K_taxo] sigma_g;
    vector<lower=0>[K_taxo] eta_g; // group mean
    vector<lower=0>[K_role] eta_r; // role mean
    vector<lower=0>[K_role] sigma_r;
}

model {
    // count model
    gamma[] ~ gamma(2, 0.01);
    // --- priors ------
    for(k in 1:K_x) {
        mu_alpha[,k] ~ normal(0, 2);
        sigma_alpha[,k] ~ inv_gamma(2,1);
        for(s in 1:K_site) {
            alpha[,s,k] ~ normal(mu_alpha[,k], sigma_alpha[,k]);
        }
    }

    // --- data model ----
    for(i in 1:N_obs) {
        real log_lambda = dot_product(W[i,1:(K_x-3)], to_vector(alpha[role[i], site[i],1:(K_x-3)])) + 
            log(img_vol[i]) +
            incld_chl[i] * dot_product(W[i,(K_x-2):K_x], to_vector(alpha[role[i], site[i],(K_x-2):K_x]));
            ;
        n[i] ~ neg_binomial_2_log(log_lambda, gamma[role[i]]);
    }

    // biomass model
    // this is where stan syntax sucks
    // -- prior ---
    sigma_g[] ~ exponential(1);;
    eta_r ~ normal(1e4, 1e4);
    sigma_r ~ exponential(1);;
    for(g in 1:K_taxo) { //yuck
        eta_g[g] ~ normal(eta_r[taxo_role[g]], sigma_r[taxo_role[g]]);
    }

    // -- data model ----
    for(i in 1:N_mes) {
        b[i] ~ normal(eta_g[taxo[i]], sigma_g[taxo[i]]);
    }
}

generated quantities {
    real MSE_obs = 0;
    real MSE_mod = 0;
    real dev_obs = 0;
    real dev_mod = 0;
 
    for (i in 1:N_obs) {
        real log_lambda = dot_product(W[i,:K_x-4], to_vector(alpha[role[i], site[i],:K_x-4])) +
                          log(img_vol[i]) +
                          incld_chl[i] * dot_product(W[i,K_x-3:], to_vector(alpha[role[i], site[i],K_x-3:]));
        
        // Sample posterior predictive
        int n_sim = neg_binomial_2_rng(exp(log_lambda), gamma[role[i]]);
    
        // MSE (squared error between data and posterior predictive mean)
        MSE_obs += square(n[i] - exp(log_lambda));
        MSE_mod += square(n_sim - exp(log_lambda));
    }

    for(i in 1:N_mes) {
        real b_sim = normal_rng(eta_g[taxo[i]], sigma_g[taxo[i]]);

        dev_obs += -2 * normal_lpdf(b[i] | eta_g[taxo[i]], sigma_g[taxo[i]]);
        dev_mod += -2 * normal_lpdf(b_sim | eta_g[taxo[i]], sigma_g[taxo[i]]);
    }
}
