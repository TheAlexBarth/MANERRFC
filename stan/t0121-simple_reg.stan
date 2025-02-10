data{
    int<lower=0> N;
    int<lower=0> K; // n pred
    matrix[N, K] x_covars;
    vector[N] y;
}

parameters {
   vector[K] beta;
   real sigma;
}

model {
    for(i in 1:K) {
        beta[i] ~ normal(0,15);
    }
    sigma ~ cauchy(0,5);

    y ~ normal(x_covars * beta, sigma);
}