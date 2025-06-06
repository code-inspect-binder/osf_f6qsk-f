data{
  int<lower=1> N; 
  int<lower=0,upper=1> Response[N];
  vector[N] bdur;
  vector[N] JumpSize;
  int<lower=0,upper=1> vel1[N];
  real<lower=0> prior_sd;
}
parameters {
  vector[5] beta; 
}
model {
  real mu; 

  //priors
  beta[1:4] ~ normal(0, 10);
  beta[5] ~ normal(0, prior_sd);

  //likelihood
  for (i in 1:N){
    mu = beta[1] + beta[2]*JumpSize[i] + beta[3]*vel1[i] 
         +  beta[4]*JumpSize[i]*vel1[i] + beta[5]*bdur[i];
    target += bernoulli_lpmf(Response[i] | Phi_approx(mu));
  }
}
