// negative  binomial for  eQTL with fixed genotypes. Allows for any mixture of gaussians for bj prior (eQTL effect).
 
data {
  int<lower=0> N; // number  individuals
  int<lower=0> K; // number of covariates
  int<lower=0> k; // number of Gaussians for eQTL effect prior
  int Y[N]; // total gene counts
  int g[N]; // rnsp geno for all individuals
  matrix[N,1+K] cov;
  vector[k] aveP; // mean for prior Gaussians for eQTL effect prior
  vector[k] sdP; // sd for prior Gaussians for eQTL effect prior
  vector[k] mixP; // log of mixing proportions for eQTL effect prior
  // simplex[k] expMixP; // log of mixing proportions for eQTL effect prior
}

parameters {
  vector[K] betas; // regression param
  real <lower=-10,upper=10> bj; // log fold change ASE
  // real bj; // log fold change ASE
  real <lower=1e-5,upper=1e5> phi; //overdipersion param for neg binom
}


model {
  // include transformed parameters of no interest
  vector[k] lps; // help for mixed gaussians
  real l1pebj = log1p(exp(bj)) - log(2);
  vector[N] intercept = rep_vector(1e-5, N); // the genetic effect; 0 if hom ref
  
  // Priors
  phi ~ gamma(1, 0.01);

  betas[1] ~ normal(6, 4); // stan normal is mean and sd
  for(i in 2:K) {
    betas[i] ~ cauchy(0, 2.5); //prior for the slopes following Gelman 2008   
  }

  if (k == 1) {
    bj ~ normal(aveP, sdP);
  } else {
    // mixture of gaussians for bj:
    for(i in 1:k) {
      lps[i] = normal_lpdf(bj | aveP[i], sdP[i]) + mixP[i];
      // lps[i] = normal_lpdf(bj | aveP[i], sdP[i]);
    }
    target += log_sum_exp(lps);
    // suggestion from stan team discussions that one of these should be faster,
    // but they don't seem to be right now
    // target += log_sum_exp(lps + mixP);
    // target += log_mix(expMixP, lps);
  }
  

  for (i in 1:N) { // log1p(exp(b_j)) - log(2) if het or bj if hom
    if (fabs(g[i]) == 1) {
      intercept[i] = l1pebj; // log1p(exp(b_afc)) - log(2) if he
    }
    if (g[i] == 2) {
      intercept[i] = bj; // b_afc if hom alt
    }
  }
  Y ~ neg_binomial_2_log_glm(cov[, 2:cols(cov)], intercept, betas, phi);
}
