// negative binomial and  ASE eQTL with unkwown rsnp genotype but fixed fsnps genotypes allowing haplotype error and reference bias correction with random intercept model. Version 2 correcting likelihood so for each individual and each genotype only those hap compatible with g are considered. Allows for any mixture of gaussians for bj prior (eQTL effect).
			      
data {
  int<lower=0> N; // number  individuals with NB info
  int<lower=0> G; // number of total genotypes for all individuals NB
  int<lower=0> K; // number of covariates
  int<lower=0> k; // number of Gaussians for eQTL effect prior
  int Y[N]; // total gene counts
  int sNB[N]; //  number of possible genotypes NB for each individual
  vector[G] gNB; // each geno NB
  vector[G] pNB; // prob for each geno NB
  matrix[N,1+K] cov; 
  vector[k] aveP; // mean for prior Gaussians for eQTL effect prior
  vector[k] sdP; // sd for prior Gaussians for eQTL effect prior
  vector[k] mixP; // log of mixing proportions for eQTL effect prior
}

transformed data {
  vector[G] log_pNB = log(pNB);
  vector[G] abs_gNB = fabs(gNB);
}

parameters {
  vector[K] betas; // regression param
  real<lower=-20,upper=20> bj; // log fold change ASE
  real<lower=1e-5,upper=1e5> phi; //overdipersion param for neg binom
}

model {
  int pos; // to advance through NB terms (1-G)
  vector[N] lmu1; // help to construct linear pred
  real lmu; // linear predictor log scale
  vector[G] ltmp; //  log NB likelihood
  real ebj; // reduce computation
  vector[k] lps; // help for mixed gaussians
  real l1pebj = log1p(ebj) - log(2);

  //priors
  phi ~ gamma(1, 0.1);
  betas[1] ~ normal(6, 4); // stan normal is mean and sd
  for(i in 2:K){
    betas[i] ~ cauchy(0, 2.5);//prior for the slopes following Gelman 2008
  }

  // mixture of gaussians for bj:
  for(i in 1:k){
    lps[i] = normal_lpdf(bj | aveP[i], sdP[i]) + mixP[i];
  }
  target += log_sum_exp(lps);
  // local variables and transformed parameters of no interest
  pos = 1;  // to advance on NB terms
  
  ebj = exp(bj);
  lmu1 = cov[, 2:cols(cov)] * betas;
  for(i in 1:N) { // lmu for each individual default to GT=0
    
    for (r in pos:(pos+sNB[i]-1)) { // genotype level, Gi=g
      // print("gNB = ", gNB[r]," r = ", r);
      lmu = lmu1[i];
      lmu = abs_gNB[r] == 1 ? lmu + l1pebj : lmu;
      lmu = gNB[r] == 2 ? lmu + bj : lmu;

      ltmp[r] = neg_binomial_2_log_lpmf(Y[i] | lmu, phi) + log_pNB[r];
    }
    target += log_sum_exp(ltmp[pos:(pos+sNB[i]-1)]);
    pos += sNB[i];
  }
}
