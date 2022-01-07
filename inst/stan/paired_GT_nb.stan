
// negative binomial and  ASE eQTL with unkwown  genotypes  allowing haplotype error, allowing for interaction term between 2 treatments in paired design, with or without covariates and ref bias correction version 2. Updated likelihood and mixed prior.
			      
data {
  int<lower=0> N; // number  individuals with NB info
  int<lower=0> K; // number of covariates
  int<lower=0> k; // number of Gaussians for eQTL effect prior
  int Y[N, 2]; // total gene counts for each treatment
  int g[N]; // rnsp geno for all individuals
  matrix[N, 1+K] cov;
  vector[k] aveP; // mean for prior Gaussians for eQTL effect prior
  vector[k] sdP; // sd for prior Gaussians for eQTL effect prior
  vector[k] mixP; // log of mixing proportions for eQTL effect prior
}

parameters {
  real anorm; // mean expression normal tissue
  real apso; // mean expression pso tissues
  real ba; // log average-fold change ASE
  real bd; // log difference-fold change ASE
  real<lower=1e-5,upper=1e5> phi; //overdipersion param for neg binom
  vector[K-1] betas; // regression parameters 
  vector[N] ui; //random term NB
  real<lower=1e-5,upper=1e5> sdnb; //sd for random term in nb
}

transformed parameters {
  real bp; // parameter of interest for psoriasis
  real bn; // parameter of interest for normal skin
  
  bp = ba + bd;
  bn = ba - bd;
}

model {
  // prior temp vars
  vector[k] lpsa; // help for mixed gaussians for ba
  vector[k] lpsd; // help for mixed gaussians for bd
  // likelihood temp vars
  vector[N] intercept[2];
  vector[2] l1pebj;
  vector[2] bj;
  vector[K] betas_all[2];
  betas_all[1] = append_row(anorm, betas);
  betas_all[2] = append_row(apso, betas);
  bj[1] = bp;
  bj[2] = bn;
  l1pebj = log1p(exp(bj) - log(2));    

  //priors
  phi ~ gamma(1, 0.1); //  based on stan code example
  
  anorm ~ normal(6, 4); // mean expression normal skin, stan normal is mean and sd
  apso ~ normal(6, 4); // mean expression pso skin
  for (i in 1:(K-1)) {
    betas[i] ~ cauchy(0, 2.5);//prior for the covariates slopes following Gelman 2008
  }

  // random effect priors
  sdnb ~ cauchy(0, 1); 
  ui ~ normal(0, sdnb);

  // mixture of gaussians for ba and bd:
  for(i in 1:k) {
    lpsa[i] = normal_lpdf(ba | aveP[i], sdP[i]) + mixP[i];
    lpsd[i] = normal_lpdf(bd | aveP[i], sdP[i]) + mixP[i];
  }
  target += log_sum_exp(lpsa);
  target += log_sum_exp(lpsd);

  for (t in 1:2) {
    intercept[t] = ui;
    for (i in 1:N) { // log1p(exp(b_j)) - log(2) if het or bj if hom
      if (fabs(g[i]) == 1) {
        intercept[t][i] = l1pebj[t];
      }
      if (g[i] == 2) {
        intercept[t][i] = bj[t];
      }
    }
    Y[, t] ~ neg_binomial_2_log_glm(
      cov[, 2:cols(cov)],
      intercept[t],
      betas_all[t],
      phi
    );
  }
}
