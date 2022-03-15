
// negative binomial and  ASE eQTL with  genotypes  allowing haplotype error, allowing for interaction term between 2 treatments in paired design, with or without covariates. Updated likelihood and mixed prior.
			      
data {
  int<lower=0> N; // number  individuals with NB info
  int<lower=0> A; // number of total individuals ASE info
  int<lower=0> L; // length of vectors with n counts,  p(H) and ai0 
  int<lower=0> K; // number of covariates
  int<lower=0> k; // number of Gaussians for eQTL effect prior
  int Y[N,2]; // total gene counts for each treatment
  int g[N]; // rnsp geno for all individuals
  int s[A]; // number of haplotypes per individual
  int gase[A]; // genotype ASE individuals
  int m[A,2]; // total ase counts for each treatment, entry of 0 if an individual doesnt have ASE ounts
  int n[L,2]; //n counts for ASE ind
  vector[L] pH; //p(H) for ASE ind
  matrix[N,1+K] cov;
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
  real<lower=1e-5,upper=1e5> theta; //the overdispersion parameter for beta binomial
  vector[K-1] betas; // regression parameters
  vector[N] ui; //random term NB
  vector[A] uasei; //random term ASE
  real<lower=1e-5> sdnb; //sd for random term in nb
  real<lower=1e-5> sdase; // sd for random ase term
}

transformed parameters {
  real bp; // parameter of interest for psoriasis
  real bn; // parameter of interest for normal skin

  bp = ba + bd;
  bn = ba - bd;
}

model {
  int pos; // to advance through NB terms (1-G)
  int posl; // to advance through ASE terms (1-L)  

  real p; // ase proportion
  vector[L] ltmp;// log BB likelihood
  real esum; // reduce computation inverse logit (rai0 + bp/bn)
  real esum0; // allelic imbalance proportion under the null
  vector[k] lpsa; // help for mixed gaussians for ba
  vector[k] lpsd; // help for mixed gaussians for bd
  
  vector[N] intercept[2]; // NB intercept
  vector[2] bj; // bp and bn
  vector[2] l1pebj; // transformed bj
  vector[K] betas_all[2]; // combination of treat-specific intercept & global betas
  betas_all[1] = append_row(anorm, betas);
  betas_all[2] = append_row(apso, betas);
  bj[1] = bp;
  bj[2] = bn;
  l1pebj = log1p(exp(bj) - log(2));

  //priors
  theta ~ gamma(1, 0.1); //  based on stan code example
  phi ~ gamma(1, 0.1);
  anorm ~ normal(6, 4); // mean expression normal skin, stan normal is mean and sd
  apso ~ normal(6, 4); // mean expression pso skin
  for(i in 1:(K - 1)){
    betas[i] ~ cauchy(0, 2.5);//prior for the covariates slopes following Gelman 2008
  }

  ui ~ normal(0, sdnb);
  uasei  ~ normal(0, sdase);

  sdnb ~ cauchy(0, 1);
  sdase ~ cauchy(0, 1);

  // mixture of gaussians for ba and bd:
  if (k == 1) {
    ba ~ normal(aveP, sdP);
    bd ~ normal(aveP, sdP);
  } else {

    // mixture of gaussians for bj:
    for(i in 1:k) {
      lpsa[i] = normal_lpdf(ba | aveP[i], sdP[i]) + mixP[i];
      lpsd[i] = normal_lpdf(bd | aveP[i], sdP[i]) + mixP[i];
    }
    target += log_sum_exp(lpsa);
    target += log_sum_exp(lpsd);
  }

  // transformed parameters of no interest
  pos = 1; // to advance on ASE terms
  /* ase = rep_vector(0,Max);  // initialize ase vector to 0s to collect ase termns for each hap pair compatible with Gi=g *\/ */
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

  pos=1;
    
  for (i in 1:A) { //ASE
    for (t in 1:2) {
      if (m[i, t] > 0) {
        esum0 = uasei[i];
        if (t==1) { // first treament	
          esum = inv_logit(esum0 + bp);
        } else {
          esum = inv_logit(esum0 + bn);
        }			   	  
        p = gase[i] == 1 ? esum : inv_logit(esum0);
        p = gase[i] == -1 ? 1-esum : p; //hap swap
        for (r in pos:(pos+s[i]-1)) {
          ltmp[r] = beta_binomial_lpmf(n[r,t] | m[i,t], p * theta, (1 - p) * theta) + log(pH[r]);
        }
        target += log_sum_exp(ltmp[pos:(pos + s[i] - 1)]);
      }
    }
    pos=pos+s[i];
  }
}	     
