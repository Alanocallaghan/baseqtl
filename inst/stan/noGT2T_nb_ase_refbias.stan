// negative binomial and  ASE eQTL with unkwown rsnp genotype but fixed fsnps genotypes allowing haplotype error, allowing for interaction term normal psoriasis, with or without covariates and ref bias correction version 2. Updated likelihood and mixed prior.
			      
data {
  int<lower=0> N; // number  individuals with NB info
  int<lower=0> G; // number of total genotypes for all individuals NB
  int<lower=0> A; // number of individuals ASE info
  int<lower=0> L; // length of vectors with n counts, gase p(H) and ai0 
  int<lower=0> K; // number of covariates
  int<lower=0> k; // number of Gaussians for eQTL effect prior
  int Y[N]; // total gene counts
  int sNB[N]; //  number of possible genotypes NB for each individual
  vector[G] gNB; // each geno NB
  vector[G] pNB; // prob for each geno NB
  int gase[L]; // genotype rsnp ASE individuals
  int m[A]; // total ase counts
  int n[L]; //n counts for ASE ind
  vector[L] pH; //p(H) for ASE ind
  vector[L] ai0; // allelic imbalance estimate for each haplotype for each sample
  vector[L] sdai0; // standard deviation for allelic imbalance estimate for each haplotype for each sample, log scale
  //int s[A]; // number of haplotypes per individual
  matrix[N,1+K] cov;
  int ASEi[N,2]; // index to link NB with ASE, first col is 1 if the individual has NB and ASE info, 0 otherwise. Second col gives index of ASE individual to relate NB with ASE
  int h2g[G]; // number of haps per genotype for ASE inds, 0 when ASE is not available
  vector[k] aveP; // mean for prior Gaussians for eQTL effect prior
  vector[k] sdP; // sd for prior Gaussians for eQTL effect prior
  vector[k] mixP; // log of mixing proportions for eQTL effect prior
  int I[N]; //indicator for skin all samples: 1=pso, 0=normal
  //int IA[A]; // indicator for skin in ASE samples: 1=pso, 0=normal
  
}

transformed data {
  int Max; // maximun number of elements in h2g
  vector[G] abs_gNB = fabs(gNB);
  vector[G] log_pNB = log(pNB);
  vector[L] log_pH = log(pH);
  Max = max(h2g);
}


parameters {
  real anorm; // mean expression normal tissue
  real apso; // mean expression pso tissues
  real ba; // log average-fold change ASE
  real bd; // log difference-fold change ASE
  real<lower=1e-5,upper=1e5> phi; //overdipersion param for neg binom
  real<lower=1e-5,upper=1e5> theta; //the overdispersion parameter for beta binomial
  vector[K-1] betas; // regression parameters 
  vector[L] rai0; // random intercept AI
}

transformed parameters {
  real bp; // parameter of interest for psoriasis
  real bn; // parameter of interest for normal skin
    
  bp = ba + bd;
  bn = ba -bd;
}

model {
  int pos; // to advance through NB terms (1-G)
  int posl; // to advance through ASE terms (1-L)
  vector[N] lmu1; // help to construct linear pred
  vector[G] lmu; // linear predictor log scale
  vector[G] ltmp; //  log NB likelihood

  real p; // ase proportion
  vector[Max] ase; //beta-binom terms
  real sAse; // sums beta-binom terms for haplotypes compatible with Gi=g
  real esum; // reduce computation inverse logit (rai0 + bp/bn)
  vector[L] esum0; // allelic imbalance proportion under the null
  vector[k] lpsa; // help for mixed gaussians for ba
  vector[k] lpsd; // help for mixed gaussians for bd
  
  /* vector[L] p; // ase proportion */
  /* vector[L] ase; //beta-binom term */
  vector[K] betasN; //regression parameters for normal skin
  vector[K] betasP; //regression parameters for pso skin
  vector[N] lmuN; //help linear pred normal skin
  vector[N] lmuP; //help linear pred pso skin
  /* real esum; // reduce computation inverse logit (la0 + bj) */
  
  //priors
  theta ~ gamma(1, 0.1); //  based on stan code example
  phi ~ gamma(1, 0.1);
  anorm ~ normal(6, 4); // mean expression normal skin, stan normal is mean and sd
  apso ~ normal(6, 4); // mean expression pso skin
  for(i in 1:(K-1)) {
    betas[i] ~ cauchy(0, 2.5);//prior for the covariates slopes following Gelman 2008
  }

  // allelic imbalance priors
  for(i in 1:L) {
    rai0[i] ~ normal(ai0[i], sdai0[i]);
  }

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
  pos = 1; // to advance on NB terms
  posl = 1; // to advance on ASE terms
  ase = rep_vector(0,Max);  // initialize ase vector to 0s to collect ase termns for each hap pair compatible with Gi=g */

  betasN = append_row(anorm, betas); //betas for normal inds
  betasP = append_row(apso, betas); // betas for pso inds
  lmuN = cov[, 2:cols(cov)] * betasN; // will be used for normal inds (based on indicator)
  lmuP = cov[, 2:cols(cov)] * betasP; // for pso inds
  
  esum0 = inv_logit(rai0);
  
  for (i in 1:N) { // lmu for each individual
    // check skin first
    if (I[i] == 1) { // psoriasis

      for (r in pos:(pos+sNB[i]-1)){ // then genotype

        lmu[r] = lmuP[i]; // G = 0
        lmu[r] = abs_gNB[r] == 1 ? lmu[r] + log1p(exp(bp)) - log(2) : lmu[r];
        lmu[r] = gNB[r] == 2 ? lmu[r] + bp : lmu[r];
        ltmp[r] = neg_binomial_2_log_lpmf(
          Y[i] | exp(lmu[r]), phi
        ) + log_pNB[r];

        if (ASEi[i,1] == 1) {  // ASE info
          for (x in 1:h2g[r]) {  // look at the haps compatibles with Gi=g
            esum = inv_logit(rai0[posl] + bp);
            p = gase[posl]==1 ? esum: esum0[posl];
            p = gase[posl]==-1 ? 1-esum : p;  // haplotype swap
            
            ase[x] = beta_binomial_lpmf(
              n[posl] | m[ASEi[i, 2]], p * theta, (1 - p) * theta
            ) + log_pH[posl];
            posl += 1;
          }
          sAse = log_sum_exp(ase[1:h2g[r]]);
          target +=  log_sum_exp(ltmp[r], sAse);
        }
      }
      if(ASEi[i,1] == 0) { // NO ASE, only NB terms for this ind
        target += log_sum_exp(ltmp[pos:(pos+sNB[i]-1)]);
      }

      pos += sNB[i];
      
    } else {
      for (r in pos:(pos+sNB[i]-1)) { //normal skin
        
        lmu[r] = lmuN[i]; // G = 0
        lmu[r] = fabs(gNB[r])==1 ? lmu[r] + log1p(exp(bn)) - log(2) : lmu[r];
        lmu[r] = gNB[r]==2 ? lmu[r] + bn : lmu[r];
        ltmp[r] = neg_binomial_2_log_lpmf(
          Y[i] | lmu[r], phi
        ) + log_pNB[r];

        if (ASEi[i,1] == 1) {  // ASE info
    
          for (x in 1:h2g[r]) {  // look at the haps compatibles with Gi=g
            
            esum = inv_logit(rai0[posl] + bn);
            p = gase[posl]==1 ? esum : esum0[posl];
            p = gase[posl]==-1 ? 1-esum : p;  // haplotype swap
            ase[x] = beta_binomial_lpmf(
              n[posl] | m[ASEi[i, 2]], p * theta, (1 - p) * theta
            ) + log_pH[posl];
            
            posl += 1;
          }
          sAse = log_sum_exp(ase[1:h2g[r]]);	 
          target += log_sum_exp(ltmp[r], sAse);
        }
      }
      
      if (ASEi[i, 1] == 0) { // NO ASE, only NB terms for this ind
        target += log_sum_exp(ltmp[pos:(pos+sNB[i]-1)]);
      }

      pos += sNB[i];
    }   
  }
}
