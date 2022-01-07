
// negative binomial and  ASE eQTL with  genotypes  allowing haplotype error, allowing for interaction term between 2 treatments in paired design, with or without covariates and ref bias correction version 2. Updated likelihood and mixed prior.
			      
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
  int m[A,2]; // total ase counts for each treatment, entry of 0 ifan individual doesnt have ASE ounts
  int n[L,2]; //n counts for ASE ind
  vector[L] pH; //p(H) for ASE ind
  real ai0[L,2]; // allelic imbalance estimate for each haplotype for each sample
  real sdai0[L,2]; // standard deviation for allelic imbalance estimate for each haplotype for each sample, log scale
  //int s[A]; // number of haplotypes per individual
  matrix[N,1+K] cov;
  // int ASEi[N,4]; // index to link NB with ASE, first col is 1 if the individual has NB and ASE info, 0 otherwise. Second col gives index of ASE individual to relate NB with ASE. Same order for second treatment.
  //int h2g[G]; // number of haps per genotype for ASE inds, 0 when ASE is not available
  vector[k] aveP; // mean for prior Gaussians for eQTL effect prior
  vector[k] sdP; // sd for prior Gaussians for eQTL effect prior
  vector[k] mixP; // log of mixing proportions for eQTL effect prior
  //int I[N]; //indicator for skin all samples: 1=pso, 0=normal
  //int IA[A]; // indicator for skin in ASE samples: 1=pso, 0=normal
  
}

transformed data {
  vector [L] logpH = log(pH);
}

/* transformed data { */
/*   int Max; // maximun number of elements in h2g */
/*   Max = max(h2g); */
/* } */


parameters {
  real anorm; // mean expression normal tissue
  real apso; // mean expression pso tissues
  real ba; // log average-fold change ASE
  real bd; // log difference-fold change ASE
  real<lower=1e-5,upper=1e5> phi; //overdipersion param for neg binom
  real<lower=1e-5,upper=1e5> theta; //the overdispersion parameter for beta binomial
  vector[K-1] betas; // regression parameters 
  matrix[L,2] rai0; // random intercept AI, easy way to access array columns 
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
  vector[L] esum; // reduce computation inverse logit (rai0 + bp/bn)
  vector[L] esum0; // allelic imbalance proportion under the null
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
  for(i in 1:(K-1)){
    betas[i] ~ cauchy(0, 2.5);//prior for the covariates slopes following Gelman 2008
  }

  // allelic imbalance priors per treatment
  for(i in 1:L) {
    for(t in 1:2){
      rai0[i, t] ~ normal(ai0[i, t], sdai0[i, t]);
    }
  }
  ui ~ normal(0, sdnb);
  uasei  ~ normal(0, sdase);

  sdnb ~ cauchy(0,1);

  sdase ~ cauchy(0,1);

  // mixture of gaussians for ba and bd:
  for(i in 1:k){
    lpsa[i] = normal_lpdf(ba | aveP[i], sdP[i]) + mixP[i];
    lpsd[i] = normal_lpdf(bd | aveP[i], sdP[i]) + mixP[i];
  }
  target += log_sum_exp(lpsa);
  target += log_sum_exp(lpsd);

  // transformed parameters of no interest
  pos = 1; // to advance on ASE terms
  /* ase = rep_vector(0,Max);  // initialize ase vector to 0s to collect ase termns for each hap pair compatible with Gi=g *\/ */

  //esum0 = inv_logit(rai0);
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
        esum0 = rai0[, t] + uasei[i];
        if (t == 1) { // first treament	
          esum = inv_logit(esum0 + bp);
        } else {
          esum = inv_logit(esum0 + bn);
        }	
        for (r in pos:(pos+s[i]-1)) {

          p = gase[i]==1 ? esum[r] : inv_logit(esum0[r]);
          p = gase[i]==-1 ? 1-esum[r] : p; //hap swap
              
          ltmp[r] = beta_binomial_lpmf(n[r,t] | m[i,t], p * theta, (1-p) * theta) + logpH[r];
        }
        target += log_sum_exp(ltmp[pos:(pos+s[i]-1)]);
      }
    }
    pos = pos + s[i];
  }
}	     
