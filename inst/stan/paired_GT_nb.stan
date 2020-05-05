
// negative binomial and  ASE eQTL with unkwown  genotypes  allowing haplotype error, allowing for interaction term between 2 treatments in paired design, with or without covariates and ref bias correction version 2. Updated likelihood and mixed prior.
			      
data {
  int<lower=0> N; // number  individuals with NB info
  int<lower=0> K; // number of covariates
  int<lower=0> k; // number of Gaussians for eQTL effect prior
  int Y[N,2]; // total gene counts for each treatment
  int g[N]; // rnsp geno for all individuals
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
  real<lower=0> phi; //overdipersion param for neg binom
  vector[K-1] betas; // regression parameters 
  vector[N] ui; //random term NB
  real<lower=0> sdnb; //sd for random term in nb
  
  }

transformed parameters {
  real bp; // parameter of interest for psoriasis
  real bn; // parameter of interest for normal skin
    
  bp = ba + bd;
  bn = ba -bd;
    
}

model {
  real lmu; // help with linear predictor log scale
  vector[k] lpsa; // help for mixed gaussians for ba
  vector[k] lpsd; // help for mixed gaussians for bd
  vector[K] betasN; //regression parameters for normal skin
  vector[K] betasP; //regression parameters for pso skin
  vector[N] lmuN; //help linear pred normal skin
  vector[N] lmuP; //help linear pred pso skin
 
  //priors
  phi ~ gamma(1,0.1); //  based on stan code example
  
  anorm ~ normal(6,4); // mean expression normal skin, stan normal is mean and sd
  apso ~ normal(6,4); // mean expression pso skin
  for(i in 1:(K-1)){
    betas[i] ~ cauchy(0,2.5);//prior for the covariates slopes following Gelman 2008
  }
 
  // random effect priors
  ui ~ normal(0, sdnb);
  sdnb ~ cauchy(0,1); 

  // mixture of gaussians for ba and bd:
  for(i in 1:k){
    lpsa[i] = normal_lpdf(ba | aveP[i], sdP[i]) + mixP[i];
    lpsd[i] = normal_lpdf(bd | aveP[i], sdP[i]) + mixP[i];
  }
  target += log_sum_exp(lpsa);
  target += log_sum_exp(lpsd);

  betasN=append_row(anorm, betas); //betas for normal inds
  betasP=append_row(apso, betas); // betas for pso inds
  lmuN=cov[,2:cols(cov)]*betasN; // will be used for normal inds (based on indicator)
  lmuP=cov[,2:cols(cov)]*betasP; // for pso inds
  
  for(i in 1:N){ //  for each individual

    // go by treatment
    for(t in 1:2){
   
      if(t == 1){ //first treatment: bp
      
      lmu = lmuP[i] + ui[i]; // G = 0
      lmu = fabs(g[i])==1 ? lmu  + log1p(1+exp(bp))-log(2) : lmu;
      lmu = g[i]==2 ? lmu + bp : lmu;
    

    } else { // second treatment
      
      lmu = lmuN[i] + ui[i]; // G = 0
      lmu = fabs(g[i])==1 ? lmu  + log1p(1+exp(bn))-log(2) : lmu;
      lmu = g[i]==2 ? lmu + bn : lmu;
   
    }
     target += neg_binomial_2_lpmf(Y[i] | exp(lmu),phi);
    }
  }
}
    

   
	     
		   
    
