data{/////////////////////////////////////////////////////////////////////
    int<lower=1> S;    // number of samples (nrow)
    array[S] int<lower=1> Species;    // index of species, each of which will have a different value for p11 and p10
    int<lower=1> Nspecies;    // number of species, each of which will have a different value for p11 and p10
    array[S] int<lower=1> L;   // index of locations or species/site combinations, each of which will have a different value psi
    int<lower=1> Nloc;   // number of locations or species/site combinations, each of which will have a different value psi
    array[S] int<lower=1> K;   // number of replicates per site (ncol)
    array[S] int<lower=0> N; // number of detections among these replicates
    array[S] int z;   // integer flag to help estimate psi parameter
}
parameters{/////////////////////////////////////////////////////////////////////
array[Nloc] real<lower=0,upper=1> psi;  //commonness parameter
array[Nspecies] real<lower=0,upper=1> p11; //true positive detection rate
array[Nspecies] real<lower=0,upper=1> p10; //false positive detection rate
}
transformed parameters{/////////////////////////////////////////////////////////////////////
}
model{/////////////////////////////////////////////////////////////////////
array[S] real p;

for (i in 1:S){
z[i] ~ bernoulli(psi[L[i]]);
p[i] = z[i]*p11[Species[i]] + (1-z[i])*p10[Species[i]];
N[i] ~ binomial(K[i], p[i]);
}; 

//priors
psi ~ beta(2,2); 
p11 ~ beta(2,2); 
p10 ~ beta(1,10);
}
generated quantities{
array[S] real<lower=0,upper=1> Occupancy_prob;    //after inferring parameters above, now calculate occupancy probability for each observation. Equation from Lahoz-Monfort et al. 2015

for (i in 1:S){
Occupancy_prob[i]  = (psi[L[i]]*(p11[Species[i]]^N[i])*(1-p11[Species[i]])^(K[i]-N[i])) 
/ ((psi[L[i]]*(p11[Species[i]]^N[i])*(1-p11[Species[i]])^(K[i]-N[i])) 
+ (((1-psi[L[i]])*(p10[Species[i]]^N[i]))*((1-p10[Species[i]])^(K[i]-N[i])))
);
}
}


