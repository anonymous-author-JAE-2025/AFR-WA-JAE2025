model_t <- nimbleCode({
  # WORK ONLY IF ALL INDIVIDUALS START at AGE 1
  # -------------------------------------------------
  # Transition Parameters:
  # phi: survival probability
  # psi: breeding probability 
  # rho: breeding success probability 
  #
  # Detection Parameters:
  # pPB et pB: recapture probability
  # alp: known success probability
  # -------------------------------------------------
  # States:
  # 1 alive PB
  # 2 alive SB                                                                                                                                                    
  # 3 alive FB
  # 4 alive PSB
  # 5 alive PFB
  # 6 alive NB
  # 7 dead
  
  # Observations (y):  
  # 1 not detected
  # 2 detected as a NB
  # 3 detected as a SB
  # 4 detected as a FB
  # 5 detected as a B
  # -------------------------------------------------
  
  # for loops, are used:
  # ac for loops over age class
  # a = age
  # t = year
  # c = cohort
  # s = state
  # i = individual
  
  ##-------------------------------------------------
  ## 1. Define the priors for the parameters
  ##-------------------------------------------------
  
  #Mean Recapture probabilities
  for (ac in 1:Nclass_pPB){pPB[ac] ~ dunif(0,1)}
  for (s in 1:Nclass_pB){pB[s] ~ dunif(0,1)}
  for (s in 1:2){alp[s] ~ dunif(0,1)}
  
  #Mean Survival and reproductive transitions parameters
  for (s in 1:Nclass_phiB){
    Mu.phiB[s] ~ dnorm(0,sd=2)
  }
  for (s in 1:Nclass_psiB){
    Mu.psiB[s] ~ dnorm(0,sd=2)
  }
  for (s in 1:Nclass_rhoB){
    Mu.rhoB[s] ~ dnorm(0,sd=2)
  }
  # 
  #for(ac in 1:(Nclass_psiPB-1)){#age from 4 to 15
  for(ac in 1:(Nclass_psiPB)){#age from 4 to 15
    Mu.psiPB[ac] ~ dnorm(0,sd=2)
    Mu.rhoPB[ac] ~ dnorm(0,sd=2)
  }
  for(ac in 1:Nclass_phiPB){
    Mu.phiPB[ac]  ~ dnorm(0.9,0.6)
  }
  
  b.phi.PB ~ dnorm(0,sd=2)
  b.psi.PB ~ dnorm(0,sd=2)
  b.rho.PB ~ dnorm(0,sd=2)
  
  sigma.eps.phi ~ dunif(0, 10)
  sigma.eps.psi ~ dunif(0, 10)
  sigma.eps.rho ~ dunif(0, 10)
  
  for(t in 1:n_year) {
    eps.phi[t] ~ dnorm(0, sd = sigma.eps.phi)
    eps.psi[t] ~ dnorm(0, sd = sigma.eps.psi)
    eps.rho[t] ~ dnorm(0, sd = sigma.eps.rho)
  }
  
  ##-------------------------------------------------
  ## 2. Derived parameters
  ##-------------------------------------------------
  # Derived survival and transition probabilities
  for (t in 1:n_year){
    for(ac in 1:Nclass_phiPB){
      logit(phiPB[ac,t]) <- Mu.phiPB[ac] + b.phi.PB * COV.phi[t] + eps.phi[t]
    }
    for(ac in 1:Nclass_psiPB) {
      logit(psiPB[ac,t]) <- Mu.psiPB[ac] + b.psi.PB * COV.psi[t] + eps.psi[t]
    }
    for(ac in 1:Nclass_psiPB) {
      logit(rhoPB[ac,t]) <- Mu.rhoPB[ac] + b.rho.PB * COV.rho[t] + eps.rho[t]
    }
  } 
  
  for (s in 1:Nclass_psiB){
    logit(psiB[s]) <- Mu.psiB[s]
  }
  
  for (s in 1:Nclass_phiB){
    logit(phiB[s]) <- Mu.phiB[s]
  }
  for (s in 1:Nclass_rhoB){
    logit(rhoB[s]) <- Mu.rhoB[s]
  }
  
  
  # Detection matrix
  for (ac in 1:Nclass_pPB){
    omega[1, 1, ac] <- 1-pPB[ac]
    omega[1, 2, ac] <- pPB[ac]
    for (s in c(2,3)) {
      omega[s, 5, ac] <- (1-alp[s-1])*pB[s-1]
      omega[s, 1, ac] <- 1-pB[s-1]
    }
    omega[2, 3, ac] <- alp[1]*pB[1]
    omega[3, 4, ac] <- alp[2]*pB[2]
    for (s in 4:6) {
      omega[s, 1, ac] <- 1-pB[s-1]
      omega[s, 2, ac] <- pB[s-1]
    }
  }
  
  # Transition matrix
  for(c in 1:n_cohorts){
    for(t in (c):(n_year)){
      gamma[1,1,c,t] <- phiPB[AgeclassphiPB[(t-c+1)],t] * (1 - psiPB[AgeclasspsiPB[(t-c+1)],t])
      gamma[1,2,c,t] <- phiPB[AgeclassphiPB[(t-c+1)],t] * psiPB[AgeclasspsiPB[(t-c+1)],t] * rhoPB[AgeclasspsiPB[(t-c+1)],t]
      gamma[1,3,c,t] <- phiPB[AgeclassphiPB[(t-c+1)],t] * psiPB[AgeclasspsiPB[(t-c+1)],t] * (1-rhoPB[AgeclasspsiPB[(t-c+1)],t])
      gamma[1,7,c,t] <- 1-phiPB[AgeclassphiPB[(t-c+1)],t]
      
      gamma[2,2,c,t] <- phiB[1] * psiB[1] * rhoB[1]
      gamma[2,3,c,t] <- phiB[1] *  psiB[1] * (1-rhoB[1])
      gamma[2,4,c,t] <- phiB[1]   * (1- psiB[1] )
      gamma[2,7,c,t] <-  1-phiB[1]
      gamma[3,2,c,t] <- phiB[2] * psiB[2] * rhoB[1]
      gamma[3,3,c,t] <- phiB[2] *  psiB[2] * (1-rhoB[1])
      gamma[3,5,c,t] <- phiB[2]  * (1- psiB[2])
      gamma[3,7,c,t] <-  1-phiB[2]
      
      for (s in 4:6){
        gamma[s,2,c,t] <- phiB[s-3] * psiB[s-1] * rhoB[1]
        gamma[s,3,c,t] <- phiB[s-3] *  psiB[s-1] * (1-rhoB[1])
        gamma[s,6,c,t] <- phiB[s-3] * (1- psiB[s-1])
        gamma[s,7,c,t] <- 1-phiB[s-3]
      }
      
    }
  }
  
  #-------------------------------------------------
  # 3. The likelihoods
  #-------------------------------------------------
  for(i in 1:n_ind){
    y[i,first[i]:n_year] ~ dDHMMmy(init = init[1:7], 
                                   probObs = omega[1:7,1:5,1:Nclass_pPB],
                                   probTrans = gamma[1:7, 1:7,first[i],first[i]:(n_year)],
                                   NUMi = NUM[i],
                                   agei = agePB[first[i],first[i]:(n_year)],
                                   len = n_year-first[i]+1, checkRowSums = 1)
  }
  
})
