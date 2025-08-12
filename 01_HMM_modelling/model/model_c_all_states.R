model_c_all_states <- nimbleCode({
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

  sigma.eps.phiPB ~ dunif(0, 10)
  sigma.eps.psiPB ~ dunif(0, 10)
  sigma.eps.rhoPB ~ dunif(0, 10)
  
  for(c in 1:n_cohorts) {
    eps.phiPB[c] ~ dnorm(0, sd = sigma.eps.phiPB)
    eps.psiPB[c] ~ dnorm(0, sd = sigma.eps.psiPB)
    eps.rhoPB[c] ~ dnorm(0, sd = sigma.eps.rhoPB)
  }
  
  for(s in 1:Nclass_phiB) {
    sigma.eps.phiB[s] ~ dunif(0, 10)
    for(c in 1:n_cohorts) {
      eps.phiB[s,c] ~ dnorm(0, sd = sigma.eps.phiB[s])
    }
  }
  
 for(s in 1:Nclass_psiB) {
  sigma.eps.psiB[s] ~ dunif(0, 10)
  for(c in 1:n_cohorts) {
    eps.psiB[s,c] ~ dnorm(0, sd = sigma.eps.psiB[s])
  }
 }
  
  ##-------------------------------------------------
  ## 2. Derived parameters
  ##-------------------------------------------------
  # Derived survival and transition probabilities
  for (c in 1:n_cohorts){
    for(ac in 1:Nclass_phiPB){
      logit(phiPB[ac,c]) <- Mu.phiPB[ac] + eps.phiPB[c]
    }
    for(ac in 1:Nclass_psiPB) {
      logit(psiPB[ac,c]) <- Mu.psiPB[ac] + eps.psiPB[c]
    }
    for(ac in 1:Nclass_psiPB) {
      logit(rhoPB[ac,c]) <- Mu.rhoPB[ac] + eps.rhoPB[c]
    }
   
  
    for (s in 1:Nclass_psiB){
      logit(psiB[s,c]) <- Mu.psiB[s] + eps.psiB[s,c]
    }
  
    for (s in 1:Nclass_phiB){
      logit(phiB[s,c]) <- Mu.phiB[s] + eps.phiB[s,c]
    }
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
      gamma[1,1,c,t] <- phiPB[AgeclassphiPB[(t-c+1)],c] * (1 - psiPB[AgeclasspsiPB[(t-c+1)],cpsi[c]])
      gamma[1,2,c,t] <- phiPB[AgeclassphiPB[(t-c+1)],c] * psiPB[AgeclasspsiPB[(t-c+1)],cpsi[c]] * rhoPB[AgeclasspsiPB[(t-c+1)],cpsi[c]]
      gamma[1,3,c,t] <- phiPB[AgeclassphiPB[(t-c+1)],c] * psiPB[AgeclasspsiPB[(t-c+1)],cpsi[c]] * (1-rhoPB[AgeclasspsiPB[(t-c+1)],cpsi[c]])
      gamma[1,7,c,t] <- 1-phiPB[AgeclassphiPB[(t-c+1)],c]
      
      gamma[2,2,c,t] <- phiB[1,c] * psiB[1,c] * rhoB[1]
      gamma[2,3,c,t] <- phiB[1,c] *  psiB[1,c] * (1-rhoB[1])
      gamma[2,4,c,t] <- phiB[1,c]   * (1- psiB[1,c] )
      gamma[2,7,c,t] <-  1-phiB[1,c]
      gamma[3,2,c,t] <- phiB[2,c] * psiB[2,c] * rhoB[1]
      gamma[3,3,c,t] <- phiB[2,c] *  psiB[2,c] * (1-rhoB[1])
      gamma[3,5,c,t] <- phiB[2,c]  * (1- psiB[2,c])
      gamma[3,7,c,t] <-  1-phiB[2,c]
      
      for (s in 4:6){
        gamma[s,2,c,t] <- phiB[s-3,c] * psiB[s-1,c] * rhoB[1]
        gamma[s,3,c,t] <- phiB[s-3,c] *  psiB[s-1,c] * (1-rhoB[1])
        gamma[s,6,c,t] <- phiB[s-3,c] * (1- psiB[s-1,c])
        gamma[s,7,c,t] <- 1-phiB[s-3,c]
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
