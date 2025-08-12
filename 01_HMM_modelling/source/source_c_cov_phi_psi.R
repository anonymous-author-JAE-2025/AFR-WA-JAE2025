

# Pass data and values to nimble ------------------------------------------------
my.constants <- list(first = first, 
                     NUM = NUM,
                     init = c(1,rep(0,n_states-1)),
                     n_year = n_year, 
                     cpsi = cpsi,
                     n_cohorts = n_cohorts,
                     n_ind = nrow(y),
                     agePB = agePB,
                     Nclass_pPB = Nclass_pPB,
                     Nclass_phiPB = Nclass_phiPB, AgeclassphiPB  = AgeclassphiPB,
                     Nclass_psiPB = Nclass_psiPB, AgeclasspsiPB  = AgeclasspsiPB, 
                     Nclass_rhoB =  Nclass_rhoB, Nclass_pB =  Nclass_pB,
                     Nclass_phiB =  Nclass_phiB, Nclass_psiB =  Nclass_psiB
)

# --- Initial states --- #
y[y==0]=1

omegaper = array(0,c(n_states,n_event,Nclass_pPB))
omegaper[,1,] =1
gammacoh = array(0,dim=c(n_states,n_states,n_cohorts,(ncol(y))))
gammacoh [n_states,n_states,,] = 1
initial.values_trend <- list( Mu.phiPB = rnorm(Nclass_phiPB, 0.9, 0.6),
                              Mu.psiPB = rnorm(Nclass_psiPB, 0, 2), # Why in the model Floriane puts a -1 ?
                              Mu.rhoPB = rnorm(Nclass_psiPB, 0, 2), # Why in the model Floriane puts a -1 ?
                              Mu.phiB = rnorm(Nclass_phiB, 0, 2),
                              Mu.psiB = rnorm(Nclass_psiB, 0, 2), 
                              Mu.rhoB = rnorm(Nclass_rhoB, 0, 2),
                              b1cpsiPB = 0, b1cphiPB = 0,
                              b2cpsiPB = 0, b2cphiPB = 0,
                              b3cpsiPB = 0, b3cphiPB = 0,
                              pPB = runif(Nclass_pPB, 0,1), 
                              pB= runif(Nclass_pB, 0,1),
                              alp = runif(2, 0,1),
                              omega = omegaper, gamma = gammacoh, # These ones need to be initialized here as squares that do not change are filled in the model
                              psiPB = array(0,c(Nclass_psiPB, n_year)), # These ones need to be initialized here as squares that do not change are filled in the model
                              rhoPB = rep(0,Nclass_psiPB) # These ones need to be initialized here as squares that do not change are filled in the model
)

parameters.to.save <- c("Mu.phiPB", "Mu.psiPB", "Mu.rhoPB",
                        "Mu.phiB", "Mu.psiB", "Mu.rhoB",
                        "pPB","pB","alp", 
                        "b1cpsiPB", "b1cphiPB")

# Bundle Data for Nimble
my.data <- list(y = y,
                COVpsi1 = COVpsi1,
                COVpsi2 = COVpsi2,
                COVpsi3 = COVpsi3,
                COVphi1 = COVphi1,
                COVphi2 = COVphi2,
                COVphi3 = COVphi3)

#Consider increasing number of iterations, burnins and thins
n.iter <- 5000
n.burnin <- 2500
n.chains <- 3
nthin = 1

# Run the nimble model --------------------------------------------------
print(Sys.time())
Rmodel <- nimbleModel(code = model, 
                      constants = my.constants,
                      data = my.data,              
                      inits = initial.values_trend, check=T)
print(Sys.time())

## configure MCMC
conf<- configureMCMC(Rmodel,  monitors = parameters.to.save,thin=nthin, enableWAIC =T)

if (!is.null(param_rm)) {
  conf$removeSamplers(param_rm)
}

if (is.null(param_rm)) {
  #conf$removeSamplers(param_rm)
}

## build MCMC
Rmcmc <- buildMCMC(conf)
## compile model and MCMC
Cmodel <- compileNimble(Rmodel,showCompilerOutput = T)
Cmcmc <- compileNimble(Rmcmc, project=Rmodel)
print(Sys.time())
samplesList <- runMCMC(Cmcmc, niter=n.iter, nburnin = n.burnin, nchains=n.chains, progressBar=T, summary = F, WAIC=T)
print(Sys.time())

# Make results summary & save results ---------------------------------------
m = (n.iter-n.burnin)/nthin
rb = array(NA,dim=c(m,n.chains,ncol(samplesList$samples$chain1)))
rb[,1,] = samplesList$samples$chain1
rb[,2,] = samplesList$samples$chain2
rb[,3,] = samplesList$samples$chain3
dimnames(rb)[[3]] <- colnames(samplesList$samples$chain1)
#output_dd = sum_nim(rb, na.rm = "TRUE")
output_dd = sum_nim(rb)
ilogit(output_dd)

waic <- samplesList$WAIC

save(output_dd, rb, waic, file = paste("outputs/", m_, "_",  sex, ".Rdata", sep = ""))