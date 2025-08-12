# Data loading and reshape -----------------------------------------------------
# Here we keep all individuals to avoid selecting only individuals that have 
# reach a stage where they could be sexed. This would make phiPB close to 1!!!
CMR_df <- read_csv(glue("{DataDir}/data_cmr/GA_CRO_AnalyseCroisee_updated_to_2020.csv")) %>%
  mutate(across(where(is.character), ~na_if(., "NA"))) %>%
  filter(StatutBaguage == "P")

### Assign Sex and select only males
CMR_df$Sex <- NA

for (i in 1:nrow(CMR_df)) {
  # Assign the good sex given genetic 
  if (!is.na(CMR_df$SexeGenetique[i]) & CMR_df$SexeGenetique[i] == "F") {
    CMR_df$Sex[i] <- "F"
  }
  if (!is.na(CMR_df$SexeGenetique[i]) & CMR_df$SexeGenetique[i] == "M") {
    CMR_df$Sex[i] <- "M"
  }
  # Now if no genetic has been done but visual assignement has been made, assign the state
  if (is.na(CMR_df$Sex[i]) & !is.na(CMR_df$SexeVisuel[i]) & CMR_df$SexeVisuel[i] == "F") {
    CMR_df$Sex[i] <- "F"
  }
  if (is.na(CMR_df$Sex[i]) & !is.na(CMR_df$SexeVisuel[i]) & CMR_df$SexeVisuel[i] == "M") {
    CMR_df$Sex[i] <- "M"
  }
  # Now assign randomly a sex if there is no sex available. Since it will allow to better estimate survival. 
  if (is.na(CMR_df$Sex[i]) & is.na(CMR_df$SexeVisuel[i]) & is.na(CMR_df$SexeVisuel[i])) {
    s <- c("F", "M")
    CMR_df$Sex[i] <- sample(s, 1, prob = c(0.5, 0.5))
  }
}

y_tmp <- CMR_df %>% 
  filter(CycleBaguage >= 1970 & CycleBaguage <= 2020) %>% 
  filter(Sex == sex) %>% 
  select("CycleBaguage", "1970":"2020") 
y_tmp <- y_tmp[,1:52] %>% as.matrix

head(y_tmp)

# breeding success codes	
# C	control/ present
# R	breeding= lay egg
# RP	successful breeder= raise successfully chick
# RM	failed breeder stage unknown: egg did not hatch or if the chick died
# RMO	failed breeder stage egg
# RMP	failed breeder stage chick rearing
# NR	nonbreeder
# Pr?Repro	pre-breeding season
# R/1O	breeding stage egg
# R/1P	breeding stage chick-rearing
# 1O	egg

### Pre-Breeder
y_tmp <- replace(y_tmp, y_tmp == "POU", 2)

### Breeder
y_tmp <- replace(y_tmp, y_tmp == "R/10", 5)
y_tmp <- replace(y_tmp, y_tmp == "R/1O", 5)
y_tmp <- replace(y_tmp, y_tmp == "R", 5)
y_tmp <- replace(y_tmp, y_tmp == "10", 5)
y_tmp <- replace(y_tmp, y_tmp == "R/1P", 5)
y_tmp <- replace(y_tmp, y_tmp == "1P", 5)

### Successful Breeder
y_tmp <- replace(y_tmp, y_tmp == "RP", 3)

### Failed Breeder
y_tmp <- replace(y_tmp, y_tmp == "RM", 4)
y_tmp <- replace(y_tmp, y_tmp == "RMP", 4)
y_tmp <- replace(y_tmp, y_tmp == "RMO", 4)

### Non Breeder
y_tmp <- replace(y_tmp, y_tmp == "NR", 2)
y_tmp <- replace(y_tmp, y_tmp == "PréRepro", 2) 
y_tmp <- replace(y_tmp, y_tmp == "C", 2) 
y_tmp <- replace(y_tmp, y_tmp == "Controleur", 2) 
# y <- matrix(as.numeric(y_tmp[8000:8100,]),    # Convert to numeric matrix
#             ncol = ncol(y_tmp), dimnames = dimnames(y_tmp))

y <- matrix(as.numeric(y_tmp),    # Convert to numeric matrix
            ncol = ncol(y_tmp), dimnames = dimnames(y_tmp))

#Keep unique life histories and number of individuals per history
y <- y %>% as_tibble()%>%
  group_by_all() %>% count

NUM = y$n #number of individuals per unique life history

y <- y[,2:52] %>% as.matrix #remove NUM from y

y[is.na(y)] <- 0

table(y)
head(y)

# Filter the individuals who are only 0's
NUM <- NUM[rowSums(y) > 0]
y <- y[rowSums(y) > 0,]

# Filter the individuals who reproduce before 6. 
get.first <- function(x) min(which(x != 0))
first <- apply(y, 1, get.first)

# Obtain the age of each monitored individual if they are alive each year. 
# get.first.B <- function(x) min(which(x == 4 | x == 5))
# first.B <- apply(y, 1, get.first.B)
# 
# diff <- first.B - first
# remove <- which(diff < 6)
# 
# y <- y[-remove,]

# Model for true data ------------------------------------------------
# -------------------------------------------------
# States (z):
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

n_cohorts = max(first)
n_year = ncol(y)
n_states = 7
n_event = 5

## MODEL Constraint ##
#Age classes are defined. They can be modified here
#WARNING: if state classes are included or modified, all transitions matrices need to be modified accordingly

## --- Survival --- ##
# phiPB: Survival probability of birds that are in a Pre-Breeding State (ex: Chicks and Juveniles)
#phiPB has: 4 age classes: 1:2, 3:8, 9:12, 13+
AgeclassphiPB = c(1,1,2,2,2,2,2,2,3,3,3,3,rep(4,n_year-12))
Nclass_phiPB = length(unique(AgeclassphiPB))

# phiB: Survival probability of breeders birds
#varies with 3 state classes: SB & PSB / FB & PFB /NB
#Warning: if state classes are modified, all transitions matrices need to be modified accordingly
Nclass_phiB = 3

## --- Probability to reproduce  --- ##
# psiPB: Breeding probability of birds that are in a Pre-Breeding State 
#psiPB has age effect in factor from 6 to >=10 
AgeclasspsiPB = c(1,1,1,1,1,2:6, rep(7,n_year-8))
Nclass_psiPB = length(unique(AgeclasspsiPB))

# psiB: Breeding probability of birds that have already reproduced
# varies with 5 states state : SB / PSB / FB / PFB /NB
#Warning: if state classes are modified, all transitions matrices need to be modified accordingly
Nclass_psiB = 5

## --- Breeding success --- ##
# rhoPB: Breeding success probability of birds that are in a Pre-Breeding State 
#same age effect as psi PB

# rhoB: Breeding success probability of birds  that have already reproduced
#constant
#Warning: if state classes are modified, all transitions matrices need to be modified accordingly
Nclass_rhoB = 1

# ---
# Detection probabilities
# ---
# pPB Detection probability for birds that are in a Pre-Breeding State
#pPB class 1-2, 3:5, 6_10, 11+
AgeclasspPB = c(1,1,2,2,2,3,4,5,6,7, rep(8,n_year-10))
Nclass_pPB = length(unique(AgeclasspPB))

# pB: Detection probability for mature birds
# varies with 5 states state : SB / PSB / FB / PFB /NB
#Warning: if state classes are modified, all detection matrices need to be modified accordingly
Nclass_pB = 5

#Create age variables
#real_age is true age per cohort
#agePB is age class for immature probability of detection over time per cohort
agePB = real_age = array(0, c(n_cohorts, n_year))
for (c in 1:n_cohorts){
  real_age[c,c:n_year] = 0:(n_year-c)
  agePB[c,(c+1):n_year] = AgeclasspPB[real_age[c,(c+1):n_year]]}

# COHORT effect
## With a limit for psi, since I want to estimate the effect of the cohort 2005 at maximum on psi.
cpsi <- c(seq(1970,2004,1), rep(2005, ncol(y_tmp)-16)) - 1970 + 1
CPSI <- 36
cohort.psi <- seq(1,n_cohorts,1)
cohort.psi[cohort.psi > CPSI ] <- CPSI
# cohort.psi <- rep(0, n_cohorts)
cohort.phi <- seq(1,n_cohorts,1)
# cohort.phi <- rep(0, n_cohorts)

# Pass data and values to nimble ------------------------------------------------
my.constants <- list(first = first, 
                     NUM = NUM,
                     init = c(1,rep(0,n_states-1)),
                     n_year = n_year, 
                     cpsi = cpsi,
                     n_ind = nrow(y),
                     agePB = agePB,
                     Nclass_pPB = Nclass_pPB,
                     Nclass_phiPB = Nclass_phiPB, AgeclassphiPB  = AgeclassphiPB,
                     Nclass_psiPB = Nclass_psiPB, AgeclasspsiPB  = AgeclasspsiPB, 
                     Nclass_rhoB =  Nclass_rhoB, Nclass_pB =  Nclass_pB,
                     Nclass_phiB =  Nclass_phiB, Nclass_psiB =  Nclass_psiB # ,
                     # trend.phi = trend.phi,
                     # trend.psi = trend.psi
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
                              b1psiPB = 0, 
                              b1phiPB = 0,
                              b2psiPB = 0, 
                              b2phiPB = 0,
                              pPB = runif(Nclass_pPB, 0, 1), 
                              pB= runif(Nclass_pB, 0, 1),
                              alp = runif(2, 0,1),
                              omega = omegaper, gamma = gammacoh, # These ones need to be initialized here as squares that do not change are filled in the model
                              psiPB = array(0,c(Nclass_psiPB, n_year)), # These ones need to be initialized here as squares that do not change are filled in the model
                              rhoPB = rep(0,Nclass_psiPB) # These ones need to be initialized here as squares that do not change are filled in the model
)

parameters.to.save <- c("Mu.phiPB", "Mu.psiPB", "Mu.rhoPB",
                        "Mu.phiB", "Mu.psiB", "Mu.rhoB",
                        "pPB","pB","alp", 
                        "b1psiPB", "b1phiPB",
                        "b2psiPB", "b2phiPB")

# Bundle Data for Nimble
my.data <- list(y = y,
                cohort.psi = cohort.psi,
                cohort.phi = cohort.phi
                )

# Consider increasing number of iterations, burnins and thins
n.iter <- 100000
n.burnin <- 50000
n.chains <- 3
nthin <- 100

# Run the nimble model --------------------------------------------------
print(Sys.time())
Rmodel <- nimbleModel(code = model, 
                      constants = my.constants,
                      data = my.data,              
                      inits = initial.values_trend, check=T)
print(Sys.time())

## configure MCMC
conf <- configureMCMC(Rmodel, monitors = parameters.to.save, thin = nthin, enableWAIC = T)

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
