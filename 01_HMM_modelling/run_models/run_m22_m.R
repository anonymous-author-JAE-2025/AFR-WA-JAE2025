##--------------------------------------------------------------------------------------------------------
## SCRIPT : Run wandering Albatross demographic model on real data
##
## Authors : Floriane Plard
## Last update : 2023-15-19
## R version 4.0.4 (2021-02-15) -- "Lost Library Book"
## Copyright (C) 2021 The R Foundation for Statistical Computing
## Platform: x86_64-w64-mingw32/x64 (64-bit)
##--------------------------------------------------------------------------------------------------------

# Remove all objects from the environment
rm(list = ls())

# Load the relevant packages
lapply(c("tidyverse", "dplyr", "tidyr", "nimble", "glue", "reshape2", "readxl"), library, character.only = TRUE)

# Define directory paths for data, functions, output, and model
WorkDir <- getwd() # Working directory
DataDir <- paste(WorkDir, "data", sep = "/") # Data folder
OutDir <- paste(WorkDir, "output", sep = "/") # Output folder
ModDir <- paste(WorkDir, "model", sep = "/") # Models folder
FunDir <- paste(WorkDir, "function", sep = "/") # Models folder
SourceDir <- paste(WorkDir, "source", sep = "/") # Models folder
source(glue("{FunDir}/functions.R"))

source(glue("{ModDir}/model_trends.R"))
model <- model_t
m_ <- "model_22"
sex <- "M"

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
y_tmp <- y_tmp[,1:52]%>%as.matrix

head(y_tmp)

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

# Covariate dataset reshape

### Divorce
divorce_f <- read_excel("data/data_cov/divorce_f.xlsx") %>% 
  mutate(year = round(as.numeric(x)),
         n = round(as.numeric(y)),
         sex = rep("Female")) %>% 
  select(year, n, sex)

divorce_m <- read_excel("data/data_cov/divorce_m.xlsx") %>% 
  mutate(year = round(as.numeric(x)),
         n = round(as.numeric(y)),
         sex = rep("Male")) %>% 
  select(year, n, sex)

divorce <- rbind(divorce_m, divorce_f)

### Widowhood
widowhood_f <- read_excel("data/data_cov/widowhood_f.xlsx") %>% 
  mutate(year = round(as.numeric(x)),
         n = round(as.numeric(y)),
         sex = rep("Female")) %>% 
  select(year, n, sex)

widowhood_m <- read_excel("data/data_cov/widowhood_m.xlsx") %>% 
  mutate(year = round(as.numeric(x)),
         n = round(as.numeric(y)),
         sex = rep("Male")) %>% 
  select(year, n, sex)

widowhood <- rbind(widowhood_m, widowhood_f)

### Available
available_f <- data.frame(year = divorce_f$year,
                          n = divorce_f$n + widowhood_f$n)

available_m <- data.frame(year = divorce_m$year,
                          n = divorce_m$n + widowhood_m$n)

available <- data.frame(year = divorce_m$year,
                        n = divorce_f$n + widowhood_f$n + divorce_m$n + widowhood_m$n)

## Load n_individuals
n_couples <- read_excel("data/data_cov/waa_crz_monitoring_231003.xlsx")
colnames(n_couples) <- c("year", "nbre_cples")

n_couples$year <- as.numeric(n_couples$year)
n_couples$nbre_cples <- as.numeric(n_couples$nbre_cples)

#fit linear regression model using data frame
mod <- lm(nbre_cples ~ year, data = n_couples)

#interpolate y value based on x value of 13
y_new <- approx(n_couples$year, n_couples$nbre_cples, xout= c(1960,1961,1962,1963,1964,1965,1966,1967,1970,1971,1972,1973,1974,1978,1979,1980))
y_new <- as.data.frame(y_new)
colnames(y_new) <- c("year", "nbre_cples")

n_couples_inter <- n_couples
n_couples_inter$interpol <- rep("no")

for (i in 1:nrow(n_couples_inter)) {
  for (j in 1:nrow(y_new)) {
    if(is.na(n_couples_inter$nbre_cples[i]) & n_couples_inter$year[i] == y_new$year[j]) {
      n_couples_inter$nbre_cples[i] <- y_new$nbre_cples[j]
      n_couples_inter$interpol[i] <- "yes"
    }
  }
}

n_couples_inter$nbre_cples <- round(n_couples_inter$nbre_cples)

### Obtain ratio of available individuals/total couples
available <- merge(available, n_couples_inter, by = "year") %>% 
  filter(year > 1969)

available$ratio <- available$n / available$nbre_cples
available$std_ratio <- (available$ratio - mean(available$ratio))/sd(available$ratio)

COV.phi <- rep(0, length(seq(1970,2020,1)))
COV.psi <- c(available$std_ratio, rep(available$std_ratio[43], (n_cohorts-41)))
COV.rho <- rep(0, length(seq(1970,2020,1)))

### Remove parameter from sampler
param_rm <- c("b.phi.PB", "b.rho.PB")

source(glue("{SourceDir}/source_t.R"))