# Clear environment and load required packages
rm(list = ls())
library(readr)
library(ggplot2)
library(R.matlab)

### Cohort effects -----------------------------------------------------------
# Function to process survival data with random effects for one sex
process_survival <- function(rb, sex_label) {
  # Get indices for fixed and random effects
  phiPB_indices <- which(dimnames(rb)[[3]] %in% paste0("Mu.phiPB[", 1:4, "]"))
  eps_phi_indices <- which(dimnames(rb)[[3]] %in% paste0("eps.phi[", 1:49, "]"))
  
  # Randomly sample 100 iterations
  set.seed(123)
  random_iterations <- sample(1:250, 100)
  
  # Create age classes vector and ages
  age_classes <- c(1,1,2,2,2,3,3,3,3,4,4,4,4,4,4)
  actual_ages <- 1:15
  cohorts <- 1:49
  
  # Initialize dataframe with all combinations
  survival_df <- expand.grid(
    iteration = 1:100,
    age = actual_ages,
    cohort = cohorts,
    sex = sex_label
  )
  survival_df$age_class <- age_classes[survival_df$age]
  
  # Calculate probabilities
  survival_probs <- numeric(nrow(survival_df))
  
  for(i in 1:100) {
    iter <- random_iterations[i]
    
    # Get fixed effects for this iteration
    phi_values <- sapply(1:4, function(p) {
      mean(rb[iter, 1:3, phiPB_indices[p]])
    })
    
    # Get random effects for this iteration
    eps_values <- sapply(1:49, function(p) {
      mean(rb[iter, 1:3, eps_phi_indices[p]])
    })
    
    # Calculate probabilities for all combinations of age and cohort
    for(c in 1:49) {
      idx <- which(survival_df$iteration == i & survival_df$cohort == c)
      fixed_effect <- phi_values[age_classes]
      random_effect <- eps_values[c]
      survival_probs[idx] <- plogis(fixed_effect + random_effect)
    }
  }
  
  survival_df$survival_prob <- survival_probs
  return(survival_df)
}

# Function to process breeding data with random effects for one sex
process_breeding <- function(rb, sex_label) {
  # Get indices for fixed and random effects
  psiPB_indices <- which(dimnames(rb)[[3]] %in% paste0("Mu.psiPB[", 1:7, "]"))
  eps_psi_indices <- which(dimnames(rb)[[3]] %in% paste0("eps.psi[", 1:49, "]"))
  
  # Randomly sample 100 iterations (using same seed for consistency)
  set.seed(123)
  random_iterations <- sample(1:250, 100)
  
  # Create age classes vector and ages
  age_classes <- c(1,1,1,1,1,2,3,4,5,6,7,7,7,7,7)
  actual_ages <- 1:15
  cohorts <- 1:49
  
  # Initialize dataframe with all combinations
  breeding_df <- expand.grid(
    iteration = 1:100,
    age = actual_ages,
    cohort = cohorts,
    sex = sex_label
  )
  breeding_df$age_class <- age_classes[breeding_df$age]
  
  # Calculate probabilities
  breeding_probs <- numeric(nrow(breeding_df))
  
  for(i in 1:100) {
    iter <- random_iterations[i]
    
    # Get fixed effects for this iteration
    psi_values <- sapply(1:7, function(p) {
      mean(rb[iter, 1:3, psiPB_indices[p]])
    })
    
    # Get random effects for this iteration
    eps_values <- sapply(1:49, function(p) {
      mean(rb[iter, 1:3, eps_psi_indices[p]])
    })
    
    # Calculate probabilities for all combinations of age and cohort
    for(c in 1:49) {
      idx <- which(breeding_df$iteration == i & breeding_df$cohort == c)
      fixed_effect <- psi_values[age_classes]
      random_effect <- eps_values[c]
      breeding_probs[idx] <- plogis(fixed_effect + random_effect)
    }
  }
  
  breeding_df$breeding_prob <- breeding_probs
  return(breeding_df)
}

# Process female model
load("outputs/model_estimates/model_2_F.Rdata")
female_survival_df <- process_survival(rb, "F")
female_breeding_df <- process_breeding(rb, "F")

# Process male model
load("outputs/model_estimates/model_2_M.Rdata")
male_survival_df <- process_survival(rb, "M")
male_breeding_df <- process_breeding(rb, "M")

# Combine the dataframes
survival_probabilities <- rbind(female_survival_df, male_survival_df)
breeding_probabilities <- rbind(female_breeding_df, male_breeding_df)

# Order the dataframes nicely
survival_probabilities <- survival_probabilities[order(survival_probabilities$sex, 
                                                       survival_probabilities$iteration,
                                                       survival_probabilities$cohort,
                                                       survival_probabilities$age),]
breeding_probabilities <- breeding_probabilities[order(breeding_probabilities$sex,
                                                       breeding_probabilities$iteration,
                                                       breeding_probabilities$cohort,
                                                       breeding_probabilities$age),]

# Display first few rows of each to verify
head(survival_probabilities)
head(breeding_probabilities)

# For survival
survival_array <- array(NA, dim = c(15, 49, 100, 2))
for(s in 1:2) {
  sex_label <- c("F", "M")[s]
  df_subset <- survival_probabilities[survival_probabilities$sex == sex_label, ]
  
  for(i in 1:100) {
    for(c in 1:49) {
      subset <- df_subset[df_subset$iteration == i & df_subset$cohort == c, ]
      survival_array[, c, i, s] <- subset$survival_prob[order(subset$age)]
    }
  }
}

# For breeding
breeding_array <- array(NA, dim = c(15, 49, 100, 2))
for(s in 1:2) {
  sex_label <- c("F", "M")[s]
  df_subset <- breeding_probabilities[breeding_probabilities$sex == sex_label, ]
  
  for(i in 1:100) {
    for(c in 1:49) {
      subset <- df_subset[df_subset$iteration == i & df_subset$cohort == c, ]
      breeding_array[, c, i, s] <- subset$breeding_prob[order(subset$age)]
    }
  }
}

# Create a list with metadata and arrays
matlab_data <- list(
  survival = survival_array,
  breeding = breeding_array,
  dimensions = list(
    dim1 = "age (1-15)",
    dim2 = "cohort (1-49)",
    dim3 = "iteration (1-100)",
    dim4 = "sex (1=F, 2=M)"
  )
)

# Export to MATLAB .mat file
writeMat("demographic_parameters_cohorts.mat",
         survival = survival_array,
         breeding = breeding_array)

# Optionally, also save as RData for R backup
save(survival_probabilities, breeding_probabilities, file = "demographic_parameters_cohorts.RData")

### Trend effects -----------------------------------------------------------

load("outputs/model_estimates/model_10_F.RData")
dimnames(rb)
