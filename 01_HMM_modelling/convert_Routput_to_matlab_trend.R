# Clear environment and load required packages
rm(list = ls())
library(readr)
library(ggplot2)
library(R.matlab)

# Function to process survival data with trend effects
process_survival_trend <- function(rb, sex_label) {
  # Get indices for parameters
  phiPB_indices <- which(dimnames(rb)[[3]] %in% paste0("Mu.phiPB[", 1:4, "]"))
  b1phi_index <- which(dimnames(rb)[[3]] == "b1phiPB")
  b2phi_index <- which(dimnames(rb)[[3]] == "b2phiPB")
  
  # Randomly sample 100 iterations
  set.seed(123)
  random_iterations <- sample(1:250, 100)
  
  # Create age classes vector and ages
  age_classes <- c(1,1,2,2,2,3,3,3,3,4,4,4,4,4,4)
  actual_ages <- 1:15
  years <- 1:49
  
  # Initialize dataframe with all combinations
  survival_df <- expand.grid(
    iteration = 1:100,
    age = actual_ages,
    year = years,
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
    
    # Get trend coefficients (note: for males, b2 will be 0)
    b1 <- mean(rb[iter, 1:3, b1phi_index])
    b2 <- mean(rb[iter, 1:3, b2phi_index])
    
    # Calculate probabilities for all combinations of age and year
    for(y in years) {
      idx <- which(survival_df$iteration == i & survival_df$year == y)
      fixed_effect <- phi_values[age_classes]
      trend_effect <- b1 * y + b2 * y^2  # b2 will be 0 for males
      survival_probs[idx] <- plogis(fixed_effect + trend_effect)
    }
  }
  
  survival_df$survival_prob <- survival_probs
  return(survival_df)
}

# Function to process breeding data with trend effects
process_breeding_trend <- function(rb, sex_label) {
  # Get indices for parameters
  psiPB_indices <- which(dimnames(rb)[[3]] %in% paste0("Mu.psiPB[", 1:7, "]"))
  b1psi_index <- which(dimnames(rb)[[3]] == "b1psiPB")
  b2psi_index <- which(dimnames(rb)[[3]] == "b2psiPB")
  
  # Randomly sample 100 iterations
  set.seed(123)
  random_iterations <- sample(1:250, 100)
  
  # Create age classes vector and ages
  age_classes <- c(1,1,1,1,1,2,3,4,5,6,7,7,7,7,7)
  actual_ages <- 1:15
  years <- 1:49
  
  # Initialize dataframe with all combinations
  breeding_df <- expand.grid(
    iteration = 1:100,
    age = actual_ages,
    year = years,
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
    
    # Get trend coefficients (note: for males, b2 will be 0)
    b1 <- mean(rb[iter, 1:3, b1psi_index])
    b2 <- mean(rb[iter, 1:3, b2psi_index])
    
    # Calculate probabilities for all combinations of age and year
    for(y in years) {
      idx <- which(breeding_df$iteration == i & breeding_df$year == y)
      fixed_effect <- psi_values[age_classes]
      trend_effect <- b1 * y + b2 * y^2  # b2 will be 0 for males
      breeding_probs[idx] <- plogis(fixed_effect + trend_effect)
    }
  }
  
  breeding_df$breeding_prob <- breeding_probs
  return(breeding_df)
}

# Process female model
load("outputs/model_estimates/model_10_F.Rdata")
female_survival_df <- process_survival_trend(rb, "F")
female_breeding_df <- process_breeding_trend(rb, "F")

# Process male model
load("outputs/model_estimates/model_4_M.Rdata")
male_survival_df <- process_survival_trend(rb, "M")
male_breeding_df <- process_breeding_trend(rb, "M")

# Combine dataframes
survival_probabilities <- rbind(female_survival_df, male_survival_df)
breeding_probabilities <- rbind(female_breeding_df, male_breeding_df)

# Convert to arrays for MATLAB export
# Reshape data into 4D arrays: age x year x iteration x sex
survival_array <- array(NA, dim = c(15, 49, 100, 2))  # 2 for both sexes
breeding_array <- array(NA, dim = c(15, 49, 100, 2))

# Fill arrays - Females (sex index 1)
for(i in 1:100) {
  for(y in 1:49) {
    subset_surv <- female_survival_df[female_survival_df$iteration == i & 
                                        female_survival_df$year == y, ]
    survival_array[, y, i, 1] <- subset_surv$survival_prob[order(subset_surv$age)]
    
    subset_breed <- female_breeding_df[female_breeding_df$iteration == i & 
                                         female_breeding_df$year == y, ]
    breeding_array[, y, i, 1] <- subset_breed$breeding_prob[order(subset_breed$age)]
  }
}

# Fill arrays - Males (sex index 2)
for(i in 1:100) {
  for(y in 1:49) {
    subset_surv <- male_survival_df[male_survival_df$iteration == i & 
                                      male_survival_df$year == y, ]
    survival_array[, y, i, 2] <- subset_surv$survival_prob[order(subset_surv$age)]
    
    subset_breed <- male_breeding_df[male_breeding_df$iteration == i & 
                                       male_breeding_df$year == y, ]
    breeding_array[, y, i, 2] <- subset_breed$breeding_prob[order(subset_breed$age)]
  }
}

# Export to MATLAB
writeMat("demographic_parameters_trend.mat",
         survival = survival_array,
         breeding = breeding_array)

# Optional: save as RData
save(survival_probabilities, breeding_probabilities, 
     file = "demographic_parameters_trend.RData")