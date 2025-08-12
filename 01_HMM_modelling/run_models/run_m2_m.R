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
lapply(c("tidyverse", "dplyr", "tidyr", "nimble", "glue", "reshape2"), library, character.only = TRUE)

# Define directory paths for data, functions, output, and model
WorkDir <- getwd() # Working directory
DataDir <- paste(WorkDir, "data", sep = "/") # Data folder
OutDir <- paste(WorkDir, "output", sep = "/") # Output folder
ModDir <- paste(WorkDir, "model", sep = "/") # Models folder
FunDir <- paste(WorkDir, "function", sep = "/") # Models folder
SourceDir <- paste(WorkDir, "source", sep = "/") # Models folder
source(glue("{FunDir}/functions.R"))

source(glue("{ModDir}/model_trends.R"))
model <- model_trends
m_ <- "model_2"
sex <- "M"

trend.phi <- "none"
trend.psi <- "none"

param_rm <- c("b1phiPB", "b2phiPB", "b1psiPB", "b2psiPB")

source(glue("{SourceDir}/source_trends.R"))
