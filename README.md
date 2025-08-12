# AFR-WA-JAE2025

# 🕊️ Albatross Demographic Analysis Framework

## 🌊 Overview

This repository contains the complete analytical framework and reproducible code for Hidden Markov Model (HMM) demographic analysis of albatross populations and Age at First Reproduction (AFR) calculations using Markov chain methods.

Our integrated approach combines state-of-the-art HMM modeling to capture complex demographic processes with subsequent AFR estimation through Markov chain Monte Carlo methods, providing comprehensive insights into albatross life history parameters and population dynamics.

## 🧬 Background

This research framework implements a two-stage analytical approach:
1. **Hidden Markov Models (HMM)** to model complex demographic states and transitions
2. **Markov Chain Monte Carlo (MCMC)** methods to derive Age at First Reproduction estimates from HMM outputs

This methodology is particularly suited for long-lived seabirds like albatrosses, where traditional demographic approaches may miss critical life history transitions and state-dependent processes.

## 📊 Key Features

- **🔍 Comprehensive HMM modeling** with trend and random effects
- **⚡ Cluster-compatible execution** for computationally intensive analyses
- **👥 Sex-specific analyses** for males and females separately
- **📈 Covariate integration** at birth, recruitment (t and t-1)
- **🔄 Seamless R-to-MATLAB workflow** for specialized AFR calculations
- **🎯 Age at First Reproduction estimation** using advanced Markov chain methods

## 🗂️ Repository Structure

```
albatross_demographic_analysis/
├── 📁 01_HMM_modelling/              # Hidden Markov Model analysis
│   ├── 📁 data/                     # Input demographic data
│   ├── 📁 function/                 # Custom functions for HMM
│   │   └── functions.R              # Core HMM functions
│   ├── 📁 model/                    # HMM model definitions
│   │   ├── model_c.R               # Birth effects model
│   │   ├── model_null.R            # Null model
│   │   ├── model_t_minus_1.R       # Recruitment effects (t-1)
│   │   ├── model_t.R               # Recruitment effects (t)
│   │   └── model_trends.R          # Temporal trend models
│   ├── 📁 outputs/                  # HMM model results
│   ├── 📁 run_models/               # Model execution scripts
│   │   ├── run_m2_f.R              # Model M2 for females
│   │   ├── run_m2_m.R              # Model M2 for males
│   │   ├── run_m10_f.R             # Model M10 for females
│   │   ├── run_m10_m.R             # Model M10 for males
│   │   ├── run_m22_f.R             # Model M22 for females
│   │   ├── run_m22_m.R             # Model M22 for males
│   │   ├── run_m23_f.R             # Model M23 for females
│   │   ├── run_m23_m.R             # Model M23 for males
│   │   ├── run_m30_f.R             # Model M30 for females
│   │   ├── run_m30_m.R             # Model M30 for males
│   │   ├── run_m36_f.R             # Model M36 for females
│   │   └── run_m36_m.R             # Model M36 for males
│   ├── 📁 source/                   # Core execution engine
│   │   ├── source_c.R              # Birth effects execution
│   │   ├── source_null.R           # Null model execution
│   │   ├── source_t.R              # Recruitment effects execution
│   │   └── source_trends.R         # Temporal trends execution
│   ├── 📄 01_HMM_modelling.Rproj    # R Project file
│   ├── 📄 convert_Routput_to_matlab_random_effects.R # R to MATLAB converter (random)
│   └── 📄 convert_Routput_to_matlab_trend.R         # R to MATLAB converter (trend)
├── 📁 02_AFR_markov_chains/         # Age at First Reproduction analysis
│   ├── 📁 data/                     # Processed HMM outputs
│   ├── 📁 functions/                # AFR calculation functions
│   ├── 📁 outputs/                  # Final AFR results
│   ├── 📁 resources/                # Supporting materials
│   ├── 📄 construct_F_matrix.m      # Fecundity matrix construction
│   ├── 📄 construct_U_matrix.m      # Survival matrix construction
│   ├── 📄 Info_init.m               # Initialization parameters
│   ├── 📄 interval_stats.m          # Interval statistics calculation
│   ├── 📄 invlogit.m                # Inverse logit transformation
│   ├── 📄 logit.m                   # Logit transformation
│   ├── 📄 longevity_stats.m         # Longevity statistics
│   ├── 📄 main_afr_cohort.m         # Main AFR analysis (cohort)
│   ├── 📄 main_afr_trend.m          # Main AFR analysis (trend)
│   ├── 📄 repro_stats.m             # Reproductive statistics
│   ├── 📄 workspace_afr_cohort.mat  # MATLAB workspace (cohort)
│   └── 📄 workspace_afr_trend.mat   # MATLAB workspace (trend)
└── 📄 README.md                     # This file
```

## 🔧 Prerequisites

### Required Software
- **R** (≥ 4.3.1) with RStan
- **MATLAB** (≥ R2020a) with Statistics Toolbox
- **Nimble** (≥ 1.3.0) for Bayesian modeling
- **HPC cluster access** (recommended for large analyses)

### Required R Packages
```r
# Core HMM and Bayesian packages
install.packages(c("nimble", "loo", "ggplot2", "dplyr"))

# Additional modeling packages
install.packages(c("tidyr", "gridExtra", "cowplot", "R.matlab"))
```

### Required MATLAB Toolboxes
- Statistics and Machine Learning Toolbox
- Optimization Toolbox (recommended)

## 🚀 Getting Started - Step by Step

### Step 1: HMM Model Execution 🔍

Navigate to the HMM modeling directory:
```r
# Set working directory
setwd("01_HMM_modelling/")
```

### Step 2: Run Example Models 📊
```r
# Execute Model M2 (trend models)
source("run_models/run_m2_f.R")  # Females
source("run_models/run_m2_m.R")  # Males

# Execute Model M10 (random effects models)  
source("run_models/run_m10_f.R") # Females
source("run_models/run_m10_m.R") # Males

# Additional model examples
source("run_models/run_m22_f.R") # Model M22 females
source("run_models/run_m22_m.R") # Model M22 males
source("run_models/run_m30_f.R") # Model M30 females
source("run_models/run_m30_m.R") # Model M30 males
```

### Step 3: Model Execution Framework 📈
```r
# The source/ folder contains the computational engine
# Each source file corresponds to a specific model type:

# Birth covariate models
source("source/source_c.R")                   # Basic birth effects
source("source/source_c_cov_psi.R")           # Birth covariates (psi)
source("source/source_c_cov_phi_psi.R")       # Birth covariates (phi & psi)

# Recruitment covariate models  
source("source/source_t.R")                   # Basic recruitment effects
source("source/source_t_cov_psi.R")           # Recruitment covariates (psi)
source("source/source_t_cov_phi.R")           # Recruitment covariates (phi)
source("source/source_t_cov_phi_psi.R")       # Recruitment covariates (phi & psi)

# Cohort and trend models
source("source/source_cohort.R")              # Basic cohort model
source("source/source_trends.R")              # Temporal trend model
source("source/source_null.R")                # Null model
```

### Step 4: Convert Outputs for MATLAB 🔄
```r
# Convert R outputs to MATLAB format
source("convert_Routput_to_matlab_random_effects.R")
source("convert_Routput_to_matlab_trend.R")
```

### Step 5: AFR Markov Chain Analysis 🎯
```matlab
% Switch to MATLAB environment
cd('02_AFR_markov_chains/')

% Run main AFR calculations
run('main_afr_cohort.m')  % Cohort-based AFR analysis
run('main_afr_trend.m')   % Trend-based AFR analysis

% Supporting calculations
run('construct_F_matrix.m')  % Build fecundity matrices
run('construct_U_matrix.m')  % Build survival matrices
run('longevity_stats.m')     % Calculate longevity statistics
run('repro_stats.m')         % Calculate reproductive statistics
```

## 📋 Detailed Workflow

### 🔄 HMM Modeling Process (`01_HMM_modelling/`)

**Core Architecture:**
- **`model/` folder**: Contains all HMM model definitions
- **`run_models/` folder**: Provides executable examples for different model configurations
- **`source/` folder**: Houses the computational engine for model execution
- **Root scripts**: Handle R-to-MATLAB conversion for downstream analysis

**Model Examples:**
- **M2**: Trend model implementation (available for both sexes)
- **M10**: Random effects modeling approach (available for both sexes)
- **M22, M23, M30, M36**: Additional model configurations with specific parameter sets

**Model Categories:**
- **Birth effect models**: `model_c.R`, `model_c_all_states.R`
- **Recruitment effect models**: `model_t.R` with temporal lags (t-1, t-2, t-3)
- **Covariate models**: Separate models for phi and psi parameters with birth (c) and recruitment (t) effects
- **Cohort models**: `model_cohort.R`, `model_cohort_random.R`, `model_cohort_quadra_random.R`
- **Trend models**: `model_trends.R` for temporal pattern analysis
- **Null model**: `model_null.R` as baseline comparison

**Sex-Specific Analysis:**
- Separate modeling pipelines for males and females
- Accounts for sex-specific demographic patterns and life history differences

### 📈 AFR Calculation Process (`02_AFR_markov_chains/`)

**MATLAB-Based Analysis:**
- **Matrix construction**: `construct_F_matrix.m` and `construct_U_matrix.m` for demographic matrices
- **Main AFR calculations**: `main_afr_cohort.m` and `main_afr_trend.m` for different analytical approaches
- **Statistical analysis**: `longevity_stats.m`, `repro_stats.m`, and `interval_stats.m` for demographic metrics
- **Utility functions**: `logit.m`, `invlogit.m` for data transformations
- **Initialization**: `Info_init.m` for parameter setup

**Key Components:**
- **Data processing**: Converts and structures HMM outputs from R
- **Markov chain simulation**: Estimates AFR distributions using demographic matrices
- **Results compilation**: Generates comprehensive demographic statistics and AFR estimates
- **Workspace management**: Pre-configured workspaces (`workspace_afr_cohort.mat`, `workspace_afr_trend.mat`)

## 🔬 Cluster Computing Integration

This framework is designed for high-performance computing environments:

### Cluster Execution Strategy
```r
# The source/ folder contains the computational engine
# The run_models/ scripts serve as lightweight job dispatchers
# This separation enables efficient cluster resource utilization

# Example: Model M2 execution
source("run_models/run_m2_f.R")  # Calls -> source("source/source_trends.R")
source("run_models/run_m2_m.R")  # Calls -> source("source/source_trends.R")

# The source files do the heavy computational work
# The run_models files handle model-specific parameters and sex separation
```

### Batch Processing
- Submit multiple model configurations simultaneously
- Parallel processing of sex-specific analyses  
- Efficient resource allocation for large parameter spaces

## 🔧 Adapting This Framework

### For Different Species
1. **Modify life history parameters** in model definitions
2. **Adjust covariate structures** based on species-specific ecology
3. **Update AFR calculation parameters** for species longevity patterns

### For Different Study Systems
1. **Customize HMM state definitions** for your demographic processes
2. **Modify covariate integration** based on available environmental data
3. **Adjust Markov chain parameters** for your data structure

## 📚 Associated Publication

**Citation:**
> [Author List]. (2025). [Title]. *Journal Name*, [Volume(Issue)], [Pages]. 
> DOI: [DOI will be added upon publication]

## 📊 Data Availability

- **Demographic data**: Available upon request
- **Model outputs**: Included in repository structure
- **AFR calculations**: Generated through provided scripts

## 👥 Authors & Contact

**Lead Author:** ****
📧 Email: ****
🏛️ Affiliation: ****

## 🙏 Acknowledgments

- **Field research teams** for data collection
- **HPC facility** for computational resources
- **Funding agencies** for project support
- **Collaborating institutions** for data sharing

## ⚖️ License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🚨 Research Implications

This analytical framework provides:
- **Advanced demographic modeling** for long-lived seabirds
- **Robust AFR estimation** critical for population projections
- **Flexible covariate integration** for environmental impact assessment
- **Scalable computational approach** for large-scale demographic studies

---

🕊️ *"Understanding albatross demography requires sophisticated methods that can capture the complexity of long-lived seabird life histories."*

**Last updated:** August 2025  
**Repository maintained by:** ****
