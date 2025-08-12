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
│   ├── 📁 model/                    # Stan/R model definitions
│   ├── 📁 outputs/                  # HMM model results
│   ├── 📁 run_models/               # Model execution examples
│   │   ├── M2_example.R            # Trend model example
│   │   ├── M10_example.R           # Random effects model example
│   │   ├── covariates_birth_c.R    # Birth covariates (c)
│   │   ├── covariates_recruit_t.R  # Recruitment covariates (t)
│   │   ├── covariates_recruit_t1.R # Recruitment covariates (t-1)
│   │   ├── males_analysis.R        # Male-specific models
│   │   └── females_analysis.R      # Female-specific models
│   ├── 📁 source/                   # Core model execution functions
│   └── 📄 *.R                      # Root conversion scripts
├── 📁 02_AFR_markov_chains/         # Age at First Reproduction analysis
│   ├── 📁 data/                     # Processed HMM outputs
│   ├── 📁 functions/                # AFR calculation functions
│   ├── 📁 outputs/                  # Final AFR results
│   ├── 📁 resources/                # Supporting materials
│   ├── 📄 *.m                      # MATLAB analysis scripts
│   └── 📄 *.mat                    # MATLAB workspace files
└── 📄 README.md                    # This file
```

## 🔧 Prerequisites

### Required Software
- **R** (≥ 4.3.1) with RStan
- **MATLAB** (≥ R2020a) with Statistics Toolbox
- **Stan** (≥ 2.21.8) for Bayesian modeling
- **HPC cluster access** (recommended for large analyses)

### Required R Packages
```r
# Core HMM and Bayesian packages
install.packages(c("rstan", "loo", "ggplot2", "dplyr"))

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

# Load core functions
source("source/hmm_core_functions.R")
```

### Step 2: Run Example Models 📊
```r
# Execute trend model example
source("run_models/M2_example.R")

# Execute random effects model example  
source("run_models/M10_example.R")

# Run sex-specific analyses
source("run_models/males_analysis.R")
source("run_models/females_analysis.R")
```

### Step 3: Covariate Integration 📈
```r
# Birth covariates
source("run_models/covariates_birth_c.R")

# Recruitment covariates (current year)
source("run_models/covariates_recruit_t.R")

# Recruitment covariates (previous year)
source("run_models/covariates_recruit_t1.R")
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

% Run main AFR calculation
run('main_afr_cohort.m')

% Execute trend analysis
run('main_afr_trend.m')
```

## 📋 Detailed Workflow

### 🔄 HMM Modeling Process (`01_HMM_modelling/`)

**Core Architecture:**
- **`model/` folder**: Contains all HMM model definitions
- **`run_models/` folder**: Provides executable examples for different model configurations
- **`source/` folder**: Houses the computational engine for model execution
- **Root scripts**: Handle R-to-MATLAB conversion for downstream analysis

**Model Examples:**
- **M2**: Demonstrates trend model implementation
- **M10**: Showcases random effects modeling approach

**Covariate Integration:**
- **Birth effects (c)**: Environmental conditions at hatching
- **Recruitment effects (t)**: Current year recruitment conditions  
- **Recruitment effects (t-1)**: Previous year recruitment conditions

**Sex-Specific Analysis:**
- Separate modeling pipelines for males and females
- Accounts for sex-specific demographic patterns and life history differences

### 📈 AFR Calculation Process (`02_AFR_markov_chains/`)

**MATLAB-Based Analysis:**
- Utilizes HMM outputs from Stage 1
- Implements advanced Markov chain methods
- Generates robust AFR estimates with uncertainty quantification

**Key Components:**
- **Data processing**: Converts and structures HMM outputs
- **Markov chain simulation**: Estimates AFR distributions
- **Results compilation**: Generates publication-ready outputs

## 🔬 Cluster Computing Integration

This framework is designed for high-performance computing environments:

### Cluster Execution Strategy
```r
# The source/ folder contains the computational core
# The run_models/ scripts serve as lightweight job dispatchers
# This separation enables efficient cluster resource utilization
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

**Lead Author:** [Your Name]  
📧 Email: [your.email@institution.edu]  
🏛️ Affiliation: [Your Institution]

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

**Last updated:** January 2025  
**Repository maintained by:** [Your Name]
