# Another Look at Inference After Prediction

This repository contains the code and results for the paper *"Another Look at Inference After Prediction"* by Jessica Gronsbell,  Jianhui Gao, Yaqi Shi, Zachary R. McCaw, and David Cheng. You can find the preprint [here](https://arxiv.org/abs/2411.19908).

## Overview

Our paper investigates the statistical efficiency of prediction-based (PB) inference methods. We analyze and compare several PB inference approaches, including the prediction-powered inference (PPI) method from [Angelopoulos et al. (2023)](https://www.science.org/doi/10.1126/science.adi6000) and the [Chen and Chen (2000)](https://www.jstor.org/stable/2680690) estimator with theoretical and numerical evaluations.

The repository includes:
- Implementation of PB inference methods discussed in the paper.
- Simulations and analyses used to generate the results in the paper.
- Code for reproduce our UK Biobank Analysis, but access to UK Biobank is required as the data cannot be released.

To run the simulations, the required packages include dplyr, broom, and parallel.

## Simple Example

``` r
# Set the working directory
# Update this path based on your local setup
setwd('~/PBInference/Scripts')

# Load dependencies
library(dplyr)
library(broom)
library(parallel)

# Load functions
source('estimators.R')
source('data_gen.R')

# Set up parameters
n_train <- 10000
n_test <- 1000
n_val <- 9000

# Scenario details are described in the paper
scenario <- "1a"

# Set seed
set.seed(2025)

# Generate data
sim_dat <- data_gen(n_train, n_test, n_val, scenario)

# Calculate coefficient and standard errors
ppi <- predpowinf(sim_dat)[, c("ppi_beta", "ppi_se")]
ppi_full <- predpowinf_full(sim_dat)[, c("ppi_full_beta", "ppi_full_se")]
cc <- chen_chen(sim_dat)[, c("cc_beta", "cc_se")]
```

## Acknowledgments

This repository builds on the code from [PredictionBasedInference](https://github.com/keshav-motwani/PredictionBasedInference) repository. We thank [Motwani and Witten (2023)](https://www.jmlr.org/papers/volume24/23-0896/23-0896.pdf) for publicly sharing their code, which allowed us to further investigate those methods in this manuscript. 


## Contact

For questions or issues, please contact [Yaqi Shi](mailto:yaqi.shi@mail.utoronto.com) or open an issue on this repository.

