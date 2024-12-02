# Another Look at Inference After Prediction

This repository contains the code and results for the paper *"Another Look at Inference After Prediction"* by Jessica Gronsbell,  Jianhui Gao, Yaqi Shi, Zachary R. McCaw, and David Cheng. You can find the preprint [here](https://arxiv.org/abs/2411.19908).

## Overview

Our paper investigates the statistical efficiency of prediction-based (PB) inference methods. We analyze and compare several approaches, including the prediction-powered inference (PPI) method from Angelopoulos et al. (2023a) and the Chen and Chen (2000) estimator with theoretical and numerical evaluations.

The repository includes:
- Implementation of PB inference methods discussed in the paper.
- Simulations and analyses used to generate the results in the paper.

## Acknowledgments

This repository builds on the code from [PredictionBasedInference](https://github.com/keshav-motwani/PredictionBasedInference) repository. We thank [Motwani and Witten (2023)](https://www.jmlr.org/papers/volume24/23-0896/23-0896.pdf) for publicly sharing their code, which allowed us to further investigate those methods in this manuscript. 

## Example

``` r
# Set the working directory
# Update this path based on your local setup
setwd('~/PBInference/Scripts')

# Load dependencies
library(dplyr)
library(gam)
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
ppi <- predpowinf(sim_dat)[2:3]
ppi_full <- predpowinf_full(sim_dat)[2:3]
cc <- chen_chen(sim_dat)[2:3]
```

### Results

| Method         | Estimate | Standard Error |
|----------------|----------|----------------|
| $\text{PPI}$   | 2.346244 | 0.03989220     |
| $\text{PPI}_a$ | 2.354460 | 0.03929328     |
| $\text{CC}$    | 2.354395 | 0.03929419     |



## Contact

For questions or issues, please contact [Yaqi Shi](mailto:yaqi.shi@mail.utoronto.com) or open an issue on this repository.

