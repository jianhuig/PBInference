# Another Look at Inference After Prediction

This repository contains the code and results for the paper *"Another Look at Inference After Prediction"* by Jessica Gronsbell,  Jianhui Gao, Yaqi Shi, Zachary R. McCaw, and David Cheng. You can find the preprint here: [arXiv link](https://arxiv.org/abs/2411.19908).

## Overview

This paper revisits prediction-based (PB) inference methods and investigates their statistical efficiency. We analyze and compare several approaches, including the prediction-powered inference (PPI) method from Angelopoulos et al. (2023a) and the Chen and Chen (2000) estimator with theoretical and numerical evaluations.

The repository includes:
- Implementation of PB inference methods discussed in the paper.
- Simulations and analyses used to generate the results in the manuscript.

## Acknowledgments

This repository builds on the code from [PredictionBasedInference](https://github.com/keshav-motwani/PredictionBasedInference) repository. We thank [Motwani and Witten (2023)](https://www.jmlr.org/papers/volume24/23-0896/23-0896.pdf) for publicly sharing their code, which allowed us to further investigate those methods in this manuscript. 

## Useful Examples

```r
# Parameter Setup
n_train <- 10000
n_tot <- 10000
n_tests <- 1000
scnes <- c("1a")
n_sim <- 1

# Data Generation
sim_dat <- data_gen(n_train, n_test, n_val, n_sim)

# Calculate coefficient and SE
ppi_df <- predpowinf(sim_dat, n_sim)
ppi_full_df <- predpowinf_full(sim_dat, n_sim)
cc_df <- chen_chen(sim_dat, n_sim)
```

## Contact

For questions or issues, please contact [Yaqi Shi](mailto:yaqi.shi@mail.utoronto.com) or open an issue on this repository.

