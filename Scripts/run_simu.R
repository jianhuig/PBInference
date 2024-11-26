
# ------------------------ Load Packages ------------------------ #

library(dplyr)
library(tidyr)
library(broom)
library(gam)
library(parallel)
library(pspa)


# ------------------------ Source Scripts ------------------------ #

source("data_gen.R")
source("estimators.R")


# ------------------------ Simulation Setup ------------------------ #

# Set up parameter
n_sim <- 2000
n_train <- 10000
n_tot <- 10000
n_tests <- c(1000, 2000, 3000, 4000, 5000)
n_cores <- detectCores()
scnes <- c("1a")

# Set seed for simulations
set.seed(2023)

# Results names
methods <- c("naive", "val*", "observed", 
             "predpowinf", "predpowinf-full", 
             "chen-chen", "sur", "chen-chen-unlab", "pspa")
beta_names <- c("nc_beta", "truth_beta", "observed_beta", 
                "ppi_beta", "ppi_full_beta", 
                "cc_beta", "sur_beta", "cc_unlab_beta", "pspa_beta")
se_names <- c("nc_sd", "truth_sd", "observed_sd", 
              "ppi_se", "ppi_full_se", 
              "cc_se", "sur_se", "cc_unlab_se", "pspa_se")


# ------------------------ Run Simulation ------------------------ #


for (i in 1:length(scnes)) {
  
  sce <- scnes[i]
  
  for (j in 1:length(n_tests)) {
    
    n_test <- n_tests[j]
    n_val <- n_tot - n_test
    
    # Data generation
    sim_dat_tv_tot <- data_gen(n_train, n_test, n_val, n_sim)
    
    # Run simulations
    sim_function <- function(k, sim_dat_tv_tot) {
      
      sim_dat_tv <- filter(sim_dat_tv_tot, sim == k)
      
      list(
        R2 = mean(sim_dat_tv$R2),
        truth_nc_df = truth_and_nc(sim_dat_tv, k),
        ppi_df = predpowinf(sim_dat_tv, k),
        ppi_full_df = predpowinf_full(sim_dat_tv, k),
        cc_df = chen_chen(sim_dat_tv, k),
        sur_df = sur(sim_dat_tv, k),
        cc_unlab_df = chen_chen_unlab(sim_dat_tv, k),
        pspa_df = pspa(sim_dat_tv, k)
      )
    }
    
    # Run the function in parallel
    results <- mclapply(1:n_sim, sim_function, sim_dat_tv_tot, mc.cores = n_cores)
    
    # Combine the results
    R2 <- do.call(rbind, lapply(results, `[[`, "R2"))
    truth_nc_df <- do.call(rbind, lapply(results, `[[`, "truth_nc_df"))
    ppi_df <- do.call(rbind, lapply(results, `[[`, "ppi_df"))
    ppi_full_df <- do.call(rbind, lapply(results, `[[`, "ppi_full_df"))
    cc_df <- do.call(rbind, lapply(results, `[[`, "cc_df"))
    sur_df <- do.call(rbind, lapply(results, `[[`, "sur_df"))
    cc_unlab_df <- do.call(rbind, lapply(results, `[[`, "cc_unlab_df"))
    pspa_df <- do.call(rbind, lapply(results, `[[`, "pspa_df"))
    
    
    df <- cbind(truth_nc_df, ppi_df, ppi_full_df, cc_df, sur_df, cc_unlab_df, pspa_df, R2)
    
    estimates <- df[, beta_names]
    colnames(estimates) <- methods
    
    reported_ses <- df[, se_names]
    colnames(reported_ses) <- methods
    
    result <- list()
    
    for (l in 1:length(methods)) {
      method <- methods[l]
      
      estimate <- estimates[, method]
      reported_se <- reported_ses[, method]
      z_stat <- estimate / reported_se
      p_value <- 2 * (1 - pnorm(abs(z_stat)))
      power <- ifelse((p_value < 0.05), 1, 0)
      cilength <- 2 * 1.96 * reported_se
      R2 <- R2
      
      result <- c(result, list(data.frame(
        n_train = n_train, n_test = n_test, n_val = n_val, R2 = R2,
        estimate = estimate, reported_se = reported_se,
        method = method, z_stat = z_stat, p_value = p_value,cilength = cilength,
        power = power, sce = sce
      )))
    }
    
    dir.create("results")
    file <- paste0("results/main_sim_power_n_test_", n_test, "_scenario_", sce, ".rds")
    saveRDS(do.call(rbind, result), file)
  }
  
}



