
# Implementation of the PPI method
predpowinf <- function(sim_dat_tv, 
                       n_sim) {
  
  # Wang et al. (2020) refer to the labeled data set as \test" data, and to the unlabeled data as \validation" data.
  
  # Identify the labeled and unlabeled set
  test_data <- filter(sim_dat_tv, set == "testing")
  test_data$diff <- test_data$pred - test_data$y
  val_data <- filter(sim_dat_tv, set == "validation")
  
  # Set up the placeholder
  df <- c()
  
  # Prediction powered estimation
  val_pred_x <- lm(pred ~ x1, data = val_data)
  
  rectifier <- lm(diff ~ x1, data = test_data)
  
  theta_hat_pp <- coef(val_pred_x) - coef(rectifier)
  
  # Calculate se
  X_tilde <- model.matrix(~ val_data$x1)
  
  X <- model.matrix(~ test_data$x1)
  
  N <- nrow(X_tilde)
  
  # Calculate Sigma tilde
  Sigma_tilde <- 1 / N * crossprod(X_tilde)
  
  # Calculate M tilde
  M_tilde <- matrix(0, nrow = ncol(X_tilde), ncol = ncol(X_tilde))
  
  for (j in 1:N) {
    M_tilde <- M_tilde + (val_data$pred[j] - crossprod(X_tilde[j, ], coef(val_pred_x))[1, 1])^2 * tcrossprod(X_tilde[j, ])
  }
  
  M_tilde <- M_tilde / N
  
  # Calculate V tilde
  Sigma_tilde_inv <- solve(Sigma_tilde)
  V_tilde <- Sigma_tilde_inv %*% M_tilde %*% Sigma_tilde_inv
  
  # Calculate Sigma
  n <- nrow(X)
  Sigma <- 1 / n * crossprod(X)
  
  # Calculate M
  M <- matrix(0, nrow = ncol(X), ncol = ncol(X))
  for (j in 1:n) {
    M <- M + (test_data$pred[j] - test_data$y[j] - crossprod(X[j, ], coef(rectifier))[1, 1])^2 * tcrossprod(X[j, ])
  }
  M <- M / n
  
  # Calculate V
  Sigma_inv <- solve(Sigma)
  V <- Sigma_inv %*% M %*% Sigma_inv
  
  # Calculate se
  se <- sqrt((V[2, 2] / n) + (V_tilde[2, 2] / N))
  
  
  # Partial correlation
  labeled_data <- filter(sim_dat_tv, set == "testing")
  unlabeled_data <- filter(sim_dat_tv, set == "validation")
  
  X_lab <- labeled_data[, "x1"]
  X_unlab <- unlabeled_data[, "x1"]
  Y_lab <- labeled_data[, "y"]
  
  pred_lab <- labeled_data[, "pred"]
  pred_unlab <- unlabeled_data[, "pred"]
  
  model_beta_lab <- lm(Y_lab ~ X_lab)
  model_gamma_lab <- lm(pred_lab ~ X_lab)
  model_gamma_unlab <- lm(pred_unlab ~ X_unlab)
  
  res_beta <- residuals(model_beta_lab)
  res_gamma_lab <- residuals(model_gamma_lab)
  res_gamma_unlab <- residuals(model_gamma_unlab)
  res_gamma_adj <- pred_lab - cbind(1, X_lab) %*% coef(model_gamma_unlab)
  
  ppi_cov <- cov(res_beta, res_gamma_lab)
  ppi_cov_adj <- cov(res_beta, res_gamma_adj)
  ppi_vbeta <- var(res_beta)
  ppi_vgammalab <- var(res_gamma_lab)
  ppi_vgammaunlab <- var(res_gamma_unlab)
  
  # results
  df <- data.frame(
    sim = n_sim,
    ppi_beta = theta_hat_pp[2],
    ppi_se = se,
    ppi_cov = ppi_cov,
    ppi_cov_adj = ppi_cov_adj,
    ppi_vbeta = ppi_vbeta,
    ppi_vgammalab = ppi_vgammalab,
    ppi_vgammaunlab = ppi_vgammaunlab
  )
  
  return(df)
}


# Implementation of PPI method using full data
predpowinf_full <- function(sim_dat_tv, 
                            n_sim) {
  
  # Wang et al. (2020) refer to the labeled data set as \test" data, and to the unlabeled data as \validation" data.
  
  # Identify the labeled set
  test_data <- filter(sim_dat_tv, set == "testing")
  test_data$diff <- test_data$pred - test_data$y
  
  # Set up the placeholder
  df <- c()
  
  # Prediction powered estimation
  val_pred_x <- lm(pred ~ x1, data = sim_dat_tv)
  
  rectifier <- lm(diff ~ x1, data = test_data)
  
  theta_hat_pp <- coef(val_pred_x) - coef(rectifier)
  
  # Calculate se
  X_tilde <- model.matrix(~ sim_dat_tv$x1)
  
  X <- model.matrix(~ test_data$x1)
  
  N <- nrow(X_tilde)
  
  # Calculate Sigma tilde
  Sigma_tilde <- 1 / N * crossprod(X_tilde)
  
  # Calculate M tilde
  M_tilde <- matrix(0, nrow = ncol(X_tilde), ncol = ncol(X_tilde))
  
  for (j in 1:N) {
    M_tilde <- M_tilde + (sim_dat_tv$pred[j] - crossprod(X_tilde[j, ], coef(val_pred_x))[1, 1])^2 * tcrossprod(X_tilde[j, ])
  }
  
  M_tilde <- M_tilde / N
  
  # Calculate V tilde
  Sigma_tilde_inv <- solve(Sigma_tilde)
  V_tilde <- Sigma_tilde_inv %*% M_tilde %*% Sigma_tilde_inv
  
  # Calculate Sigma
  n <- nrow(X)
  Sigma <- 1 / n * crossprod(X)
  
  # Calculate M
  M <- matrix(0, nrow = ncol(X), ncol = ncol(X))
  for (j in 1:n) {
    M <- M + (test_data$pred[j] - test_data$y[j] - crossprod(X[j, ], coef(rectifier))[1, 1])^2 * tcrossprod(X[j, ])
  }
  M <- M / n
  
  # Calculate V
  Sigma_inv <- solve(Sigma)
  V <- Sigma_inv %*% M %*% Sigma_inv
  
  # Calculate se
  se <- sqrt((V[2, 2] / n) + (V_tilde[2, 2] / N))
  
  # Partial correlation
  labeled_data <- filter(sim_dat_tv, set == "testing")
  
  X_lab <- labeled_data[, "x1"]
  X_all <- sim_dat_tv[, "x1"]
  Y_lab <- labeled_data[, "y"]
  
  pred_lab <- labeled_data[, "pred"]
  pred_all <- sim_dat_tv[, "pred"]
  
  model_beta_lab <- lm(Y_lab ~ X_lab)
  model_gamma_lab <- lm(pred_lab ~ X_lab)
  model_gamma_all <- lm(pred_all ~ X_all)
  
  res_beta <- residuals(model_beta_lab)
  res_gamma_lab <- residuals(model_gamma_lab)
  res_gamma_all <- residuals(model_gamma_all)
  res_gamma_adj <- pred_lab - cbind(1, X_lab) %*% coef(model_gamma_all)
  
  ppi_full_cov <- cov(res_beta, res_gamma_lab)
  ppi_full_cov_adj <- cov(res_beta, res_gamma_adj)
  ppi_full_vbeta <- var(res_beta)
  ppi_full_vgammalab <- var(res_gamma_lab)
  ppi_full_vgammaall <- var(res_gamma_all)
  
  # results
  df <- data.frame(
    sim = n_sim,
    ppi_full_beta = theta_hat_pp[2],
    ppi_full_se = se,
    ppi_full_cov = ppi_full_cov,
    ppi_full_cov_adj = ppi_full_cov_adj,
    ppi_full_vbeta = ppi_full_vbeta,
    ppi_full_vgammalab = ppi_full_vgammalab,
    ppi_full_vgammaall = ppi_full_vgammaall
  )
  
  return(df)
}


# Implementation of naive and classic method with calculation for true beta
truth_and_nc <- function(sim_dat_tv, n_sim) {
  
  val_data <- filter(sim_dat_tv, set == "validation")
  observed_data <- filter(sim_dat_tv, set == "testing")
  
  df <- c()
  
  truth_model <- lm(y ~ x1, data = sim_dat_tv) # val_data
  nc_model <- lm(pred ~ x1, data = sim_dat_tv) # val_data
  observed_model <- lm(y ~ x1, data = observed_data)
  
  truth_beta <- truth_model %>%
    tidy() %>%
    filter(term == "x1") %>%
    pull(estimate)
  truth_sd <- truth_model %>%
    tidy() %>%
    filter(term == "x1") %>%
    pull(std.error)
  truth_t <- truth_model %>%
    tidy() %>%
    filter(term == "x1") %>%
    pull(statistic)
  truth_p <- truth_model %>%
    tidy() %>%
    filter(term == "x1") %>%
    pull(p.value)
  
  
  nc_beta <- nc_model %>%
    tidy() %>%
    filter(term == "x1") %>%
    pull(estimate)
  nc_sd <- nc_model %>%
    tidy() %>%
    filter(term == "x1") %>%
    pull(std.error)
  nc_t <- nc_model %>%
    tidy() %>%
    filter(term == "x1") %>%
    pull(statistic)
  nc_p <- nc_model %>%
    tidy() %>%
    filter(term == "x1") %>%
    pull(p.value)
  
  
  observed_beta <- observed_model %>%
    tidy() %>%
    filter(term == "x1") %>%
    pull(estimate)
  observed_sd <- observed_model %>%
    tidy() %>%
    filter(term == "x1") %>%
    pull(std.error)
  observed_t <- observed_model %>%
    tidy() %>%
    filter(term == "x1") %>%
    pull(statistic)
  observed_p <- observed_model %>%
    tidy() %>%
    filter(term == "x1") %>%
    pull(p.value)
  
  # results
  df <- data.frame(
    sim = n_sim,
    truth_beta = truth_beta,
    truth_sd = truth_sd,
    truth_t = truth_t,
    truth_p = truth_p,
    nc_beta = nc_beta,
    nc_sd = nc_sd,
    nc_t = nc_t,
    nc_p = nc_p,
    observed_beta = observed_beta,
    observed_sd = observed_sd,
    observed_t = observed_t,
    observed_p = observed_p
  )
  
  return(df)
}


# Implementation of the Chen and Chen method
chen_chen <- function(sim_dat_tv, 
                      n_sim) {
  
  # Identify the validation (unlabeled) and the testing (labeled) data
  labeled_data <- filter(sim_dat_tv, set == "testing")
  unlabeled_data <- filter(sim_dat_tv, set == "validation")
  all_data <- sim_dat_tv
  
  n <- nrow(labeled_data)
  N <- nrow(all_data)
  
  X_lab <- labeled_data[, "x1"]
  X_all <- all_data[, "x1"]
  
  Y_lab <- labeled_data[, "y"]
  
  pred_lab <- labeled_data[, "pred"]
  pred_all <- all_data[, "pred"]
  
  model_beta_lab <- lm(Y_lab ~ X_lab)
  model_gamma_lab <- lm(pred_lab ~ X_lab)
  model_gamma_all <- lm(pred_all ~ X_all)
  
  beta_lab <- coef(model_beta_lab)
  gamma_lab <- coef(model_gamma_lab)
  gamma_all <- coef(model_gamma_all)
  
  X_lab_int <- cbind(1, X_lab)
  
  D1 <- -crossprod(X_lab_int, X_lab_int) / n
  D2 <- -crossprod(X_lab_int, X_lab_int) / n
  
  S2 <- X_lab_int * residuals(model_beta_lab)
  S_tilde2 <- X_lab_int * residuals(model_gamma_lab)
  
  C11 <- crossprod(S2, S2) / n
  C12 <- crossprod(S2, S_tilde2) / n
  C22 <- crossprod(S_tilde2, S_tilde2) / n
  
  D1_inv <- solve(D1)
  C22_inv <- solve(C22)
  
  D1C12 <- D1_inv %*% C12
  D1C12C22 <- D1C12 %*% C22_inv
  
  theta_hat_cc <- as.vector(beta_lab - D1C12C22 %*% D2 %*% (gamma_lab - gamma_all))
  
  Omega <- (D1_inv %*% C11 %*% D1_inv - (1 - n / N) * D1C12C22 %*% t(D1C12)) / n
  
  se_beta_bar <- sqrt(diag(Omega))
  
  # results
  df <- data.frame(
    sim = n_sim,
    cc_beta = theta_hat_cc[2],
    cc_se = se_beta_bar[2]
  )
  
  return(df)
}


# Implementation of the Seemingly Unrelated Regression (SUR) method
sur <- function(sim_dat_tv, 
                n_sim) {
  
  # Identify the validation (unlabeled) and the testing (labeled) data
  labeled_data <- filter(sim_dat_tv, set == "testing")
  unlabeled_data <- filter(sim_dat_tv, set == "validation")
  all_data <- sim_dat_tv
  
  n <- nrow(labeled_data)
  N <- nrow(all_data)
  
  X_lab <- labeled_data[, "x1"]
  X_all <- all_data[, "x1"]
  
  Y_lab <- labeled_data[, "y"]
  
  pred_lab <- labeled_data[, "pred"]
  pred_all <- all_data[, "pred"]
  
  model_beta_lab <- lm(Y_lab ~ X_lab)
  model_gamma_lab <- lm(pred_lab ~ X_lab)
  model_gamma_all <- lm(pred_all ~ X_all)
  model_cor_lab_all <- lm(Y_lab ~ X_lab + pred_lab)
  
  beta_lab <- coef(model_beta_lab)
  gamma_lab <- coef(model_gamma_lab)
  gamma_all <- coef(model_gamma_all)
  
  theta_tilde <- sum(residuals(model_beta_lab) * residuals(model_gamma_lab)) / sum(residuals(model_gamma_lab)^2)
  theta_hat_sur <- as.vector(beta_lab - theta_tilde %*% (gamma_lab - gamma_all))
  
  # Calculate the variance
  sigma_1_sqr <- (summary(model_beta_lab)$sigma)^2
  sigma_1_2_sqr <- (summary(model_cor_lab_all)$sigma)^2
  
  X_lab_int <- cbind(1, X_lab)
  X_all_int <- cbind(1, X_all)
  
  x_lab_prod_inv <- solve(t(X_lab_int) %*% X_lab_int)
  x_all_prod_inv <- solve(t(X_all_int) %*% X_all_int)
  
  se <- sqrt(diag(sigma_1_2_sqr * ((n - 1 - 1) / (n - 1 - 2)) * (x_lab_prod_inv - x_all_prod_inv) + sigma_1_sqr * x_all_prod_inv))
  
  # results
  df <- data.frame(
    sim = n_sim,
    sur_beta = theta_hat_sur[2],
    sur_se = se[2]
  )
  
  return(df)
}


# Implementation of the Chen and Chen method with modification of using unlabeled data
chen_chen_unlab <- function(sim_dat_tv, 
                            n_sim) {
  
  # Identify the validation (unlabeled) and the testing (labeled) data
  labeled_data <- filter(sim_dat_tv, set == "testing")
  unlabeled_data <- filter(sim_dat_tv, set == "validation")
  all_data <- unlabeled_data
  
  n <- nrow(labeled_data)
  N <- nrow(all_data)
  
  X_lab <- labeled_data[, "x1"]
  X_all <- all_data[, "x1"]
  
  Y_lab <- labeled_data[, "y"]
  
  pred_lab <- labeled_data[, "pred"]
  pred_all <- all_data[, "pred"]
  
  model_beta_lab <- lm(Y_lab ~ X_lab)
  model_gamma_lab <- lm(pred_lab ~ X_lab)
  model_gamma_all <- lm(pred_all ~ X_all)
  
  beta_lab <- coef(model_beta_lab)
  gamma_lab <- coef(model_gamma_lab)
  gamma_all <- coef(model_gamma_all)
  
  X_lab_int <- cbind(1, X_lab)
  
  D1 <- -crossprod(X_lab_int, X_lab_int) / n
  D2 <- -crossprod(X_lab_int, X_lab_int) / n
  
  S2 <- X_lab_int * residuals(model_beta_lab)
  S_tilde2 <- X_lab_int * residuals(model_gamma_lab)
  
  C11 <- crossprod(S2, S2) / n
  C12 <- crossprod(S2, S_tilde2) / n
  C22 <- crossprod(S_tilde2, S_tilde2) / n
  
  D1_inv <- solve(D1)
  C22_inv <- solve(C22)
  
  D1C12 <- D1_inv %*% C12
  D1C12C22 <- D1C12 %*% C22_inv
  
  
  theta_hat_cc <- as.vector(beta_lab - D1C12C22 %*% D2 %*% (gamma_lab - gamma_all))
  
  Omega <- (D1_inv %*% C11 %*% D1_inv - (1 - n / N) * D1C12C22 %*% t(D1C12)) / n
  
  se_beta_bar <- sqrt(diag(Omega))
  
  
  # results
  df <- data.frame(
    sim = n_sim,
    cc_unlab_beta = theta_hat_cc[2],
    cc_unlab_se = se_beta_bar[2]
  )
  
  return(df)
}


# Implementation of the PSPA method
pspa <- function(sim_dat_tv, 
                 n_sim){
  
  # Data extraction
  X_lab <- sim_dat_tv %>% dplyr::filter(set == "testing") %>% pull(x1)
  X_unlab <- sim_dat_tv %>% dplyr::filter(set == "validation") %>% pull(x1)
  
  Y_lab <- sim_dat_tv %>% dplyr::filter(set == "testing") %>% pull(y)
  
  Yhat_lab <- sim_dat_tv %>% dplyr::filter(set == "testing") %>% pull(pred)
  Yhat_unlab <- sim_dat_tv %>% dplyr::filter(set == "validation") %>% pull(pred)
  
  pspa_df <- pspa_y(X_lab, X_unlab, Y_lab, Yhat_lab, Yhat_unlab, alpha = 0.05, method = "ols")[2, 1:2]
  
  # results
  df <- data.frame(
    sim = n_sim,
    pspa_beta = pspa_df$Estimate,
    pspa_se = pspa_df$Std.Error
  )
  
  return(df)
  
}


# Helper function for check the coverage
coverage <- function(true, 
                     est, 
                     se) {
  
  ((true >= est - 1.96 * se) & (true <= est + 1.96 * se))
  
}



















