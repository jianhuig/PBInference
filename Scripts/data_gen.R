

# Calculate the true value of y
true_y <- function(df, 
                   sce) {
  
  if (sce == "1a") {
    
    df <- df %>% mutate(y = x1 + x2 + x3 + 0.5 * x4 + 0.5 * x5 + 0.25 * x6 + e_g)
    
  } else if (sce == "1b") {
    
    df <- df %>% mutate(y = x1 + x2 + x3 + 0.5 * x4 + 0.5 * x5 + 0.25 * x6 + e_g)
    
  } else if (sce == "2a") { 
    
    df <- df %>% mutate(y = x1 + 0.2 * x2 + 0.1 * x3 + 0.1 * x4 + 0.1 * x5 + 0.1 * x6 +
                          0.1 * x7 + 0.1 * x8 + 0.1 * x9 + 0.1 * x10 + e_g)
  } else if (sce == "2b") {
    
    df <- df %>% mutate(y = 2 * x1 + 0.25 * x2 + 0.25 * x3 + 0.25 * x4 + 0.25 * x1 * x4 + e_g)
    
  } else {
    
    df <- df %>% mutate(y = x1 + e_g)
    print("Scenario DNE")
    
  }
  
  return(df$y)
}

# Fit the prediction model to generate yhat
prediction_model <- function(train_data, 
                             sce) {
  
  if (sce == "1a") {
    
    pred_model <- gam(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11, data = train_data)
    
  } else if (sce == "1b") {
    
    pred_model <- gam(y ~ x1 + X6 + x1 * x6, data = train_data)
    
  } else if (sce == "2a") {
    
    pred_model <- gam(y ~ x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10, data = train_data)
    
  } else if (sce == "2b") {
    
    pred_model <- gam(y ~ x2 + x3 + x4 + x2 * x3 + x2 * x4 + x3 * x4, data = train_data)
    
  } else {
    
    pred_model <- gam(y ~ x1, data = train_data)
    print("Scenario DNE")
    
  }
  
  return(pred_model)
}

# Generate the data for simulation
data_gen <- function(n_train, 
                     n_test, 
                     n_val, 
                     n_sim) {
  
  # Set up mu and sigma
  correlation <- 0.4
  sigma <- matrix(correlation, nrow = 11, ncol = 11)
  diag(sigma) <- 1
  mu <- rep(0, 11)
  
  # Generate training data
  train_data <- as.data.frame(MASS::mvrnorm(n_train, mu = mu, Sigma = sigma))
  
  train_data <- train_data %>%
    rename_with(~ paste0("x", 1:11)) %>%
    mutate(
      e_g = rnorm(n_train, mean = 0, sd = 1),
      set = "training",
      sim = 0
    )
  
  train_data$y <- true_y(train_data, sce)
  
  # Set the model for generating true y
  pred_model <- prediction_model(train_data, sce)
  
  # Generate test and validation data
  generate_sim_data <- function(i) {
    n <- n_test + n_val
    
    # Generate test and validation data
    data <- as.data.frame(MASS::mvrnorm(n, mu = mu, Sigma = sigma))
    
    data <- data %>%
      rename_with(~ paste0("x", 1:11)) %>%
      mutate(
        e_g = rnorm(n, mean = 0, sd = 1),
        set = rep(c("testing", "validation"), c(n_test, n_val)),
        sim = i
      )
    
    data$y <- true_y(data, sce)
    
    # Generate prediction
    pred <- predict(pred_model, data)
    data$pred <- pred
    
    # Calculate R2
    SST <- sum((data$y - mean(data$y))^2)
    SSR <- sum((data$y - pred)^2)
    data$R2 <- 1 - (SSR / SST)
    
    return(data)
  }
  
  # Run the function in parallel
  sim_dat_list <- mclapply(1:n_sim, generate_sim_data, mc.cores = n_cores)
  
  # Combine the list of data frames into a single data frame
  test_val_data <- do.call(rbind, sim_dat_list)
  
  return(test_val_data)
}



















