library(dplyr)
library(ipd)
source("chenchen_JG.R")
source("INT.R")

# DEXA
dexa <- read.table("DEXA.tab", header = TRUE)
ancestry <- read.table("Subpopulation.tab", header = TRUE)
dexa <- dexa %>% left_join(ancestry)
dexa$f.23248.2.0 <- dexa$f.23248.2.0 / 1000

dat <- readRDS('full_data.rds')
dat <- as.data.frame(dat)
colnames(dat) <- c("f.eid", "BMD", "PA", "age", "smoke", "drink", "driving", "computer", "tv", "sex")
dat <- dat[!dat$smoke == -3,]
dat <- dat[!dat$drink == -3,]
dat <- dat[!dat$driving == -3 & !dat$driving == -1,]
dat$driving[dat$driving == -10] <- 0.5
dat <- dat[!dat$computer == -3 & !dat$computer == -1,]
dat$computer[dat$computer == -10] <- 0.5
dat <- dat[!dat$tv == -3 & !dat$tv == -1,]
dat$tv[dat$tv == -10] <- 0.5
dat$SB <- dat$driving + dat$computer + dat$tv
dat$PA <- factor(dat$PA)

dexa <- dexa %>% left_join(dat, by = c("f.eid" = "f.eid"))

pheno_id <- "f.23248.2.0" 
cov_id <- c("f.21022.0.0", "f.22001.0.0", "f.23104.0.0", "PA", "smoke", "SB", "drink", "f.50.0.0", "f.23107.0.0", "f.23108.0.0")


model <- c("strong")
if(model == "strong"){
  W <- NULL
} else {
  W <- c("f.21022.0.0", "f.22001.0.0", "f.23104.0.0") # age, sex, BMI
}

ancestry_id <- c(3, 3001, 3002, 3003)

# South Asian subpopulation
test <- dexa %>%
  filter(f.21000.0.0 %in% ancestry_id) %>% 
  select(all_of(c(pheno_id, cov_id))) %>% 
  tidyr::drop_na(all_of(cov_id))

train <- dexa %>%
  filter(f.21000.0.0 %in% c(1, 1001, 1002, 1003)) %>%
  filter(!is.na(f.23248.2.0)) %>% 
  select(all_of(c(pheno_id, cov_id))) %>% 
  tidyr::drop_na(all_of(cov_id)) %>%
  select(-W)

# Train the model
model <- ranger::ranger(f.23248.2.0 ~ ., data = train, num.trees = 1000)

train$yhat <- predict(model, train)$predictions
summary(lm(f.23248.2.0 ~ yhat, data = train))

test$yhat <- predict(model, test)$predictions
summary(lm(f.23248.2.0 ~ yhat, data = test))

# Inverse normal transformation
test <- INT(test, "f.23248.2.0")
test <- INT(test, "yhat")

test$f.23248.2.0 <- test$f.23248.2.0_int
test$yhat <- test$yhat_int

test$set_label <- ifelse(is.na(test$f.23248.2.0), "unlabeled", "labeled")

formula <- f.23248.2.0 - yhat ~   f.21022.0.0 + f.22001.0.0 + f.23104.0.0

ppi <- ipd(formula, method = "ppi", model = "ols", data = test, label = "set_label") 
summary(ppi)

fake_dat <- test %>% filter(set_label == "labeled") %>% mutate(set_label = "unlabeled")
ppia <- ipd(formula, method = "ppi", model = "ols", data = rbind(test, fake_dat), label = "set_label")
summary(ppia)

cc <- chen_chen(Y = test %>% filter(set_label == "labeled") %>% pull(f.23248.2.0),
                Xlab = test %>% filter(set_label == "labeled") %>% select(f.21022.0.0, f.22001.0.0, f.23104.0.0),
                flab = test %>% filter(set_label == "labeled") %>% pull(yhat),
                fall = test %>% pull(yhat),
                Xall = test %>% select(f.21022.0.0, f.22001.0.0, f.23104.0.0))

result <- c()
result <- rbind(result, t(summary(std)$coefficients[, c("Estimate", "Std. Error")]))
result <- rbind(result, result[nrow(result),] / result[2,])
result <- rbind(result, t(summary(ppi)$coefficients[, c("Estimate", "Std.Error")]))
result <- rbind(result, result[nrow(result),] / result[2,])
result <- rbind(result, t(summary(ppia)$coefficients[, c("Estimate", "Std.Error")]))
result <- rbind(result, result[nrow(result),] / result[2,])
result <- rbind(result, t(cc))
result <- rbind(result, result[nrow(result),] / result[2,])

rownames(result) <- rep(c("Estimate", "SE", "SE Ratio"), 4)

result <- as.data.frame(result) %>% tibble::rownames_to_column("rowname")
result$rowname <- rep(c("Estimate", "SE", "SE Ratio"), 4)

result <- result %>%
  mutate(
    across(
      where(is.numeric), 
      ~ case_when(
        rowname == "Estimate" ~ formatC(., format = "f", digits = 3),
        rowname == "SE" ~ formatC(., format = "f", digits = 5),
        rowname == "SE Ratio" ~ formatC(., format = "f", digits = 3),
        TRUE ~ as.character(.)
      )
    )
  )

result