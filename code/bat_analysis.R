##required libraries
library(data.table)
library(dplyr)
library(magrittr)
library(ranger)
library(grf)

## required helper
source('code/helper/Qini_calculation_and_plot.R')
source('code/helper/TMSE_calculation.R')

## set seed
set.seed(121122)

## choose target
target <- 'converted'
prob_W <- 0.75

## read Hillstrom data
data <- fread('data/df_bat_preproc.csv')

## delete menth closing campaign and encode treatment
names(data)[which(names(data)== 'treatmentGroup')] <- 'W'
y_data <- data[, target, with = F]
data$Y <- data[, target, with = F]
data <- data[,-c('converted', 'checkoutAmount'), with = F]
data <- as.data.table(model.matrix(Y~-1+., data))
data$Y <- y_data

## split in training and test set
idx_train <- sample(1:nrow(data), round(0.9*nrow(data)), replace = F)
dt_train <- data[idx_train,]
dt_test <- data[-idx_train,]

## build random forest model and adjust target
rf <-  ranger(Y ~ ., data = dt_train[,-c('W'), with = F], min.node.size = 1000, num.trees = 2000)
rf0 <- ranger(Y ~ ., data = dt_train[W==0,-c('W'), with = F], min.node.size = 1000, num.trees = 2000)
rf1 <- ranger(Y ~ ., data = dt_train[W==1,-c('W'), with = F], min.node.size = 1000, num.trees = 2000)
pred_rf <- predict(rf, dt_test)$predictions
pred_rf0 <- predict(rf0, dt_test)$predictions
pred_rf1 <- predict(rf1, dt_test)$predictions
dt_test$Y_dr <- dt_test$Y-(prob_W*pred_rf0+(1-prob_W)*pred_rf1)
dt_test$Y_cond <- dt_test$Y-pred_rf
dt_test$Y_uc <- dt_test$Y-mean(dt_train$Y)

## calculate causal random forest
cf <- causal_forest(X = dt_train[,-c('W', 'Y'), with = F], Y = dt_train$Y, W = dt_train$W, W.hat = prob_W, min.node.size = 100)
pred_cf <- predict(cf, dt_test[,-c('W', 'Y_dr', 'Y_cond', 'Y_uc', 'Y'), with = F])$predictions

## plot Qini curve
plot_qini_curve(p_preds = pred_cf, p_dt_test = dt_test) 

## calculate transformed mean squared error and store results
Wp <- (dt_test$W-prob_W)/(prob_W*(1-prob_W))
options(scipen=999)
res_orig <- dt_test %$% calculate_Delta_TMSE(p_TY = Wp*Y, p_pred1 = pred_cf, p_pred2 = 0, p_round = 9)
res_cond <- dt_test %$% calculate_Delta_TMSE(p_TY = Wp*Y_cond, p_pred1 = pred_cf, p_pred2 = 0, p_round = 9)
res_uc <- dt_test %$% calculate_Delta_TMSE(p_TY = Wp*Y_uc, p_pred1 = pred_cf, p_pred2 = 0, p_round = 9)
res_dr <- dt_test %$% calculate_Delta_TMSE(p_TY = Wp*Y_dr, p_pred1 = pred_cf, p_pred2 = 0, p_round = 9)

TMSE_res <- data.frame(orig = c(as.character(res_orig$tmse), as.character(res_orig$CI), as.character(res_orig$width_CI)),
                       uc = c(as.character(res_uc$tmse), as.character(res_uc$CI), as.character(res_uc$width_CI)),
                       cond = c(as.character(res_cond$tmse), as.character(res_cond$CI), as.character(res_cond$width_CI)),
                       dr = c(as.character(res_dr$tmse), as.character(res_dr$CI), as.character(res_dr$width_CI)))
row.names(TMSE_res) <- c('TMSE', 'CI', 'CI_width')
                         
