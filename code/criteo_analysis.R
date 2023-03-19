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
target <- 'conversion'
prob_W <- 0.85

## read Criteo data
data <- fread('data/criteo-uplift-v2.1.csv')

## delete menth closing campaign and encode treatment
names(data)[which(names(data)== 'treatment')] <- 'W'
y_data <- data[, target, with = F]
data$Y <- data[, target, with = F]
data <- data[,-c('conversion', 'visit', 'exposure'), with = F]
data <- as.data.table(model.matrix(Y~-1+., data))
data$Y <- y_data

## split in training and test set
idx_train <- sample(1:nrow(data), round(0.1*nrow(data)), replace = F)
dt_train <- data[idx_train,]
dt_test <- data[-idx_train,]

## build random forest models and causal random forest
rf <-  ranger(Y ~ ., data = dt_train[,-c('W'), with = F], min.node.size = 1000, num.trees = 500)
rf0 <- ranger(Y ~ ., data = dt_train[W==0,-c('W'), with = F], min.node.size = 1000, num.trees = 500)
rf1 <- ranger(Y ~ ., data = dt_train[W==1,-c('W'), with = F], min.node.size = 1000, num.trees = 500)
idx_crf <- sample(1:nrow(dt_train),100000)
cf <- causal_forest(X = dt_train[idx_crf,-c('W', 'Y'), with = F], Y = dt_train$Y[idx_crf], W = dt_train$W[idx_crf], min.node.size = 1000, W.hat = prob_W, num.trees = 500)

## make predictions on the test data (split data in 10 parts for computational reasons)
pred_rf <- c()
pred_rf0 <- c()
pred_rf1 <- c()
pred_cf <- c()
for(i_fold in 1:12){
  pred_rf <- c(pred_rf, predict(rf, dt_test[(1:1000000)+(i_fold-1)*1000000,])$predictions)
  pred_rf0 <- c(pred_rf0, predict(rf0, dt_test[(1:1000000)+(i_fold-1)*1000000,])$predictions)
  pred_rf1 <- c(pred_rf1, predict(rf1, dt_test[(1:1000000)+(i_fold-1)*1000000,])$predictions)
  pred_cf <- c(pred_cf, predict(cf, dt_test[(1:1000000)+(i_fold-1)*1000000,-c('W', 'Y'), with = F])$predictions)
}
pred_rf <- c(pred_rf, predict(rf, dt_test[12000001:nrow(dt_test),])$predictions)
pred_rf0 <- c(pred_rf0, predict(rf0, dt_test[12000001:nrow(dt_test),])$predictions)
pred_rf1 <- c(pred_rf1, predict(rf1, dt_test[12000001:nrow(dt_test),])$predictions)
pred_cf <- c(pred_cf, predict(cf, dt_test[12000001:nrow(dt_test),-c('W', 'Y'), with = F])$predictions)

## save predictions
fwrite(data.table(pred_rf, pred_rf0, pred_rf1, pred_cf), 'simulation_results/criteo/preds_test.csv')

## add adjusted outcomes to test set
dt_test$Y_dr <- dt_test$Y-(prob_W*pred_rf0+(1-prob_W)*pred_rf1)
dt_test$Y_cond <- dt_test$Y-pred_rf
dt_test$Y_uc <- dt_test$Y-mean(dt_train$Y)

## plot Qini curve
plot_qini_curve(p_preds = pred_cf, p_dt_test = dt_test) 

## calculate transformed mean squared error and store results
Wp <- (dt_test$W-prob_W)/(prob_W*(1-prob_W))
options(scipen=999)
res_orig <- dt_test %$% calculate_Delta_TMSE(p_TY = Wp*Y, p_pred1 = pred_cf, p_pred2 = 0, p_round = 6)
res_cond <- dt_test %$% calculate_Delta_TMSE(p_TY = Wp*Y_cond, p_pred1 = pred_cf, p_pred2 = 0, p_round = 6)
res_uc <- dt_test %$% calculate_Delta_TMSE(p_TY = Wp*Y_uc, p_pred1 = pred_cf, p_pred2 = 0, p_round = 6)
res_dr <- dt_test %$% calculate_Delta_TMSE(p_TY = Wp*Y_dr, p_pred1 = pred_cf, p_pred2 = 0, p_round = 6)

TMSE_res <- data.frame(orig = c(as.character(res_orig$tmse), as.character(res_orig$CI), as.character(res_orig$width_CI)),
                       uc = c(as.character(res_uc$tmse), as.character(res_uc$CI), as.character(res_uc$width_CI)),
                       cond = c(as.character(res_cond$tmse), as.character(res_cond$CI), as.character(res_cond$width_CI)),
                       dr = c(as.character(res_dr$tmse), as.character(res_dr$CI), as.character(res_dr$width_CI)))
row.names(TMSE_res) <- c('TMSE', 'CI', 'CI_width')
