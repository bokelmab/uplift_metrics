## required libraries
library(ranger)
library(data.table)
library(dplyr)
library(magrittr)
library(pROC) ## AUC
library(grf)
library(uplift)
library(tools4uplift)
library(MASS)

## seed
set.seed(031022)

## simulate training data
data <- generate_causal_data(10000, 6, dgp = 'nw2', sigma.noise = 0.5)
dt_train <- data.table(Y = data$Y, X1 = data$X[,1], X2 = data$X[,2], X3 = data$X[,3], 
                       X4 = data$X[,4], X5 = data$X[,5], X6 = data$X[,6], W = data$W, tau = data$tau, 
                       m = data$m, e = data$e)

## train causal random forest
cf <- causal_forest(X = dt_train[,c('X1', 'X2', 'X3', 'X4', 'X5', 'X6'), with = F], Y = dt_train$Y, W = dt_train$W,
                    num.trees=1000, sample.fraction=0.5,  
                    honesty=T, tune.parameters='none', honesty.fraction=0.5, ci.group.size=2, compute.oob.predictions=T)

## function to calculate Qini curve
calculate_qini <- function(p_data, p_theor){
  
  ate_deciles <- c()
  for(i_dec in 1:10){
    data_dec <- p_data[pred_cf > quantile(p_data$pred_cf, 1-i_dec/10),]
    if(p_theor){
      ate_new <- mean(data_dec[W==1,]$tau)*sum(data_dec$W)
    }else{
      ate_new <- sum(data_dec[W==1,]$Y)-sum(data_dec[W==0,]$Y)*sum(data_dec$W)/sum(1-data_dec$W)
    }
    
    ate_deciles <- c(ate_deciles, ate_new)
  }
  return(ate_deciles)
}



## simulate test data
sim_theo <- NULL
sim_emp <- NULL
for(i_sim in 1:100){
  data <- generate_causal_data(5000, 6, dgp = 'nw2', sigma.noise = 0.5)
  dt_test <- data.table(Y = data$Y, X1 = data$X[,1], X2 = data$X[,2], X3 = data$X[,3], 
                        X4 = data$X[,4], X5 = data$X[,5], X6 = data$X[,6], W = data$W, tau = data$tau, 
                        m = data$m, e = data$e)
  dt_test$pred_cf <- predict(cf, newdata = dt_test[,c('X1', 'X2', 'X3', 'X4', 'X5', 'X6'), with = F])$predictions
  emp_qini <- calculate_qini(dt_test, F)
  theo_qini <- calculate_qini(dt_test, T)
  sim_theo %<>% rbind(theo_qini)
  sim_emp %<>% rbind(emp_qini)
}
ave_theo <- apply(sim_theo, 2, mean)

## plot results
layout(matrix(c(1,1,2,3,4,4), ncol=2, byrow=TRUE), heights = c(0.1,0.8, 0.1))

## title
par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5,0.5,"Qini curve and theoretical performance measure",cex=2,font=1.5)

par(mar=c(3, 3, 1, 2))
## plot single estimated Qini curve
plot((0:10)/10*2500, c(0, ave_theo), type = 'l', ylim = c(min(sim_emp), max(sim_emp)), lwd = 3, ylab = 'incremental gain', xlab = expression('N'[W]*'(s)'),
     mgp=c(1.5,0.5,0))
lines((0:10)/10*2500, c(0, sim_emp[1,]), type = 'l', col = 'red', lwd = 3)

## plot confidence intervall
plot((0:10)/10*2500, c(0, ave_theo), type = 'l', ylim = c(min(sim_emp), max(sim_emp)), lwd = 3, ylab = 'incremental gain', xlab = expression('N'[W]*'(s)'),
     mgp=c(1.5,0.5,0))
for(i_emp in 1:nrow(sim_emp)){
  lines((0:10)/10*2500, c(0, sim_emp[i_emp,]), type = 'l', col = 'red', lwd = 1)
}
lines((0:10)/10*2500, c(0, apply(sim_emp, 2, quantile, probs = 0.975)), col = 'black', lwd = 3, lty = 3)
lines((0:10)/10*2500, c(0, apply(sim_emp, 2, quantile, probs = 0.025)), col = 'black', lwd = 3, lty = 3)
lines(lines((0:10)/10*2500, c(0, ave_theo), type = 'l', lwd = 3))

##legend
par(mar=c(1, 1, 0.5, 0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c('black', 'red', 'black')
legend(x = "bottom",inset = 0,
       legend = c(expression('ATE'[s]%.%'N'[W]), 'Qini', '95% CI'), 
       col=plot_colors, lty = c(1,1,3), lwd=5, cex=1, horiz = TRUE)


