## required libraries
library(data.table)
library(magrittr)

## required functions
source('code/helper/Qini_calculation_and_plot.R')

############################################################################################################################
################### Sample size reduction for TMSE #########################################################################
var_red_TMSE <- readRDS('simulation_results/summary_var_red/var_red_TMSE.RDS')

setting <- row.names(var_red_TMSE)[1]
dgp_list <- c("aw2", "nw2")
sigma_noise_list <- c(0.5, 1, 2)

################## conditional mean adjustment ############################################################################
tmse_results <- NULL
for(i_dgp in 1:length(dgp_list)){
  
  for(i_noise in 1:length(sigma_noise_list)){
    
    path_sim_data <- paste0('simulation_results/', dgp_list[i_dgp], '/noise_', sigma_noise_list[i_noise])
    var_red_TMSE_sel <- var_red_TMSE[row.names(var_red_TMSE) == paste0(dgp_list[i_dgp], sigma_noise_list[i_noise]), 'cond']
    
    tmse_cond <- c()
    tmse_orig <- c()
    tmse_cond_red <- c()
    for(i_list_data in 1:10){
      
      list_data <- readRDS(paste0(path_sim_data, '/list_data', i_list_data, '.RDS'))
      for(i_sim in 1:length(list_data)){
        
        ## get reduced data set
        dt_test <- list_data[[i_sim]]
        dt_red <- dt_test[1:ceiling((1-var_red_TMSE_sel/100)*nrow(dt_test)),]
        
        ## calculate TMSE with adjusted outcomes
        calc_mse <- function(p_Y, p_W, p_P, p_pred){
          return(mean((p_W*p_Y/p_P - (1-p_W)*p_Y/(1-p_P) - p_pred)^2))
        }
        diff_cond_red <- (dt_red %$% calc_mse(Y_cond, W, e, 0)) - (dt_red %$% calc_mse(Y_cond, W, e, pred_cf))
        diff_cond <- (dt_test %$% calc_mse(Y_cond, W, e, 0)) - (dt_test %$% calc_mse(Y_cond, W, e, pred_cf))
        diff_orig <- (dt_test %$% calc_mse(Y, W, e, 0)) - (dt_test %$% calc_mse(Y, W, e, pred_cf))
        tmse_cond_red <- c(tmse_cond_red, diff_cond_red)
        tmse_cond <- c(tmse_cond, diff_cond)
        tmse_orig <- c(tmse_orig, diff_orig)
      }
    }
    
    tmse_results %<>% rbind(data.frame(mean_cond_red = mean(tmse_cond_red), mean_cond = mean(tmse_cond), mean_orig = mean(tmse_orig),
                                       sd_cond_red = sd(tmse_cond_red), sd_cond = sd(tmse_cond), sd_orig = sd(tmse_orig), frac = (1-var_red_TMSE_sel/100)))
    
    
  }
}

## plot 
layout(matrix(c(1,2,3), ncol=1, byrow=TRUE), heights = c(0.03,0.34, 0.04))

## title
par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5,0.5,expression('Test set size reduction MSE'[W]),cex=2, font=2)
par(mar=c(4, 4, 1, 4))
#plot.new()
tmse_results %$% plot(frac, (sd_cond_red/sd_orig)^2, ylim = c(0,1.1), xlim = c(0,1), 
                      ylab = expression('Var['*Delta*'MSE'[W]*'(cond)]/Var['*Delta*'MSE'[W]*'(orig)]'), xlab = 'fraction of test set', mgp=c(2,0.8,0),
                      main = '', cex.axis = 1.3, cex.lab = 1.3)
tmse_results %$% points(rep(1, nrow(tmse_results)), (sd_cond/sd_orig)^2)
tmse_results[1,] %$% lines(c(frac, 1), c((sd_cond_red/sd_orig)^2, (sd_cond/sd_orig)^2), col = 'red', lwd = 2)
tmse_results[2,] %$% lines(c(frac, 1), c((sd_cond_red/sd_orig)^2, (sd_cond/sd_orig)^2), col = 'blue', lwd = 2)
tmse_results[3,] %$% lines(c(frac, 1), c((sd_cond_red/sd_orig)^2, (sd_cond/sd_orig)^2), col = 'green', lwd = 2)
tmse_results[4,] %$% lines(c(frac, 1), c((sd_cond_red/sd_orig)^2, (sd_cond/sd_orig)^2), col = 'orange', lwd = 2)
tmse_results[5,] %$% lines(c(frac, 1), c((sd_cond_red/sd_orig)^2, (sd_cond/sd_orig)^2), col = 'yellow', lwd = 2)
tmse_results[6,] %$% lines(c(frac, 1), c((sd_cond_red/sd_orig)^2, (sd_cond/sd_orig)^2), col = 'black', lwd = 2)

##legend
par(mar=c(1, 1, 0.5, 0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c('red', 'blue', 'green', 'orange', 'yellow', 'black')
legend(x = "bottom",inset = 0,
       legend = c('aw0.5', 'aw1', 'aw2', 'nw0.5', 'nw1', 'nw2'), 
       col=plot_colors, lwd=5, cex=1, horiz = TRUE)


row.names(tmse_results) <- paste0(rep(dgp_list,c(3,3)), sigma_noise_list)
saveRDS(tmse_results, 'simulation_results/summary_var_red/tmse_ss_reduction.RDS')

################## doubly-robust adjustment ############################################################################
tmse_results <- NULL
for(i_dgp in 1:length(dgp_list)){
  
  for(i_noise in 1:length(sigma_noise_list)){
    
    path_sim_data <- paste0('simulation_results/', dgp_list[i_dgp], '/noise_', sigma_noise_list[i_noise])
    var_red_TMSE_sel <- var_red_TMSE[row.names(var_red_TMSE) == paste0(dgp_list[i_dgp], sigma_noise_list[i_noise]), 'dr']
    
    tmse_dr <- c()
    tmse_orig <- c()
    tmse_dr_red <- c()
    for(i_list_data in 1:10){
      
      list_data <- readRDS(paste0(path_sim_data, '/list_data', i_list_data, '.RDS'))
      for(i_sim in 1:length(list_data)){
        
        ## get reduced data set
        dt_test <- list_data[[i_sim]]
        dt_red <- dt_test[1:ceiling((1-var_red_TMSE_sel/100)*nrow(dt_test)),]
        
        ## calculate TMSE with adjusted outcomes
        calc_mse <- function(p_Y, p_W, p_P, p_pred){
          return(mean((p_W*p_Y/p_P - (1-p_W)*p_Y/(1-p_P) - p_pred)^2))
        }
        diff_dr_red <- (dt_red %$% calc_mse(Y_dr, W, e, 0)) - (dt_red %$% calc_mse(Y_dr, W, e, pred_cf))
        diff_dr <- (dt_test %$% calc_mse(Y_dr, W, e, 0)) - (dt_test %$% calc_mse(Y_dr, W, e, pred_cf))
        diff_orig <- (dt_test %$% calc_mse(Y, W, e, 0)) - (dt_test %$% calc_mse(Y, W, e, pred_cf))
        tmse_dr_red <- c(tmse_dr_red, diff_dr_red)
        tmse_dr <- c(tmse_dr, diff_dr)
        tmse_orig <- c(tmse_orig, diff_orig)
      }
    }
    
    tmse_results %<>% rbind(data.frame(mean_dr_red = mean(tmse_dr_red), mean_dr = mean(tmse_dr), mean_orig = mean(tmse_orig),
                                       sd_dr_red = sd(tmse_dr_red), sd_dr = sd(tmse_dr), sd_orig = sd(tmse_orig), frac = (1-var_red_TMSE_sel/100)))
    
    
  }
}

## plot 
layout(matrix(c(1,2,3), ncol=1, byrow=TRUE), heights = c(0.03,0.34, 0.04))

## title
par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5,0.5,expression('Test set size reduction MSE'[W]),cex=2, font=2)
par(mar=c(4, 4, 1, 4))
#plot.new()
tmse_results %$% plot(frac, (sd_dr_red/sd_orig)^2, ylim = c(0,1.1), xlim = c(0,1), 
                      ylab = expression('Var['*Delta*'MSE'[W]*'(dr)]/Var['*Delta*'MSE'[W]*'(orig)]'), xlab = 'fraction of test set', mgp=c(2,0.8,0),
                      main = '', cex.axis = 1.3, cex.lab = 1.3)
tmse_results %$% points(rep(1, nrow(tmse_results)), (sd_dr/sd_orig)^2)
tmse_results[1,] %$% lines(c(frac, 1), c((sd_dr_red/sd_orig)^2, (sd_dr/sd_orig)^2), col = 'red', lwd = 2)
tmse_results[2,] %$% lines(c(frac, 1), c((sd_dr_red/sd_orig)^2, (sd_dr/sd_orig)^2), col = 'blue', lwd = 2)
tmse_results[3,] %$% lines(c(frac, 1), c((sd_dr_red/sd_orig)^2, (sd_dr/sd_orig)^2), col = 'green', lwd = 2)
tmse_results[4,] %$% lines(c(frac, 1), c((sd_dr_red/sd_orig)^2, (sd_dr/sd_orig)^2), col = 'orange', lwd = 2)
tmse_results[5,] %$% lines(c(frac, 1), c((sd_dr_red/sd_orig)^2, (sd_dr/sd_orig)^2), col = 'yellow', lwd = 2)
tmse_results[6,] %$% lines(c(frac, 1), c((sd_dr_red/sd_orig)^2, (sd_dr/sd_orig)^2), col = 'black', lwd = 2)

##legend
par(mar=c(1, 1, 0.5, 0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c('red', 'blue', 'green', 'orange', 'yellow', 'black')
legend(x = "bottom",inset = 0,
       legend = c('aw0.5', 'aw1', 'aw2', 'nw0.5', 'nw1', 'nw2'), 
       col=plot_colors, lwd=5, cex=1, horiz = TRUE)


row.names(tmse_results) <- paste0(rep(dgp_list,c(3,3)), sigma_noise_list)
saveRDS(tmse_results, 'simulation_results/summary_var_red/tmse_ss_reduction_dr.RDS')

################## unconditional mean adjustment ############################################################################
tmse_results <- NULL
for(i_dgp in 1:length(dgp_list)){
  
  for(i_noise in 1:length(sigma_noise_list)){
    
    path_sim_data <- paste0('simulation_results/', dgp_list[i_dgp], '/noise_', sigma_noise_list[i_noise])
    var_red_TMSE_sel <- var_red_TMSE[row.names(var_red_TMSE) == paste0(dgp_list[i_dgp], sigma_noise_list[i_noise]), 'uc']
    
    tmse_uc <- c()
    tmse_orig <- c()
    tmse_uc_red <- c()
    for(i_list_data in 1:10){
      
      list_data <- readRDS(paste0(path_sim_data, '/list_data', i_list_data, '.RDS'))
      for(i_sim in 1:length(list_data)){
        
        ## get reduced data set
        dt_test <- list_data[[i_sim]]
        dt_red <- dt_test[1:ceiling((1-var_red_TMSE_sel/100)*nrow(dt_test)),]
        
        ## calculate TMSE with adjusted outcomes
        calc_mse <- function(p_Y, p_W, p_P, p_pred){
          return(mean((p_W*p_Y/p_P - (1-p_W)*p_Y/(1-p_P) - p_pred)^2))
        }
        diff_uc_red <- (dt_red %$% calc_mse(Y_uc, W, e, 0)) - (dt_red %$% calc_mse(Y_uc, W, e, pred_cf))
        diff_uc <- (dt_test %$% calc_mse(Y_uc, W, e, 0)) - (dt_test %$% calc_mse(Y_uc, W, e, pred_cf))
        diff_orig <- (dt_test %$% calc_mse(Y, W, e, 0)) - (dt_test %$% calc_mse(Y, W, e, pred_cf))
        tmse_uc_red <- c(tmse_uc_red, diff_uc_red)
        tmse_uc <- c(tmse_uc, diff_uc)
        tmse_orig <- c(tmse_orig, diff_orig)
      }
    }
    
    tmse_results %<>% rbind(data.frame(mean_uc_red = mean(tmse_uc_red), mean_uc = mean(tmse_uc), mean_orig = mean(tmse_orig),
                                       sd_uc_red = sd(tmse_uc_red), sd_uc = sd(tmse_uc), sd_orig = sd(tmse_orig), frac = (1-var_red_TMSE_sel/100)))
    
    
  }
}

## plot 
layout(matrix(c(1,2,3), ncol=1, byrow=TRUE), heights = c(0.03,0.34, 0.04))

## title
par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5,0.5,expression('Test set size reduction MSE'[W]),cex=2, font=2)
par(mar=c(4, 4, 1, 4))
#plot.new()
tmse_results %$% plot(frac, (sd_uc_red/sd_orig)^2, ylim = c(0,1.1), xlim = c(0,1), 
                      ylab = expression('Var['*Delta*'MSE'[W]*'(uc)]/Var['*Delta*'MSE'[W]*'(orig)]'), xlab = 'fraction of test set', mgp=c(2,0.8,0),
                      main = '', cex.axis = 1.3, cex.lab = 1.3)
tmse_results %$% points(rep(1, nrow(tmse_results)), (sd_uc/sd_orig)^2)
tmse_results[1,] %$% lines(c(frac, 1), c((sd_uc_red/sd_orig)^2, (sd_uc/sd_orig)^2), col = 'red', lwd = 2)
tmse_results[2,] %$% lines(c(frac, 1), c((sd_uc_red/sd_orig)^2, (sd_uc/sd_orig)^2), col = 'blue', lwd = 2)
tmse_results[3,] %$% lines(c(frac, 1), c((sd_uc_red/sd_orig)^2, (sd_uc/sd_orig)^2), col = 'green', lwd = 2)
tmse_results[4,] %$% lines(c(frac, 1), c((sd_uc_red/sd_orig)^2, (sd_uc/sd_orig)^2), col = 'orange', lwd = 2)
tmse_results[5,] %$% lines(c(frac, 1), c((sd_uc_red/sd_orig)^2, (sd_uc/sd_orig)^2), col = 'yellow', lwd = 2)
tmse_results[6,] %$% lines(c(frac, 1), c((sd_uc_red/sd_orig)^2, (sd_uc/sd_orig)^2), col = 'black', lwd = 2)

##legend
par(mar=c(1, 1, 0.5, 0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c('red', 'blue', 'green', 'orange', 'yellow', 'black')
legend(x = "bottom",inset = 0,
       legend = c('aw0.5', 'aw1', 'aw2', 'nw0.5', 'nw1', 'nw2'), 
       col=plot_colors, lwd=5, cex=1, horiz = TRUE)


row.names(tmse_results) <- paste0(rep(dgp_list,c(3,3)), sigma_noise_list)
saveRDS(tmse_results, 'simulation_results/summary_var_red/tmse_ss_reduction_uc.RDS')


#############################################################################################################
################################# sample size reduction for Qini curve ######################################
var_red_qini <- readRDS('simulation_results/summary_var_red/var_red_qini.RDS')

################################ conditional mean adjustment ###############################################
for(i_dgp in 1:length(dgp_list)){
  
  for(i_noise in 1:length(sigma_noise_list)){
    
    path_sim_data <- paste0('simulation_results/', dgp_list[i_dgp], '/noise_', sigma_noise_list[i_noise])
    var_red_qini_sel <- var_red_qini[row.names(var_red_qini) == paste0(dgp_list[i_dgp], sigma_noise_list[i_noise]), 'cond_min']
    qini_orig <- NULL
    qini_cond <- NULL
    for(i_list_data in 1:10){
      
      list_data <- readRDS(paste0(path_sim_data, '/list_data', i_list_data, '.RDS'))
      
      for(i_sim in 1:length(list_data)){
        
        ## get reduced data set
        dt_test <- list_data[[i_sim]]
        dt_red <- dt_test[1:ceiling((1-var_red_qini_sel/100)*nrow(dt_test)),]
        
        qini_curve <- c()
        qini_curve_cond <- c()
        for(i_dec in 9:0){
          
          ## Qini curve and standard deviation with original target 
          res_qini <- dt_test %$% calc_CATE(p_Y = Y, p_W = W, p_pred = pred_cf, p_cut = i_dec/10)
          qini_curve <- c(qini_curve, res_qini$qini_mean)
          
          ## Qini curve and standard deviation with conditional mean adjusted target
          res_qini_cond <- dt_red %$% calc_CATE(p_Y = Y_cond, p_W = W, p_pred = pred_cf, p_cut = i_dec/10)
          res_qini_cond$qini_mean <- res_qini_cond$qini_mean * nrow(dt_test)/nrow(dt_red)
          qini_curve_cond <- c(qini_curve_cond, res_qini_cond$qini_mean)
        }
        
        ## add results per decile to previous simulation runs
        qini_orig %<>% rbind(data.frame(first = qini_curve[1], second = qini_curve[2], third = qini_curve[3], fourth = qini_curve[4], fifth = qini_curve[5],
                                        sixth = qini_curve[6], seventh = qini_curve[7], eigth = qini_curve[8], nineth = qini_curve[9], tenth = qini_curve[10]))
        qini_cond %<>% rbind(data.frame(first = qini_curve_cond[1], second = qini_curve_cond[2], third = qini_curve_cond[3], fourth = qini_curve_cond[4], fifth = qini_curve_cond[5],
                                        sixth = qini_curve_cond[6], seventh = qini_curve_cond[7], eigth = qini_curve_cond[8], nineth = qini_curve_cond[9], tenth = qini_curve_cond[10]))
      }
      
      
    }
    qini_orig_mean <- apply(qini_orig, 2, mean)
    upper_orig <- apply(qini_orig, 2, quantile, p = 0.975)
    lower_orig <- apply(qini_orig, 2, quantile, p = 0.025)
    qini_cond_mean <- apply(qini_cond, 2, mean)
    upper_cond <- apply(qini_cond, 2, quantile, p = 0.975)
    lower_cond <- apply(qini_cond, 2, quantile, p = 0.025)
    saveRDS(data.frame(qini_orig_mean = qini_orig_mean, upper_orig = upper_orig, lower_orig = lower_orig,
                       qini_cond_mean = qini_cond_mean, upper_cond = upper_cond, lower_cond = lower_cond), paste0(path_sim_data, '/ss_red_qini.RDS'))
    
  }
}

plot_qini <- function(p_dgp , p_noise){
  
  ## get results
  ss_red_qini <- readRDS(paste0('simulation_results/', p_dgp, '/noise_', p_noise, '/ss_red_qini.RDS'))
  
  ## plot qini curve
  ss_red_qini %$% plot((1:10)/10*2500, qini_orig_mean, ylim = c(min(lower_orig, lower_cond), max(upper_orig, upper_cond)), type = 'l',
       xlab = '', ylab = '', lwd = 2, yaxt='n', xaxt='n')
  ss_red_qini %$% lines((1:10)/10*2500, lower_orig, lty = 3, lwd = 2)
  ss_red_qini %$% lines((1:10)/10*2500, upper_orig, lty = 3, lwd = 2)
  ss_red_qini %$% lines((1:10)/10*2500, qini_cond_mean, lty = 1, lwd = 2, col = 'green')
  ss_red_qini %$% lines((1:10)/10*2500, lower_cond, lty = 3, lwd = 2, col = 'green')
  ss_red_qini %$% lines((1:10)/10*2500, upper_cond, lty = 3, lwd = 2, col = 'green')
  
}

### Plot Qini results
layout(matrix(c(1,1,1,2,2,2,3,4,5,6,6,6,7,8,9,10,10,10), ncol=3, byrow=TRUE), heights = c(0.03,0.025,0.15,0.015,0.15,0.04))

## title
par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5,0.5,"Test set size reduction Qini",cex=2,font=1)

## aw 
par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5,0.5,"Setting aw",cex=1.5,font=1.5)
par(mar=c(1.2, 0.4, 0, 0.4), mgp=c(3, 0.3, 0))
plot_qini(p_dgp = 'aw2', p_noise = 0.5)
plot_qini(p_dgp = 'aw2', p_noise = 1)
plot_qini(p_dgp = 'aw2', p_noise = 2)

## nw 
par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5,0.5,"Setting nw",cex=1.5,font=1.5)
par(mar=c(1.2, 0.4, 0, 0.4), mgp=c(3, 0.3, 0))
plot_qini(p_dgp = 'nw2', p_noise = 0.5)
plot_qini(p_dgp = 'nw2', p_noise = 1)
plot_qini(p_dgp = 'nw2', p_noise = 2)

##legend
par(mar=c(1, 1, 0.5, 0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c('black', 'green')
legend(x = "bottom",inset = 0,
       legend = c('orig', 'cond'), 
       col=plot_colors, lwd=5, cex=1, horiz = TRUE)


################################ doubly-robust adjustment ###############################################
for(i_dgp in 1:length(dgp_list)){
  
  for(i_noise in 1:length(sigma_noise_list)){
    
    path_sim_data <- paste0('simulation_results/', dgp_list[i_dgp], '/noise_', sigma_noise_list[i_noise])
    var_red_qini_sel <- var_red_qini[row.names(var_red_qini) == paste0(dgp_list[i_dgp], sigma_noise_list[i_noise]), 'dr_min']
    qini_orig <- NULL
    qini_dr <- NULL
    for(i_list_data in 1:10){
      
      list_data <- readRDS(paste0(path_sim_data, '/list_data', i_list_data, '.RDS'))
      
      for(i_sim in 1:length(list_data)){
        
        ## get reduced data set
        dt_test <- list_data[[i_sim]]
        dt_red <- dt_test[1:ceiling((1-var_red_qini_sel/100)*nrow(dt_test)),]
        
        qini_curve <- c()
        qini_curve_dr <- c()
        for(i_dec in 9:0){
          
          ## Qini curve and standard deviation with original target 
          res_qini <- dt_test %$% calc_CATE(p_Y = Y, p_W = W, p_pred = pred_cf, p_cut = i_dec/10)
          qini_curve <- c(qini_curve, res_qini$qini_mean)
          
          ## Qini curve and standard deviation with doubly-robust adjusted target
          res_qini_dr <- dt_red %$% calc_CATE(p_Y = Y_dr, p_W = W, p_pred = pred_cf, p_cut = i_dec/10)
          res_qini_dr$qini_mean <- res_qini_dr$qini_mean * nrow(dt_test)/nrow(dt_red)
          qini_curve_dr <- c(qini_curve_dr, res_qini_dr$qini_mean)
        }
        
        ## add results per decile to previous simulation runs
        qini_orig %<>% rbind(data.frame(first = qini_curve[1], second = qini_curve[2], third = qini_curve[3], fourth = qini_curve[4], fifth = qini_curve[5],
                                        sixth = qini_curve[6], seventh = qini_curve[7], eigth = qini_curve[8], nineth = qini_curve[9], tenth = qini_curve[10]))
        qini_dr %<>% rbind(data.frame(first = qini_curve_dr[1], second = qini_curve_dr[2], third = qini_curve_dr[3], fourth = qini_curve_dr[4], fifth = qini_curve_dr[5],
                                        sixth = qini_curve_dr[6], seventh = qini_curve_dr[7], eigth = qini_curve_dr[8], nineth = qini_curve_dr[9], tenth = qini_curve_dr[10]))
      }
      
      
    }
    qini_orig_mean <- apply(qini_orig, 2, mean)
    upper_orig <- apply(qini_orig, 2, quantile, p = 0.975)
    lower_orig <- apply(qini_orig, 2, quantile, p = 0.025)
    qini_dr_mean <- apply(qini_dr, 2, mean)
    upper_dr <- apply(qini_dr, 2, quantile, p = 0.975)
    lower_dr <- apply(qini_dr, 2, quantile, p = 0.025)
    saveRDS(data.frame(qini_orig_mean = qini_orig_mean, upper_orig = upper_orig, lower_orig = lower_orig,
                       qini_dr_mean = qini_dr_mean, upper_dr = upper_dr, lower_dr = lower_dr), paste0(path_sim_data, '/ss_red_qini.RDS'))
    
  }
}

plot_qini <- function(p_dgp , p_noise){
  
  ## get results
  ss_red_qini <- readRDS(paste0('simulation_results/', p_dgp, '/noise_', p_noise, '/ss_red_qini.RDS'))
  
  ## plot qini curve
  ss_red_qini %$% plot((1:10)/10*2500, qini_orig_mean, ylim = c(min(lower_orig, lower_dr), max(upper_orig, upper_dr)), type = 'l',
                       xlab = '', ylab = '', lwd = 2, yaxt='n', xaxt='n')
  ss_red_qini %$% lines((1:10)/10*2500, lower_orig, lty = 3, lwd = 2)
  ss_red_qini %$% lines((1:10)/10*2500, upper_orig, lty = 3, lwd = 2)
  ss_red_qini %$% lines((1:10)/10*2500, qini_dr_mean, lty = 1, lwd = 2, col = 'green')
  ss_red_qini %$% lines((1:10)/10*2500, lower_dr, lty = 3, lwd = 2, col = 'green')
  ss_red_qini %$% lines((1:10)/10*2500, upper_dr, lty = 3, lwd = 2, col = 'green')
  
}

### Plot Qini results
layout(matrix(c(1,1,1,2,2,2,3,4,5,6,6,6,7,8,9,10,10,10), ncol=3, byrow=TRUE), heights = c(0.03,0.025,0.15,0.015,0.15,0.04))

## title
par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5,0.5,"Test set size reduction Qini",cex=2,font=1)

## aw 
par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5,0.5,"Setting aw",cex=1.5,font=1.5)
par(mar=c(1.2, 0.4, 0, 0.4), mgp=c(3, 0.3, 0))
plot_qini(p_dgp = 'aw2', p_noise = 0.5)
plot_qini(p_dgp = 'aw2', p_noise = 1)
plot_qini(p_dgp = 'aw2', p_noise = 2)

## nw 
par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5,0.5,"Setting nw",cex=1.5,font=1.5)
par(mar=c(1.2, 0.4, 0, 0.4), mgp=c(3, 0.3, 0))
plot_qini(p_dgp = 'nw2', p_noise = 0.5)
plot_qini(p_dgp = 'nw2', p_noise = 1)
plot_qini(p_dgp = 'nw2', p_noise = 2)

##legend
par(mar=c(1, 1, 0.5, 0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c('black', 'green')
legend(x = "bottom",inset = 0,
       legend = c('orig', 'dr'), 
       col=plot_colors, lwd=5, cex=1, horiz = TRUE)









