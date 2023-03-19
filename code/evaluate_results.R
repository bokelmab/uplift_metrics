## distribution of CATE10% estimation
plot_dist_TEest <- function(p_dgp , p_noise, p_main = NULL){
  
  ## get results
  results_simulation <- readRDS(paste0('simulation_results/', p_dgp, '/noise_', p_noise, '/simulation_results.RDS'))
  
  ## plot Qini curve
  data_orig <- data.table(first = results_simulation$qini_orig$first, second = results_simulation$qini_orig$second, third = results_simulation$qini_orig$third,
                          fourth = results_simulation$qini_orig$fourth, fifth = results_simulation$qini_orig$fifth, sixth = results_simulation$qini_orig$sixth,
                          seventh = results_simulation$qini_orig$seventh, eith = results_simulation$qini_orig$eith, nineth = results_simulation$qini_orig$nineth,
                          tenth = results_simulation$qini_orig$tenth)
  data_cond <- data.table(first = results_simulation$qini_cond$first, second = results_simulation$qini_cond$second, third = results_simulation$qini_cond$third,
                          fourth = results_simulation$qini_cond$fourth, fifth = results_simulation$qini_cond$fifth, sixth = results_simulation$qini_cond$sixth,
                          seventh = results_simulation$qini_cond$seventh, eith = results_simulation$qini_cond$eith, nineth = results_simulation$qini_cond$nineth,
                          tenth = results_simulation$qini_cond$tenth)
  data_dr <- data.table(first = results_simulation$qini_dr$first, second = results_simulation$qini_dr$second, third = results_simulation$qini_dr$third,
                          fourth = results_simulation$qini_dr$fourth, fifth = results_simulation$qini_dr$fifth, sixth = results_simulation$qini_dr$sixth,
                          seventh = results_simulation$qini_dr$seventh, eith = results_simulation$qini_dr$eith, nineth = results_simulation$qini_dr$nineth,
                          tenth = results_simulation$qini_dr$tenth)
  qini_orig <- apply(data_orig, 2, mean)
  upper_orig <- apply(data_orig, 2, quantile, p = 0.975)
  lower_orig <- apply(data_orig, 2, quantile, p = 0.025)
  qini_cond <- apply(data_cond, 2, mean)
  upper_cond <- apply(data_cond, 2, quantile, p = 0.975)
  lower_cond <- apply(data_cond, 2, quantile, p = 0.025)
  qini_dr <- apply(data_dr, 2, mean)
  upper_dr <- apply(data_dr, 2, quantile, p = 0.975)
  lower_dr <- apply(data_dr, 2, quantile, p = 0.025)
  plot((1:10)/10*2500, qini_orig, ylim = c(min(lower_orig, lower_cond, lower_dr), max(upper_orig, upper_cond, upper_dr)), type = 'l',
       xlab = '', ylab = '', lwd = 2, yaxt='n')
  lines((1:10)/10*2500, lower_orig, lty = 3, lwd = 2)
  lines((1:10)/10*2500, upper_orig, lty = 3, lwd = 2)
  lines((1:10)/10*2500, qini_cond, lty = 1, lwd = 2, col = 'green')
  lines((1:10)/10*2500, lower_cond, lty = 3, lwd = 2, col = 'green')
  lines((1:10)/10*2500, upper_cond, lty = 3, lwd = 2, col = 'green')
  lines((1:10)/10*2500, qini_dr, lty = 1, lwd = 2, col = 'orange')
  lines((1:10)/10*2500, lower_dr, lty = 3, lwd = 2, col = 'orange')
  lines((1:10)/10*2500, upper_dr, lty = 3, lwd = 2, col = 'orange')
  
  ## calculate mean variance reduction
  sd_curve <- apply(data_orig, 2, sd)
  sd_curve_dr <- apply(data_dr, 2, sd)
  sd_curve_cond <- apply(data_cond, 2, sd)
  perc_red_var_dr <- (1-(sd_curve_dr/sd_curve)^2)*100
  perc_red_var_cond <- (1-(sd_curve_cond/sd_curve)^2)*100
  
  ## provide standard errors
  return(data.frame(cond_mean = round(mean(perc_red_var_cond),3), cond_min = round(min(perc_red_var_cond),3), cond_max = round(max(perc_red_var_cond),3),
             dr_mean = round(mean(perc_red_var_dr),3), dr_min = round(min(perc_red_var_dr),3), dr_max = round(max(perc_red_var_dr),3),
             row.names = paste0(p_dgp, p_noise)))
}

### Plot Qini results
var_red_qini <- NULL
layout(matrix(c(1,1,1,2,3,4,5,5,5,6,7,8,9,9,9), ncol=3, byrow=TRUE), heights = c(0.025,0.15,0.025,0.15,0.06))

## aw 
par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5,0.5,"Setting aw",cex=1.5,font=1.5)
par(mar=c(1.2, 0.4, 0, 0.4), mgp=c(3, 0.3, 0))
var_red_qini %<>% rbind(plot_dist_TEest(p_dgp = 'aw2', p_noise = 0.5, 'yes'))
var_red_qini %<>% rbind(plot_dist_TEest(p_dgp = 'aw2', p_noise = 1, 'yes'))
var_red_qini %<>% rbind(plot_dist_TEest(p_dgp = 'aw2', p_noise = 2, 'yes'))

## nw 
par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5,0.5,"Setting nw",cex=1.5,font=1.5)
par(mar=c(1.2, 0.4, 0, 0.4), mgp=c(3, 0.3, 0))
var_red_qini %<>% rbind(plot_dist_TEest(p_dgp = 'nw2', p_noise = 0.5, 'yes'))
var_red_qini %<>% rbind(plot_dist_TEest(p_dgp = 'nw2', p_noise = 1, 'yes'))
var_red_qini %<>% rbind(plot_dist_TEest(p_dgp = 'nw2', p_noise = 2, 'yes'))

##legend
par(mar=c(1, 1, 0.5, 0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c('black', 'green', 'orange')
legend(x = "bottom",inset = 0,
       legend = c('orig', 'cond', 'dr'), 
       col=plot_colors, lwd=5, cex=1, horiz = TRUE)

## distr
plot_mse_metric <- function(p_comp_name, p_dgp , p_noise, p_main = NULL){
  
  ## get results
  results_simulation <- readRDS(paste0('simulation_results/', p_dgp, '/noise_', p_noise, '/simulation_results.RDS'))
  
  mse_orig <- results_simulation[['mse_orig']]
  mse_uc <- results_simulation[['mse_uc']]
  mse_cond <- results_simulation[['mse_cond']]
  mse_true <- results_simulation[['mse_true']]
  mse_dr <- results_simulation[['mse_dr']]
  mse_comp_orig <- results_simulation[[paste0('mse_orig_', p_comp_name)]]
  mse_comp_uc <- results_simulation[[paste0('mse_uc_', p_comp_name)]]
  mse_comp_cond <- results_simulation[[paste0('mse_cond_', p_comp_name)]]
  mse_comp_true <- results_simulation[[paste0('mse_true_', p_comp_name)]]
  mse_comp_dr <- results_simulation[[paste0('mse_dr_', p_comp_name)]]
  dens_mse_orig <- density(mse_orig-mse_comp_orig)
  dens_mse_uc <- density(mse_uc-mse_comp_uc)
  dens_mse_cond <- density(mse_cond-mse_comp_cond)
  dens_mse_dr <- density(mse_dr-mse_comp_dr)
  
  main <- ifelse(is.null(p_main), paste0('MSE model vs. ', p_comp_name, ' (', dgp, ', sig=', sigma_noise, ', var_exp=', round(mean(results_simulation$var_exp),3),')'),
                 paste0('var_exp=', round(mean(results_simulation$var_exp),3)))
  plot(dens_mse_orig, main = '', 
       ylim = c(min(c(dens_mse_orig$y,dens_mse_cond$y, dens_mse_dr$y)), max(c(dens_mse_orig$y,dens_mse_cond$y, dens_mse_dr$y))), 
       xlab = '', ylab = '', lwd = 2, yaxt='n')
  lines(dens_mse_cond, col = 'green', lwd = 2)
  lines(dens_mse_uc, col = 'blue', lwd = 2)
  lines(dens_mse_dr, col = 'orange', lwd = 2)
  abline(v = mean(mse_true-mse_comp_true), col = 'red', lwd = 2)
  abline(v = 0, col = 'black', lty = 3, lwd = 2)
  
  ## provide standard errors
  return(data.frame(orig = round(mean(mse_orig > mse_comp_orig),3), cond = round(mean(mse_cond > mse_comp_cond),3),
                    uncond = round(mean(mse_uc > mse_comp_uc),3), dr = round(mean(mse_dr > mse_comp_dr),3), row.names = paste0(p_dgp, '_', p_comp_name, '_',p_noise)))
  
}

##### MSE graphic for aw2 #########################
res_MSE <- NULL
layout(matrix(c(1,1,1,2,3,4,5,5,5,6,7,8,9,9,9,10,11,12,13,13,13), ncol=3, byrow=TRUE), heights = c(0.025,0.15,0.025,0.15,0.025,0.15,0.06))

## optimal estimator
par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5,0.5,"CF against optimal",cex=1.5,font=1.5)
par(mar=c(1.2, 0.4, 0, 0.4), mgp=c(3, 0.3, 0))
res_MSE %<>% rbind(plot_mse_metric('opt', p_dgp = 'aw2', p_noise = 0.5, 'yes'))
par(mar=c(1.2, 0.4, 0, 0.4), mgp=c(3, 0.3, 0))
res_MSE %<>% rbind(plot_mse_metric('opt', p_dgp = 'aw2', p_noise = 1, 'yes'))
par(mar=c(1.2, 0.4, 0, 0.4), mgp=c(3, 0.3, 0))
res_MSE %<>% rbind(plot_mse_metric('opt', p_dgp = 'aw2', p_noise = 2, 'yes'))

## trivial estimator
par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5,0.5,"CF against trivial",cex=1.5,font=1.5)
par(mar=c(1.2, 0.4, 0, 0.4), mgp=c(3, 0.3, 0))
res_MSE %<>% rbind(plot_mse_metric('triv', p_dgp = 'aw2', p_noise = 0.5, 'yes'))
par(mar=c(1.2, 0.4, 0, 0.4), mgp=c(3, 0.3, 0))
res_MSE %<>% rbind(plot_mse_metric('triv', p_dgp = 'aw2', p_noise = 1, 'yes'))
par(mar=c(1.2, 0.4, 0, 0.4), mgp=c(3, 0.3, 0))
res_MSE %<>% rbind(plot_mse_metric('triv', p_dgp = 'aw2', p_noise = 2, 'yes'))

## worse estimator
## trivial estimator
par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5,0.5,"CF against worse",cex=1.5,font=1.5)
par(mar=c(1.2, 0.4, 0, 0.4), mgp=c(3, 0.3, 0))
res_MSE %<>% rbind(plot_mse_metric('worse', p_dgp = 'aw2', p_noise = 0.5, 'yes'))
par(mar=c(1.2, 0.4, 0, 0.4), mgp=c(3, 0.3, 0))
res_MSE %<>% rbind(plot_mse_metric('worse', p_dgp = 'aw2', p_noise = 1, 'yes'))
par(mar=c(1.2, 0.4, 0, 0.4), mgp=c(3, 0.3, 0))
res_MSE %<>% rbind(plot_mse_metric('worse', p_dgp = 'aw2', p_noise = 2, 'yes'))

##legend
par(mar=c(1, 1, 0.5, 0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c('black', 'green', 'blue', 'orange')
legend(x = "bottom",inset = 0,
       legend = c("orig", "cond", "uncond","dr"), 
       col=plot_colors, lwd=5, cex=1, horiz = TRUE)

##### MSE graphic for nw2 #########################
layout(matrix(c(1,1,1,2,3,4,5,5,5,6,7,8,9,9,9,10,11,12,13,13,13), ncol=3, byrow=TRUE), heights = c(0.025,0.15,0.025,0.15,0.025,0.15,0.06))

## optimal estimator
par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5,0.5,"CF against optimal",cex=1.5,font=1.5)
par(mar=c(1.2, 0.4, 0, 0.4), mgp=c(3, 0.3, 0))
res_MSE %<>% rbind(plot_mse_metric('opt', p_dgp = 'nw2', p_noise = 0.5, 'yes'))
par(mar=c(1.2, 0.4, 0, 0.4), mgp=c(3, 0.3, 0))
res_MSE %<>% rbind(plot_mse_metric('opt', p_dgp = 'nw2', p_noise = 1, 'yes'))
par(mar=c(1.2, 0.4, 0, 0.4), mgp=c(3, 0.3, 0))
res_MSE %<>% rbind(plot_mse_metric('opt', p_dgp = 'nw2', p_noise = 2, 'yes'))

## trivial estimator
par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5,0.5,"CF against trivial",cex=1.5,font=1.5)
par(mar=c(1.2, 0.4, 0, 0.4), mgp=c(3, 0.3, 0))
res_MSE %<>% rbind(plot_mse_metric('triv', p_dgp = 'nw2', p_noise = 0.5, 'yes'))
par(mar=c(1.2, 0.4, 0, 0.4), mgp=c(3, 0.3, 0))
res_MSE %<>% rbind(plot_mse_metric('triv', p_dgp = 'nw2', p_noise = 1, 'yes'))
par(mar=c(1.2, 0.4, 0, 0.4), mgp=c(3, 0.3, 0))
res_MSE %<>% rbind(plot_mse_metric('triv', p_dgp = 'nw2', p_noise = 2, 'yes'))

## worse estimator
## trivial estimator
par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5,0.5,"CF against worse",cex=1.5,font=1.5)
par(mar=c(1.2, 0.4, 0, 0.4), mgp=c(3, 0.3, 0))
res_MSE %<>% rbind(plot_mse_metric('worse', p_dgp = 'nw2', p_noise = 0.5, 'yes'))
par(mar=c(1.2, 0.4, 0, 0.4), mgp=c(3, 0.3, 0))
res_MSE %<>% rbind(plot_mse_metric('worse', p_dgp = 'nw2', p_noise = 1, 'yes'))
par(mar=c(1.2, 0.4, 0, 0.4), mgp=c(3, 0.3, 0))
res_MSE %<>% rbind(plot_mse_metric('worse', p_dgp = 'nw2', p_noise = 2, 'yes'))

##legend
par(mar=c(1, 1, 0.5, 0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c('black', 'green', 'blue', 'orange')
legend(x = "bottom",inset = 0,
       legend = c("orig", "cond", "uncond","dr"), 
       col=plot_colors, lwd=5, cex=1, horiz = TRUE)

### variance reduction for TMSE difference between CF and triv
settings <- data.frame(setting = c('aw2', 'aw2', 'aw2', 'nw2', 'nw2', 'nw2'), sigma = c(0.5, 1, 2, 0.5, 1, 2)) 
var_red_TMSE <- NULL
for(i_set in 1:nrow(settings)){
  
  ## get results
  results_simulation <- readRDS(paste0('simulation_results/', settings$setting[i_set], '/noise_', settings$sigma[i_set], '/simulation_results.RDS'))
  mse_orig <- results_simulation[['mse_orig']]
  mse_uc <- results_simulation[['mse_uc']]
  mse_cond <- results_simulation[['mse_cond']]
  mse_true <- results_simulation[['mse_true']]
  mse_dr <- results_simulation[['mse_dr']]
  mse_comp_orig <- results_simulation[['mse_orig_triv']]
  mse_comp_uc <- results_simulation[['mse_uc_triv']]
  mse_comp_cond <- results_simulation[['mse_cond_triv']]
  mse_comp_true <- results_simulation[['mse_true_triv']]
  mse_comp_dr <- results_simulation[['mse_dr_triv']]
  
  var_red_TMSE %<>% rbind(data.frame(uc = 1-var(mse_uc-mse_comp_uc)/var(mse_orig-mse_comp_orig), 
                            cond = 1-var(mse_cond-mse_comp_cond)/var(mse_orig-mse_comp_orig), 
                            dr = 1-var(mse_dr-mse_comp_dr)/var(mse_orig-mse_comp_orig), row.names = paste0(settings$setting[i_set], settings$sigma[i_set])))
}
var_red_TMSE <- round(var_red_TMSE*100,1)


## prepare table with MSE results
## change to 1-P[false decision] where CF performs worse than compteting estimator
res_MSE[substr(row.names(res_MSE),5,7)=='opt', c('orig', 'cond', 'uncond', 'dr')] <- 1 - res_MSE[substr(row.names(res_MSE),5,7)=='opt', c('orig', 'cond', 'uncond', 'dr')]
res_MSE[row.names(res_MSE) == 'nw2_triv_2', c('orig', 'cond', 'uncond', 'dr')] <- 1 - res_MSE[row.names(res_MSE) == 'nw2_triv_2', c('orig', 'cond', 'uncond', 'dr')]
res_MSE <- res_MSE * 100
res_MSE <- res_MSE[, c("orig", "uncond", "cond", "dr")]

## save results
path_results <- 'simulation_results/summary_var_red'
dir.create(path_results)
saveRDS(var_red_qini, paste0(path_results, '/var_red_qini.RDS'))
saveRDS(var_red_TMSE, paste0(path_results, '/var_red_TMSE.RDS'))
