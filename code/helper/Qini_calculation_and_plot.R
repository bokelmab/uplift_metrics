## calculate Qini curve
calc_CATE <- function(p_Y, p_W, p_pred, p_cut){
  y_treat <- p_Y[p_pred > quantile(p_pred, p_cut) & p_W == 1]
  y_contr <- p_Y[p_pred > quantile(p_pred, p_cut) & p_W == 0]
  qini_mean <- (mean(y_treat)-mean(y_contr))*length(y_treat)
  qini_sd <- sqrt(var(y_treat)/length(y_treat)+var(y_contr)/length(y_contr))*length(y_treat)
  return(list(qini_mean = qini_mean, qini_sd = qini_sd))
}

plot_qini_curve <- function(p_preds, p_dt_test){
  
  qini_curve <- c()
  sd_curve <- c()
  qini_curve_dr <- c()
  sd_curve_dr <- c()
  qini_curve_cond <- c()
  sd_curve_cond <- c()
  for(i_dec in 9:0){
    
    ## Qini curve and standard deviation with original target 
    res_qini <- p_dt_test %$% calc_CATE(p_Y = Y, p_W = W, p_pred = p_preds, p_cut = i_dec/10)
    qini_curve <- c(qini_curve, res_qini$qini_mean)
    sd_curve <- c(sd_curve, res_qini$qini_sd)
    
    ## Qini curve and standard deviation with dr-adjusted target
    res_qini_dr <- p_dt_test %$% calc_CATE(p_Y = Y_dr, p_W = W, p_pred = p_preds, p_cut = i_dec/10)
    qini_curve_dr <- c(qini_curve_dr, res_qini_dr$qini_mean)
    sd_curve_dr <- c(sd_curve_dr, res_qini_dr$qini_sd)
    
    ## Qini curve and standard deviation with conditional mean adjusted target
    res_qini_cond <- p_dt_test %$% calc_CATE(p_Y = Y_cond, p_W = W, p_pred = p_preds, p_cut = i_dec/10)
    qini_curve_cond <- c(qini_curve_cond, res_qini_cond$qini_mean)
    sd_curve_cond <- c(sd_curve_cond, res_qini_cond$qini_sd)
  }
  
  ## plot Qini curve with CI
  layout(matrix(c(1,2,3,4), ncol=2, byrow=TRUE), heights = c(0.9, 0.1))
  par(mar=c(3, 3, 3, 1), mgp=c(1.5, 0.4, 0))
  upper <- qini_curve + 1.96*sd_curve
  lower <- qini_curve - 1.96*sd_curve
  upper_dr <- qini_curve_dr + 1.96*sd_curve_dr
  lower_dr <- qini_curve_dr - 1.96*sd_curve_dr
  upper_cond <- qini_curve_cond + 1.96*sd_curve_cond
  lower_cond <- qini_curve_cond - 1.96*sd_curve_cond
  plot((0:10)/10 * sum(p_dt_test$W), c(0, qini_curve), type = 'l', main = 'Qini curve', ylim = c(min(c(lower,0)), max(upper)), ylab = 'incremental gain', xlab = expression('N'[W]*'(s)'), lwd = 3)
  lines((0:10)/10 * sum(p_dt_test$W), c(0,lower), lty = 3, col = 'black', lwd = 3)
  lines((0:10)/10 * sum(p_dt_test$W), c(0,upper), lty = 3, col = 'black', lwd = 3)
  lines((0:10)/10 * sum(p_dt_test$W), c(0, qini_curve_dr), col = 'orange', lwd = 3)
  lines((0:10)/10 * sum(p_dt_test$W), c(0,lower_dr), lty = 3, col = 'orange', lwd = 3)
  lines((0:10)/10 * sum(p_dt_test$W), c(0,upper_dr), lty = 3, col = 'orange', lwd = 3)
  lines((0:10)/10 * sum(p_dt_test$W), c(0, qini_curve_cond), col = 'green', lwd = 3)
  lines((0:10)/10 * sum(p_dt_test$W), c(0,lower_cond), lty = 3, col = 'green', lwd = 3)
  lines((0:10)/10 * sum(p_dt_test$W), c(0,upper_cond), lty = 3, col = 'green', lwd = 3)
  lines(c(0,1) * sum(p_dt_test$W),c(0, qini_curve[length(qini_curve)]), col = 'black', lty = 2, lwd = 3)
  lines(c(0,1) * sum(p_dt_test$W),c(0, qini_curve_dr[length(qini_curve_dr)]), col = 'orange', lty = 2, lwd = 3)
  lines(c(0,1) * sum(p_dt_test$W),c(0, qini_curve_cond[length(qini_curve_cond)]), col = 'green', lty = 2, lwd = 3)
  
  ## variance reduction per decile
  perc_red_var_dr <- (1-c(1, (sd_curve_dr/sd_curve)^2))*100
  perc_red_var_cond <- (1-c(1, (sd_curve_cond/sd_curve)^2))*100
  par(mar=c(3, 2, 3, 2))
  plot((0:10)/10 * sum(p_dt_test$W), perc_red_var_dr, type = 'l', lwd = 3, main = '% reduction of variance', xlab = expression('N'[W]*'(s)'), ylab = '', col = 'orange', 
       ylim = c(min(c(perc_red_var_dr, perc_red_var_cond)), max(c(perc_red_var_dr, perc_red_var_cond))))
  lines((0:10)/10 * sum(p_dt_test$W), perc_red_var_cond, type = 'l', lwd = 3, xlab = expression('N'[W]*'(s)'), ylab = '', col = 'green')
  res_var_red <- rbind(perc_red_var_dr, perc_red_var_cond)
  row.names(res_var_red) <- c('dr', 'cond')
  print(res_var_red)
  
  ##legend
  par(mar=c(1, 1, 0.5, 0))
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  plot_colors <- c('black', 'green', 'orange')
  legend(x = "bottom",inset = 0,
         legend = c('orig', 'cond', 'dr'), 
         col=plot_colors, lty = c(1,1,1), lwd=5, cex=1, horiz = TRUE, seg.len = 0.5, text.width = 0.07, bty = 'n')
  plot.new()
  
  
}