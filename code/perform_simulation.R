perform_simulation <- function(p_dgp, p_sigma_noise, p_ntrain, p_ntest, p_nsim, p_seed, p_path_noise){
  
  ## simulate data
  set.seed(p_seed)
  data <- generate_causal_data(p_ntrain, 6, dgp = p_dgp, sigma.noise = p_sigma_noise)
  dt_train <- data.table(Y = data$Y, X1 = data$X[,1], X2 = data$X[,2], X3 = data$X[,3], 
                         X4 = data$X[,4], X5 = data$X[,5], X6 = data$X[,6], W = data$W, tau = data$tau, 
                         m = data$m, e = data$e)
  
  ## train models
  cf <- causal_forest(X = dt_train[,c('X1', 'X2', 'X3', 'X4', 'X5', 'X6'), with = F], Y = dt_train$Y, W = dt_train$W,
                      num.trees=1000, sample.fraction=0.5,  
                      honesty=T, tune.parameters='none', honesty.fraction=0.5, ci.group.size=2, compute.oob.predictions=T,
                      seed = p_seed)
  rf <- ranger(Y ~ X1 + X2 + X3 + X4 + X5 + X6, data = dt_train, min.node.size = 20, seed = p_seed)
  rf0 <- ranger(Y ~ X1 + X2 + X3 + X4 + X5 + X6, data = dt_train[W == 0,], min.node.size = 20, seed = p_seed)
  rf1 <- ranger(Y ~ X1 + X2 + X3 + X4 + X5 + X6, data = dt_train[W == 1,], min.node.size = 20, seed = p_seed)
  
  ## evaluate according to MSE
  calc_mse <- function(p_Y, p_W, p_P, p_pred){
    return(mean((p_W*p_Y/p_P - (1-p_W)*p_Y/(1-p_P) - p_pred)^2))
  }
  
  ## simulate distribution of metric
  mse_orig <- c()
  mse_uc <- c()
  mse_cond <- c()
  mse_orig_triv <- c()
  mse_uc_triv <- c()
  mse_cond_triv <- c()
  mse_orig_opt <- c()
  mse_uc_opt <- c()
  mse_cond_opt <- c()
  mse_orig_worse <- c()
  mse_uc_worse <- c()
  mse_cond_worse <- c()
  mse_true <- c()
  mse_true_triv <- c()
  mse_true_opt <- c()
  mse_true_worse <- c()
  mse_dr <- c()
  mse_dr_triv <- c()
  mse_dr_opt <- c()
  mse_dr_worse <- c()
  qini_orig <- list(first = c(), second = c(), third = c(), fourth = c(), fifth = c(), sixth = c(), seventh = c(), 
                    eith = c(), nineth = c(), tenth = c())
  qini_cond <- list(first = c(), second = c(), third = c(), fourth = c(), fifth = c(), sixth = c(), seventh = c(), 
                    eith = c(), nineth = c(), tenth = c())
  qini_dr <- list(first = c(), second = c(), third = c(), fourth = c(), fifth = c(), sixth = c(), seventh = c(), 
                    eith = c(), nineth = c(), tenth = c())
  var_exp <- c()
  var_Wp_orig <- c()
  var_Wp_uc <- c()
  var_Wp_cond <- c()
  var_Wp_dr <- c()
  list_data <- list()
  set.seed(p_seed)
  for(i in 1:p_nsim){
    
    ## generate test data for iteration
    data <- generate_causal_data(p_ntest, 6, dgp = p_dgp, sigma.noise = p_sigma_noise) ## aw2, ai1, ai2, nw2
    dt_test <- data.table(Y = data$Y, X1 = data$X[,1], X2 = data$X[,2], X3 = data$X[,3], 
                          X4 = data$X[,4], X5 = data$X[,5], X6 = data$X[,6], W = data$W, tau = data$tau, 
                          m = data$m, e = data$e)
    
    ## make predictions
    dt_test$pred_cf <- predict(cf, newdata = dt_test[,c('X1', 'X2', 'X3', 'X4', 'X5', 'X6'), with = F])$predictions
    dt_test$pred_rf <- predict(rf, data = dt_test)$predictions
    dt_test$pred_dr <- (1-dt_test$e)*predict(rf1, data = dt_test)$predictions + dt_test$e*predict(rf0, data = dt_test)$predictions
    dt_test[, Y_cond := Y-pred_rf]
    dt_test[, Y_uc:= Y-mean(dt_train$Y)]
    dt_test[, Y_dr := Y-pred_dr]
    
    ## save data without features for further analysis (effect of sample size reduction)
    if(i %% 1000 == 0){
      list_data[[1000]] <- dt_test[, c('Y', 'W', 'e', 'Y_cond', 'Y_uc', 'Y_dr', 'pred_cf')]
      saveRDS(list_data, paste0(p_path_noise, '/list_data', i %/% 1000, '.RDS'))
      print('Saving 1000 test data sets')
      list_data <- list()
    }else{
      list_data[[i %% 1000]] <- dt_test[, c('Y', 'W', 'e', 'Y_cond', 'Y_uc', 'Y_dr', 'pred_cf')]
    }
    
    ## percentage of variance explained by RF
    var_exp <- c(var_exp, round(1-var(dt_test$Y_cond)/var(dt_test$Y), 2))
    
    ## Variance of transformed outcome
    dt_test[, Wp := (W-e)*(e*(1-e))]
    var_Wp_orig <- c(var_Wp_orig, var(dt_test$Wp*dt_test$Y))
    var_Wp_uc <- c(var_Wp_uc, var(dt_test$Wp*dt_test$Y_uc))
    var_Wp_cond <- c(var_Wp_cond, var(dt_test$Wp*dt_test$Y_cond))
    var_Wp_dr <- c(var_Wp_dr, var(dt_test$Wp*dt_test$Y_dr))
    
    ## calculate mse
    mse_orig <- c(mse_orig, calc_mse(p_Y = dt_test$Y, p_W = dt_test$W, p_P = dt_test$e, p_pred = dt_test$pred_cf))
    mse_uc <- c(mse_uc, calc_mse(p_Y = dt_test$Y_uc, p_W = dt_test$W, p_P = dt_test$e, p_pred = dt_test$pred_cf))
    mse_cond <- c(mse_cond, calc_mse(p_Y = dt_test$Y_cond, p_W = dt_test$W, p_P = dt_test$e, p_pred = dt_test$pred_cf))
    mse_true <- c(mse_true, mean((dt_test$tau-dt_test$pred_cf)^2))
    mse_dr <- c(mse_dr, calc_mse(p_Y = dt_test$Y_dr, p_W = dt_test$W, p_P = dt_test$e, p_pred = dt_test$pred_cf))
    
    ## trivial CATE prediction
    mse_orig_triv <- c(mse_orig_triv, calc_mse(p_Y = dt_test$Y, p_W = dt_test$W, p_P = dt_test$e, p_pred = 0))
    mse_uc_triv <- c(mse_uc_triv, calc_mse(p_Y = dt_test$Y_uc, p_W = dt_test$W, p_P = dt_test$e, p_pred = 0))
    mse_cond_triv <- c(mse_cond_triv, calc_mse(p_Y = dt_test$Y_cond, p_W = dt_test$W, p_P = dt_test$e, p_pred = 0))
    mse_true_triv <- c(mse_true_triv, mean((dt_test$tau-0)^2))
    mse_dr_triv <- c(mse_dr_triv, calc_mse(p_Y = dt_test$Y_dr, p_W = dt_test$W, p_P = dt_test$e, p_pred = 0))
    
    ## optimal CATE prediction
    mse_orig_opt <- c(mse_orig_opt, calc_mse(p_Y = dt_test$Y, p_W = dt_test$W, p_P = dt_test$e, p_pred = dt_test$tau))
    mse_uc_opt <- c(mse_uc_opt, calc_mse(p_Y = dt_test$Y_uc, p_W = dt_test$W, p_P = dt_test$e, p_pred = dt_test$tau))
    mse_cond_opt <- c(mse_cond_opt, calc_mse(p_Y = dt_test$Y_cond, p_W = dt_test$W, p_P = dt_test$e, p_pred = dt_test$tau))
    mse_true_opt <- c(mse_true_opt, mean((dt_test$tau-dt_test$tau)^2))
    mse_dr_opt <- c(mse_dr_opt, calc_mse(p_Y = dt_test$Y_dr, p_W = dt_test$W, p_P = dt_test$e, p_pred = dt_test$tau))
    
    ## worse CATE prediction
    pred_worse <- dt_test$pred_cf + rnorm(n = nrow(dt_test), sd = 0.1*sd(dt_test$pred_cf))
    mse_orig_worse <- c(mse_orig_worse, calc_mse(p_Y = dt_test$Y, p_W = dt_test$W, p_P = dt_test$e, p_pred = pred_worse))
    mse_uc_worse <- c(mse_uc_worse, calc_mse(p_Y = dt_test$Y_uc, p_W = dt_test$W, p_P = dt_test$e, p_pred = pred_worse))
    mse_cond_worse <- c(mse_cond_worse, calc_mse(p_Y = dt_test$Y_cond, p_W = dt_test$W, p_P = dt_test$e, p_pred = pred_worse))
    mse_true_worse <- c(mse_true_worse, mean((dt_test$tau-pred_worse)^2))
    mse_dr_worse <- c(mse_dr_worse, calc_mse(p_Y = dt_test$Y_dr, p_W = dt_test$W, p_P = dt_test$e, p_pred = pred_worse))
    
    ## Qini curve - CATE estimation at each decile
    for(i_dec in 9:0){
      
      ## Qini curve with original outcome 
      results_qini_orig <- dt_test %$% calc_CATE(p_Y = Y, p_W = W, p_pred = pred_cf, p_cut = i_dec/10)
      qini_orig[[10-i_dec]] <- c(qini_orig[[10-i_dec]], results_qini_orig$qini_mean)
      
      ## Qini curve with conditional mean adjusted outcome  
      results_qini_cond <- dt_test %$% calc_CATE(p_Y = Y_cond, p_W = W, p_pred = pred_cf, p_cut = i_dec/10)
      qini_cond[[10-i_dec]] <- c(qini_cond[[10-i_dec]], results_qini_cond$qini_mean)
      
      ## Qini curve with doubly-robust adjusted outcome  
      results_qini_dr <- dt_test %$% calc_CATE(p_Y = Y_dr, p_W = W, p_pred = pred_cf, p_cut = i_dec/10)
      qini_dr[[10-i_dec]] <- c(qini_dr[[10-i_dec]], results_qini_dr$qini_mean)
      
    }
    
    ## time
    if(i%% 100 == 0){
      print(i)
    }
  }
  
  ## return results
  return(list(var_exp = var_exp, mse_orig = mse_orig, mse_uc = mse_uc, mse_cond = mse_cond, mse_orig_triv = mse_orig_triv, mse_uc_triv = mse_uc_triv, mse_cond_triv = mse_cond_triv,
              mse_orig_opt = mse_orig_opt, mse_uc_opt = mse_uc_opt, mse_cond_opt = mse_cond_opt, mse_orig_worse = mse_orig_worse, mse_uc_worse = mse_uc_worse, mse_cond_worse = mse_cond_worse,
              mse_true = mse_true, mse_true_triv = mse_true_triv, mse_true_opt = mse_true_opt, mse_true_worse = mse_true_worse,
              mse_dr = mse_dr, mse_dr_triv = mse_dr_triv, mse_dr_opt = mse_dr_opt, mse_dr_worse = mse_dr_worse,
              qini_orig = qini_orig, qini_cond = qini_cond, qini_dr = qini_dr, var_Wp_orig = var_Wp_orig, var_Wp_uc = var_Wp_uc, var_Wp_cond = var_Wp_cond, var_Wp_dr = var_Wp_dr
              ))
  
}