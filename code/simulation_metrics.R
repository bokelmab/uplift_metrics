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

## required functions
source('code/perform_simulation.R')
source('code/helper/Qini_calculation_and_plot.R')

## simulation parameters
n_train <- 10000
n_test <- 5000
n_sim <- 10000
dgp_list <- c("aw2", "nw2")
sigma_noise_list <- c(0.5, 1, 2)
seed <- 7323

for(i_dgp in 1:length(dgp_list)){
  
  dgp <- dgp_list[i_dgp]
  path_dgp <- paste0('simulation_results/', dgp)
  dir.create(path_dgp)
  
  for(i_sigma in 1:length(sigma_noise_list)){
    
    sigma_noise <- sigma_noise_list[i_sigma]
    path_noise <- paste0(path_dgp, '/noise_' , sigma_noise)
    dir.create(path_noise)
    results_simulation <- perform_simulation(p_dgp = dgp, p_sigma = sigma_noise, p_ntrain = n_train, p_ntest = n_test, p_nsim = n_sim, 
                                             p_seed = seed, p_path_noise = path_noise)
    
    ## save results
    saveRDS(results_simulation, paste0(path_noise, '/simulation_results.RDS'))
    
  }
  
}




