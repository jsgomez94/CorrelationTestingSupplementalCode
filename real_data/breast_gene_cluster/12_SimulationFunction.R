#################################################
#################################################
#################################################
##
## Important function: FullSimulation()
##
#################################################
#################################################
#################################################

full_simulation <- function(run_id) {
  
  #################################################
  #################################################
  df_cleanskew <- read.csv(
  file = paste0("04_Breast_CleanGeneExpression.csv"),
  header = TRUE,
  row.names = 1)

  #################################################
  #################################################
  n_cleanskew <- dim(df_cleanskew)[1]
  
  ## p_cleanskew <- 20 ## Small trial before full simulation.
  p_cleanskew <- dim(df_cleanskew)[2]
  
  ## range_ns <- seq(from = 10, to = 30, 10) ## Small trial before full simulation.
  range_ns <- seq(from = 25, to = 300, 25)
  

  # pvals_ours                  <- matrix(0, ncol = length(range_ns), nrow = nsim)
  # pvals_cor_CaiZhang2016      <- matrix(0, ncol = length(range_ns), nrow = nsim)
  # pvals_cor_ZhengEtAl2019_max <- matrix(0, ncol = length(range_ns), nrow = nsim)
  # pvals_cor_ZhengEtAl2019_sum <- matrix(0, ncol = length(range_ns), nrow = nsim)

  pvals_ours                  <- rep(0, length(range_ns))
  pvals_cor_CaiZhang2016      <- rep(0, length(range_ns))
  pvals_cor_ZhengEtAl2019_max <- rep(0, length(range_ns))
  pvals_cor_ZhengEtAl2019_sum <- rep(0, length(range_ns))

  for (n_index in 1:length(range_ns)) {
    
    time_start <- Sys.time()
    
    n_val <- range_ns[n_index]
    n1 <- n_val
    n2 <- n_val
  
    subsample <- sample(1:n_cleanskew, size = n1 + n2, replace = FALSE)

    df_data1 <- as.matrix(df_cleanskew[subsample[1:n1], 1:p_cleanskew])
    df_data2 <- as.matrix(df_cleanskew[subsample[n1 + (1:n2)], 1:p_cleanskew])

    eta_output <- xeta_yeta_calculation_old(df_data1, df_data2)
  
    outputs_ours <- zhang_fun(df_data1, df_data2, 1000, NULL)
    pvals_ours[n_index] <- outputs_ours
  
    outputs_cor_CaiZhang2016 <- cor_CaiZhang2016_fun(df_data1, df_data2, eta_output, FALSE)
    pvals_cor_CaiZhang2016[n_index] <- outputs_cor_CaiZhang2016[2]
  
    outputs_cor_ZhengEtAl2019 <- cor_ZhengEtAl2019_fun(df_data1, df_data2, eta_output)
    pvals_cor_ZhengEtAl2019_max[n_index] <- outputs_cor_ZhengEtAl2019[4]
    pvals_cor_ZhengEtAl2019_sum[n_index] <- outputs_cor_ZhengEtAl2019[2]
    
    time_end <- Sys.time()
    time_simrun <- difftime(time1 = time_end, time2 = time_start, units = "m")
    print(paste0(" / n = ", n1 + n2, ": ", time_simrun))
      
    rm(df_data1, df_data2, eta_output) 
  }

  results_full <- NULL
  results_full <- rbind(
    pvals_cor_CaiZhang2016,
    pvals_cor_ZhengEtAl2019_sum,
    pvals_cor_ZhengEtAl2019_max,
    pvals_ours)
  results_full <- cbind(
    c(
        "Cai & Zhang (2016)",
        "Zheng et al. (2019) sum",
        "Zheng et al. (2019) max",
        "Permutation Test"),
    results_full)

  colnames(results_full) <- c("Method", paste0("n1 = ", range_ns))

  print(results_full)

  write.csv(x = results_full, file = paste0("result_data/result_data_",run_id,".csv"))

}

