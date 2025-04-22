#################################################
#################################################
#################################################
##
## Important function: FullSimulation()
##
#################################################
#################################################
#################################################

full_simulation <- function(id_task, run_id, run_type) {
  
  wd <- getwd()
  writepath <- ifelse(run_type == 1, "/mvt5_comp/txt_exps/","/mvt5_comp/txt_full/")
  nsim      <- ifelse(run_type == 1, 2, 20)

  for(ind_sim in 0:(nsim - 1)) {
    args <- create_parameters(id_task)
    p <- args$p
    sp <- args$sp
    model <- args$model
    nnn <- args$nnn

    split_vals <- 1
    n1 <-  split_vals * nnn
    n2 <-  (2 - split_vals) * nnn

    sigma0 <- do.call(paste0("m", as.character(model)), list(p))
    sigma_list <- mat_pair(sigma0, p, sp)

    sigma1 <- sigma_list[[1]]
    sigma2 <- sigma_list[[2]]


    ## Simulation with T-distribution
    x1 <- rmvt(n1 + 1, sigma = sigma1, df = 5)
    x2 <- rmvt(n2 + 1, sigma = sigma2, df = 5)

    x1 <- x1[1:n1, ]
    x2 <- x2[1:n2, ]

    ## Covariance testing:
    cov_LiChen2012_pv <- equalCovs(x1, x2, n1, n2)[2]
    cov_CaiLiuXia2013_pv <- cov_CaiLiuXia2013_fun(x1, x2)[2]
    cov_ZhengEtAl2020_pv <- cov_ZhengEtAl2020_fun(x1, x2)
    
    ## Correlation testing:
    cor_CaiZhang2016_pv <- cor_CaiZhang2016_fun(x1, x2)[2]
    cor_ZhengEtAl2019_pv <- cor_ZhengEtAl2019_fun(x1, x2)[c(2, 4)]
    
    ## Our proposals:
    zhang_pv <- zhang_fun(x1, x2, 500)

    cov_zhang_pvB1000 <- c(zhang_pv[4:9])#, zhang_trunc_pv[-(1:2)])
    cor_zhang_pvB1000 <- c(zhang_pv[1:3])#, zhang_trunc_pv[1:2])
    
    cov_zhang_pvB200 <- c(zhang_pv[9 + 4:9])#, zhang_trunc_pv[-(1:2)])
    cor_zhang_pvB200 <- c(zhang_pv[9 + 1:3])#, zhang_trunc_pv[1:2])

    cov_zhang_pvB50 <- c(zhang_pv[18 + 4:9])#, zhang_trunc_pv[-(1:2)])
    cor_zhang_pvB50 <- c(zhang_pv[18 + 1:3])#, zhang_trunc_pv[1:2])
    

    mre <- c(
      n1, #1
      n2, #2
      p, #3
      model, #4
      sp, #5
      max(abs(sigma1 - sigma2)), #6
      sum(abs(sigma1 - sigma2)), #7
      sum((sigma1 - sigma2)^2), #8
      sum((sigma1 - sigma2) != 0), #9
      max(abs(cov2cor(sigma1) - cov2cor(sigma2))), #10
      sum(abs(cov2cor(sigma1) - cov2cor(sigma2))), #11
      sum((cov2cor(sigma1) - cov2cor(sigma2))^2), #12
      sum((cov2cor(sigma1) - cov2cor(sigma2)) != 0), #13
      cov_LiChen2012_pv, #14
      cov_CaiLiuXia2013_pv, #15
      cov_ZhengEtAl2020_pv, #16
      cov_zhang_pvB1000, # Without thresholding: 17,18,19,20,21,22
      cov_zhang_pvB200, # Without thresholding: 23,24,25,26,27,28
      cov_zhang_pvB50, # Without thresholding: 29,30,31,32,33,34
      cor_CaiZhang2016_pv, #35
      cor_ZhengEtAl2019_pv, #36, 37
      cor_zhang_pvB1000, # Without thresholding: 38, 39, 40
      cor_zhang_pvB200, # Without thresholding: 41, 42, 43
      cor_zhang_pvB50 # Without thresholding: 44, 45, 46
    )
    mre <- round(mre, 3)
    print(length(mre))

    print(mre)
    if (sp == 0) {
      filename <- paste0(
        id_task, "_",
        nsim * run_id + ind_sim, "_m", model,
        "_n1_", n1,
        "_n2_", n2,
        "_p", p,
        "sizetxt")
    } else {
      filename <- paste0(
        id_task, "_",
        nsim * run_id + ind_sim, "_m", model,
        "_n1_", n1,
        "_n2_", n2,
        "_p", p,
        "powertxt")
    }
    print(paste(
      "The file name is:",
      filename))
    write.table(
      t(mre),
      file = paste0(wd, writepath, filename),
      row.names = FALSE,
      col.names = FALSE)
  }
  
  return(1)

}