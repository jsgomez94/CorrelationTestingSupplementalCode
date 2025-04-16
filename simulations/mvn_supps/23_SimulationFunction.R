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
  writepath <- ifelse(run_type == 1, "/mvn_supps/txt_exps/","/mvn_supps/txt_full/")
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

    x1 <- rmvnorm(n1 + 1, rep(0, p), sigma1)
    x2 <- rmvnorm(n2 + 1, rep(0, p), sigma2)

    x1 <- x1[1:n1, ]
    x2 <- x2[1:n2, ]

    ## Our proposals:
    cor_zhang_pv <- zhang_fun(x1, x2, 200)
    
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
      cor_zhang_pv
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