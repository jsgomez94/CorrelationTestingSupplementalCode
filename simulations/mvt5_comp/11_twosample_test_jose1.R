##############################################################################################################################################
######################### The Proposed Testing Procedure:
##############################################################################################################################################

zhang_fun <- function(X1, X2, B = 1000) {
  
  ###############################
  ###############################
  ## Variance-stabilizing transformation:
  H_fun <- function(m) {
    return(log((1 + m) / (1 - m)) / 2)
  }

  ###############################
  ###############################
  ## Value Initialization
  p <- ncol(X1)
  Ip <- diag(rep(1, p))
  Ipp <- rep(1, p) %*% t(rep(1, p))
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  n <- n1 + n2
  X <- rbind(X1, X2)
  
  cov.hat1 <- cov(X1) * (n1 - 1) / n1
  cov.hat2 <- cov(X2) * (n2 - 1) / n2
  
  cor.hat1 <- cor(X1)
  cor.hat2 <- cor(X2)
  
  sd.vec1 <- (diag(cov.hat1))^0.5
  sd.vec2 <- (diag(cov.hat2))^0.5
  var.vec1 <- (diag(cov.hat1))
  var.vec2 <- (diag(cov.hat2))
  var.mat1 <- sd.vec1 %*% t(sd.vec1)
  var.mat2 <- sd.vec2 %*% t(sd.vec2)
  
  ###############################
  ###############################
  ## Calculation of correlation statistics:
  
  w1 <- 1 / (abs(Ipp - (cor.hat1^2 - Ip)))
  w2 <- 1 / (abs(Ipp - (cor.hat2^2 - Ip)))
  transf1 <- H_fun(cor.hat1 - diag(diag(cor.hat1))) + diag(diag(cor.hat1))
  transf2 <- H_fun(cor.hat2 - diag(diag(cor.hat2))) + diag(diag(cor.hat2))

  cor.dif1 <- abs(cor.hat1 * w1  - cor.hat2 * w2)
  cor.dif2 <- abs(transf1 - transf2)
  cor.dif3 <- abs(cor.hat1  - cor.hat2)

  stat_q_cor1 <- sum(cor.dif1^2)
  stat_q_cor2 <- sum(cor.dif2^2)
  stat_q_cor3 <- sum(cor.dif3^2)
  
  ###############################
  ###############################
  ## Calculation of the covariance statistics:
  cov.dif1 <- abs(cor.hat1 * w1 * var.mat1 - cor.hat2 * w2 * var.mat2)
  cov.dif2 <- abs(cor.hat1 * w1 * var.mat1 / var.mat2  - cor.hat2 * w2 * var.mat2 / var.mat1)
  cov.dif3 <- abs(cor.hat1 * w1 * var.mat2 / var.mat1  - cor.hat2 * w2 * var.mat1 / var.mat2)
  
  cov.transf1 <- abs(transf1 * var.mat1 - transf2 * var.mat2)
  cov.transf2 <- abs(transf1 * var.mat1 / var.mat2  - transf2 * var.mat2 / var.mat1)
  cov.transf3 <- abs(transf1 * var.mat2 / var.mat1  - transf2 * var.mat1 / var.mat2)


  stat_dif_cov1 <- sum(cov.dif1^2)
  stat_dif_cov2 <- sum(cov.dif2^2)
  stat_dif_cov3 <- sum(cov.dif3^2)

  stat_transf_cov1 <- sum(cov.transf1^2)
  stat_transf_cov2 <- sum(cov.transf2^2)
  stat_transf_cov3 <- sum(cov.transf3^2)
  
  
  
  ###############################
  ###############################
  ## Resampling calculations:
  bt.stat_q_cor1 <- rep(0, B)
  bt.stat_q_cor2 <- rep(0, B)
  bt.stat_q_cor3 <- rep(0, B)
  
  bt.stat_dif_cov1 <- rep(0, B)
  bt.stat_dif_cov2 <- rep(0, B)
  bt.stat_dif_cov3 <- rep(0, B)

  bt.stat_transf_cov1 <- rep(0, B)
  bt.stat_transf_cov2 <- rep(0, B)
  bt.stat_transf_cov3 <- rep(0, B)
  
  for (j in 1:B) {
    b <- sample(1:n, n, replace = FALSE)
    b1 <- b[1:n1]
    b2 <- b[(n1 + 1):n]
    
    Xb1 <- X[b1, ]
    Xb2 <- X[b2, ]
    
    bt.cov.hat1 <- cov(Xb1) * (n1 - 1) / n1
    bt.cov.hat2 <- cov(Xb2) * (n2 - 1) / n2
    bt.cor.hat1 <- cor(Xb1)
    bt.cor.hat2 <- cor(Xb2)
    
    bt.sd.vec1 <- (diag(bt.cov.hat1))^0.5
    bt.sd.vec2 <- (diag(bt.cov.hat2))^0.5
    
    bt.var.vec1 <- diag(bt.cov.hat1)
    bt.var.vec2 <- diag(bt.cov.hat2)
    
    bt.var.mat1 <- bt.sd.vec1 %*% t(bt.sd.vec1)
    bt.var.mat2 <- bt.sd.vec2 %*% t(bt.sd.vec2)
    
    bt.w1 <- 1 / (abs(Ipp - (bt.cor.hat1^2 - Ip)))
    bt.w2 <- 1 / (abs(Ipp - (bt.cor.hat2^2 - Ip)))
    bt.transf1 <- H_fun(bt.cor.hat1 - diag(diag(bt.cor.hat1))) + diag(diag(bt.cor.hat1))
    bt.transf2 <- H_fun(bt.cor.hat2 - diag(diag(bt.cor.hat2))) + diag(diag(bt.cor.hat2))
    
    ###############################
    ###############################
    ## Calculation of resampling correlation statistics:
    bt.cor.dif1 <- abs(bt.cor.hat1 * bt.w1 - bt.cor.hat2 * bt.w2)
    bt.cor.dif2  <- abs(bt.transf1 - bt.transf2)
    bt.cor.dif3  <- abs(bt.cor.hat1 - bt.cor.hat2)

    ###############################
    ###############################
    ## Calculation of the covariance statistics:
    
    bt.cov.dif1 <- abs(bt.cor.hat1 * bt.w1 * bt.var.mat1 - bt.cor.hat2 * bt.w2 * bt.var.mat2)
    bt.cov.dif2 <- abs(bt.cor.hat1 * bt.w1 * bt.var.mat1 / bt.var.mat2  - bt.cor.hat2 * bt.w2 * bt.var.mat2 / bt.var.mat1)
    bt.cov.dif3 <- abs(bt.cor.hat1 * bt.w1 * bt.var.mat2 / bt.var.mat1  - bt.cor.hat2 * bt.w2 * bt.var.mat1 / bt.var.mat2)
  
    bt.cov.transf1 <- abs(bt.transf1 * bt.var.mat1 - bt.transf2 * bt.var.mat2)
    bt.cov.transf2 <- abs(bt.transf1 * bt.var.mat1 / bt.var.mat2  - bt.transf2 * bt.var.mat2 / bt.var.mat1)
    bt.cov.transf3 <- abs(bt.transf1 * bt.var.mat2 / bt.var.mat1  - bt.transf2 * bt.var.mat1 / bt.var.mat2)

    bt.stat_q_cor1[j] <- sum(bt.cor.dif1^2)
    bt.stat_q_cor2[j] <- sum(bt.cor.dif2^2)
    bt.stat_q_cor3[j] <- sum(bt.cor.dif3^2)
    
    bt.stat_dif_cov1[j] <- sum(bt.cov.dif1^2)
    bt.stat_dif_cov2[j] <- sum(bt.cov.dif2^2) 
    bt.stat_dif_cov3[j] <- sum(bt.cov.dif3^2)

    bt.stat_transf_cov1[j] <- sum(bt.cov.transf1^2)
    bt.stat_transf_cov2[j] <- sum(bt.cov.transf2^2)
    bt.stat_transf_cov3[j] <- sum(bt.cov.transf3^2)
  }
  
  ##################################################################
  ##################################################################
  ##################################################################
  ## Output guide for full B = 1000.
  ## pval1: cor test with original Yongli proposal.
  ## pval2: cor test with variance stabilizing transformation H.
  ##
  ## pval3 : cov 1 with original Yongli proposal + variance factor
  ## pval4 : cov with original Yongli proposal + var ratio
  ## pval5 : cov with original Yongli proposal + (var ratio)^{-1}
  ##
  ## pval6 : cov 1 with transformation H + variance factor
  ## pval7 : cov with transformation H + var ratio
  ## pval8 : cov with transformation H + (var ratio)^{-1}

  zhang_q_cor1_pv <- mean(stat_q_cor1 < bt.stat_q_cor1)
  zhang_q_cor2_pv <- mean(stat_q_cor2 < bt.stat_q_cor2)
  zhang_q_cor3_pv <- mean(stat_q_cor3 < bt.stat_q_cor3)

  zhang_dif_cov1_pv <- mean(stat_dif_cov1 < bt.stat_dif_cov1)
  zhang_dif_cov2_pv <- mean(stat_dif_cov2 < bt.stat_dif_cov2)
  zhang_dif_cov3_pv <- mean(stat_dif_cov3 < bt.stat_dif_cov3)

  zhang_transf_cov1_pv <- mean(stat_transf_cov1 < bt.stat_transf_cov1) 
  zhang_transf_cov2_pv <- mean(stat_transf_cov2 < bt.stat_transf_cov2) 
  zhang_transf_cov3_pv <- mean(stat_transf_cov3 < bt.stat_transf_cov3)
  
  ##################################################################
  ##################################################################
  ##################################################################
  ## Output guide for full B = 200.
  ## pval1: cor test with original Yongli proposal.
  ## pval2: cor test with variance stabilizing transformation H.
  ##
  ## pval3 : cov 1 with original Yongli proposal + variance factor
  ## pval4 : cov with original Yongli proposal + var ratio
  ## pval5 : cov with original Yongli proposal + (var ratio)^{-1}
  ##
  ## pval6 : cov 1 with transformation H + variance factor
  ## pval7 : cov with transformation H + var ratio
  ## pval8 : cov with transformation H + (var ratio)^{-1}

  zhang_q_cor1_pv200 <- mean(stat_q_cor1 < bt.stat_q_cor1[1:200])
  zhang_q_cor2_pv200 <- mean(stat_q_cor2 < bt.stat_q_cor2[1:200])
  zhang_q_cor3_pv200 <- mean(stat_q_cor3 < bt.stat_q_cor3[1:200])
  
  zhang_dif_cov1_pv200 <- mean(stat_dif_cov1 < bt.stat_dif_cov1[1:200])
  zhang_dif_cov2_pv200 <- mean(stat_dif_cov2 < bt.stat_dif_cov2[1:200])
  zhang_dif_cov3_pv200 <- mean(stat_dif_cov3 < bt.stat_dif_cov3[1:200])

  zhang_transf_cov1_pv200 <- mean(stat_transf_cov1 < bt.stat_transf_cov1[1:200]) 
  zhang_transf_cov2_pv200 <- mean(stat_transf_cov2 < bt.stat_transf_cov2[1:200]) 
  zhang_transf_cov3_pv200 <- mean(stat_transf_cov3 < bt.stat_transf_cov3[1:200])

  ##################################################################
  ##################################################################
  ##################################################################
  ## Output guide for full B = 50.
  ## pval1: cor test with original Yongli proposal.
  ## pval2: cor test with variance stabilizing transformation H.
  ##
  ## pval3 : cov 1 with original Yongli proposal + variance factor
  ## pval4 : cov with original Yongli proposal + var ratio
  ## pval5 : cov with original Yongli proposal + (var ratio)^{-1}
  ##
  ## pval6 : cov 1 with transformation H + variance factor
  ## pval7 : cov with transformation H + var ratio
  ## pval8 : cov with transformation H + (var ratio)^{-1}

  zhang_q_cor1_pv50 <- mean(stat_q_cor1 < bt.stat_q_cor1[201:250])
  zhang_q_cor2_pv50 <- mean(stat_q_cor2 < bt.stat_q_cor2[201:250])
  zhang_q_cor3_pv50 <- mean(stat_q_cor3 < bt.stat_q_cor3[201:250])
  
  zhang_dif_cov1_pv50 <- mean(stat_dif_cov1 < bt.stat_dif_cov1[201:250])
  zhang_dif_cov2_pv50 <- mean(stat_dif_cov2 < bt.stat_dif_cov2[201:250])
  zhang_dif_cov3_pv50 <- mean(stat_dif_cov3 < bt.stat_dif_cov3[201:250])

  zhang_transf_cov1_pv50 <- mean(stat_transf_cov1 < bt.stat_transf_cov1[201:250]) 
  zhang_transf_cov2_pv50 <- mean(stat_transf_cov2 < bt.stat_transf_cov2[201:250]) 
  zhang_transf_cov3_pv50 <- mean(stat_transf_cov3 < bt.stat_transf_cov3[201:250])

  mre <- c(
    zhang_q_cor1_pv, zhang_q_cor2_pv, zhang_q_cor3_pv,
    zhang_dif_cov1_pv, zhang_dif_cov2_pv, zhang_dif_cov3_pv,
    zhang_transf_cov1_pv, zhang_transf_cov2_pv, zhang_transf_cov1_pv,
    
    zhang_q_cor1_pv200, zhang_q_cor2_pv200, zhang_q_cor3_pv200,
    zhang_dif_cov1_pv200, zhang_dif_cov2_pv200, zhang_dif_cov3_pv200,
    zhang_transf_cov1_pv200, zhang_transf_cov2_pv200, zhang_transf_cov1_pv200,
    
    zhang_q_cor1_pv50, zhang_q_cor2_pv50, zhang_q_cor3_pv50,
    zhang_dif_cov1_pv50, zhang_dif_cov2_pv50, zhang_dif_cov3_pv50,
    zhang_transf_cov1_pv50, zhang_transf_cov2_pv50, zhang_transf_cov1_pv50
  )
  
  return(mre)

}

##############################################################################################################################################
######################### The Proposed Testing Procedure:
##############################################################################################################################################

zhang_trunc_fun <- function(X1, X2, B) {
  
  ###############################
  ###############################
  ## Value Initialization
  p <- ncol(X1)
  Ip <- diag(rep(1, p))
  Ipp <- rep(1, p) %*% t(rep(1, p))
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  n <- n1 + n2
  
  ###############################
  ###############################
  ## Calculating thresholds:
  e1 <- rep(1, n1)
  e2 <- rep(1, n2)
  x1bar <- apply(X1, MARGIN = 2, mean)
  x2bar <- apply(X2, MARGIN = 2, mean)
  xx1 <- X1 - ((e1) %*% t(x1bar)) ## Centered x data
  xx2 <- X2 - ((e2) %*% t(x2bar)) ## Centered y data  
  mom4_xx1 <- apply(xx1, MARGIN = 2, function(x) mean(x^4))
  mom4_xx2 <- apply(xx2, MARGIN = 2, function(x) mean(x^4))
  mom4_max <- max(c(mom4_xx1, mom4_xx2))
  threshold <- 2 * sqrt(sqrt(mom4_max * n / log(p)))  ## TODO: DECIDE ON THE CHOICE OF THRESHOLD....
  
  ###############################
  ###############################
  ## Thresholding data:
  X1_thresholded <- X1 * (abs(X1) < threshold) + threshold * (2 * ((X1) >= 0) - 1) * (abs(X1) >= threshold)
  X2_thresholded <- X2 * (abs(X2) < threshold) + threshold * (2 * ((X2) >= 0) - 1) * (abs(X2) >= threshold)

  return(zhang_fun(X1_thresholded, X2_thresholded, B))

}

##############################################################################################################################################
######################### Covariance tests:
##############################################################################################################################################

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Li & Chen (2012): Two Sample Tests for High-Dimensional Covariance Matrices
###### We use implementation in package "covEqual" 
###### package inclusion found in "./1_requirements.R"

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Cai, Liu and Xia (2013)  Two-Sample Covariance Matrix Testing and Support Recovery in High-Dimensional and Sparse Settings
cov_CaiLiuXia2013_fun <- function(x, y) {
  n1<-nrow(x); n2<-nrow(y); p<-ncol(x)
  e1<-rep(1,n1); e2<-rep(1,n2)
  
  xbar = apply(x,2,mean); ybar = apply(y,2,mean);
  xcov = cov(x)*(n1-1)/n1;  ycov = cov(y)*(n2-1)/n2;
  xx = x-((e1)%*%t(xbar));  yy = y-((e2)%*%t(ybar));
  
  t1 = (t(xx^2))%*%(xx^2)/n1 -2*((t(xx))%*%(xx)/n1)*(xcov) + xcov^2;
  t2 = (t(yy^2))%*%(yy^2)/n2 -2*((t(yy))%*%(yy)/n2)*(ycov) + ycov^2;
  
  Mhat = max(((xcov-ycov)^2)/(t1/n1+t2/n2));
  cri = Mhat - 4*log(p)+log(log(p))
  pv = 1- exp(- 1/sqrt(8*pi)* exp (-cri/2))
  c(Mhat,pv)  
   
}

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Zheng, Lin, Guo and Yin (2020): Testing Homogeneity of High-Dimensional Covariance matrices.
cov_ZhengEtAl2020_fun <- function(x, y) {
  
  trace <- function(m) {
    return(sum(diag(m)))
  }

  ############ Initialization:
  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  xbar <- apply(x, MARGIN = 2, mean)
  ybar <- apply(y, MARGIN = 2, mean)
  S1 <- cov(x) ## REMEMBER: cov function uses 1/(n-1) instead of 1/n
  S2 <- cov(y)
  

  ################################ 
  ## Calculating test statistic:
  ################################

  ############ Calculating theta1, theta2 matrices:
  e1 <- rep(1, n1)
  e2 <- rep(1, n2)
  xx <- x - ((e1) %*% t(xbar)) ## Centered x data
  yy <- y - ((e2) %*% t(ybar)) ## Centered y data  
  theta1 <- (t(xx^2) %*% (xx^2) / n1) - (2 * (t(xx) %*% xx / n1) * S1) + (S1^2)
  theta2 <- (t(yy^2) %*% (yy^2) / n2) - (2 * (t(yy) %*% yy / n2) * S2) + (S2^2)

  ############ Calculating delta12
  delta12 <- (S1 - S2)^2 / ((theta1 / n1) + (theta2 / n2))

  ############ Calculating regularizing constants:
  q <- (-2) * log((-1) * sqrt(8 * pi) * log(0.985)) 
  s_n1_n2_p <- ((log(log(n1 / 2 + n2 / 2)) - 1)^2 / 4 + 1) * (4 * log(p) - log(log(p))) + q
  K0 <- p^2
  w12 <- 1

  ############ Test statistic:
  T_21 <- w12 * trace((S1 - S2) %*% (S1 - S2))
  T_22 <- K0 * 1 * (max(delta12) > s_n1_n2_p)

  T2 <- T_21 + T_22
  
  ################################ 
  ## Calculating null distribution:
  ################################
  
  ############ muhat21:
  muhat21 <- (n1^2 - n1 - 1) * (trace(S1))^2 / (n1 * (n1 - 1)^2)            ## xx part
  muhat21 <- muhat21 + (n2^2 - n2 - 1) * (trace(S2))^2 / (n2 * (n2 - 1)^2)  ## yy part

  ############ muhat2:
  ## muhat2 <- (trace((xx^2) %*% t(xx^2)) - (2 * trace((xx) %*% t(xx)) * trace(S1)) + n1 * (trace(S1)^2)) / ((n1 - 2)^2)
  ## muhat2 <- muhat2 - n1 * (trace(S1 %*% S1) - trace(S1)^2 / (n1 - 2)) / ((n1 + 2)^2)

  ## muhat2 <- (trace((yy^2) %*% t(yy^2)) - (2 * trace((yy) %*% t(yy)) * trace(S2)) + n2 * (trace(S2)^2)) / ((n2 - 2)^2)
  ## muhat2 <- muhat2 - n2 * (trace(S2 %*% S2) - trace(S2)^2 / (n2 - 2)) / ((n2 + 2)^2)

  prod_x <- rep(0, n1)
  for(l in 1:n1) {
    prod_x[l] <- sum(xx[l, ]* xx[l, ])
    
  }
  muhat2_hand <- sum( (prod_x - trace(S1))^2 ) / ((n1 - 2)^2)
  prod_y <- rep(0, n2)
  for(l in 1:n2) {
    prod_y[l] <- sum(yy[l, ] * yy[l, ])
  }
  muhat2_hand <- muhat2_hand + sum( (prod_y - trace(S2))^2 ) / ((n2 - 2)^2)
  muhat2_hand <- muhat2_hand - n1 * (trace(S1 %*% S1) - trace(S1)^2 / (n1 - 2) ) / ((n1 + 2)^2)
  muhat2_hand <- muhat2_hand - n2 * (trace(S2 %*% S2) - trace(S2)^2 / (n2 - 2) ) / ((n2 + 2)^2)


  ############ sigmahat2:
  S  <- ((n1 - 1) * S1 + (n2 - 1) * S2) / (n1 + n2 - 2)
  trsigma2hat <- trace(S %*% S) - (trace(S))^2 / (n1 + n2 - 2)
  sigma2hat2 <- 4 * (1 / (n1 - 1) + 1 / (n2 - 1))^2 * (trsigma2hat)^2
  sigmahat2 <- sqrt(sigma2hat2)

  ################################ 
  ## Calculating p-value:
  ################################

  norm_stat <- (T2 - muhat21 - muhat2_hand) / sigmahat2
  pval <- 1 - pnorm(norm_stat, mean = 0, sd = 1)
  return(pval)

}


##############################################################################################################################################
######################### Correlation tests:
##############################################################################################################################################

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Cai and Zhang (2016)  Inference for High Dimensional Differential Correlation Matrices. (original)
cor_CaiZhang2016_fun <- function(x,y) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  e1 <- rep(1,n1)
  e2 <- rep(1,n2)
  
  xbar <- apply(x, 2, mean)
  ybar <- apply(y, 2, mean)
  xcov <- cov(x) * (n1-1)/n1
  ycov <- cov(y) * (n2 - 1) / n2;
  xcor <- cor(x)
  ycor <- cor(y)
  sig2x <- diag(xcov)
  sig2y <- diag(ycov)
  xx <- x - ((e1) %*% t(xbar))
  yy <- y - ((e2) %*% t(ybar))
  
  #xa<-array(0, dim=c(p,p, n1))
  #ya<-array(0, dim=c(p,p, n2))
  
  #for (k in 1: n1){
  #  xa[,,k]<-(xx[k,]%*%t(xx[k,])/sqrt(sig2x%*%t(sig2x))-xcor/2*(replicate(p, xx[k,]^2/sig2x)+t(replicate(p, xx[k,]^2/sig2x))))^2
  #}
  
  #for (k in 1: n2){
  #  ya[,,k]<-(yy[k,]%*%t(yy[k,])/sqrt(sig2y%*%t(sig2y))-ycor/2*(replicate(p, yy[k,]^2/sig2y)+t(replicate(p, yy[k,]^2/sig2y))))^2
  #}
  
  #xeta = apply(xa,c(1,2),mean)
  #yeta = apply(ya,c(1,2),mean)
  
  xeta <- matrix(0, p, p)
  yeta <- matrix(0, p, p)
  for (i in 1:p){
    for (j in 1:p) {
      zx<-rep(0,n1)
      for (k in 1:n1) {
        zx[k] <-  xx[k,i] * xx[k,j] / sqrt(sig2x[i] * sig2x[j]) - xcor[i,j] / 2 * (xx[k,i]^2 / sig2x[i] + xx[k,j]^2 / sig2x[j]) ## OLD VERSION
        # zx[k] <-  xx[k, i] * xx[k, j]  - xcov[i, j]
      }
      xeta[i, j] <- sum(zx^2) / n1 ## OLD VERSION
      # xeta[i, j] <- sum(zx^2) / (n1 * sig2x[i] * sig2x[j])
      
      zy<-rep(0, n2)
      for (k in 1:n2) {
        zy[k] <- yy[k, i]*yy[k, j] / sqrt(sig2y[i] * sig2y[j]) - ycor[i, j] / 2 * (yy[k,i]^2 / sig2y[i] + yy[k,j]^2/sig2y[j]) ## OLD VERSION
        # zy[k] <-  yy[k, i] * yy[k, j]  - ycov[i, j]
      }
      yeta[i, j] = sum(zy^2) / n2; ## OLD VERSION
      # yeta[i, j] <- sum(zy^2) / (n2 * sig2y[i] * sig2y[j])
    }
  }
  TIJ<-(xcor - ycor)^2 / (xeta / n1 + yeta / n2)
  Mhat = max(TIJ[row(TIJ)<col(TIJ)])
  #cri1 = 4*log(p)-log(log(p))-2*log (-(8*pi)^(1/2)*log (1-alpha))
  cri <- Mhat - 4 * log(p) + log(log(p))
  pv = 1 - exp(- 1 / sqrt(8 * pi) * exp (-cri / 2))
  c(Mhat,pv)
}


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Zheng, Cheng, Guo and Zhu (2019): Test for High-Dimensional Correlation Matrices
cor_ZhengEtAl2019_fun <- function(x, y) { 
  
  tr <- function(x){sum(diag(x))}
  
  ###############################
  ###############################
  ## Initialization:
  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  e1 <- rep(1,n1)
  e2 <- rep(1,n2)
  
  xbar <- apply(x, 2, mean)
  ybar <- apply(y, 2, mean)
  xcov <- cov(x) * (n1 - 1) / n1
  ycov <- cov(y) * (n2 - 1) / n2
  xcor <- cor(x)
  ycor <- cor(y)
  sig2x <- diag(xcov)
  sig2y <- diag(ycov)
  xx <- x - ((e1) %*% t(xbar))
  yy <- y - ((e2) %*% t(ybar))
  
  ###############################
  ###############################
  ## Calculation of Theta:
  xeta <- matrix(0, p, p)
  yeta <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      zx <- rep(0, n1)
      for (k in 1:n1) {
        zx[k] <- xx[k, i] * xx[k, j] / sqrt(sig2x[i] * sig2x[j]) - xcor[i, j] / 2 * (xx[k, i]^2 / sig2x[i] + xx[k, j]^2 / sig2x[j])
      }
      xeta[i, j] <- sum(zx^2) / n1
      zy<-rep(0, n2)
      for (k in 1:n2) {
        zy[k] <- yy[k,i] * yy[k, j] / sqrt(sig2y[i] * sig2y[j]) - ycor[i, j] / 2 * (yy[k, i]^2 / sig2y[i] + yy[k, j]^2 / sig2y[j])
      }
      yeta[i, j] <- sum(zy^2) / n2
    }
  }
  
  ###############################
  ###############################
  ## Calculating Mn2 statistic:
  TIJ <- (xcor - ycor)^2 / (xeta / n1 + yeta / n2)
  
  Tn2 <- max(TIJ[row(TIJ)<col(TIJ)])
  Ln2 <- tr((xcor-ycor)%*%(xcor-ycor))
  
  ## First test statistic:
  alpha <- 0.05
  q_val <- 1 - alpha / 2
  u0_prime <- (-2) * log( -log(q_val) / sqrt(8 * pi))
  s_n1_n2_p <- (4 + (log(log(n1 + n2)) - 1)^2) * (log(p) - 0.25 * log(log(p))) + u0_prime
  c02 <- p * p
  Mn2 <- Ln2 + c02 * (Tn2 > s_n1_n2_p)

  ###############################
  ###############################
  ## Calculating Mn2_prime statistic:
  calpha2 <- qnorm(q_val)
  c02_prime <- calpha2 / u0_prime
  output_jz2 <- jz2_fun(xx,yy)
  mu12 <- output_jz2[1]
  a0_hat <- output_jz2[2]
  Mn2_prime <- max(
    (Ln2 - mu12) / (2 * (1 / (n1 - 1) + 1 / (n2 - 1)) * p * a0_hat),
    c02_prime * (Tn2 - 4 * log(p) + log(log(p))))
  
  ###############################
  ###############################
  ## Conclusion:
  stat <- (Mn2 - mu12) / (2 * (1 / (n1 - 1) + 1 / (n2 - 1)) * p * a0_hat)
  pval <- 1 - pnorm(stat, mean = 0, sd = 1)
  
  stat_prime <- Mn2_prime
  pval_prime <- 2 * (1 - pnorm(stat_prime, mean = 0, sd = 1))

  output <- c(stat, pval, stat_prime, pval_prime)
  return(output)
  
}

####################################
####################################
## Auxiliary function: calculate the mean for two sample test
jz2_fun <- function(xx, yy) {
  
  tr <- function(x){sum(diag(x))}

  n1 <- nrow(xx)
  n2 <- nrow(yy)
  p <- ncol(xx)
  meany <- apply(yy, 2, mean)
  meanx <- apply(xx, 2, mean)
  Ip <- diag(rep(1, p))
  S1 <- cor(xx)
  Sig1 = cov(xx)
  S2 <- cor(yy)
  Sig2<- cov(yy)

  ###############################
  ###############################
  ## Kurtosis beta1
  Zhat1 <- rep(0, n1)
  for (j in 1:n1) {
    Zhat1[j] <- t(xx[j, ] - meanx) %*% (xx[j, ] - meanx) 
  }
  Vhat1 <- 1 / (n1 - 1) * sum((Zhat1 - mean(Zhat1))^2)
  betax <- (n1 * Vhat1 - 2 * ((n1 - 1) * tr(Sig1 %*% Sig1) - (tr(Sig1))^2)) / (n1 * tr(Sig1^2))
  
  ###############################
  ###############################
  ## Kurtosis beta2
  Zhat2 <- rep(0, n2);
  for (j in 1:n2) {
    Zhat2[j] <- t(yy[j, ] - meany) %*% (yy[j, ] - meany) 
  }
  Vhat2 <- 1 / (n2 - 1) * sum((Zhat2 - mean(Zhat2))^2)
  betay <- (n2 * Vhat2 - 2 * ((n2 - 1) * tr(Sig2 %*% Sig2) - (tr(Sig2))^2)) / (n2 * tr(Sig2^2))
 
  ###############################
  ###############################
  ## Pre-calculations XXX:
  ss5x <- rep(0,n1)
  ss4x <- rep(0,n1)
  ss3x <- rep(0,n1)

  ttx <- diag(1 / sqrt(diag(Sig1)))
  xxx <- ttx %*% (t(xx) - matrix(rep(meanx, n1), p, n1))
  for (j in 1:n1) {
    xxj <- xxx[,j]
    tempx <- xxj %*% t(xxj)
    tempx1 <- diag(diag(tempx) - rep(1, p))
   
    ss3x[j] <- tr(S2 %*% tempx %*% tempx1)
    ss4x[j] <- tr(S2 %*% (S1 - tempx / (n1 - 1)) %*% tempx1 %*% tempx1)
    ss5x[j] <- tr((S1 - tempx / (n1 - 1)) %*% tempx1 %*% S2 %*% tempx1)
  
  }

  ###############################
  ###############################
  ## Pre-calculations YYY:
  ss5y <- rep(0, n2)
  ss4y <- rep(0, n2)
  ss3y <- rep(0, n2)

  tty <- diag(1 / sqrt(diag(Sig2)))
  yyy <- tty %*% (t(yy) - matrix(rep(meany, n2), p, n2))
  for (j in 1:n2) {
    yyj <- yyy[,j]
    tempy <- yyj %*% t(yyj)
    tempy1 <- diag(diag(tempy) - rep(1, p));
   
    ss3y[j] <- tr(S1 %*% tempy %*% tempy1)
    ss4y[j] <- tr(S1 %*% (S2 - tempy / (n2 - 1)) %*% tempy1 %*% tempy1)
    ss5y[j] <- tr((S2 - tempy / (n2 - 1)) %*% tempy1 %*% S1 %*% tempy1)
  }

  ###############################
  ###############################  
  ## Calculation of A0
  w1 <- n1 / (n1 + n2)
  w2 <- 1 - w1
  a1hat <- (tr(S1 %*% S1) - (n1^2 - n1 - 1) * p^2 / (n1 * (n1 - 1)^2) - betax * n1 * p / ((n1 - 1)^2)) * (n1 - 1)^2 / ((n1^2 - n1 + 2) * p)
  a2hat <- (tr(S2 %*% S2) - (n2^2 - n2 - 1) * p^2 / (n2 * (n2 - 1)^2) - betay * n2 * p / ((n2 - 1)^2)) * (n2 - 1)^2 / ((n2^2 - n2 + 2) * p)
  a0hat <- w1 * a1hat + w2 * a2hat

  ###############################  
  ###############################
  ## Calculation of B0
  tt31 <- 1 / (p * betax) * (sum(ss3x) / n1 - 2 * p * a0hat)
  tt32 <- 1 / (p * betay) * (sum(ss3y) / n1 - 2 * p * a0hat)
  b0hat <- w1 * tt31 + w2 * tt32  ##b0hat


  ###############################  
  ###############################
  ## Calculation of C0
  tt41 <- 1 / (p * betax) * (sum(ss4x) / n1 - 2 * p * a0hat)
  tt42 <- 1 / (p * betay) * (sum(ss4y) / n2 - 2 * p * a0hat)
  c0hat <- w1 * tt41 + w2 * tt42  ##c0hat
  
  ###############################  
  ###############################
  ## Calculation of D0
  d1hat <- 1 / p * sum((ss5x)) / (n1)
  d2hat <- 1 / p * sum((ss5y)) / (n2)
  
  ###############################  
  ###############################
  ## Mean:
  
  cc1 <- (n1^2 - n1 - 1) * p^2 / (n1 * (n1 - 1)^2) + (n2^2 - n2 - 1) * p^2 / (n2 * (n2 - 1)^2)
  cc2 <- -(p * (2 * n1^2 + n1 + 1) / ((n1 - 1)^3) + p * (2 * n2^2 + n2 + 1) / ((n2 - 1)^3)) * a0hat
  cc3 <- betax * n1 * p / ((n1 - 1)^2) + betay * n2 * p / ((n2 - 1)^2)
  
  cc4 <- -(2 * p * n1 * betax / ((n1 - 1)^2) + 2 * p * n2 * betay / ((n2 - 1)^2)) * b0hat
  cc5 <- (betax * p * (n1^2 - 5*n1) / (2 * (n1 - 1)^3) + betay * p * (n2^2 - 5*n2) / (2 * (n2 - 1)^3)) * c0hat
  cc6 <- (p * (n1 - 3) * n1 * d1hat / (2 * (n1 - 1)^3) + p * (n2 - 3) * n2 * d2hat/ (2 * (n2 - 1)^3))
  
  muz12 <- cc1 + cc2 + cc3 + cc4 + cc5 + cc6
  return(c(muz12, a0hat))
}

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Zhou, Han, Zhang and Liu (2019): An extreme-value approach for testing the equality 
##                                of large U-statistic based correlation matrices
cor_ZhouEtAl2019_fun <- function(x, y) {
  return(1) ## TODO: program
  ## It resulted that this test does not fit our scenarios of interest.
  ## From this, we will not consider it as part of our simulations.
}