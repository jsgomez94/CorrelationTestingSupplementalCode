##########################
##########################
##Model 1
## Model 2 of chenzhangzhong(2010)
m1 <- function(p) {
  id_p <- diag(rep(1, p))
  sigma0 <- id_p

  return(sigma0)
}

##########################
##########################
##Model 2 : Model 2.1 of Zheng et al. (2019)
m2 <- function(p) {
  id_p <- diag(rep(1, p))
  sigma0 <- id_p
  d <- diag(runif(p, 0.5, 2.5))

  rho1 <- 0.25
  sigma0 <- (rho1)^abs(col(sigma0) - row(sigma0))
  sigma0 <- d^(1 / 2) %*% sigma0 %*% d^(1 / 2)

  return(sigma0)
}

##########################
##########################
## Original Model 2 of Cai et al (2013)
## (Two-sample covariance matrix testing and
## support recovery in high-dimensional and sparse settings)
m3 <- function(p) {
  id_p <- diag(rep(1, p))
  sigma0 <- id_p

  d <- diag(runif(p, 0.5, 2.5))
  sigma0 <- (0.5)^abs(col(sigma0) - row(sigma0))
  sigma0 <- d^(1 / 2) %*% sigma0 %*% d^(1 / 2)

  return(sigma0)
}

##########################
##########################
## Original Model 4 of Cai et al (2013)
## (Two-sample covariance matrix testing and
## support recovery in high-dimensional and sparse settings)
#m4 <- function(p) {
#  id_p <- diag(rep(1, p))
#  sigma0 <- id_p
#
#  d <- diag(runif(p, 1, 5))
#  sign <- (-1)^(col(sigma0) + row(sigma0))
#  magnitude <- (0.4)^(abs(col(sigma0) - row(sigma0))^(1 / 10))
#  sigma0 <- d %*% (sign * magnitude) %*% d
#
#  return(sigma0)
#}

##########################
##########################
## Original Model 1 of Cai et al (2013)
## (Two-sample covariance matrix testing and
## support recovery in high-dimensional and sparse settings)
m4 <- function(p) {
  id_p <- diag(rep(1, p))
  sigma0 <- id_p
  d <- diag(runif(p, 0.5, 2.5))
  for (k in 1:floor(p / 5)) {
    sigma0[
      (row(sigma0) != col(sigma0)) &
      (row(sigma0) >= 5 * (k - 1) + 1) &
      (row(sigma0) <= 5 * k) &
      (col(sigma0) >= 5 * (k - 1) + 1) &
      (col(sigma0) <= 5 * k)
    ] <- 0.5
  }
  sigma0 <- d^(1 / 2) %*% sigma0 %*% d^(1 / 2)

  return(sigma0)
}

##########################
##########################
## Original Model 3 of Cai et al (2013)
## Two-sample covariance matrix testing and
## support recovery in high-dimensional and sparse settings)
m5 <- function(p) {

  id_p <- diag(rep(1, p))
  sigma0 <- id_p

  sigma0[col(sigma0) < row(sigma0)] <- rbinom(p * (p - 1) / 2, 1, 0.05) * 0.5
  sigma0 <- t(sigma0) + sigma0
  diag(sigma0) <- 1
  d <- diag(runif(p, 0.5, 2.5))
  delta <- abs(min(eigen(sigma0)$values)) + 0.05
  sigma0 <- d^(1 / 2) %*% (sigma0 + delta * id_p) %*% d^(1 / 2) / (1 + delta)

  return(sigma0)
}


####################################################
####################################################
####################################################
####################################################

## Dense diff <- 2
## Sparse diff <- 1
## No diff <- 0
mat_pair <- function(sigma0, p, sp) {
  id_p <- diag(rep(1, p))
  mat1_p <- rep(1, p) %*% t(rep(1, p))

  if (sp == 0) {
    sigma1 <- sigma0
    sigma2 <- sigma0
  }
  if (sp == 1) {
    u_mat <- matrix(0, p, p)
    uu_mat <- matrix(1:p^2, p, p)
    u_mat[
      matrix(
        uu_mat %in% sample(uu_mat[upper.tri(uu_mat)], 4, replace = FALSE), p, p
      )
    #V2: ] <- runif(4, 0, 4) * max(diag(sigma0)) ## Original (model2)
    #V3: ] <- runif(4, 0, 1) * max(diag(sigma0)) ## More challenging by reducing entry size. size 4 -> size 1.
    #V4: ] <- runif(4, 0, 2) * max(diag(sigma0)) ## less challenging size 1 -> size 2.
    #V5: ] <- runif(4, 0, 1.5) * max(diag(sigma0)) ## more challenging size 2 -> size 1.5.
    #V5: ] <- runif(4, 0, 1.25) * max(diag(sigma0)) ## even more challenging size 1.5 -> size 1.25.
    ] <- runif(4, 0, 1.25) * max(diag(sigma0)) ## 
    u_mat <- u_mat + t(u_mat)

    delta <- abs(
      min(c(eigen(sigma0)$values, eigen(sigma0 + u_mat)$values))
      ) + 0.05
    sigma1 <- sigma0 + delta * id_p
    sigma2 <- sigma1 + u_mat
  }
  if (sp == 2) {
    #V2: u_mat <- 0.05 * (mat1_p - id_p) ## More challenging by reducing entry size. 0.05 -> 0.01.
    #V3: u_mat <- 0.01 * (mat1_p - id_p) ## More challenging by reducing entry size. 0.05 -> 0.01.
    #V4: u_mat <- 0.03 * (mat1_p - id_p) ## More challenging by reducing entry size. 0.05 -> 0.03.
    u_mat <- 0.025 * (mat1_p - id_p) ## More challenging by reducing entry size. 0.05 -> 0.03.

    delta <- abs(
      min(c(eigen(sigma0)$values, eigen(sigma0 + u_mat)$values))
      ) + 0.05
    sigma1 <- sigma0 + delta * id_p
    sigma2 <- sigma1 + u_mat
  }

  return(list(sigma1, sigma2))
}



####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
## Unused functions.

##########################
##########################
##Model 7 : Model 2.2 of Zheng et al. (2019)
m7 <- function(p, sp) {
  id_p <- diag(rep(1, p))
  mat1_p <- rep(1, p) %*% t(rep(1, p))
  sigma2 <- sigma1 <- id_p
  rho1 <- 0.5
  sigma1 <- (rho1)^abs(col(sigma1) - row(sigma1))
  if (sp == 0) {
    sigma2 <- sigma1
  } else {
    sigma2 <- sigma1 + 0.05 * (mat1_p - id_p)
  }

  list(sigma1, sigma2)
}


##########################
##########################
##Model 8 : Model 2.2 of Zheng et al. (2019)
m8 <- function(p, sp) {
  id_p <- diag(rep(1, p))
  mat1_p <- rep(1, p) %*% t(rep(1, p))
  sigma2 <- sigma1 <- id_p
  rho1 <- 0.5
  sigma1 <- (rho1)^abs(col(sigma1) - row(sigma1))
  if (sp == 0) {
    sigma2 <- sigma1
  } else {
    sigma2 <- sigma1 + 0.08 * (mat1_p - id_p)
  }

  list(sigma1, sigma2)
}


##########################
##########################
##Model 9 : Model 2.3 of Zheng et al. (2019)
m9 <- function(p, sp) {

  id_p <- diag(rep(1, p))
  sigma2 <- sigma1 <- id_p
  rho1 <- 0.05
  sigma1 <- (rho1)^abs(col(sigma1) - row(sigma1))
  if (sp == 0) {
    sigma2 <- sigma1
  } else {
    sigma2 <- sigma1 +
      exp(0.008 * p) /
      (1 + exp(0.008 * p)) *
      (id_p[, 2] %*% t(id_p[, 1]) + id_p[, 1] %*% t(id_p[, 2]))
    }

  list(sigma1, sigma2)
}

##########################
##########################
##Model 10 : Model 2.3 of Zheng et al. (2019)
m10 <- function(p, sp) {
  id_p <- diag(rep(1, p))
  sigma2 <- sigma1 <- id_p
  rho1 <- 0.1
  sigma1 <- (rho1)^abs(col(sigma1) - row(sigma1))
  if (sp == 0) {
    sigma2 <- sigma1
  } else {
    sigma2 <- sigma1 +
      exp(0.008 * p) /
      (1 + exp(0.008 * p)) *
      (id_p[, 2] %*% t(id_p[, 1]) + id_p[, 1] %*% t(id_p[, 2]))
  }

  list(sigma1, sigma2)
}