# Load packages
packages <- c('dplyr','glmnet','ncvreg','Rcpp','doParallel')
for (pkg in packages){
  if (!requireNamespace(pkg, quietly = TRUE)){
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

#-------------------------------------------------------------------------------
# bar_ccd
#-------------------------------------------------------------------------------
# Estimates a BAR model using Cyclical Coordinate Descent (CCD) in C++
# sourceCpp('tsBAR/src/bar.cpp');

#-------------------------------------------------------------------------------
#                                 CMP_DGP
#-------------------------------------------------------------------------------
# Cubadda-Margaritella-Prifti DGP (CMP_DGP)
CMP_DGP <- function(T_, K, r, rho, sigma = 0.1, cov.mat = 'Toeplitz', eig.max=0.75, eig.min=0.25, burnin=50) {
  T_new <- T_ + burnin
  if (r == 1) {
    alternate_sign <- rep(c(1, -1), length.out = K)
    d <- alternate_sign * 0.5
    A <- diag(d)
  } else {
    alternate_sign <- rep(c(1, -1), length.out = r)
    d  <- alternate_sign * seq(eig.max, eig.min, length.out = r)
    D  <- diag(d)
    tt <- seq(1, r, 1)
    C  <- matrix(0, r, r)
    for (i in 1:r){
      C[,i] = cos(2*i*tt*pi/r) + sin(2*i*tt*pi/r)
    }
    A <- kronecker(diag(K/r), solve(C)) %*%
      kronecker(diag(K/r), D) %*%
      kronecker(diag(K/r), C)
  }

  eps <- cov.toeplitz(K = K, T_ = T_new, sigma = sigma, rho = rho)$errors

  # dataset
  X     <- matrix(nrow = ncol(A), ncol = T_new)
  X[,1] <- eps[,1]
  for (t in 2:T_new) {
    X[,t] <- A %*% X[,t-1] + eps[,t]
  }

  dataset <- t(X[,burnin:T_new])
  return(list('dataset' = dataset, 'A' = A))
}

#-------------------------------------------------------------------------------
# cov.toeplitz
#-------------------------------------------------------------------------------
# Generates a Toeplitz covariance matrix
cov.toeplitz <- function(T_, K, sigma, rho) {

  # generate covariance matrix
    m.ind <- matrix(rep(seq(1,K),K),K,K)
    if(rho == 0) {
      corr <- diag(rep(1,K))
    }else{
      corr <- exp(log(rho) * abs(m.ind - t(m.ind)))                             # correlation matrix
    }
    if(length(sigma)==1){
      sigma <- rep(sigma,K)                                                     # standard deviations
    }
    S <- diag(sigma)%*%corr%*%diag(sigma)                                       # covariance matrix

  # generate errors and starting values
  C   <- t(chol(S))
  e   <- matrix(rnorm(K*T_),K,T_)
  eps <- C%*%e

  # output
  list(cov = S, errors = eps)
}

#-------------------------------------------------------------------------------
#                                 create_lags
#-------------------------------------------------------------------------------
# Creates lags
create_lags <- function(x, p = 1, include.original = TRUE, trim = TRUE) {
  x  <- as.matrix(x)
  n  <- nrow(x)
  K  <- ncol(x)
  lx <- matrix(0, nrow = n, ncol = (p + include.original) * K)
  if (is.null(colnames(x))) {
    c.names <- rep('', K)
  } else {
    c.names <- colnames(x)
  }
  colnames(lx) <- rep(c.names, p + include.original)
  for (i in (1 - include.original):p) {
    cols               <- K * (i - 1 + include.original) + 1:K
    lx[(1+i):n, cols]  <- x[1:(n-i), ]
    colnames(lx)[cols] <- paste(c.names, ' l', i, sep = '')
  }
  return(lx[(1 + trim * p):n, , drop = FALSE])
}

#-------------------------------------------------------------------------------
# ic.glmnet
#-------------------------------------------------------------------------------
# Elastic-net estimate with information criteria tuning the lambda parameter
ic.glmnet <- function(x, y, crit = c('bic','ebic','aic','aicc','hqc'), ...) {
  if (is.matrix(x) == FALSE) {
    x <- as.matrix(x)
  }
  if (is.vector(y) == FALSE) {
    y <- as.vector(y)
  }
  crit      <- match.arg(crit)
  n         <- length(y)

  # model estimation
  model     <- glmnet::glmnet(x = x, y = y, ...)
  coef      <- coef(model)
  lambda    <- model$lambda
  df        <- model$df
  dev_ratio <- model[['dev.ratio']]
  null_dev  <- model$nulldev

  # fitted values, residuals, mean squared error, Sum of squared errors
  yhat      <- cbind(1,x) %*% coef
  res       <- (y- yhat)
  mse       <- colMeans(res^2)
  sse       <- colSums(res^2)

  # information critera: bic, ebic, aic, aicc, hqc.
  nvar      <- df                                                                # number of active variables
  bic       <- n * log(mse) + nvar * log(n)                                      # bayesian information criterion
  ebic      <- n * log(mse) + nvar * log(n) + (2 * nvar * 0.5 * log(ncol(x)))    # extended bayesian information criterion
  aic       <- n * log(mse) + 2 * nvar                                           # Akaike information criterion
  aicc      <- aic + (2 * nvar * (nvar+1)) / (n-nvar-1)                          # corrected Akaike information criterion
  hqc       <- n * log(mse) + 2 * nvar * log(log(n))                             # Hannan-Quinn information criterion

  sst       <- (n - 1) * var(y)
  r2        <- 1 - (sse / sst)                                                  # R square
  adjr2     <- 1 - ((1 - r2) * (n - 1) / (nrow(x) - nvar - 1))                  # Adjusted R square

  crit      <- switch(crit, bic = bic, ebic = ebic, aic = aic, aicc = aicc, hqc = hqc)

  selected  <- min(which(crit == min(crit)))                                    # lambda position minimizing the information criteria

  ic        <- c(bic = bic[selected], ebic = ebic[selected], aic = aic[selected],
                 aicc = aicc[selected], hqc = hqc[selected])

  result    <- list(coefficients = coef[,selected],
                    ic            = ic,
                    lambda        = lambda[selected],
                    lambda.index  = selected,
                    nvar          = nvar[selected],
                    lambda.max    = lambda[1],
                    dev.ratio     = dev_ratio[selected],
                    null.dev      = null_dev,
                    glmnet        = model,
                    residuals     = res[,selected],
                    fitted.values = yhat[,selected],
                    ic.range      = crit,
                    call          = match.call()
  )

  class(result) <- 'glmnet'

  return(result)
}

#-------------------------------------------------------------------------------
# ic.scad
#-------------------------------------------------------------------------------
# Estimates a SCAD model with information criteria tuning the lambda parameter
ic.scad <- function(y, x, crit = c('bic','ebic','aic','aicc','hqc'), alpha = 1,
                    nlambda = 100, family = 'gaussian', epsilon = 1e-4,
                    thresh = 1e-9, max_itt = 1e6, ...) {
  if (is.matrix(x) == FALSE) {
    x <- as.matrix(x)
  }
  if (is.vector(y) == FALSE) {
    y <- as.vector(y)
  }
  crit   <- match.arg(crit)
  n      <- length(y)

  # model estimation
  model  <- ncvreg::ncvreg(X = x, y = y, family = family, alpha = alpha,
                           penalty = 'SCAD', ...)
  coef   <- model$beta
  lambda <- model$lambda

  # fitted values, residuals, mean squared error, Sum of squared errors
  yhat   <- cbind(1,x) %*% coef
  res    <- (y - yhat)
  mse    <- colMeans(res^2)
  sse    <- colSums(res^2)

  # information criteria: bic, ebic, aic, aicc, hqc.
  nvar   <- colSums(coef[-1,] != 0)                                             # number of active variables
  bic    <- n * log(mse) + nvar * log(n)                                        # bayesian information criterion
  ebic   <- n * log(mse) + nvar * log(n) + (2 * nvar * 0.5 * log(ncol(x)))      # extended bayesian information criterion
  aic    <- n * log(mse) + 2 * nvar                                             # Akaike information criterion
  aicc   <- aic + (2 * nvar * (nvar+1)) / (n-nvar-1)                            # corrected Akaike information criterion
  hqc    <- n * log(mse) + 2 * nvar * log(log(n))                               # Hannan-Quinn information criterion

  sst      <- (n - 1) * var(y)
  r2       <- 1 - (sse / sst)                                                   # R square
  adjr2    <- 1 - ((1 - r2) * (n - 1) / (nrow(x) - nvar - 1))                   # Adjusted R square

  crit     <- switch(crit, bic = bic,ebic = ebic, aic = aic, aicc = aicc, hqc = hqc)

  selected <- min(which(crit == min(crit)))                                  # lambda position minimizing the information criteria

  ic       <- c(bic = bic[selected], ebic = ebic[selected], aic = aic[selected],
                aicc = aicc[selected], hqc = hqc[selected])

  result   <- list(coefficients  = coef[,selected],
                   ic            = ic,
                   lambda        = lambda[selected],
                   lambda.index  = selected,
                   nvar          = nvar[selected],
                   lambda.max    = lambda[1],
                   SCAD          = model,
                   residuals     = res[,selected],
                   fitted.values = yhat[,selected],
                   ic.range      = crit,
                   call          = match.call()
                   )

  class(result) <- 'SCAD'

  return(result)
}

#-------------------------------------------------------------------------------
# var.lasso
#-------------------------------------------------------------------------------
# LASSO/Adaptive LASSO applied to a VAR system of equations
var.lasso <- function(y, x, K, p = 1, alpha = 1, adaptive = c('none','lasso','ridge'),
                      crit = c('bic','ebic','aic','aicc','hqc'), nlambda = 100,
                      standardize = TRUE, intercept = TRUE, ...) {

  adaptive <- match.arg(adaptive)
  crit <- match.arg(crit)

  # 1st step
  if(adaptive == 'lasso'){
    eq.coef <- matrix(nrow = K*p, ncol = K)
    for(i in 1:K){
      eq <- ic.glmnet(x = x, y = y[,i], crit = crit, alpha = 1,
                      standardize = standardize, intercept = intercept)
      eq.coef[,i] <- eq$coefficients[-1]
    }
    # adaptive weights
    ada.w	<- abs(eq.coef)^(-1)
  } else if (adaptive == 'ridge'){
    eq.coef <- matrix(nrow = K*p, ncol = K)
    for(i in 1:K){
      eq <- ic.glmnet(x = x, y = y[,i], crit = crit, alpha = 0,
                      standardize = standardize, intercept = intercept)
      eq.coef[,i] <- eq$coefficients[-1]
    }
    # adaptive weights
    ada.w	<- abs(eq.coef)^(-1)
  } else {
    ada.w <- NULL
  }

  # 2nd step
  alasso.ic   <- list()
  alasso.est	<- list('x' = x, 'y' = y, 'coefficients' = NULL, 'ada.w' = ada.w,
                     'rss' = NULL, 'lambda' = NULL, 'lambda.max' = NULL,
                     'lambda.index' = NULL
  )
  for(i in 1:K){
    # Plain LASSO/Ridge
    if (is.null(ada.w)) {
      alasso <- ic.glmnet(x = x, y = y[,i], crit = crit, alpha = alpha,
                             standardize = standardize, intercept = intercept, ...)

      # collect estimates
      alasso.ic[[i]] <- list('ic'           = crit,
                             'ic.all'       = alasso$ic,
                             'coefficients' = alasso$coefficients,
                             'rss'          = sum(alasso$residuals^2),
                             'lambda'       = alasso$lambda,
                             'lambda.max'   = alasso$lambda.max,
                             'lambda.index' = alasso$lambda.index
      )
    } else {
      # Adaptive LASSO/Ridge
      if (sum(ada.w[,i] == Inf) == length(ada.w[,i])) {
        warning('Adaptive weights all equal to Inf, all candidated excluded in the first step. Adaptive Lasso not estimated.')
        alasso              <- list()
        alasso$rss		      <- sum((y[,i] - mean(y[,i]))^2)
        alasso$coefficients <- rbind(mean(y[,i]), matrix(0,nrow = ncol(x), ncol = 1))
        alasso$lambda       <- 0
        alasso$ic           <- c(0)
        alasso$lambda.max   <- NA
        alasso$lambda.index <- NA
      } else{
        wpen                     <- ada.w[,i]
        wpen[which(wpen == Inf)] <- 0                                           # variables with Inf weights are manually removed
        alasso <- ic.glmnet(x = x, y = y[,i], crit = crit, alpha = alpha,
                               standardize = standardize, intercept = intercept,
                               exclude = which(ada.w[,i] == Inf),                   # exclude variables with Inf weights
                               penalty.factor = wpen, ...)
      }

      # collect estimates
      alasso.ic[[i]] <- list('ic'           = crit,
                             'ic.all'       = alasso$ic,
                             'coefficients' = alasso$coefficients,
                             'rss'          = sum(alasso$residuals^2),
                             'lambda'       = alasso$lambda,
                             'lambda.max'   = alasso$lambda.max,
                             'lambda.index' = alasso$lambda.index
      )
    }
  }

  # collect results
  for (i in 1:K) {
    alasso.est$coefficients <- cbind(alasso.est$coefficients, alasso.ic[[i]]$coefficients)
    alasso.est$rss          <- cbind(alasso.est$rss, alasso.ic[[i]]$rss)
    alasso.est$lambda	      <- cbind(alasso.est$lambda, alasso.ic[[i]]$lambda)
    alasso.est$lambda.max	  <- cbind(alasso.est$lambda.max, alasso.ic[[i]]$lambda.max)
    alasso.est$lambda.index	<- cbind(alasso.est$lambda.index, alasso.ic[[i]]$lambda.index)
  }

  return(alasso.est)
}

#-------------------------------------------------------------------------------
# var.scad
#-------------------------------------------------------------------------------
# SCAD applied to a VAR system of equations
var.scad <- function(y, x, K, p = 1, crit = c('bic','ebic','aic','aicc','hqc'),
                     alpha = 1, nlambda = 100, family = 'gaussian',
                     epsilon = 1e-4, thresh = 1e-9, max_itt = 1e6,...) {
  scad.ic <- list()
  scad.est	<- list('x' = x, 'y' = y, 'coefficients' = NULL, 'rss' = NULL,
                   'lambda' = NULL, 'lambda.max' = NULL, 'lambda.index' = NULL)
  for (i in 1:K) {
    model_est    <- ic.scad(x, y = y[,i], crit = crit, alpha = alpha, ...)
    scad.ic[[i]] <- list('ic'           = crit,
                         'ic.all'       = model_est$ic,
                         'coefficients' = model_est$coefficients,
                         'rss'          = sum(model_est$residuals^2),
                         'lambda'       = model_est$lambda,
                         'lambda.max'   = model_est$lambda.max,
                         'lambda.index' = model_est$lambda.index
                         )
  }

  # collect results
  for (i in 1:K) {
    scad.est$coefficients <- cbind(scad.est$coefficients, scad.ic[[i]]$coefficients)
    scad.est$rss          <- cbind(scad.est$rss, scad.ic[[i]]$rss)
    scad.est$lambda	      <- cbind(scad.est$lambda, scad.ic[[i]]$lambda)
    scad.est$lambda.max	  <- cbind(scad.est$lambda.max, scad.ic[[i]]$lambda.max)
    scad.est$lambda.index	<- cbind(scad.est$lambda.index, scad.ic[[i]]$lambda.index)
  }

  return(scad.est)
}
