# Load packages
packages <- c('dplyr','readr','tidyr','stringr','zoo','remotes','glmnet','ncvreg',
              'doParallel','xtable','ggplot2')
for (pkg in packages){
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

#-------------------------------------------------------------------------------
# bar_ccd_var() / bar_ccd_reg()
#-------------------------------------------------------------------------------
# Estimates a BAR model using Cyclical Coordinate Descent (CCD) in C++
if (!requireNamespace('tsBAR', quietly = TRUE)) {
  remotes::install_github('orestprifti/tsBAR')
}
library('tsBAR', character.only = TRUE)

#-------------------------------------------------------------------------------
# cmp_dgp()
#-------------------------------------------------------------------------------
# Cubadda-Margaritella-Prifti DGP (cmp_dgp)
cmp_dgp <- function(T_, K, r, rho, sigma = 0.1, cov.mat = 'Toeplitz', eig.max=0.75, eig.min=0.25, burnin=50) {
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
# cov.toeplitz()
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
# create_lags()
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
# ic.glmnet()
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
# ic.scad()
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
# var.lasso()
#-------------------------------------------------------------------------------
# LASSO/Adaptive LASSO applied to a VAR system of equations
var.lasso <- function(y, x, p = 1, alpha = 1, adaptive = c('none','lasso','ridge'),
                      crit = c('bic','ebic','aic','aicc','hqc'), nlambda = 100,
                      standardize = TRUE, intercept = TRUE, ...) {

  adaptive <- match.arg(adaptive)
  crit <- match.arg(crit)

  K <- ncol(x)
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
    if (is.null(ada.w)) {                                                       # plain LASSO/Ridge
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
    } else {                                                                    # adaptive LASSO/Ridge
      if (sum(ada.w[,i] == Inf) == length(ada.w[,i])) {                        
        warning('Adaptive weights all equal to Inf, all candidated excluded in the first step. Adaptive Lasso not estimated.')
        alasso              <- list()
        alasso$rss		      <- sum((y[,i] - mean(y[,i]))^2)
        alasso$coefficients <- rbind(mean(y[,i]), matrix(0, nrow = ncol(x), ncol = 1))
        alasso$lambda       <- 0
        alasso$ic           <- c(0)
        alasso$lambda.max   <- NA
        alasso$lambda.index <- NA
      } else{
        wpen                     <- ada.w[,i]
        wpen[which(wpen == Inf)] <- 0                                           # variables with Inf weights are manually removed
        alasso <- ic.glmnet(x = x, y = y[,i], crit = crit, alpha = alpha,
                               standardize = standardize, intercept = intercept,
                               exclude = which(ada.w[,i] == Inf),               # exclude variables with Inf weights
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
# univar.lasso()
#-------------------------------------------------------------------------------
# LASSO/Adaptive LASSO applied to a single equation
univar.lasso <- function(y, x, p = 1, alpha = 1, adaptive = c('none','lasso','ridge'),
                      crit = c('bic','ebic','aic','aicc','hqc'), nlambda = 100,
                      standardize = TRUE, intercept = TRUE, ...) {
  
  adaptive <- match.arg(adaptive)
  crit <- match.arg(crit)
  
  # 1st step
  if(adaptive == 'lasso'){
    eq <- ic.glmnet(x = x, y = y, crit = crit, alpha = 1,
                    standardize = standardize, intercept = intercept)
    eq.coef <- eq$coefficients[-1]
    ada.w	<- abs(eq.coef)^(-1)                                                  # adaptive weights
  } else if (adaptive == 'ridge'){
    eq <- ic.glmnet(x = x, y = y, crit = crit, alpha = 0,
                    standardize = standardize, intercept = intercept)
    eq.coef <- eq$coefficients[-1]
    ada.w	<- abs(eq.coef)^(-1)                                                  # adaptive weights
  } else {
    ada.w <- NULL
  }
  
  # 2nd step
  if (is.null(ada.w)) {
    alasso <- ic.glmnet(x = x, y = y, crit = crit, alpha = alpha,               # plain LASSO/Ridge
                        standardize = standardize, intercept = intercept, ...)
  } else {                                                                      # adaptive LASSO/Ridge
    if (sum(ada.w == Inf) == length(ada.w)) {
      warning('Adaptive weights all equal to Inf, all candidated excluded in the first step. Adaptive Lasso not estimated.')
      alasso              <- list()
      alasso$rss		      <- sum((y - mean(y))^2)
      alasso$coefficients <- rbind(mean(y), matrix(0, nrow = ncol(x), ncol = 1))
      alasso$lambda       <- 0
      alasso$ic           <- c(0)
      alasso$lambda.max   <- NA
      alasso$lambda.index <- NA
    } else{
      wpen                     <- ada.w
      wpen[which(wpen == Inf)] <- 0                                             # variables with Inf weights are manually removed
      alasso <- ic.glmnet(x = x, y = y, crit = crit, alpha = alpha,
                          standardize = standardize, intercept = intercept,
                          exclude = which(ada.w == Inf),                        # exclude variables with Inf weights
                          penalty.factor = wpen, ...)
    }
  }
  
  # collect results
  alasso.est	<- list('x' = x, 'y' = y, 'coefficients' = alasso$coefficients, 'ada.w' = ada.w,
                     'rss' = alasso$rss, 'lambda' = alasso$lambda, 'lambda.max' = alasso$lambda.max,
                     'lambda.index' = alasso$lambda.index
                     )

  return(alasso.est)
}

#-------------------------------------------------------------------------------
# var.scad()
#-------------------------------------------------------------------------------
# SCAD applied to a VAR system of equations
var.scad <- function(y, x, p = 1, crit = c('bic','ebic','aic','aicc','hqc'),
                     alpha = 1, nlambda = 100, family = 'gaussian',
                     epsilon = 1e-4, thresh = 1e-9, max_itt = 1e6,...) {
  K <- ncol(x)
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

#-------------------------------------------------------------------------------
# transform_fredmd()
#-------------------------------------------------------------------------------
# Do transformations to FRED-MD using the t-code
transform_fredmd <- function(file = '', date_start = NULL, date_end = NULL, transform = TRUE, skip = 2) {
  # error checking
  if (!is.logical(transform))
    stop("'transform' must be logical.")
  if ((class(date_start) != 'Date') && (!is.null(date_start)))
    stop("'date_start' must be Date or NULL.")
  if ((class(date_end) != 'Date') && (!is.null(date_end)))
    stop("'date_end' must be Date or NULL.")
  
  if (class(date_start) == 'Date') {
    if (as.numeric(format(date_start, '%d')) != 1)
      stop("'date_start' must be Date whose day is 1.")
    if (date_start < as.Date('1959-01-01'))
      stop("'date_start' must be later than 1959-01-01.")
  }
  
  if (class(date_end) == 'Date') {
    if (as.numeric(format(date_end, '%d')) != 1)
      stop("'date_end' must be Date whose day is 1.")
  }
  
  # prepare raw data
  rawdata <- readr::read_csv(file, col_names = FALSE, col_types = cols(X1 = col_date(format = '%m/%d/%Y')),
                             skip = skip)
  
  rawdata <- as.data.frame(rawdata)
  row_to_remove = c()
  for (row in (nrow(rawdata)-20):nrow(rawdata)){
    if(!any(is.finite(unlist(rawdata[row, ])))){
      row_to_remove = c(row_to_remove,row)
    }
  }
  if(length(row_to_remove)>0){
    rawdata = rawdata[-row_to_remove,]
  }
  
  attrdata          <- utils::read.csv(file, header = FALSE, nrows = skip)
  header            <- c('date', unlist(attrdata[1,2:ncol(attrdata)]))
  colnames(rawdata) <- header
  
  # store transformation codes as tcode
  tcode <- unlist(attrdata[skip,2:ncol(attrdata)])
  
  # subfunction transxf: data transformation based on tcodes
  transxf <- function(x, tcode) {
    n <- length(x)                                                              # number of observations (including missing values)
    small <- 1e-06                                                              # value close to zero
    
    # allocate output variable
    y  <- rep(NA, n)
    y1 <- rep(NA, n)
    
    # transformation
    if (tcode == 1) {
      y <- x                                                                    # case 1: Level (i.e. no transformation): x(t)
    } else if (tcode == 2) {
      y[2:n] <- (x[2:n] - x[1:(n - 1)])                                         # case 2: First difference: x(t)-x(t-1)
    } else if (tcode == 3) {
      y[3:n] <- (x[3:n] - 2 * x[2:(n - 1)] + x[1:(n - 2)])                      # case 3: Second difference: (x(t)-x(t-1))-(x(t-1)-x(t-2))
    } else if (tcode == 4) {
      if (min(x, na.rm = TRUE) > small) y <- log(x)                             # case 4: Natural log: ln(x)
    } else if (tcode == 5) {                                                    # case 5: First difference of natural log: ln(x)-ln(x-1)
      if (min(x, na.rm = TRUE) > small) {
        x      <- log(x)
        y[2:n] <- (x[2:n] - x[1:(n - 1)])
      }
    } else if (tcode == 6) {                                                    # case 6: Second difference of natural log:
      if (min(x, na.rm = TRUE) > small) {                                       # (ln(x)-ln(x-1))-(ln(x-1)-ln(x-2))
        x      <- log(x)
        y[3:n] <- (x[3:n] - 2 * x[2:(n - 1)] + x[1:(n - 2)])
      }
    } else if (tcode == 7) {                                                    # case 7: First difference of percent change:
      y1[2:n] <- (x[2:n] - x[1:(n - 1)])/x[1:(n - 1)]                           # (x(t)/x(t-1)-1)-(x(t-1)/x(t-2)-1)
      y[3:n]  <- (y1[3:n] - y1[2:(n - 1)])
    }
    return(y)
  }
  
  # fill single NA values by the mean of the previous and following observations
  for (j in 2:ncol(rawdata)){
    for (i in 2:nrow(rawdata)){
      if (is.na(rawdata[i,j]) & !is.na(rawdata[i-1,j]) & !is.na(rawdata[i+1,j])){
        rawdata[i,j] <- (as.numeric(rawdata[i-1,j]) + as.numeric(rawdata[i+1,j])) / 2
      }
    }
  }
  
  # transform data
  if (transform) {
    N          <- ncol(rawdata)
    data       <- rawdata
    data[,2:N] <- NA
    # perform transformation using subfunction transxf()
    for (i in 2:N) {
      temp <- transxf(rawdata[, i], tcode[i - 1])
      data[, i] <- temp
    }
  } else {
    data <- rawdata
  }
  
  # Null case of date_start and date_end
  if (is.null(date_start))
    date_start <- as.Date('1959-01-01')
  if (is.null(date_end))
    date_end <- data[, 1][nrow(data)]
  
  # subset data
  index_start <- which.max(data[,1] == date_start)
  index_end   <- which.max(data[,1] == date_end)
  
  outdata <- data[index_start:index_end,]
  class(outdata) <- c('data.frame','fredmd')
  return(outdata)
}

#-------------------------------------------------------------------------------
# MCSprocedure()
#-------------------------------------------------------------------------------
setClass("SSM", representation(Info = "list", Bootstrap = "list",
                               show = "matrix"))
setMethod("show", "SSM", function(object) {
  cat(paste("\n------------------------------------------"))
  cat(paste("\n-          Superior Set of Models        -"))
  cat(paste("\n------------------------------------------\n"))
  show(object@show)
  cat(paste("\nDetails"))
  cat(paste("\n------------------------------------------\n"))
  cat(paste("\nNumber of eliminated models\t:\t"))
  cat(object@Info$n_elim)
  cat(paste("\nStatistic\t:\t"))
  cat(object@Info$statistic)
  cat(paste("\nElapsed Time\t:\t"))
  print(object@Info$elapsed.time)
})

boot.block <- function(x, v, n, k) {
  startIndexes = sample(1:(n - k), v + 1)
  blocks = do.call(c, lapply(startIndexes, function(p, k) p +
                               seq_len(k), k = k))
  return(blocks[1:n])
}

d_b_i_mean.fun <- function(x, M, model.names) {
  export = sapply(1:M, function(l) sum(x[paste(model.names[l],
                                               model.names[-l], sep = ".")])/(M - 1))
  return(export)
}

LossVaR <- function(realized, evaluated, which = "asymmetricLoss",
                    type = "normal", delta = 25, tau) {
  if (delta <= 0) {
    stop("delta must be greater than zero")
  }
  if (which != "asymmetricLoss") {
    stop("Not valid choice of \"which\". Valid choices are \"asymmetricLoss\"")
  }
  if (which == "asymmetricLoss" & all(type != c("normal", "differentiable"))) {
    stop("Not valid choice of \"type\". Valid choices are \"normal\" and \"differentiable\" ")
  }
  
  if (tau < 0 || tau > 1) {
    stop("tau must be in (0,1)")
  }
  
  if (is.matrix(evaluated)) {
    m = ncol(evaluated)
  } else {
    m = 1
  }
  if (!is.matrix(evaluated)) {
    evaluated = as.matrix(evaluated)
  }
  if (!is.matrix(realized)) {
    realized = as.matrix(realized)
  }
  if (ncol(realized) > 1) {
    stop("realized must of one column")
  }
  if (nrow(realized) != nrow(evaluated)) {
    stop("the number of realized and evaluated observation must be the same")
  }
  
  if (which == "asymmetricLoss") {
    if (type == "normal") {
      dummies = matrix(0, nrow = length(realized), ncol = m)
      for (j in 1:m) {
        dummies[realized < evaluated[, j], j] = 1
      }
      loss = (tau - dummies) * (realized - evaluated)
    }
    if (type == "differentiable") {
      loss = (tau - (1 + exp(delta * (realized - evaluated)))^(-1)) *
        (realized - evaluated)
    }
  }
  return(loss)
}

LossVol <- function(realized, evaluated, which = "SE1") {
  
  if (!which %in% c("SE1", "SE2", "QLIKE", "R2LOG", "AE1",
                    "AE2")) {
    stop("Not valid choice of \"which\". Valid choices are \"SE1\",\"SE2\",\"QLIKE\",\"R2LOG\",\"AE1\",\"AE2\" ")
  }
  
  if (is.matrix(evaluated)) {
    m = ncol(evaluated)
  } else {
    m = 1
  }
  if (!is.matrix(evaluated)) {
    evaluated = as.matrix(evaluated)
  }
  if (!is.matrix(realized)) {
    realized = as.matrix(realized)
  }
  
  if (ncol(realized) > 1) {
    stop("realized must of one column")
  }
  if (nrow(realized) != nrow(evaluated)) {
    stop("the number of realized and evaluated observation must be the same")
  }
  
  if (which == "SE1") {
    loss = (evaluated - realized)^2
  }
  if (which == "SE2") {
    loss = (evaluated^2 - realized^2)^2
  }
  if (which == "QLIKE") {
    loss = log(evaluated^2) + realized^2 * evaluated^(-2)
  }
  if (which == "R2LOG") {
    loss = (log(realized^2 * evaluated^-2))^2
  }
  if (which == "AE1") {
    loss = abs(evaluated - realized)
  }
  if (which == "AE2") {
    loss = abs(evaluated^2 - realized^2)
  }
  return(loss)
}

LossLevel <- function(realized, evaluated, which = "SE") {
  
  if (all(which != c("SE", "AE")))
    stop("Not valid choice of \"which\". Valid choices are \"SE\",\"AE\" ")
  
  if (is.matrix(evaluated)) {
    m = ncol(evaluated)
  } else {
    m = 1
  }
  if (!is.matrix(evaluated)) {
    evaluated = as.matrix(evaluated)
  }
  if (!is.matrix(realized)) {
    realized = as.matrix(realized)
  }
  if (ncol(realized) > 1)
    stop("realized must of one column")
  if (nrow(realized) != nrow(evaluated))
    stop("the number of realized and evaluated observation must be the same")
  
  if (which == "SE") {
    loss = (evaluated - realized)^2
  }
  if (which == "AE") {
    loss = abs(evaluated - realized)
  }
  return(loss)
}

MCSprocedure <- function(Loss, alpha = 0.15, B = 5000, cl = NULL,
                         ram.allocation = TRUE, statistic = "Tmax", k = NULL, min.k = 3,
                         verbose = TRUE) {
  
  time.start = Sys.time()
  Loss = as.matrix(Loss)
  M_start = ncol(Loss)
  
  # no points in colnames(Loss)
  colnames(Loss) = gsub(".", "_", colnames(Loss), fixed = T)
  if (is.null(colnames(Loss))) {
    colnames(Loss) = paste("model", 1:ncol(Loss), sep = "_")
  }
  if (is.null(cl)) {
    cl_number = 0
  }
  if (!is.null(cl)) {
    max.cores = length(cl)
  }
  
  B = round(B)
  if (B < 1000) {
    cat(paste("Warning: B is small"))
  }
  
  if (any(is.na(Loss))) {
    stop("NAs in Loss are not allowed")
  }
  if (any(abs(Loss) == Inf)) {
    stop("Inf in Loss are not allowed")
  }
  
  N = nrow(Loss)
  
  repeat {
    
    M = ncol(Loss)
    model.names = colnames(Loss)
    
    col.names.d = do.call(c,
                          lapply(1:M, function(x) {
                            paste(model.names[x], model.names[-x], sep = ".")
                          }))
    
    d = do.call(cbind,
                lapply(1:M, function(x) {
                  Loss[, x] - Loss[, -x]
                }))
    
    colnames(d) = col.names.d
    
    d_ij_mean = colMeans(d)
    
    d_i_mean = sapply(1:M, function(x) {
      sum(d_ij_mean[paste(model.names[x], model.names[-x], sep = ".")])/(M - 1)
    })
    
    names(d_i_mean) = model.names
    block.length = NULL
    foo = expand.grid(1:M, 1:M)[, 2:1]
    foo = foo[foo[, 1] != foo[, 2], ]
    index = col.names.d[foo[, 1] < foo[, 2]]
    
    if (ncol(d) > 2) {
      
      d_cut = as.list(as.data.frame(d[, index]))
      
      if (!is.null(cl)) {
        k = max(na.omit(as.numeric(parSapply(cl, d_cut, function(x) {
          try(ar(x)$order, silent = TRUE)
        }))))
      }
      
      if (is.null(cl)) {
        k = max(na.omit(as.numeric(sapply(d_cut, function(x) {
          try(ar(x)$order, silent = TRUE)
        }))))
      }
      
    } else {
      k = ar(d[, 1])$order
    }
    
    if (k < min.k) {
      k = min.k
    }
    
    v = ceiling(nrow(d)/k)
    n = nrow(d)
    
    if (!is.null(cl)) {
      indexes_b = parLapply(cl, 1:B, boot.block, v = v, n = n, k = k)
    }
    if (is.null(cl)) {
      indexes_b = lapply(1:B, boot.block, v = v, n = n, k = k)
    }
    
    if (!is.null(cl)) {
      if (ram.allocation) {
        weigth.in.cl = (as.numeric(object.size(d))/1e+06) *
          3.2
        
        used.memory = sum(sapply(ls(), function(x) object.size(get(x))))/1e+06
        
        cl_number = min(as.numeric(ceiling((memory.limit() * 0.8 - used.memory)/weigth.in.cl) - 1), max.cores)
        
      } else {
        cl_number = max.cores
      }
      
    } else {
      cl_number = 1
    }
    
    d = as.data.frame(d)
    
    if (cl_number > 1) {
      d_ij_avg_resampled = parLapply(cl[1:cl_number], indexes_b,
                                     function(x, d, N) {
                                       colSums((d[x, ,drop = FALSE]))/N
                                     }, d = d, N = N)
    } else {
      d_ij_avg_resampled = lapply(indexes_b, function(x, d, N) {
        colSums((d[x, ,drop = FALSE]))/N
      }, d = d, N = N)
    }
    
    if (!is.null(cl)) {
      d_b_i_mean = do.call(rbind, parLapply(cl, d_ij_avg_resampled, d_b_i_mean.fun, M = M, model.names = model.names))
    }
    if (is.null(cl)) {
      d_b_i_mean = do.call(rbind, lapply(d_ij_avg_resampled, d_b_i_mean.fun, M = M, model.names = model.names))
    }
    
    d_b_i_var = colSums(t(apply(d_b_i_mean, 1, function(z) (z - d_i_mean)^2)))/B
    names(d_b_i_var) = model.names
    d_ij_avg_resampled = do.call(rbind, d_ij_avg_resampled)
    d_var_ij = matrix(colSums(t(apply(d_ij_avg_resampled, 1, function(z) (z - d_ij_mean)^2)))/B, nrow = 1)
    colnames(d_var_ij) = names(d_ij_mean)
    
    TR = max((abs(d_ij_mean))/(d_var_ij^0.5))
    TM = max(d_i_mean/(d_b_i_var^0.5))
    
    Tb_R = sapply(1:B, function(i) {
      max(abs(d_ij_avg_resampled[i, ] - d_ij_mean)/(d_var_ij^0.5))
    })
    
    Tb_M = sapply(1:B, function(i) {
      max((d_b_i_mean[i, ] - d_i_mean)/(d_b_i_var^0.5))
    })
    
    Pr = length(which(TR < Tb_R))/B
    Pm = length(which(TM < Tb_M))/B
    
    # elimination rules
    v_i_M = d_i_mean/(d_b_i_var^0.5)
    v_i_M_order = order(v_i_M, decreasing = TRUE)
    TR_all = ((d_ij_mean))/(d_var_ij^0.5)
    v_i_R = sapply(1:M, function(x) {
      max(TR_all[1, paste(model.names[x], model.names[-x], sep = ".")])
    })
    
    names(v_i_R) = model.names
    v_i_R_order = order(v_i_R, decreasing = TRUE)
    Pm_h0_mk = NULL
    
    for (i in 1:M) {
      if (i == 1) {
        model_temp = model.names
      } else {
        model_temp = model.names[-v_i_M_order[1:(i - 1)]]
      }
      TM_temp = max(d_i_mean[model_temp]/(d_b_i_var[model_temp]^0.5))
      Pm_h0_mk = c(Pm_h0_mk, length(which(TM_temp < Tb_M))/B)
    }
    
    mcs_Pm = sapply(1:M, function(x) max(Pm_h0_mk[1:x]))
    
    # For TR statistic
    combine.names = names(d_ij_mean)
    Pr_h0_mk = NULL
    
    for (i in 1:M) {
      remove = NULL
      if (i == 1) {
        model_temp = combine.names
      } else {
        
        remove = do.call(c, lapply(1:(i - 1), function(x) {
          which(gsub(model.names[v_i_R_order[x]], "", combine.names) != combine.names)
        }))
        
        model_temp = combine.names[-remove]
      }
      if (i < M) {
        TR_temp = max((abs(d_ij_mean[model_temp]))/(d_var_ij[1, model_temp]^0.5))
        Pr_h0_mk = c(Pr_h0_mk, length(which(TR_temp < Tb_R))/B)
      } else {
        Pr_h0_mk = c(Pr_h0_mk, 1)
      }
    }
    
    mcs_Pr = sapply(1:M, function(x) max(Pr_h0_mk[1:x]))
    matrix_show = matrix(NA, M, 7, dimnames = list(model.names,
                                                   c("Rank_M", "v_M", "MCS_M", "Rank_R", "v_R", "MCS_R",
                                                     "Loss")))
    
    matrix_show[, "v_M"] = v_i_M
    matrix_show[names(sort(v_i_M)), "Rank_M"] = 1:M
    matrix_show[model.names[v_i_M_order], "MCS_M"] = mcs_Pm
    matrix_show[, "v_R"] = v_i_R
    if (length(names(sort(v_i_R))) == M) {
      matrix_show[names(sort(v_i_R)), "Rank_R"] = 1:M
    }
    matrix_show[model.names[v_i_R_order], "MCS_R"] = mcs_Pr
    matrix_show[, "Loss"] = colMeans(Loss)
    
    rm(list = "indexes_b")
    
    if (!is.null(cl))
      clusterEvalQ(cl, {
        rm(list = ls())
        gc()
      })
    
    gc()
    
    if (statistic == "Tmax") {
      p2test = Pm
    }
    if (statistic == "TR") {
      p2test = Pr
    }
    
    if (p2test > alpha | all(d_var_ij == 0)) {
      if (verbose) {
        cat(paste("\n###########################################################################################################################\n"))
        cat(paste("Superior Set Model created\t:\n"))
        show(matrix_show)
        cat(paste("p-value\t:\n"))
        print(p2test)
        cat(paste("\n###########################################################################################################################"))
      }
      break
    } else {
      if (statistic == "Tmax")
        eliminate = which(v_i_M == max(v_i_M))
      if (statistic == "TR")
        eliminate = which(v_i_R == max(v_i_R))
      
      if (verbose)
        cat(paste("\nModel", model.names[eliminate],
                  "eliminated", Sys.time()))
      Loss = as.matrix(Loss[, -eliminate])
      colnames(Loss) = model.names[-eliminate]
      if (ncol(Loss) == 1) {
        if (verbose) {
          cat(paste("\n###########################################################################################################################\n"))
          cat(paste("Superior Set Model created\t:\n"))
          matrix_show = matrix(matrix_show[-eliminate, ], nrow = 1, dimnames = list(colnames(Loss),
                                                                                    colnames(matrix_show)))
          show(matrix_show)
          cat(paste("p-value\t:\n"))
          show(p2test)
          cat(paste("\n###########################################################################################################################"))
        }
        break
      }
    }
  }
  
  elapsed.time = Sys.time() - time.start
  n_elim = M_start - nrow(matrix_show)
  out = new("SSM", show = matrix_show, Info = list(model.names = rownames(matrix_show),
                                                   elapsed.time = elapsed.time, statistic = statistic, n_elim = n_elim,
                                                   mcs_pvalue = p2test, alpha = alpha, B = B, k = k), Bootstrap = list(TR = list(Stat = TR,
                                                                                                                                 BootDist = Tb_R), Tmax = list(Stat = TM, BootDist = Tb_M)))
  return(out)
}

#-------------------------------------------------------------------------------
# mcs.test()
#-------------------------------------------------------------------------------
mcs.test <- function(models, alpha, statistic){ 
  loss.level <- data.frame(models[[1]]^2)
  for(mod in 2:length(models)){
    loss.level <- data.frame(loss.level, models[[mod]]^2)
  }
  colnames(loss.level) <- names(models)
  invisible(
    capture.output(
      MCS <- MCSprocedure(Loss = loss.level,
                          alpha = alpha,
                          B = 5000,
                          statistic = statistic) %>% 
        suppressWarnings()
      )
    )
  
  return(MCS)
}