rm(list=ls())
packages <- c('lubridate','dplyr','glmnet','MASS','pracma','zoo','Matrix','foreach','doParallel','parallel','Rcpp','tsBARrcpp',
              'tsBARLSrcpp','magic','matrixStats','OpenMx')
for (pkg in packages){
  if (!requireNamespace(pkg, quietly=T)){
    install.packages(pkg, dependencies=T)
  }
  library(pkg, character.only=T)
}

################################################################################
#                                FUNCTIONS
################################################################################
#-------------------------------------------------------------------------------
#                               ic.glmnet
#-------------------------------------------------------------------------------
# to estimate lasso with info criteria tuning of the lambda parameter
ic.glmnet = function (x, y, crit=c('bic','ebic','aic','aicc','hqc'),...)
{
  if (is.matrix(x) == FALSE) {
    x = as.matrix(x)
  }
  if (is.vector(y) == FALSE) {
    y = as.vector(y)
  }
  crit=match.arg(crit)
  n=length(y)
  model = glmnet(x = x, y = y, ...)
  coef = coef(model)
  lambda = model$lambda
  df = model$df
  dev_ratio = model[['dev.ratio']]
  null_dev= model$nulldev
  
  yhat=cbind(1,x)%*%coef  #x*beta_hat
  residuals = (y- yhat)
  mse = colMeans(residuals^2)
  sse = colSums(residuals^2)
  
  #defining criteria
  nvar = df #+ 1
  bic = n*log(mse)+nvar*log(n)
  ebic= n*log(mse)+nvar*log(n)+(2*nvar*0.5*log(ncol(x)))
  aic = n*log(mse)+2*nvar
  aicc = aic+(2*nvar*(nvar+1))/(n-nvar-1)
  hqc = n*log(mse)+2*nvar*log(log(n))
  
  sst = (n-1)*var(y)
  r2 = 1 - (sse/sst) #r square
  adjr2 = (1 - (1 - r2) * (n - 1)/(nrow(x) - nvar - 1)) #adj r square
  
  crit=switch(crit,bic=bic,ebic=ebic,aic=aic,aicc=aicc,hqc=hqc)
  
  selected=best.model = which(crit == min(crit))
  
  ic=c(bic=bic[selected],ebic=ebic[selected],aic=aic[selected],aicc=aicc[selected],hqc=hqc[selected])
  
  result=list(coefficients=coef[,selected],ic=ic,lambda=lambda[selected],lambda.index=selected,
              nvar=nvar[selected],lambda.max=lambda[1],dev.ratio=dev_ratio[selected],
              null.dev=null_dev,glmnet=model,residuals=residuals[,selected],
              fitted.values=yhat[,selected],ic.range=crit,call=match.call())
  
  class(result)='ic.glmnet'
  return(result)
}

#-------------------------------------------------------------------------------
#                           var.lasso_univariate
#-------------------------------------------------------------------------------
var.lasso_univariate <- function(y, x, np, alpha, adaptive=c('none','ols','lasso','group','ridge'), 
                                 standardize=F, intercept=F, crit=c('bic','ebic','aic','aicc','hqc'), nlambda,...) {
  
  adaptive <- match.arg(adaptive)
  crit <- match.arg(crit)
  
  #=============================================================================
  # Adaptive weights
  #=============================================================================
  if(adaptive=='lasso'){
    # GLM regression
    eq <- ic.glmnet(x=x, y=y, crit=crit, alpha=1, standardize=standardize, intercept=intercept)
    eq.coef <- eq$coefficients[-1]
    # Weights for adaptive lasso
    ada.w	<- abs(eq.coef)^(-1) 
  } else if (adaptive=='ridge'){
    # GLM regression
    eq <- ic.glmnet(x=x, y=y, crit=crit, alpha=0, standardize=standardize, intercept=intercept)
    eq.coef <- eq$coefficients[-1]
    # Weights for adaptive lasso
    ada.w	<- abs(eq.coef)^(-1)
  } else {
    ada.w <- NULL
  }
  
  #=============================================================================
  # Final Lasso/Ridge
  #=============================================================================
  lasso.est	<- list('x'=x, 'y'=y, 'coefficients'=NULL, 'ada.w'=ada.w ,'rss'=NULL, 'lambda'=NULL, 'lambda.max'=NULL, 'lambda.index'=NULL)	
  # Simple Lasso/Ridge
  #---------------------------------------------------------------------------
  if(is.null(ada.w)){
    lasso <- ic.glmnet(x=x, y=y, crit=crit, alpha=alpha,...)
    # Collect estimates
    lasso.ic <- list('ic'=crit,'ic.all'=lasso$ic,'coefficients'=lasso$coefficients,'rss'=sum(lasso$residuals^2), 
                     'lambda'=lasso$lambda, 'lambda.max'=lasso$lambda.max, 'lambda.index'=lasso$lambda.index)
  } else {
    #---------------------------------------------------------------------------
    # Adaptive Lasso/Ridge
    #---------------------------------------------------------------------------
    # in case all the adaptive weights are equal to zero (empty model) 
    if(sum(ada.w==Inf)==length(ada.w)){
      warning('Adaptive weights all equal to Inf, all candidated excluded in the first step. Adaptive Lasso not estimated.')
      lasso <- list()
      lasso$rss		<- sum((y-mean(y))^2)
      lasso$coefficients	<- rbind(mean(y), matrix(0,nrow=ncol(x),ncol=1))
      lasso$lambda <- 0
      lasso$ic <- c(0)
      lasso$lambda.max <- NA
      lasso$lambda.index <- NA
    } else{
      wpen <- ada.w
      wpen[which(wpen==Inf)] <- 0 # variables with Inf weights are manually removed
      lasso <- ic.glmnet(x=x, y=y, crit=crit, alpha=alpha, exclude=which(ada.w==Inf), penalty.factor=wpen, standardize=standardize, intercept=intercept) # exclude variables with Inf weights
    }
    # Collect estimates of Adaptive Lasso/Ridge
    lasso.ic <- list('ic'=crit,'ic.all'=lasso$ic,'coefficients'=lasso$coefficients,'rss'=sum(lasso$residuals^2), 'lambda'=lasso$lambda, 
                     'lambda.max'=lasso$lambda.max, 'lambda.index'=lasso$lambda.index)
  }
  
  # Collect results of interest
  lasso.est$coefficients <- cbind(lasso.est$coefficients, lasso.ic$coefficients)
  lasso.est$rss <- cbind(lasso.est$rss, lasso.ic$rss)
  lasso.est$lambda	<- cbind(lasso.est$lambda, lasso.ic$lambda)
  lasso.est$lambda.max	<- cbind(lasso.est$lambda.max, lasso.ic$lambda.max)
  lasso.est$lambda.index	<- cbind(lasso.est$lambda.index, lasso.ic$lambda.index)
  
  return(lasso.est)
}

#-------------------------------------------------------------------------------
#                                 create_lags
#-------------------------------------------------------------------------------
# to create lags
create_lags <- function(y, p=1, include.original=TRUE, trim=TRUE) {
  x <- as.matrix(y)
  n <- nrow(x)
  k <- ncol(x)
  lx <- matrix(0, nrow = n, ncol = (p + include.original) * k)
  if (is.null(colnames(x))) {
    c.names <- rep('', k)
  } else {
    c.names <- colnames(x)
  }
  colnames(lx) <- rep(c.names, p + include.original)
  for (i in (1 - include.original):p) {
    cols <- k * (i - 1 + include.original) + 1:k
    lx[(1+i):n, cols] <- x[1:(n-i), ]
    colnames(lx)[cols] <- paste(c.names, ' l', i, sep = '')
  }
  return(lx[(1 + trim * p):n, , drop = FALSE])
}

#-------------------------------------------------------------------------------
#                                first_row
#-------------------------------------------------------------------------------
# Give the first row with non-NA values in all columns
first_row <- function(x) {
  x <- as.data.frame(x)
  n <- last_row(x)
  for (i in 1:n) {
    sub <- x[i, ]
    if (all(complete.cases(sub))) {
      return(i)
    }
  }
  return(NA)
}

#-------------------------------------------------------------------------------
#                                last_row
#-------------------------------------------------------------------------------
# Give the last row with non-NA values in all columns
last_row <- function(x){
  row <- nrow(x)
  for (i in 1:length(x)) {
    a <- which(!is.na(x[,i]))
    b <- max(a)
    row <- min(b, row)
  }
  return(row)
}

#-------------------------------------------------------------------------------
#                                   ic.scad
#-------------------------------------------------------------------------------
ic.scad <- function(y, x, crit=c('bic','ebic','aic','aicc','hqc'), alpha=1, nlambda=100, family='gaussian', epsilon=1e-4, thresh=1e-9, max_itt=1e6,...) {
  if (is.matrix(x) == FALSE) {
    x = as.matrix(x)
  }
  if (is.vector(y) == FALSE) {
    y = as.vector(y)
  }
  crit=match.arg(crit)
  n=length(y)
  model <- ncvreg(X=x, y=y, family=family, alpha=alpha, penalty="SCAD", ...)
  coef <- model$beta
  lambda <- model$lambda
  
  yhat=cbind(1,x)%*%coef  #x*beta_hat
  residuals = (y - yhat)
  mse = colMeans(residuals^2)
  sse = colSums(residuals^2)
  
  #defining criteria
  nvar = colSums(coef[-1,] != 0) #+ 1
  bic = n*log(mse)+nvar*log(n)
  ebic= n*log(mse)+nvar*log(n)+(2*nvar*0.5*log(ncol(x)))
  aic = n*log(mse)+2*nvar
  aicc = aic+(2*nvar*(nvar+1))/(n-nvar-1)
  hqc = n*log(mse)+2*nvar*log(log(n))
  
  sst = (n-1)*var(y)
  r2 = 1 - (sse/sst) #r square
  adjr2 = (1 - (1 - r2) * (n - 1)/(nrow(x) - nvar - 1)) #adj r square
  
  crit=switch(crit,bic=bic,ebic=ebic,aic=aic,aicc=aicc,hqc=hqc)
  
  selected=best.model = min(which(crit == min(crit)))
  
  ic=c(bic=bic[selected],ebic=ebic[selected],aic=aic[selected],aicc=aicc[selected],hqc=hqc[selected])
  
  result=list(coefficients=coef[,selected],ic=ic,lambda=lambda[selected],lambda.index=selected,
              nvar=nvar[selected],lambda.max=lambda[1],
              SCAD=model,residuals=residuals[,selected],
              fitted.values=yhat[,selected],ic.range=crit,call=match.call())
  
  class(result)='SCAD'
  return(result)
}

#-------------------------------------------------------------------------------
#                                   var.scad
#-------------------------------------------------------------------------------
var.scad <- function(y, x, p, lags=1, crit=c('bic','ebic','aic','aicc','hqc'), alpha=1, nlambda=100, family='gaussian', epsilon=1e-4, thresh=1e-9, max_itt=1e6,...) {
  scad.ic <- list()
  scad.est	<- list('x'=x, 'y'=y, 'coefficients'=NULL,'rss'=NULL, 'lambda'=NULL, 'lambda.max'=NULL, 'lambda.index'=NULL)	
  for(i in 1:p){
    yi <- y[,i]
    scad <- ic.scad(x, y=yi, crit=crit,...)
    scad.ic[[i]] <- list('ic'=crit,'ic.all'=scad$ic,'coefficients'=scad$coefficients,'rss'=sum(scad$residuals^2), 
                         'lambda'=scad$lambda, 'lambda.max'=scad$lambda.max, 'lambda.index'=scad$lambda.index)
  }
  
  for(i in 1:p){
    scad.est$coefficients <- cbind(scad.est$coefficients, scad.ic[[i]]$coefficients)
    scad.est$rss <- cbind(scad.est$rss, scad.ic[[i]]$rss)
    scad.est$lambda	<- cbind(scad.est$lambda, scad.ic[[i]]$lambda)
    scad.est$lambda.max	<- cbind(scad.est$lambda.max, scad.ic[[i]]$lambda.max)
    scad.est$lambda.index	<- cbind(scad.est$lambda.index, scad.ic[[i]]$lambda.index)
  }
  
  return(scad.est)
}

#-------------------------------------------------------------------------------
#                             m_lasso
#-------------------------------------------------------------------------------
m_lasso <- function(x, y, lambda_glmnet, upsilon = NULL,
                    full_path = FALSE, tol_glmnet = 1e-4) {
  if (missing(x) || missing(y) || missing(lambda_glmnet)) {
    stop("Not enough input arguments.")
  }
  pq <- ncol(x)  # columns in x = pq = total number of regressors
  p <- ncol(y)   # q = autoregressive order, p = dim(output)
  if (nrow(x) != nrow(y)) {
    stop("Number of observations (rows in x and y) do not match.")
  }
  if (is.null(upsilon)) {
    upsilon <- matrix(1, p, pq) # default to unit loadings
  }
  that <- matrix(NA, p, pq) # placeholder for estimates
  if (full_path == FALSE) { # if full path *not* requested (= default)...
    for (i in 1:p) {        # use glmnet w/ single lambda
      xtilde <- sweep(x, 2, upsilon[i, ], FUN = "/")  # rescale using loadings
      fit_i <- glmnet(xtilde, y[, i], family = "gaussian",
                      standardize = FALSE, intercept = FALSE,
                      lambda = lambda_glmnet, thresh = tol_glmnet)
      that[i, ] <- coef(fit_i)[-1] / upsilon[i, ] # original scaling
    }
  } else {            # if full path requested (not default)...
    for (i in 1:p) {  # use glmnet's lambda sequence and interpolate
      xtilde <- sweep(x, 2, upsilon[i, ], FUN = "/")  # rescale using loadings
      fit_i <- glmnet(xtilde, y[, i], family = "gaussian",
                      standardize = FALSE, intercept = FALSE,
                      thresh = tol_glmnet)
      # ^-- glmnet constructs its own lambda sequence
      that_i_temp <- coef(fit_i, s = lambda_glmnet)[-1] # interpolate
      that[i, ] <- that_i_temp / upsilon[i, ]           # original scale
    }
  }
  return(that) # return as matrix
}

#-------------------------------------------------------------------------------
#                             m_refit
#-------------------------------------------------------------------------------
m_refit <- function(x, y, that) {
  p <- ncol(y)                          # num responses
  sel <- that != 0                      # flag selections
  full_rank <- logical(p)               # full rank flag (initialize FALSE)
  num_selections <- rowSums(sel)        # precompute num selections
  for (i in which(num_selections > 0)) {
    xi <- as.matrix(x[, sel[i, ]])      # get active regressors
    fit_i <- ls_sol(xi, y[, i])         # refit
    that[i, sel[i, ]] <- fit_i$sol      # overwrite (keeping zeros)
    full_rank[i] <- !fit_i$solve_failed # interpret solve failure as low rank
  }
  full_rank[num_selections == 0] <- TRUE  # default full rank for no selections
  refit <- list(that = that, full_rank = full_rank)
  return(refit)
}

#-------------------------------------------------------------------------------
#                             ls_sol
#-------------------------------------------------------------------------------
ls_sol <- function(x, y) {
  # Try to compute the least squares solution using solve()
  tryCatch({
    sol <- solve(crossprod(x), crossprod(x, y))
    solve_failed <- FALSE
    return(list(sol = sol, solve_failed = solve_failed))
  }, error = function(e) {
    # If solve fails, compute Moore-Penrose solution and flag failure
    svd_x <- svd(x)
    d <- svd_x$d  # singular values
    tol <- max(dim(x)) * .Machine$double.eps      # scale free tolerance
    d_inv <- ifelse(d > tol * max(d), 1 / d, 0)   # inverse singular values
    x_plus <- svd_x$v %*% (d_inv * t(svd_x$u))    # Moore-Penrose inverse
    sol <- x_plus %*% y                           # Moore-Penrose solution
    solve_failed <- TRUE
    return(list(sol = sol, solve_failed = solve_failed))
  })
}

#-------------------------------------------------------------------------------
#                             var.lasso.data_driven_lambda
#-------------------------------------------------------------------------------
var.lasso.data_driven_lambda <- function(x, y, post = TRUE, intercept = TRUE, c = 1.1,
                                         gamma = NULL, k = 15, upsilon = NULL, tol_ups = 1e-3, warn = TRUE,
                                         full_path = FALSE, tol_glmnet = 1e-4) {
  if (missing(x) || missing(y)) {
    stop("Not enough inputs. Provide both x and y.")
  }
  if (nrow(x) != nrow(y)) {
    stop("Number of observations (rows in x and y) do not match.")
  }
  n <- nrow(x)      # number of observations
  nx <- ncol(x)     # number of predictor variables
  ny <- ncol(y)     # number of response variables
  
  # Set default probability tolerance if not provided
  if (is.null(gamma)) {
    gamma <- 0.1 / log(max(c(n, nx)))
  }
  
  # Penalty level based on Belloni, Chernozhukov, Chen, Hansen (2012 ECTA)
  lambda_star <- 2 * c * sqrt(n) * qnorm(1 - gamma / (2 * nx * ny))
  lambda_glmnet <- lambda_star / (2 * n) # penalty for glmnet
  
  # If intercepts are requested, demean the data before proceeding
  if (intercept) {
    ybar <- colMeans(y)     # mean of response variables
    xbar <- colMeans(x)     # mean of predictor variables
    y <- sweep(y, 2, ybar)  # demean responses
    x <- sweep(x, 2, xbar)  # demean predictors
    ybar <- as.matrix(ybar) # response means as matrix
    xbar <- as.matrix(xbar) # predictor means as matrix
  }
  # Note: Means are stored to back out intercepts later
  
  # Initial penalty loadings and estimates
  ups_init <- sqrt((1 / n) * crossprod(y^2, x^2))
  
  that_init <- m_lasso(x, y, lambda_glmnet, upsilon = ups_init,
                       full_path = full_path, tol_glmnet = tol_glmnet) # var.lasso(x, y, p=p, lambda=lambda_glmnet, nlambda = 1, thresh = tol_glmnet)
  if (post) {
    refit <- m_refit(x, y, that_init) # refit initial estimates
    that_init <- refit$that
    full_rank_post_init <- refit$full_rank
  } else {
    full_rank_post_init <- NULL
  }
  if (intercept) {
    intr_init <- ybar - that_init %*% xbar # back out intercepts
  } else {
    intr_init <- NULL
  }
  
  # Update penalty loadings and estimates iteratively
  ups_refi <- array(NA, dim = c(ny, nx, k))  # refined penalty loadings
  that_refi <- array(NA, dim = c(ny, nx, k)) # refined estimates
  rel_diffs_ups <- numeric(k)                # relative changes in loadings
  if (post) {
    full_rank_post_refi <- matrix(NA, ny, k) # full rank flags
  } else {
    full_rank_post_refi <- NULL
  }
  if (intercept) {
    intr_refi <- matrix(NA, ny, k) # placeholder for intercepts
  } else {
    intr_refi <- NULL
  }
  for (l in 1:k) {
    if (l == 1) {
      that_old <- that_init # use initial estimates for first iteration
    } else {
      that_old <- that_refi[, , l - 1] # use previous estimates for others
    }
    res_old <- y - x %*% t(that_old) # calculate residuals
    ups_new <- sqrt((1 / n) * crossprod(res_old^2, x^2)) # update penalty load.
    ups_refi[, , l] <- ups_new # store updated penalty loadings
    that_new <- m_lasso(x, y, lambda_glmnet, ups_new, tol_glmnet) # new est.
    if (post) {
      refit <- m_refit(x, y, that_new) # refit new estimates
      that_new <- refit$that
      full_rank_post_refi[, l] <- refit$full_rank
    }
    that_refi[, , l] <- that_new # store new estimates
    if (intercept) {
      intr_refi[, l] <- ybar - that_new %*% xbar # back out intercepts
    }
    if (l == 1) {
      ups_old <- ups_init # compare with initial penalty loadings
    } else {
      ups_old <- ups_refi[, , l - 1] # compare with previous penalty loadings
    }
    diff_ups <- ups_new - ups_old # change in penalty loadings
    rel_diff_ups <- sqrt(sum(diff_ups^2)) /
      (sqrt(sum(ups_old^2)) + .Machine$double.eps) # relative change
    rel_diffs_ups[l] <- rel_diff_ups # store relative change
    if (rel_diff_ups <= tol_ups) {
      break # stop if change is small enough
    }
  }
  ups <- ups_refi[, , l] # set refined penalty loadings
  if (intercept) {
    intr <- intr_refi[, l] # set intercepts
  } else {
    intr <- NULL
  }
  that <- that_refi[, , l] # set estimates
  
  # Optional warning if maximum updates reached without convergence
  if (warn && l == k && rel_diff_ups > tol_ups) {
    warning("Maximum number of updates reached. Consider increasing K.")
    cat(sprintf("Relative change in penalty loadings is %3.1g percent\n",
                100 * rel_diff_ups))
  }
  
  # Return results as a list
  fit <- list(
    lambda = lambda_star,
    ups_init = ups_init,
    ups_refi = ups_refi[, , 1:l],
    ups = ups,
    intr_init = intr_init,
    intr_refi = intr_refi[, 1:l],
    intr = intr,
    that_init = that_init,
    that_refi = that_refi[, , 1:l],
    that = that,
    full_rank_post_init = full_rank_post_init,
    full_rank_post_refi = full_rank_post_refi[, 1:l],
    k_term = l,
    rel_diff_ups_term = rel_diff_ups,
    rel_diffs_ups = rel_diffs_ups[1:l]
  )
  return(fit)
}

#-------------------------------------------------------------------------------
#                             sqrt_lasso
#-------------------------------------------------------------------------------
sqrt_lasso <- function(x, y, lambda, max_iter = 100, rel_tol_norm = 1e-4) {
  # Compute necessary dimensions
  n <- nrow(x)
  p <- ncol(x)
  # Ridge matrix for the initial beta estimate
  ridge_matrix <- diag(x = lambda, nrow = p)
  beta <- solve(crossprod(x) + ridge_matrix, crossprod(x, y))
  # Precompute quantities for the algorithm
  iter <- 0
  xx <- crossprod(x) / n
  xy <- crossprod(x, y) / n
  xx_diag <- diag(xx)
  x_beta <- x %*% beta
  error <- y - x_beta
  qhat <- mean(error^2)
  # Small constant to avoid division by zero
  epsilon <- .Machine$double.eps
  while (iter < max_iter) {
    iter <- iter + 1
    beta_old <- beta
    # Coordinate descent update for each beta_j
    for (j in 1:p) {
      s0 <- xx[j, ] %*% beta - xx_diag[j] * beta[j] - xy[j]
      beta_j_old <- beta[j]
      if (n^2 < (lambda^2) / xx_diag[j]) {
        beta[j] <- 0
      } else if (s0 > (lambda / n) * sqrt(qhat)) {
        beta[j] <- ((lambda / sqrt(n^2 - lambda^2 / xx_diag[j])) *
                      sqrt(max(qhat - (s0^2 / xx_diag[j]), 0)) - s0) /
          xx_diag[j]
      } else if (s0 < -(lambda / n) * sqrt(qhat)) {
        beta[j] <- (-(lambda / sqrt(n^2 - lambda^2 / xx_diag[j])) *
                      sqrt(max(qhat - (s0^2 / xx_diag[j]), 0)) - s0) /
          xx_diag[j]
      } else {
        beta[j] <- 0
      }
      if (beta[j] != beta_j_old) { # Update mean-square error if change
        x_beta <- x_beta + x[, j] * (beta[j] - beta_j_old)
        error <- y - x_beta
        qhat <- mean(error^2)
      }
    }
    # Check for convergence using the relative change in the Euclidean norm
    norm_diff <- norm(beta - beta_old, type = "2")          # numerator
    norm_beta_old <- norm(beta_old, type = "2") + epsilon   # denominator
    if (norm_diff / norm_beta_old < rel_tol_norm) {         # check tolerance
      break                                                 # break if met
    }
  }
  return(beta)
}

#-------------------------------------------------------------------------------
#                             m_sqrt_lasso
#-------------------------------------------------------------------------------
m_sqrt_lasso <- function(x, y, lambda, upsilon = NULL,
                         max_iter = 100, rel_tol_norm = 1e-4) {
  # Check input validity
  if (missing(x) || missing(y) || missing(lambda)) {
    stop("Not enough input arguments.")
  }
  pq <- ncol(x)  # columns in x = pq = total number of regressors
  p <- ncol(y)   # q = autoregressive order, p = dim(output)
  if (nrow(x) != nrow(y)) {
    stop("Number of observations (rows in x and y) do not match.")
  }
  # If no penalty loadings provided, default to sqrt(column means of x^2)
  if (is.null(upsilon)) {
    # upsilon <- matrix(1, p, pq) # default to unit loadings
    upsilon <- sqrt(matrix(colMeans(x^2), p, pq, byrow = TRUE))
  }
  if (nrow(upsilon) != p || ncol(upsilon) != pq) {
    stop("Invalid dimensions for upsilon.")
  }
  that <- matrix(NA, p, pq) # placeholder for estimates
  for (i in 1:p) {
    # Rescale x for the i-th equation
    xtilde_i <- sweep(x, 2, upsilon[i, ], FUN = "/")
    # Perform sqrt_lasso estimation with max_iter and rel_tol_norm
    that_i_temp <- sqrt_lasso(xtilde_i, y[, i], lambda,
                              max_iter = max_iter,
                              rel_tol_norm = rel_tol_norm)
    # Store estimates back to original scale
    that[i, ] <- that_i_temp / upsilon[i, ]
  }
  return(that) # return as matrix
}

#-------------------------------------------------------------------------------
#                             var.sqrt_lasso.data_driven_lambda
#-------------------------------------------------------------------------------
var.sqrt_lasso.data_driven_lambda <- function(x, y, q = 1, post = TRUE, intercept = TRUE,
                                              c = 1.1, gamma = NULL, upsilon = NULL,
                                              max_iter = 100, rel_tol_norm = 1e-4) {
  
  # INITIALIZE
  # Unpack data and construct predictors and response
  n <- nrow(x)  # n plus q
  p <- ncol(y)       # p = dim of output
  pq <- p * q           # pq = total number of variables
  
  # Probability tolerance
  if (is.null(gamma)) {
    gamma <- 0.1 / log(max(c(n, pq)))
  }
  # Penalty level in eyes of square-root LASSO
  lambda_star_sqrt <- c * sqrt(n) * qnorm(1 - gamma / (2 * p * pq))
  # Note: Sqrt-LASSO objective cancels out the "2" in the KKT conditions
  # If intercepts requested, demean before proceeding
  if (intercept) {
    ybar <- colMeans(y)     # response means as p x 1 vector
    xbar <- colMeans(x)     # predictor means as pq x 1 vector
    y <- sweep(y, 2, ybar)  # demean responses
    x <- sweep(x, 2, xbar)  # demean predictors
  }
  # Note: Means are stored to back out intercepts later
  # If no penalty loadings provided, default to sqrt(column means of x^2)
  if (is.null(upsilon)) {
    # upsilon <- matrix(1, p, pq) # default to unit loadings
    upsilon <- sqrt(matrix(colMeans(x^2), p, pq, byrow = TRUE))
  }
  if (nrow(upsilon) != p || ncol(upsilon) != pq) {
    stop("Invalid dimensions for upsilon.")
  }
  # ESTIMATE
  # Estimate VAR model using square-root LASSO
  that <- m_sqrt_lasso(x, y, lambda = lambda_star_sqrt,
                       upsilon = upsilon,
                       max_iter = max_iter,
                       rel_tol_norm = rel_tol_norm)
  # Back out intercepts if requested
  if (intercept) {
    intr <- ybar - that %*% xbar
  } else {
    intr <- NULL
  }
  # POST-PROCESS if requested
  if (post) {
    refit <- m_refit(x, y, that)
    that_post <- refit$that
    full_rank_post <- refit$full_rank
    # Back out intercepts (based on refitted coefficients) if requested
    if (intercept) {
      intr_post <- ybar - that_post %*% xbar
    } else {
      intr_post <- NULL
    }
  } else {
    intr_post <- NULL
    that_post <- NULL
    full_rank_post <- NULL
  }
  # RETURN
  return(list(intr = intr, that = that,
              intr_post = intr_post, that_post = that_post,
              full_rank_post = full_rank_post))
}


################################################################################
#                             FRED-MD DATASET
################################################################################
dataset <- readRDS('/home/prifti/My_scratch/Code/FRED/dataset_FRED_MD.rds')

################################################################################
#                             MODEL FORECASTING
################################################################################
start_oos <- as.yearmon('1990-01')
end_oos <- as.yearmon('2025-03')
# start_oos <- as.yearqtr('1990 Q1')
# end_oos <- as.yearqtr('2025 Q1')
freq_time <- 12
max_itt <- 100000
thresh <- 1e-9

params <- expand.grid(
  variable = c('CPIAUCSL','INDPRO','UNRATE'),# c('GDPC1'),
  lags  = c(1),
  lagy = c(4),
  rollw = c(360,240,120,0),
  h     = c(1,3,6,12),
  stringsAsFactors=F)

final_results <- data.frame()
for (ix in 1:nrow(params)) {
  current <- params[ix,]
  variable <- current$variable
  h <- current$h
  rollw <- current$rollw
  lags <- current$lags
  lagy <- current$lagy
  row_trim <- max(lagy, lags)+h-1
  step_back <- rollw/freq_time
  info_parallel <- c('h','freq_time','lags','lagy','start_oos','end_oos','dataset','rollw','step_back','row_trim')
  cl <- makeCluster(48)
  registerDoParallel(cl)
  clusterExport(cl, varlist=info_parallel)
  clusterEvalQ(cl, {source('/home/prifti/My_scratch/Code/FRED/functions_EA2_parallel.R')}) 
  fcast_result <- foreach(t=seq(start_oos,end_oos,by=1/freq_time), .combine=rbind, .export=info_parallel) %dopar% {
    if(rollw==0) {
      tmp.data <- dataset[dataset$date<=t, !(names(dataset) %in% c("date"))]
    } else {
      tmp.data <- dataset[dataset$date<=t & dataset$date>=t-step_back, !(names(dataset) %in% c("date"))]
    }
    y_lags <- create_lags(tmp.data[variable], p=lagy+h-1, include.original=F, trim=F) %>% tail(-row_trim)
    x_lags <- create_lags(tmp.data[,!(names(tmp.data) %in% c(variable))], p=lags+h-1, include.original=F, trim=F) %>% tail(-row_trim)
    
    if (h>1) {
      y_lags <- y_lags[,-(1:(h-1))]
      x_lags <- x_lags[,-(1:((h-1)*(ncol(tmp.data)-1)))]
    }
    
    # Variable y
    y_lag <- y_lags %>% head(-h)
    y <- tmp.data[variable] %>% tail(-row_trim) %>% head(-h) %>% as.matrix()
    T_ <- length(y)
    ym <- mean(y)
    ys <- sqrt(sum((y - ym)^2)/(T_))
    y_pred <- tmp.data[variable] %>% tail(1) %>% as.numeric()
    Y_c <- y - ym
    Y_cs <- Y_c / sqrt(sum(Y_c^2)/T_)
    
    # Variable x
    x_lag <- x_lags %>% head(-h)
    x <- cbind(y_lag, x_lag)
    p <- ncol(x)
    xm <- apply(x, 2, function(x) mean(x))
    xs <- apply(x, 2, function(x) sqrt(sum((x-mean(x))^2)/T_))
    x_pred <- cbind(y_lags, x_lags) %>% tail(1) %>% as.matrix()
    X_c <- apply(x, 2, function(x) x - mean(x))
    X_cs <- apply(X_c, 2, function(x) x / sqrt(sum(x^2)/T_))
    
    adaLASSO_ridge <- var.lasso_univariate(y=Y_cs, x=X_cs, np=p, alpha=1, adaptive='ridge', crit='bic', standardize=F, intercept=F)$coefficients[-1] * ys / xs
    adaLASSO_lasso <- var.lasso_univariate(y=Y_cs, x=X_cs, np=p, alpha=1, adaptive='lasso', crit='bic', standardize=F, intercept=F)$coefficients[-1] * ys / xs
    LASSO <- var.lasso_univariate(y=Y_cs, x=X_cs, np=p, alpha=1, adaptive='none', crit='bic', standardize=F, intercept=F)$coefficients[-1] * ys / xs
    #BAR_model <- BAR_CCD_Rcpp(y=Y_cs, x=X_cs, flg_lambda_min=T, c_lambda_min=0.01)
    # BAR_01_model <- BAR_CCD_Rcpp(y=Y_cs, x=X_cs, flg_lambda_min=T, c_lambda_min=0.1)
    # BAR_05_model <- BAR_CCD_Rcpp(y=Y_cs, x=X_cs, flg_lambda_min=T, c_lambda_min=0.5)
    #BAR <- BAR_model$coefficients * ys / xs
    # BAR_01 <- BAR_01_model$coefficients * ys / xs
    BAR_model <- BAR_CD(y=Y_cs, x=X_cs)
    BAR <- BAR_model$beta * ys / xs    
    #BAR_LS <- BAR_LS_Rcpp(y=Y_cs, x=X_cs, flg_lambda_min=F, svd=TRUE)$coefficients * ys / xs
    
    results <- data.frame()
    models <- list(#'BAR'=BAR,
      #'BAR_LS'=BAR_LS,
      'BAR'=BAR,
      'adaLASSO_ridge'=adaLASSO_ridge, 
      'adaLASSO_lasso'=adaLASSO_lasso, 
      'LASSO'=LASSO)
    
    for (model_name in names(models)) {
      mod <- models[[model_name]]
      a0 = ym - xm%*%mod
      fcast <- a0 + x_pred%*%mod
      err <- fcast - y_pred
      RMSFE <- err^2
      MAFE <- abs(err)
      nvar <- sum(mod != 0)
      
      results <- rbind(results,
                       data.frame('variable'=variable,
                                  'h'=h,
                                  'lags'=lags,
                                  'lagy'=lagy,
                                  'rolling'=rollw,
                                  'time'=t,
                                  'model'=model_name,
                                  'nvar'=nvar,
                                  'RMSFE'=RMSFE,
                                  'MAFE'=MAFE,
                                  'err'=err,
                                  'fcast'=fcast
                       )
      )
    }
    
    # RW prediction error
    fcast <- tmp.data[variable] %>% head(-h) %>% tail(1) %>% as.numeric()
    err <- fcast - y_pred
    RMSFE <- (err)^2
    MAFE <- abs(err)
    
    results <- rbind(results,
                     data.frame('variable'=variable,
                                'h'=h,
                                'lags'=lags,
                                'lagy'=lagy,
                                'rolling'=rollw,
                                'time'=t,
                                'model'='RW',
                                'nvar'=NA,
                                'RMSFE'=RMSFE,
                                'MAFE'=MAFE,
                                'err'=err,
                                'fcast'=fcast
                     )
    )
    row.names(results) <- NULL
    return(results)
  }
  stopCluster(cl)
  # Store intermediate results
  final_results <- rbind(final_results, fcast_result)
  write.csv(final_results, file=paste0('/home/prifti/My_scratch/Code/FRED/Results/tmp_results/US_BAR3_MD_',variable,'_forecast_h',h,'_lags',lags,'_lagy',lagy,'_rolling',rollw,'_',year(start_oos),'_',year(end_oos),'_',format(Sys.Date(),'%Y_%m_%d'),'(no_T)lambda_automax_cmin0001(step_back_rollw).csv'), row.names=F)
}

# Store complete results
write.csv(final_results, file=paste0('/home/prifti/My_scratch/Code/FRED/Results/US_BAR3_MD_forecast_',year(start_oos),'_',year(end_oos),'_',format(Sys.Date(),'%Y_%m_%d'),'(no_T)lambda_automax_cmin0001(step_back_rollw).csv'), row.names=F)