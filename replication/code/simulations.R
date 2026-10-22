rm(list=ls())

#-------------------------------------------------------------------------------
# 1. FUNCTIONS
#-------------------------------------------------------------------------------
source('R/functions.R')

#-------------------------------------------------------------------------------
# 2. SETUP
#-------------------------------------------------------------------------------
nmc           <- 10                                                             # number of Monte Carlo simulations
criteria      <- 'ebic'
p             <- 1                                                              # VAR(p)
save_heatmaps <- FALSE                                                          # TRUE if you want to save coefficient for heatmaps, FALSE if you don't

all_mc_results <- data.frame()                                                  # initialize table for the results of each MC design

mc_designs <- expand.grid(
  T_  = c(50, 100, 500),
  r   = c('1', '2', '5', '10'),
  rho = c(0, 0.5),
  K   = c(20, 40, 80, 160, 320),
  stringsAsFactors = FALSE
  ) %>%
  rowwise() %>%
  mutate(r = eval(parse(text = gsub('K', as.character(K), r))))

#-------------------------------------------------------------------------------
# 3. MONTE CARLO
#-------------------------------------------------------------------------------
for(i in c(1:nrow(mc_designs))) {
  current <- mc_designs[i,]
  rho <- current$rho
  K <- current$K
  T_ <- current$T_                                                              # sample size: 50, 100, 500
  r <- current$r                                                                # sparsity level: 1, 2, 5, 10

  # Monte Carlo
  info_parallel <- c('K','T_','rho','r','criteria','p')
  cl            <- makeCluster(detectCores() - 1)                               # use all the cores leaving one free
  registerDoParallel(cl)
  clusterExport(cl, varlist = info_parallel)
  clusterEvalQ(cl, {source('R/functions.R')})                                   # upload in each core the functions
  mc_result <- foreach(imc = 1:nmc, .combine=rbind, .export = info_parallel) %dopar% {
    set.seed(imc)
    DGP <- CMP_DGP(K = K, T_ = T_+1, r = r, rho = rho)
    A <- DGP[['A']] %>% unlist() %>% as.matrix()
    DGP_data <- DGP[['dataset']] %>% unlist() %>% as.matrix()

    supp_A <- A != 0                                                            # true support recovery
    sum_A <- sum(supp_A)                                                        # number of parameters different from zero in matrix A

    dataset_lags <- create_lags(x = DGP_data, p = p)                            # create lag = 1

    # Y
    y      <- dataset_lags[,(1:K)] %>% head(-1)                                 # response variables: Y in 1:T
    ym     <- colMeans(y)                                                       # mean of each response variable in y
    y_pred <- dataset_lags[,(1:K)] %>% tail(1)                                  # Y at T+1 for evaluating forecast accuracy

    # X
    x      <- dataset_lags[,-(1:K)] %>% head(-1)                                # explanatory variables: X in 1:T
    xm     <- colMeans(x)                                                       # mean of each explanatory variable in x
    x_pred <- dataset_lags[,-(1:K)] %>% tail(1)                                 # X at T+1 for evaluating forecast accuracy

    # Broken Adaptive Ridge (BAR)
    # BAR_est <- bar_ccd(y = y, x = x, crit = criteria,
    #                    intercept = TRUE, standardize = TRUE)
    # BAR <- BAR_est$beta                                                         # coefficients without intercept

    # Adaptive LASSO with Ridge initialization
    aLASSO_ridge_est <- var.lasso(y=y, x=x, K=K, p = p, adaptive = 'ridge',
                                  crit = criteria, standardize = TRUE,
                                  intercept = TRUE)
    aLASSO_ridge <- aLASSO_ridge_est$coefficients[-1,]                          # take the coefficients discarding the intercept

    # Adaptive LASSO with LASSO initialization
    aLASSO_lasso_est <- var.lasso(y=y, x=x, K=K, p = p, adaptive = 'lasso',
                                  crit = criteria, standardize = TRUE,
                                  intercept = TRUE)
    aLASSO_lasso <- aLASSO_lasso_est$coefficients[-1,]                          # take the coefficients discarding the intercept

    # LASSO
    LASSO_est <- var.lasso(y = y, x = x, K = K, p = p, alpha = 1, adaptive = 'none',
                           crit = criteria, standardize = TRUE, intercept = TRUE)
    LASSO <- LASSO_est$coefficients[-1,]                                        # take the coefficients discarding the intercept

    # SCAD
    SCAD_est <- var.scad(y = y, x = x, K = K, crit = criteria,
                         standardize = TRUE, intercept = TRUE)
    SCAD <- SCAD_est$coefficients[-1,]                                          # take the coefficients discarding the intercept


    models <- list(#'BAR'          = BAR,
                   'aLASSO_ridge' = aLASSO_ridge,
                   'aLASSO_lasso' = aLASSO_lasso,
                   'LASSO'        = LASSO,
                   'SCAD'         = SCAD
                   )

    imc_result <- data.frame()
    for (model_name in names(models)) {
      mod   <- models[[model_name]]                                             # model coefficients
      a0    <- ym - xm %*% mod                                                  # intercept
      nvar  <- sum(mod != 0)                                                    # number of non-zero coefficients

      # ----------------------- Performance measures ---------------------------
      supp_Ahat <- matrix(mod != 0, nrow = K)                                                     # estimated support

      # 1. Correct sparsity
      correct_sparsity <- all(supp_A == supp_Ahat)

      # 2. True model included
      true_model_included <- all(supp_Ahat[supp_A])

      # 3. Fraction of relevant variables included
      frac_relevant_included <- sum(supp_Ahat[supp_A]) / sum_A

      # 4. Number of variables included
      num_variables_included <- sum(supp_Ahat) / K

      # 5. RMSE of A (estimation accuracy)
      RMSE <- sum((mod - A)^2)                                                  # root mean squared error

      # 6. RMSFE one-step ahead (forecast accuracy)
      fcast <- a0 + x_pred %*% mod                                              # point forecast
      err   <- (a0 + x_pred %*% mod) - y_pred                                   # forecast error
      RMSFE <- mean(err^2)                                                      # root mean squared forecast error

      # 7. False Positive Rate (fpr) and False Negative Rate (fnr)
      TP  <- sum( supp_A  &  supp_Ahat)                                         # true positive
      FP  <- sum(!supp_A  &  supp_Ahat)                                         # false positive
      TN  <- sum(!supp_A  & !supp_Ahat)                                         # true negative
      FN  <- sum( supp_A  & !supp_Ahat)                                         # false negative
      fpr <- FP / (FP + TN)                                                     # false positive rate
      fnr <- FN / (FN + TP)                                                     # false negative rate

      # Collect results
      imc_result <- rbind(imc_result,
                          data.frame('model'                  = model_name,
                                     'correct_sparsity'       = correct_sparsity,
                                     'true_model_included'    = true_model_included,
                                     'frac_relevant_included' = frac_relevant_included,
                                     'num_variables_included' = num_variables_included,
                                     'RMSE'                   = RMSE,
                                     'RMSFE'                  = RMSFE,
                                     'fpr'                    = fpr,
                                     'fnr'                    = fnr
                                     )
                          )

    }

    return(imc_result)
  }
  stopCluster(cl)

  mc_result <- mc_result %>%
    group_by(model) %>%
    summarise(
      across(where(~ is.numeric(.x) || is.logical(.x)), ~ mean(.x, na.rm = TRUE))
    ) %>%
    mutate(
      RMSE = sqrt(RMSE),
      RMSFE = sqrt(RMSFE)
    ) %>%
    data.frame()

  all_mc_results <- rbind(all_mc_results,
                          cbind(
                            data.frame('K' = K,'T' = T_, 'r' = r, 'rho' = rho),
                            mc_result
                            )
                          )
  # Save results for each MC design
  write.csv(mc_result,
            file = paste0('output/mc/MC_',nmc,'_',criteria,'_K',K,'_T',T_,'_r',r,'_rho',rho,'.csv'),
            row.names = FALSE)
}

# Save results of all MC designs
write.csv(all_mc_results,
          file = paste0('output/mc/MC_',nmc,'_',criteria,'_mc',mc,'.csv'),
          row.names = FALSE)
