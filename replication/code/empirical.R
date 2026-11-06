rm(list=ls())

#-------------------------------------------------------------------------------
# 1. FUNCTIONS
#-------------------------------------------------------------------------------
source(file.path('code','functions.R'))

#-------------------------------------------------------------------------------
# 2. FRED-MD
#-------------------------------------------------------------------------------
vintage   <- '2026-06'
fred_data <- transform_fredmd(file = file.path('data',                          # transform the FRED-MD
                                               paste0(vintage,'-MD.csv')), 
                              transform = TRUE)
colnames_to_rm <- c('ACOGNO','UMCSENTx','TWEXAFEGSMTHx','ANDENOx','VIXCLSx')    # variables to be removed   
fred_data      <- fred_data[,setdiff(names(fred_data), colnames_to_rm)]         # remove the variables
fred_data$date <- zoo::as.yearmon(fred_data$date)
dataset        <- fred_data[fred_data$date >= zoo::as.yearmon('1960-01'),]      # data starting from January 1960

#-------------------------------------------------------------------------------
# 3. FORECASTING
#-------------------------------------------------------------------------------
start_oos <- zoo::as.yearmon('1990-01')                                         # start out-of-sample
end_oos   <- zoo::as.yearmon('2025-12')                                         # end out-of-sample
freq_time <- 12                                                                 # time frequency: 12 for monthly
criteria  <- 'ebic'                                                             # information criterion to use

fcast_spec <- expand.grid(
  variable = 'CPIAUCSL',                                                        # target variable: inflation
  lags     = 1,                                                                 # no. of lags for X à la Medeiros et al. (2021)
  lagy     = 4,                                                                 # no. of autoregressive lags à la Medeiros et al. (2021)
  h        = c(1, 3, 6, 12),                                                    # forecast horizon
  rollw    = c(360, 0),                                                         # 360 = rolling window estimation à la Medeiros et al. (2021)
  stringsAsFactors = FALSE                                                      # 0 = expanding window estimation                                                   
  )

all_forecast <- data.frame()                                                    # initialize table for all forecasting results

cat('Forecasts running...\n')
pb <- txtProgressBar(min = 0, max = nrow(fcast_spec), style = 3,                # progress bar
                     width = max(50,nrow(fcast_spec)))
for (i in 1:nrow(fcast_spec)) {
  setTxtProgressBar(pb, i)
  
  current   <- fcast_spec[i,]
  variable  <- current$variable
  h         <- current$h
  rollw     <- current$rollw
  lags      <- current$lags
  lagy      <- current$lagy
  row_trim  <- max(lagy, lags)+h-1
  step_back <- rollw / freq_time
  info_parallel <- c('variable','h','freq_time','lags','lagy','start_oos',
                     'end_oos','dataset','rollw','step_back','row_trim')
  cl <- makeCluster(detectCores() - 1)                                          # use all the cores leaving one free
  registerDoParallel(cl)
  clusterExport(cl, varlist = info_parallel)
  clusterEvalQ(cl, {source(file.path('code','functions.R'))})                   # upload the functions in each core
  spec_forecast <- foreach(t = seq(start_oos, end_oos, by = 1 / freq_time),
                           .combine = rbind) %dopar% {
    if(rollw == 0) {                                                   
      tmp.data <- dataset[dataset$date <= t, !(names(dataset) %in% c('date'))]  # expanding window estimation
    } else {
      if (t >= zoo::as.yearmon('2001-01') & rollw == 360) {                     # à la Medeiros et al. (2021)
        step_back <- 492 / freq_time   
      }
      tmp.data <- dataset[dataset$date <= t & dataset$date >= t - step_back,    # rolling window estimation
                          !(names(dataset) %in% c('date'))]
    }
                            
    # create autoregressive lags          
    y_lags <- create_lags(tmp.data[variable],
                          p = lagy + h - 1,
                          include.original = FALSE,
                          trim = FALSE) %>% 
      tail(-row_trim)
    
    # create lags for the explanatory variables
    x_lags <- create_lags(tmp.data[,!(names(tmp.data) %in% c(variable))],
                          p = lags + h - 1,
                          include.original = FALSE,
                          trim = FALSE) %>% 
      tail(-row_trim)
    
    # remove the 
    if (h > 1) {
      y_lags <- y_lags[,-(1:(h-1))]
      x_lags <- x_lags[,-(1:((h-1) * (ncol(tmp.data)-1)))]
    }
    
    # target variable y
    y_lag  <- y_lags %>% head(-h)
    y      <- tmp.data[variable] %>% tail(-row_trim) %>% head(-h) %>% as.matrix()
    T_     <- length(y)
    ym     <- mean(y)
    ys     <- sqrt(sum((y - ym)^2) / (T_))
    y_pred <- tmp.data[variable] %>% tail(1) %>% as.numeric()
    y_c    <- y - ym
    y_cs   <- y_c / sqrt(sum(y_c^2) / T_)
    
    # explanatory variables X
    x_lag  <- x_lags %>% head(-h)
    x      <- cbind(y_lag, x_lag)
    xm     <- apply(x, 2, function(x) mean(x))
    xs     <- apply(x, 2, function(x) sqrt(sum((x-mean(x))^2)/T_))
    x_pred <- cbind(y_lags, x_lags) %>% tail(1) %>% as.matrix()
    x_c    <- apply(x, 2, function(x) x - mean(x))
    x_cs   <- apply(x_c, 2, function(x) x / sqrt(sum(x^2)/T_))
    
    # BAR
    BAR_est          <- bar_ccd_reg(y = y_cs, x = x_cs, crit = criteria)
    BAR              <- BAR_est$beta * ys / xs                                  # go back to the original scale and remove intercept 
    
    # Adaptive LASSO with Ridge initialization
    aLASSO_ridge_est <- univar.lasso(y = y_cs, x = x_cs, alpha = 1, adaptive = 'ridge',
                                     crit = criteria, standardize = FALSE, intercept = FALSE)
    aLASSO_ridge     <- aLASSO_ridge_est$coefficients[-1] * ys / xs             # go back to the original scale and remove intercept
    
    # Adaptive LASSO with LASSO initialization
    aLASSO_lasso_est <- univar.lasso(y = y_cs, x = x_cs, alpha=1, adaptive = 'lasso',
                                     crit = criteria, standardize = FALSE, intercept = FALSE)
    aLASSO_lasso     <- aLASSO_lasso_est$coefficients[-1] * ys / xs             # go back to the original scale and remove intercept
    
    # LASSO
    LASSO_est        <- univar.lasso(y = y_cs, x = x_cs, alpha = 1, adaptive = 'none',
                              crit = criteria, standardize = FALSE, intercept = FALSE)
    LASSO            <- LASSO_est$coefficients[-1] * ys / xs                    # go back to the original scale and remove intercept
  
    # SCAD
    SCAD_est         <- ic.scad(y = y_cs, x = x_cs, crit = criteria)
    SCAD             <- SCAD_est$coefficients[-1] * ys / xs                     # go back to the original scale and remove intercept

    # AR(1)
    AR1_est          <- arima(y_cs, order = c(1,0,0), include.mean = FALSE)
    AR1              <- AR1_est$coef['ar1'] %>% as.numeric()
    
    t_forecast <- data.frame()
    models <- list('BAR'          = BAR,
                   'aLASSO_ridge' = aLASSO_ridge, 
                   'aLASSO_lasso' = aLASSO_lasso, 
                   'LASSO'        = LASSO, 
                   'SCAD'         = SCAD,
                   'AR1'          = AR1
                   )
    
    for (model_name in names(models)) {
      mod   <- models[[model_name]]
      if (model_name == 'AR1') {
        nvar  <- 1
        a0    <- ym - ym * AR1
        fcast <- a0 + AR1 * x_pred[,paste(variable,'l1')] %>% as.numeric()
      } else {
        nvar  <- sum(mod != 0)
        a0    <- ym - xm %*% mod
        fcast <- a0 + x_pred %*% mod
      }
      err   <- fcast - y_pred
      RMSFE <- err^2
      MAFE  <- abs(err)
      
      t_forecast <- rbind(t_forecast,
                        data.frame('variable' = variable,
                                   'criteria' = criteria,
                                   'h'        = h,
                                   'lags'     = lags,
                                   'lagy'     = lagy,
                                   'rolling'  = rollw,
                                   'time'     = t,
                                   'model'    = model_name,
                                   'nvar'     = nvar,
                                   'RMSFE'    = RMSFE,
                                   'MAFE'     = MAFE,
                                   'err'      = err,
                                   'fcast'    = fcast
                                   )
                        )
    }
    
    row.names(t_forecast) <- NULL
    
    return(t_forecast)
  }
  stopCluster(cl)

  all_forecast <- rbind(all_forecast, spec_forecast)
  
  # save results for each forecast specification
  write.csv(spec_forecast,
            file = file.path('output','forecasting','csv',paste0('var',variable,'_h',h,'_lags',lags,'_lagy',lagy,'_roll',rollw,'_',criteria,'.csv')),
            row.names = FALSE)
}
close(pb)

# save results of all MC designs
write.csv(all_forecast,
          file = file.path('output','forecasting','csv',paste0('FINAL_FORECASTING_',criteria,'.csv')),
          row.names = FALSE)

#-------------------------------------------------------------------------------
# 3. PREDICTIVE STATISTICS & MODEL CONFIDENCE SET
#-------------------------------------------------------------------------------
tab_fcast_results <- read.csv(file.path('output','forecasting','csv',
                                        paste0('FINAL_FORECASTING_',criteria,'.csv')),
                              sep = ',',
                              dec = '.') %>%
  data.frame()

# out-of-sample periods
ls_period <- list(					   
  '1990-2019' = c(1990, 2019),
  '2020-2021' = c(2020, 2021),
  '2022-2025' = c(2022, 2025)
)

# predictive statistics: RMSFE, MAFE
pred_stats <- data.frame()
for (p in seq_along(ls_period)) {
  period_name <- names(ls_period)[p]
  
  tmp_fcast_result <- tab_fcast_results %>%
    dplyr::mutate(year = lubridate::year(zoo::as.yearmon(time))) %>%
    dplyr::select(-time) %>%
    dplyr::filter(year >= ls_period[[p]][1], year <= ls_period[[p]][2]) %>%
    dplyr::mutate(dplyr::across(-c(variable, criteria, model, h, rolling, lags, lagy), as.numeric)) %>%
    dplyr::group_by(variable, criteria, model, h, rolling, lags, lagy) %>%
    dplyr::summarise(
      RMSFE   = mean(RMSFE, na.rm = TRUE),
      MAFE    = mean(MAFE, na.rm = TRUE),
      nvar    = mean(nvar, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    dplyr::mutate(
      RMSFE  = sqrt(RMSFE),
      period = period_name
    ) %>%
    dplyr::relocate(period, .before = 1) %>%
    dplyr::arrange(variable, criteria, period, h, rolling, lags, lagy, model)
  
  pred_stats <- dplyr::bind_rows(pred_stats, tmp_fcast_result)
}

# RMSFE, MAFE ratios w.r.t. AR(1)
pred_ratio <- pred_stats %>%
  dplyr::left_join(
    pred_stats %>%
      dplyr::filter(model == 'AR1') %>%
      dplyr::select(variable, criteria, period, h, rolling, lags, lagy,
                    RMSFE_AR1 = RMSFE,
                    MAFE_AR1  = MAFE),
    by = c('variable','criteria','period','h','rolling','lags','lagy')
  ) %>%
  dplyr::mutate(
    ratio_RMSFE = RMSFE / RMSFE_AR1,
    ratio_MAFE  = MAFE  / MAFE_AR1
  ) %>%
  dplyr::filter(model != 'AR1') %>%
  dplyr::arrange(variable, criteria, match(period, c('1990-2019','2020-2021','2022-2025')),
                 h, rolling, lags, lagy, model)

# MCS specifications
mcs_spec <- expand.grid(
  alpha        = c(0.50),
  T_statistics = c('Tmax'),
  lags         = unique(pred_stats$lags),
  p            = seq_along(ls_period),
  h            = unique(pred_stats$h),
  var          = unique(pred_stats$variable),
  rolling      = unique(pred_stats$rolling),
  loss_mcs     = c('absolute','square'),
  stringsAsFactors = FALSE)

# produce SSM of the MCS
info_parallel <- c('tab_fcast_results','ls_period','mcs_spec')
cl            <- makeCluster(25-1)                                              # use all the cores leaving one free
registerDoParallel(cl)
clusterExport(cl, varlist = info_parallel)
clusterEvalQ(cl, {source(file.path('code','functions.R'))})
mcs_results <- foreach(i = 1:nrow(mcs_spec), .combine = rbind) %dopar% {
  current      <- mcs_spec[i,]
  alpha        <- current$alpha
  T_statistics <- current$T_statistics
  loss_mcs     <- current$loss_mcs
  lags         <- current$lags
  h            <- current$h
  var          <- current$var
  p            <- current$p
  rolling      <- current$rolling
  period_name  <- names(ls_period)[p]
  
  if (T_statistics == 'Tmax') {
    rank_type <- 'Rank_M'
  } else if (T_statistics == 'TR') {
    rank_type <- 'Rank_R'
  }
  
  tmp_mcs <- tab_fcast_results %>% 
    dplyr::mutate(time = zoo::as.yearmon(time)) %>%
    dplyr::mutate(year = lubridate::year(time)) %>% 
    dplyr::filter(year >= ls_period[[p]][1], year <= ls_period[[p]][2]) %>%
    dplyr::filter(rolling ==!! rolling, h ==!! h, lags ==!! lags)

  tmp_mcs <- tmp_mcs[tmp_mcs$variable == var,]
  
  if (loss_mcs == 'absolute') {
    tmp_mcs_result <- mcs.test(models = 
                                 list('BAR'         = abs(tmp_mcs[tmp_mcs$model == 'BAR', 'err']),
                                      'aLASSO_ridge'= abs(tmp_mcs[tmp_mcs$model == 'aLASSO_ridge', 'err']),
                                      'aLASSO_lasso'= abs(tmp_mcs[tmp_mcs$model == 'aLASSO_lasso', 'err']),
                                      'LASSO'       = abs(tmp_mcs[tmp_mcs$model == 'LASSO', 'err']),
                                      'SCAD'        = abs(tmp_mcs[tmp_mcs$model == 'SCAD', 'err'])
                                     ),
                               alpha     = alpha, 
                               statistic = T_statistics)
  } else if (loss_mcs == 'square') {
    tmp_mcs_result <- mcs.test(models =
                                 list('BAR'         = (tmp_mcs[tmp_mcs$model == 'BAR', 'err'])^2,
                                      'aLASSO_ridge'= (tmp_mcs[tmp_mcs$model == 'aLASSO_ridge', 'err'])^2,
                                      'aLASSO_lasso'= (tmp_mcs[tmp_mcs$model == 'aLASSO_lasso', 'err'])^2,
                                      'LASSO'       = (tmp_mcs[tmp_mcs$model == 'LASSO', 'err'])^2,
                                      'SCAD'        = (tmp_mcs[tmp_mcs$model == 'SCAD', 'err'])^2
                                      ),
                               alpha     = alpha,
                               statistic = T_statistics)
  }

  result <- data.frame('loss'         = loss_mcs,
                       'period'       = period_name,
                       'h'            = h,
                       'variable'     = var,
                       'lags'         = lags,
                       'rolling'      = rolling,
                       'T_statistics' = T_statistics,
                       'alpha'        = alpha,
                       'model'        = rownames(tmp_mcs_result@show),
                       'rank_model'   = as.vector(tmp_mcs_result@show[,rank_type]),
                       'pvalue'       = round(tmp_mcs_result@Info[['mcs_pvalue']]+1e-12,2)
                       )
  
  return(result)
}
stopCluster(cl)

# save MCS results
write.csv(mcs_results,
          file = file.path('output','forecasting','csv',paste0('MCS_',criteria,'.csv')),
          row.names = FALSE)

# add SSM of MCS to the predictive statistics ratios
pred_ratio <- pred_ratio %>%
  mutate(
    absMCS = as.numeric(
      paste('absolute', model, period, variable, criteria, rolling, lags, h) %in%
        paste(mcs_results$loss, mcs_results$model, mcs_results$period,
              mcs_results$variable, mcs_results$rolling,
              mcs_results$lags, mcs_results$h)
    ),
    sqMCS = as.numeric(
      paste('square', model, period, variable, criteria, rolling, lags, h) %in%
        paste(mcs_results$loss, mcs_results$model, mcs_results$period,
              mcs_results$variable, mcs_results$rolling,
              mcs_results$lags, mcs_results$h)
    )
  )

# LaTeX
pred_ratio_long <- pred_ratio %>%
  dplyr::select(variable, criteria, period, h, rolling, lagy, lags, model,
                starts_with('ratio_RMSFE'),
                starts_with('ratio_MAFE'),
                absMCS,
                sqMCS
  ) %>%
  tidyr::pivot_longer(
    cols         = c(starts_with('ratio_'), absMCS, sqMCS),
    names_to     = 'measure',
    names_prefix = 'ratio_',
    values_to    = 'value'
  )

pred_ratio_wide <- pred_ratio_long %>%
  tidyr::pivot_wider(
    id_cols     = c(variable, criteria, period, rolling, lagy, lags, model),
    names_from  = c(h, measure), 
    values_from = value) %>% 
  arrange(variable, criteria, match(period, c('1990-2019','2020-2021','2022-2025')),
          rolling, lags, lagy, model)

#-------------------------------------------------------------------------------
# 4. LATEX TABLES
#-------------------------------------------------------------------------------
# LaTeX table specifications
latex_spec <- expand.grid(
  variable = unique(pred_ratio_long$variable),
  rolling  = unique(pred_ratio_long$rolling),
  lags     = unique(pred_ratio_long$lags),
  stringsAsFactors = FALSE)

models_order <- c('BAR','aLASSO_ridge','aLASSO_lasso','LASSO','SCAD')
measures     <- c('RMSFE','MAFE')
horizons     <- unique(pred_ratio_long$h)
periods      <- pred_ratio_long %>% 
  distinct(period) %>%
  mutate(end = as.numeric(stringr::str_extract(period, '[0-9]{4}$'))) %>%
  arrange(end) %>%
  pull(period)

n_horizons   <- length(horizons)
n_periods    <- length(periods)
n_models     <- length(unique(pred_ratio_long$model))

# produce LaTeX tables
for (i in 1:nrow(latex_spec)) {
  current      <- latex_spec[i,]
  cur_variable <- current$variable
  cur_rolling  <- current$rolling
  cur_lags     <- current$lags
  
  tmp_tab <- pred_ratio_wide %>% 
    dplyr::filter(variable == cur_variable, rolling == cur_rolling, lags == cur_lags) %>% 
    dplyr::arrange(match(model, models_order)) %>%
    dplyr::mutate(period = factor(period, levels = unique(period))) %>%
    dplyr::select(period, everything())
  
  tmp_tab2 <- tmp_tab %>% mutate(dplyr::across(where(is.numeric), ~ sprintf('%.2f', .))) %>% 
    dplyr::mutate(period = factor(period, levels = periods)) %>%
    arrange(
      match(period, c('1990-2019','2020-2021','2022-2025')),
      match(model, models_order)
    ) %>% 
    as.data.frame()
  
  rmsfe_cols <- grep('_RMSFE$', colnames(tmp_tab2), value = TRUE)
  mafe_cols  <- grep('_MAFE$', colnames(tmp_tab2), value = TRUE)
  absmcs_cols  <- grep('_absMCS$', colnames(tmp_tab2), value = TRUE)
  absmcs_vals <- tmp_tab2[, absmcs_cols]
  sqmcs_cols  <- grep('_sqMCS$', colnames(tmp_tab2), value = TRUE)
  sqmcs_vals <- tmp_tab2[, sqmcs_cols]
  series <- sub('_RMSFE$', '', rmsfe_cols)
  extra_cols <- setdiff(colnames(tmp_tab2), c(rmsfe_cols, mafe_cols, absmcs_cols, sqmcs_cols))  
  
  # take the best RMSFE and MAFE ratios for each out-of-sample period
  for (p in periods) {
    idx <- tmp_tab2$period == p
    for (col in setdiff(names(tmp_tab2), c(extra_cols, absmcs_cols, sqmcs_cols))) {
      min_val <- min(tmp_tab2[idx, col], na.rm = TRUE)
      n_min   <- sum(tmp_tab2[idx, col] == min_val, na.rm = TRUE)
      tmp_tab2[idx, col] <- ifelse(tmp_tab2[idx, col] == min_val,
                                   if (n_min == 1) {
                                     paste0('\\underline{\\textbf{', tmp_tab2[idx, col], '}}')
                                   } else {
                                     paste0('\\textbf{', tmp_tab2[idx, col], '}')
                                   },
                                   tmp_tab2[idx, col])
    }
  }
  
  tmp_tab2$period_tex <- rep(NA_character_, nrow(tmp_tab2))
  for(p in periods){
    idx                          <- which(tmp_tab2$period == p)
    tmp_tab2$period_tex[idx[1]]  <- paste0('\\multirow{', n_models, '}{*}{', p, '}')
    tmp_tab2$period_tex[idx[-1]] <- '{}'
  }
  
  # to highlight the SSM of the MCS
  for (hh in unique(sub('_RMSFE$', '', rmsfe_cols))) {
    rmsfe_col <- paste0(hh, '_RMSFE')
    sqmcs_col <- paste0(hh, '_sqMCS')
    if (rmsfe_col %in% names(tmp_tab2) && sqmcs_col %in% names(tmp_tab2)) {
      tmp_tab2[[rmsfe_col]] <- ifelse(
        as.numeric(tmp_tab2[[sqmcs_col]]) == 1,
        paste0('\\cellcolor{gray!30}{', tmp_tab2[[rmsfe_col]], '}'),
        tmp_tab2[[rmsfe_col]]
      )
    }
    mafe_col  <- paste0(hh, '_MAFE')
    absmcs_col <- paste0(hh, '_absMCS')
    if (mafe_col %in% names(tmp_tab2) && absmcs_col %in% names(tmp_tab2)) {
      tmp_tab2[[mafe_col]] <- ifelse(
        as.numeric(tmp_tab2[[absmcs_col]]) == 1,
        paste0('\\cellcolor{cyan!10}{', tmp_tab2[[mafe_col]], '}'),
        tmp_tab2[[mafe_col]]
      )
    }
  }
  
  tmp_tab3 <- tmp_tab2 %>% dplyr::select(-all_of(sqmcs_cols), -all_of(absmcs_cols))

  header1 <- paste(
    '\\toprule\n\\multirow{2}{*}{period} & \\multirow{2}{*}{model} &',
    paste(sapply(horizons, function(h) {
      paste0('\\multicolumn{', length(measures), '}{c}{h=', h, '}')
    }), collapse = ' & '),
    '\\\\'
  )
  
  start <- 3
  header2 <- paste(
    sapply(seq_along(horizons), function(i) {
      end   <- start + length(measures) - 1
      rule  <- paste0('\\cmidrule(lr){', start, '-', end, '}')
      start <<- end + 1
      rule
    }),
    collapse = ' '
  )
  
  header3 <- paste('{} & {} &',
                   paste(rep(measures, length(horizons)), collapse = ' & '),
                   '\\\\ \n \\midrule\n'
                   )
  
  bottomrule <- paste0('\\bottomrule\n',
                       '\\multicolumn{',
                       2 + length(rep(measures, length(horizons))),
                       '}{l}{\\begin{minipage}{\\dimexpr\\columnwidth-',
                       if (ncol(tmp_tab3) > 13) {0} else {10},
                       '\\tabcolsep}\n',
                       '\\footnotesize\\textit{Notes}: the RMSFE and the MAFE are relative to the AR(1) benchmark. ',
                       'Bold and underlined entries denote the best-performing model(s) for each forecast horizon and period. ', 
                       'Grey cells indicate models belonging to the Set of Superior Models (SSM) under the 50\\% MCS procedure using the $T_{\\text{max}}$ statistic with squared errors, whereas blue cells indicate models included in the SSM when absolute errors are used.\n',
                       '\\end{minipage}\n}')
  
  addtorow         <- list()
  addtorow$command <- c(paste(c(header1, header2, header3), collapse = '\n'))
  addtorow$command <- c(addtorow$command, bottomrule) 
  addtorow$pos     <- list(0)
  addtorow$pos     <- c(addtorow$pos, nrow(tmp_tab3))
  
  file_out         <- file.path('output','forecasting','tables_tex',
                                paste0('tab_',cur_variable,'_lags',cur_lags,'_roll',cur_rolling,'.tex')
                                )
  con              <- file(file_out, open = 'wt')
  
  cat('\\begin{table}[H]\n', file = con)
  cat('\\centering\n', file = con)
  
  caption <- if (cur_rolling == 0) {
    paste0('RMSFE and MAFE ratios relative to the AR(1) using an expanding window estimation.')
  } else {
    paste0('RMSFE and MAFE ratios relative to the AR(1).')
  }
  
  cat(paste0('\\caption{',caption,'}\n'), file = con)
  cat(paste0('\\label{tab:',cur_variable,'_lags',cur_lags,'_roll',cur_rolling,'}\n'), file = con)
  cat('\\resizebox{\\textwidth}{!}{%\n', file = con)
  
  rows_with_midrule <- c()
  tab               <- tmp_tab3
  tab$period <- vapply(tab$period, function(x) {
    if (length(x) == 0 || is.null(x)) return(NA_real_)
    as.character(x)
  }, FUN.VALUE = character(1))
  for (j in seq_along(periods)) {
    idx <- which(tab$period == periods[j])
    if (length(idx) == 0) next
    if (j < length(periods)) {
      rows_with_midrule <- c(rows_with_midrule, max(idx))
    }
  }
  
  tab <- tab %>% dplyr::select(period_tex, model, everything(),
                               -c(period, variable, criteria, rolling, lags, lagy))

  tab_print        <- xtable(tab, include.rownames = FALSE)
  align(tab_print) <- c('c','c','c', rep('c',ncol(tab)-2))
  
  # modify the models' name using parenthesis
  tab_print$model <- sub('_(.*)$', '\\(\\1\\)', tab_print$model)
  
  print(tab_print,
        type                   = 'latex',
        file                   = con,
        include.rownames       = FALSE,
        include.colnames       = FALSE,
        booktabs               = TRUE,
        add.to.row             = addtorow,
        hline.after            = rows_with_midrule,
        sanitize.text.function = identity,
        comment                = FALSE,
        floating               = FALSE
        )
  
  cat('} % end resizebox\n', file = con)
  cat('\\end{table}\n', file = con)
  close(con)
}

#-------------------------------------------------------------------------------
# 5. CSSFED PLOT
#-------------------------------------------------------------------------------
# import data
input_file  <- paste0('FINAL_FORECASTING_',criteria,'.csv')
data_fcast <- read.csv(file.path('output','forecasting','csv', input_file),
                       sep = ',',
                       dec = '.', 
                       stringsAsFactors = FALSE, 
                       check.names = FALSE)

# convert time to yearmon
orig_loc <- Sys.getlocale('LC_TIME')

# if months in csv file are in English: Jan, Feb, Mar, Apr, etc.
Sys.setlocale('LC_TIME', 'C')

time_original <- trimws(as.character(data_fcast$time))

time_converted <- zoo::as.yearmon(time_original, format = '%b %Y')

# conversion control
if (anyNA(time_converted)) {
  invalid_time <- unique(time_original[is.na(time_converted)])
  stop(paste0("Impossible to convert some values of variable 'time'. Examples: ",
              paste(head(invalid_time, 10), collapse = ', ')))
}

data_fcast$time <- time_converted

# back to original local
Sys.setlocale('LC_TIME', orig_loc)

# prepare prediction statistics and rename 'aLASSO' models
pred_stats <- data_fcast %>%
  dplyr::mutate(
    dplyr::across(
      -c(variable, criteria, model, h, rolling, lags, lagy, time),
      ~ suppressWarnings(as.numeric(.x))
    )
  ) %>% 
  dplyr::mutate(
    model = dplyr::recode(
      model,
      'aLASSO_ridge' = 'aLASSO(ridge)',
      'aLASSO_lasso' = 'aLASSO(lasso)'
      )
  )

# check benchmark model
benchmark_model <- 'AR1'
if (!benchmark_model %in% pred_stats$model) {
  stop(paste0("The benchmark '",benchmark_model,"' is not present in variabile 'model'. Available models: ",
              paste(sort(unique(pred_stats$model)), collapse = ', ')))
}

# compute cumulative sum squared forecast error differential (CSSFED)
pred_cssfed <- pred_stats %>%
  inner_join(
    pred_stats %>%
      dplyr::filter(model == benchmark_model) %>%
      dplyr::select(variable, criteria, h, rolling, lags, lagy, time, err_benchmark = err),
    by = c('variable','criteria','h','rolling','lags','lagy','time')) %>%
  dplyr::filter(model != benchmark_model) %>%
  dplyr::mutate(diff_sq = err^2 - err_benchmark^2) %>%
  dplyr::arrange(variable, criteria, h, rolling, lags, lagy, model,time) %>%
  dplyr::group_by(variable, criteria, h, rolling, lags, lagy, model) %>%
  dplyr::mutate(CCSFED = cumsum(diff_sq)) %>%
  ungroup()

# graphical settings
model_order <- c('BAR','aLASSO(ridge)','aLASSO(lasso)','LASSO','SCAD')

cols <- c(
  'BAR'           = '#009E73',                                                  # green
  'aLASSO(ridge)' = '#E69F00',                                                  # orange
  'aLASSO(lasso)' = '#1F77B4',                                                  # blue
  'LASSO'         = '#D62728',                                                  # red
  'SCAD'          = '#4D4D4D'                                                   # grey
)

line_types <- c(
  'BAR'            = 'solid',
  'aLASSO(ridge)'  = 'twodash',
  'aLASSO(lasso)'  = 'dashed',
  'LASSO'          = 'dotdash',
  'SCAD'           = 'dotted'
)

# check models not included in graphical settings
models_to_plot <- setdiff(unique(pred_cssfed$model),benchmark_model)
missing_models <- setdiff(models_to_plot,names(cols))

if (length(missing_models) > 0) {
  warning(paste0("The following models don't have an assigned color and ",
                 'will not be correctly represented: ',
                 paste(missing_models, collapse = ', ')))
}

# values used in loops
rollings <- sort(unique(pred_cssfed$rolling[!is.na(pred_cssfed$rolling)]))
variables <- sort(unique(pred_cssfed$variable[!is.na(pred_cssfed$variable)]))
criteria_values <- sort(unique(pred_cssfed$criteria[!is.na(pred_cssfed$criteria)]))

# produce and save graphs
for (r in rollings) {
  for (v in variables) {
    for (crit in criteria_values) {
      tmp <- pred_cssfed %>%
        filter(variable == v, rolling == r, criteria == crit, model %in% model_order) %>%
        dplyr::mutate(
          model = factor(model, levels = model_order),
          h     = factor(h, levels = sort(unique(h)))
        ) %>%
        dplyr::arrange(h, model, time)
      
      # don't save empty figures
      if (nrow(tmp) == 0) next
      
      p <- ggplot(tmp,aes(x = time, y = CCSFED, color = model, linetype = model,
                          group = model)) +
        
        # horizontal line at zero
        geom_hline(yintercept = 0,linetype = 'dashed',linewidth = 0.45,color = 'grey45') +
        
        # model series
        geom_line(linewidth = 0.9,alpha = 0.95,na.rm = TRUE) +
        
        # horizon panels
        facet_wrap(~ h,ncol = 2,scales = 'free_y',
                   labeller = labeller(h = function(x) paste0('h = ', x))) +
        
        # colors
        scale_color_manual(values = cols,breaks = model_order,drop = FALSE) +
        
        # line-types
        scale_linetype_manual(values = line_types,breaks = model_order,drop = FALSE) +
        
        # time axis
        scale_x_yearmon(format = '%Y',n = 6,expand = expansion(mult = c(0.01, 0.02))) +
        
        scale_y_continuous(labels = function(x) {
          format(x,scientific = FALSE,trim = TRUE,digits = 8)
        }) +
        
        # labels
        labs(title = NULL,subtitle = NULL,x = NULL,y = 'Cumulative squared forecast error differential',
             color = NULL,linetype = NULL) +
        
        # legend
        guides(
          color = guide_legend(nrow = 2,byrow = TRUE,
                               override.aes = list(
                                 linewidth = 1.2,alpha = 1)),
          linetype = guide_legend(nrow = 2,byrow = TRUE)) +
        
        theme_bw(
          base_size = 13
        ) +
        
        theme(
          # background
          plot.background = element_rect(fill = 'white',color = NA),
          panel.background = element_rect(fill = 'white',color = NA),
          
          # grid
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(linewidth = 0.25,color = 'grey88'),
          panel.grid.minor = element_blank(),
          
          # panel border
          panel.border = element_rect(fill = NA,color = 'grey35',linewidth = 0.5),
          panel.spacing = unit(1,'lines'),
          
          # panel labels
          strip.background = element_rect(fill = 'grey94',color = 'grey35',linewidth = 0.5),
          
          strip.text = element_text(size = 13,face = 'bold',color = 'black',
                                    margin = margin(t = 5,r = 5,b = 5,l = 5)),
          
          # axis x and y
          axis.text.x = element_text(size = 11,color = 'black'),
          axis.text.y = element_text(size = 11,color = 'black'),
          axis.title.y = element_text(size = 12,margin = margin(r = 10)),
          axis.ticks = element_line(color = 'grey35',linewidth = 0.4),
          
          # Legenda sotto il grafico
          legend.position = 'bottom',
          legend.direction = 'horizontal',
          legend.justification = 'center',
          
          # Riquadro esterno della legenda
          legend.box.background = element_rect(fill = 'white',color = 'grey30',linewidth = 0.6),
          legend.background = element_blank(),
          legend.key = element_rect(fill = 'white',color = NA),
          legend.key.width = unit(1.8,'cm'),
          legend.key.height = unit(0.5,'cm'),
          legend.text = element_text(size = 11,color = 'black'),
          legend.spacing.x = unit(0.25,'cm'),
          legend.spacing.y = unit(0.1,'cm'),
          legend.margin = margin(t = 7,r = 12,b = 7,l = 12),
          
          plot.margin = margin(t = 8,r = 10,b = 8,l = 8)
        )
      
      safe_crit <- gsub('[^A-Za-z0-9_-]','_',as.character(crit))
      safe_variable <- gsub('[^A-Za-z0-9_-]','_',as.character(v))
      safe_rolling <- gsub('[^A-Za-z0-9_-]','_',as.character(r))
      
      fname <- paste0('plot_',safe_variable,'_',safe_crit,'_roll',safe_rolling,'_vintage',vintage,'.pdf')
      
      ggsave(filename = file.path('output','forecasting','figures',fname),
             plot     = p,
             device   = grDevices::cairo_pdf,
             width    = 12,
             height   = 7.5,
             units    = 'in')
    }
  }
}