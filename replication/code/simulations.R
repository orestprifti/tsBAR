rm(list=ls())

#-------------------------------------------------------------------------------
# 1. FUNCTIONS
#-------------------------------------------------------------------------------
source(file.path('code', 'functions.R'))

#-------------------------------------------------------------------------------
# 2. SETUP
#-------------------------------------------------------------------------------
nmc           <- 20                                                             # number of Monte Carlo simulations
criteria      <- 'ebic'                                                         # information criterion to use
p             <- 1                                                              # VAR(p)
save_heatmaps <- FALSE                                                          # TRUE if you want to save coefficient for heatmaps, FALSE if you don't

all_mc_results <- data.frame()                                                  # initialize table for the results of each MC design

mc_designs <- expand.grid(
  T_  = c(50, 100, 500),
  r   = c('1', '2', '5', '10'),
  rho = c(0), #0.5),
  K   = c(20), #40, 80, 160, 320),
  stringsAsFactors = FALSE
  ) %>%
  rowwise() %>%
  dplyr::mutate(r = eval(parse(text = gsub('K', as.character(K), r))))

#-------------------------------------------------------------------------------
# 3. MONTE CARLO SIMULATIONS
#-------------------------------------------------------------------------------
cat('Monte Carlo running...\n')
pb <- txtProgressBar(min = 0, max = nrow(mc_designs), style = 3,                # progress bar
                     width = max(50,nrow(mc_designs)))
for(i in c(1:nrow(mc_designs))) {
  setTxtProgressBar(pb, i)
  
  current <- mc_designs[i,]
  rho     <- current$rho
  K       <- current$K
  T_      <- current$T_                                                         # sample size: 50, 100, 500
  r       <- current$r                                                          # sparsity level: 1, 2, 5, 10

  # Monte Carlo
  info_parallel <- c('K','T_','rho','r','criteria','p')
  cl            <- makeCluster(11 - 1)                               # use all the cores leaving one free
  registerDoParallel(cl)
  clusterExport(cl, varlist = info_parallel)
  clusterEvalQ(cl, {source('code/functions.R')})                                # upload in each core the functions
  mc_result <- foreach(imc = 1:nmc, .combine = rbind, .export = info_parallel) %dopar% {
    set.seed(imc)
    DGP <- cmp_dgp(K = K, T_ = T_+1, r = r, rho = rho)
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
    BAR_est          <- bar_ccd_var(y = y, x = x, crit = criteria,
                                    intercept = TRUE, standardize = TRUE)
    BAR              <- BAR_est$beta                                            # coefficients without intercept

    # Adaptive LASSO with Ridge initialization
    aLASSO_ridge_est <- var.lasso(y = y, x = x, p = p, adaptive = 'ridge',
                                  crit = criteria, standardize = TRUE,
                                  intercept = TRUE)
    aLASSO_ridge     <- aLASSO_ridge_est$coefficients[-1,]                      # take the coefficients discarding the intercept

    # Adaptive LASSO with LASSO initialization
    aLASSO_lasso_est <- var.lasso(y = y, x = x, p = p, adaptive = 'lasso',
                                  crit = criteria, standardize = TRUE,
                                  intercept = TRUE)
    aLASSO_lasso     <- aLASSO_lasso_est$coefficients[-1,]                      # take the coefficients discarding the intercept

    # LASSO
    LASSO_est        <- var.lasso(y = y, x = x, p = p, alpha = 1, adaptive = 'none',
                                  crit = criteria, standardize = TRUE, intercept = TRUE)
    LASSO            <- LASSO_est$coefficients[-1,]                             # take the coefficients discarding the intercept

    # SCAD
    SCAD_est         <- var.scad(y = y, x = x, crit = criteria,
                         standardize = TRUE, intercept = TRUE)
    SCAD             <- SCAD_est$coefficients[-1,]                              # take the coefficients discarding the intercept


    models <- list('BAR'          = BAR,
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

      # 1. RMSFE one-step ahead (forecast accuracy)
      fcast <- a0 + x_pred %*% mod                                              # point forecast
      err   <- (a0 + x_pred %*% mod) - y_pred                                   # forecast error
      RMSFE <- mean(err^2)                                                      # root mean squared forecast error
      
      # 2. RMSE of A (estimation accuracy)
      RMSE <- sum((mod - A)^2)                                                  # root mean squared error
      
      # 3. Correct sparsity
      correct_sparsity <- all(supp_A == supp_Ahat)

      # 4. True model included
      true_model_included <- all(supp_Ahat[supp_A])

      # 5. Fraction of relevant variables included
      frac_relevant_included <- sum(supp_Ahat[supp_A]) / sum_A

      # 6. False Positive Rate (fpr) and False Negative Rate (fnr)
      TP  <- sum( supp_A  &  supp_Ahat)                                         # true positive
      FP  <- sum(!supp_A  &  supp_Ahat)                                         # false positive
      TN  <- sum(!supp_A  & !supp_Ahat)                                         # true negative
      FN  <- sum( supp_A  & !supp_Ahat)                                         # false negative
      fpr <- FP / (FP + TN)                                                     # false positive rate
      fnr <- FN / (FN + TP)                                                     # false negative rate
      
      # 7. Number of variables included
      num_variables_included <- sum(supp_Ahat) / K

      # Collect results
      imc_result <- rbind(imc_result,
                          data.frame('model'                  = model_name,
                                     'RMSFE'                  = RMSFE,
                                     'RMSE'                   = RMSE,
                                     'correct_sparsity'       = correct_sparsity,
                                     'true_model_included'    = true_model_included,
                                     'frac_relevant_included' = frac_relevant_included,
                                     'fpr'                    = fpr,
                                     'fnr'                    = fnr,
                                     'num_variables_included' = num_variables_included
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
    dplyr::mutate(
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
            file = file.path('output','mc','csv',paste0('K',K,'_T',T_,'_r',r,'_rho',rho,'_',criteria,'_mc',nmc,'.csv')),
            row.names = FALSE)
}
close(pb)

# Save results of all MC designs
write.csv(all_mc_results,
          file = file.path('output','mc','csv',paste0('FINAL_MC',nmc,'_',criteria,'.csv')),
          row.names = FALSE)

#-------------------------------------------------------------------------------
# 3. MONTE CARLO TABLES
#-------------------------------------------------------------------------------
tab_mc_results <- read.csv(file.path('output','mc','csv',paste0('FINAL_MC',nmc,'_',criteria,'.csv')),
                          sep = ',',
                          dec = '.'
                          ) %>% data.frame()

# check duplicates
if (all(duplicated(tab_mc_results))) {
  error('There are duplicates in the Monte Carlo results. Please check the csv file!')
}

# rename performance measures
tab_mc_results <- tab_mc_results %>%
  rename('Correct sparsity'            = correct_sparsity,
         'True model included'         = true_model_included,
         'Share of relevant variables' = frac_relevant_included,
         'False Positive Rate'         = fpr,
         'False Negative Rate'         = fnr,
         'Number of variables'         = num_variables_included
         )

models_order   <- c('BAR','aLASSO_ridge','aLASSO_lasso','LASSO','SCAD')
T_order        <- unique(tab_mc_results$T) %>% sort()
measures_order <- setdiff(colnames(tab_mc_results), c('K','T','r','rho','model'))
measures_order <- c('RMSFE','RMSE',setdiff(measures_order,c('RMSFE','RMSE','Number of variables')),'Number of variables')

measure_headers <- sapply(seq_along(measures_order), function(i) {
  prefix <- if (i == 1) '$K$' else '{}'
  paste0('\\addlinespace\n', prefix, ' & \\multicolumn{', 
         length(models_order) * length(T_order), 
         '}{c}{\\underline{\\textit{', measures_order[i], '}}}\\\\ \n')
})

measures_max <- c('Correct sparsity','True model included','Share of relevant variables')
measures_min <- c('RMSFE','RMSE','False Positive Rate','False Negative Rate')

n_models <- length(unique(tab_mc_results$model))
mc_spec   <- expand.grid(
  rho = unique(tab_mc_results$rho),
  r   = unique(tab_mc_results$r),
  stringsAsFactors = FALSE
  )

for(i in c(1:nrow(mc_spec))) {
  current <- mc_spec[i,]
  rho     <- current$rho
  r       <- current$r
  
  tab_long <-  tab_mc_results %>%
    filter(r == !!r, rho == !!rho) %>%
    dplyr::select(-c('r','rho')) %>% 
    tidyr::pivot_longer(
      cols      = all_of(measures_order),
      names_to  = 'measure',
      values_to = 'value'
    )
  
  tab_wide <- tab_long %>%
    dplyr::mutate(
      model   = factor(model, levels = models_order),
      T       = factor(T, levels = T_order),
      measure = factor(measure, levels = measures_order)
    ) %>%
    arrange(measure, K, T, model) %>%
    tidyr::unite('mod_T', model, T, sep = '_T') %>%
    tidyr::pivot_wider(names_from = mod_T, values_from = value)
  
  tab_wide <- tab_wide %>%
    group_by(measure, K) %>%
    dplyr::mutate(across(
      contains('_T'),
      .fns = \(x) {
        cur_T       <- sub('.*_T', '', cur_column())
        same_T_cols <- grep(paste0('_T', cur_T, '$'), names(cur_data()), value = TRUE)
        vals        <- as.numeric(unlist(cur_data()[, same_T_cols]))
        vals        <- vals[!is.na(vals)]
        x_num       <- as.numeric(x)
        
        if (unique(measure) == 'Number of variables') {
          return(sprintf('%.3f', x_num))
        }
        
        if (length(vals) == 0) return(sprintf('%.3f', x_num))
        
        vals_rounded <- round(vals, 3)
        x_rounded <- round(x_num, 3)
        
        if (unique(measure) %in% measures_min) {
          best <- min(vals_rounded, na.rm = TRUE)
        } else if (unique(measure) %in% measures_max) {
          best <- max(vals_rounded, na.rm = TRUE)
        } else {
          best <- NA
        }
        
        n_best          <- sum(vals_rounded == best, na.rm = TRUE)
        best_unique     <- n_best == 1
        best_not_unique <- n_best < n_models | (n_best == n_models & best != 0)
        
        ifelse(
          !is.na(best) & !is.na(x_rounded) & x_rounded == best & best_not_unique & !best_unique,
          sprintf('\\textbf{%.3f}', x_num),
          ifelse(
            !is.na(best) & !is.na(x_rounded) & x_rounded == best & best_unique,
            sprintf('\\underline{\\textbf{%.3f}}', x_num),
            sprintf('%.3f', x_num)
          )
        )
      }
    )) %>%
    ungroup()
  
  tab <- tab_wide %>%
    dplyr::select(K,
                  paste0('BAR_T',          T_order),
                  paste0('aLASSO_ridge_T', T_order),
                  paste0('aLASSO_lasso_T', T_order),
                  paste0('LASSO_T',        T_order),
                  paste0('SCAD_T',         T_order)
                  )
  
  n_variables <- unique(tab$K)
  
  header1 <- paste0(
    '& \\multicolumn{',length(T_order),'}{c}{\\multirow{2}{*}{BAR}}',
    '& \\multicolumn{',length(T_order),'}{c}{Adaptive LASSO}',
    '& \\multicolumn{',length(T_order),'}{c}{Adaptive LASSO}',
    '& \\multicolumn{',length(T_order),'}{c}{\\multirow{2}{*}{LASSO}}',
    '& \\multicolumn{',length(T_order),'}{c}{\\multirow{2}{*}{SCAD}}\\\\'
  )
  
  header2 <- paste0(
    '& \\multicolumn{',length(T_order),'}{c}{}',
    '& \\multicolumn{',length(T_order),'}{c}{(1\\textsuperscript{st} step: Ridge)}',
    '& \\multicolumn{',length(T_order),'}{c}{(1\\textsuperscript{st} step: Lasso)}',
    '& \\multicolumn{',length(T_order),'}{c}{}',
    '& \\multicolumn{',length(T_order),'}{c}{}\\\\'
  )
  
  start   <- 2
  header3 <- paste(
    sapply(seq_along(models_order), function(i) {
      end <- start + length(T_order) - 1
      rule <- paste0('\\cmidrule(lr){', start, '-', end, '}')
      start <<- end + 1
      rule}), 
    collapse = ' ')
  
  header4 <- paste('$T$ &',
                   paste(rep(T_order, length(models_order)),collapse = ' & '),
                   '\\\\ \n\\midrule \n'
                   )
  
  rows_with_headers <- seq(0, length(n_variables) * (length(measures_order)-1), length(n_variables))
  
  addtorow         <- list()
  addtorow$pos     <- list(-1)
  addtorow$command <- c(paste('\\toprule\n', header1, '\n', header2, '\n', header3, '\n', header4))
  
  addtorow$pos <- c(addtorow$pos, as.list(rows_with_headers))
  addtorow$command <- c(addtorow$command, measure_headers)
  
  rows_with_midrule <- c(0, which(diff(as.numeric(tab_wide$measure)) != 0))
  addtorow$pos <- c(addtorow$pos, as.list(rows_with_midrule))
  addtorow$command <- c(addtorow$command, rep('\\addlinespace[0.5ex]\n', length(rows_with_midrule)))
  
  # Print results in a table
  tab_print        <- xtable(tab)
  align(tab_print) <- c(rep('l',2), rep('c',ncol(tab_print)-1)) 
  
  caption  <- paste0('MC simulations: $r$=',r,', $\\rho$=',rho,', ', toupper(criteria), '.')
  label    <- paste0('tab:MC_r', r, '_rho', rho)
  
  file_out <- paste0('tab_MC_r', r, '_rho', rho, '.tex')
  con      <- file(file.path('output','mc','tables_tex',file_out), open = 'wt')
  
  cat('\\begin{table}[H]\n',              file = con)
  cat('\\centering\n',                    file = con)
  cat(paste0('\\caption{',caption,'}\n'), file = con)
  cat(paste0('\\label{',label,'}\n'),     file = con)
  cat('\\resizebox{\\textwidth}{!}{%\n',  file = con)
  
  print(
    tab_print,
    type                   = 'latex',
    file                   = con,
    include.rownames       = FALSE,
    include.colnames       = FALSE,
    add.to.row             = addtorow,
    hline.after            = c(nrow(tab_print)),
    sanitize.text.function = identity,
    booktabs               = TRUE,
    comment                = FALSE,
    floating               = FALSE
  )
  
  cat('} % end resizebox\n', file = con)
  cat('\\end{table}\n', file = con)
  close(con)
}

#-------------------------------------------------------------------------------
# 3. PLOT WIN RATE
#-------------------------------------------------------------------------------
# calculate winners
winners <- tab_mc_results %>%
  dplyr::mutate(
    RMSFE                         = round(RMSFE, 3),
    RMSE                          = round(RMSE, 3),
    `Correct sparsity`            = round(`Correct sparsity`, 3),
    `True model included`         = round(`True model included`, 3),
    `Share of relevant variables` = round(`Share of relevant variables`, 3)
  ) %>% 
  group_by(K, r, T, rho) %>%
  dplyr::mutate(
    best_RMSFE = (RMSFE == min(RMSFE, na.rm = TRUE)),
    best_RMSE  = (RMSE  == min(RMSE,  na.rm = TRUE)),
    best_Correct_Sparsity = ifelse(`Correct sparsity` == max(`Correct sparsity`, na.rm = TRUE) &
                                     max(`Correct sparsity`, na.rm = TRUE)!=0,1,0),
    best_Model_Included = ifelse(`True model included` == max(`True model included`, na.rm = TRUE) &
                                   max(`True model included`, na.rm = TRUE)!=0,1,0),
    best_Share_Relevant = (`Share of relevant variables` == max(`Share of relevant variables`, na.rm = TRUE))
  ) %>%
  ungroup()

# sum wins for each model
sum_wins <- winners %>%
  group_by(model) %>%
  summarise(
    wins_RMSE             = sum(best_RMSE, na.rm = TRUE),
    wins_RMSFE            = sum(best_RMSFE, na.rm = TRUE),
    wins_Correct_Sparsity = sum(best_Correct_Sparsity, na.rm = TRUE),
    wins_Model_Inlcuded   = sum(best_Model_Included, na.rm = TRUE),
    wins_Share_Relevant   = sum(best_Share_Relevant, na.rm = TRUE),
    .groups = 'drop'
    )

# total number of Monte Carlo designs 
tot_designs <- winners %>%
  dplyr::distinct(K, T, r, rho) %>% 
  nrow()

# calculate win rate
data_win_rate <- sum_wins %>%
  dplyr::mutate(
    share_RMSFE = wins_RMSFE / tot_designs,
    share_RMSE = wins_RMSE / tot_designs,
    share_Correct_Sparsity = wins_Correct_Sparsity / tot_designs,
    share_Model_Inlcuded = wins_Model_Inlcuded / tot_designs,
    share_Relevant = wins_Share_Relevant / tot_designs
  ) %>% 
  dplyr::mutate(
    share_RMSFE            = round(share_RMSFE, 3),
    share_RMSE             = round(share_RMSE, 3),
    share_Correct_Sparsity = round(share_Correct_Sparsity, 3),
    share_Model_Inlcuded   = round(share_Model_Inlcuded, 3),
    share_Relevant         = round(share_Relevant, 3)
  ) %>% 
  dplyr::select(c('model', 'share_RMSFE', 'share_RMSE', 'share_Correct_Sparsity',
                  'share_Model_Inlcuded', 'share_Relevant'))

# prepare data for plotting win rates
win_rate <- data_win_rate %>%
  dplyr::select(
    model,
    share_RMSFE,
    share_RMSE,
    share_Correct_Sparsity,
    share_Model_Inlcuded,
    share_Relevant
  ) %>%
  tidyr::pivot_longer(
    cols      = -model,
    names_to  = 'metric',
    values_to = 'share'
  ) %>%
  dplyr::mutate(
    share  = 100 * share,
    metric = recode(
      metric,
      share_RMSFE            = 'RMSFE',
      share_RMSE             = 'RMSE',
      share_Correct_Sparsity = 'Correct sparsity',
      share_Model_Inlcuded   = 'True model included',
      share_Relevant         = 'Relevant coefficients\nselected'
    ),
    model = recode(
      model,
      'aLASSO_ridge' = 'aLASSO(ridge)',
      'aLASSO_lasso' = 'aLASSO(lasso)'
    ),
    model = factor(
      model,
      levels = c(
        'BAR',
        'aLASSO(ridge)',
        'aLASSO(lasso)',
        'LASSO',
        'SCAD'
      )
    )
  )

# create win-rate plot
plot_win_rate <- ggplot(
  win_rate,
  aes(
    x = share,
    y = factor(
      metric,
      levels = c(
        'Relevant coefficients\nselected',
        'True model included',
        'Correct sparsity',
        'RMSE',
        'RMSFE'
      )
    ),
    fill = model
  )
) +
  geom_col(
    position = position_dodge(width = 0.75,reverse = TRUE),
    width = 0.6,
    color = 'black',
    linewidth = 0
  ) +
  geom_text(
    aes(label = sprintf('%.1f', share)),
    position = position_dodge(width = 0.75,reverse=TRUE),
    hjust = -0.2,
    size = 3.5
  ) +
  scale_fill_manual(
    values = c(
      'BAR'           = 'royalblue',
      'aLASSO(ridge)' = 'orange',
      'aLASSO(lasso)' = 'forestgreen',
      'LASSO'         = 'firebrick',
      'SCAD'          = 'purple'
    )
  ) +
  guides(
    fill = guide_legend(override.aes = list(colour = NA))
  ) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(clip = 'off') +
  labs(
    x    = 'Average win rate (%)',
    y    = NULL,
    fill = NULL
  ) +
  theme_classic(base_family = 'Times New Roman') +
  theme(
    legend.position  = 'bottom',
    legend.direction = 'horizontal',
    axis.line.x      = element_line(linewidth = 0.2),
    axis.text.y      = element_text(size = 11, margin = margin(r = 15)),
    axis.text.x      = element_text(size = 10),
    axis.title.x     = element_text(size = 12),
    axis.ticks.x     = element_blank(),
    legend.text      = element_text(size = 10),
    # remove axis-y line
    axis.line.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    # vertical lines on axis-x breaks 
    panel.grid.major.x = element_line(
      colour    = 'grey75',
      linetype  = 'dotted',
      linewidth = 0.2
    ),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(t = 0, r = 60, b = 0, l = 10)
  )

# save win-rate plot
ggsave(
  plot     = plot_win_rate,
  filename = file.path('output','mc','figures','plot_win_rate.pdf'),
  device   = grDevices::cairo_pdf,
  width    = 7,
  height   = 6,
  units    = 'in'
)