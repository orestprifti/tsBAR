# ============================================================
# Replication code: Broken Adaptive Ridge
# Run after you opened "replication.Rproj"
# ============================================================

# 1. Packages
required <- c('here','Rcpp')

missing <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing) > 0L) {
  stop(
    'Install the required packages: ',
    paste(missing, collapse = ', ')
  )
}

# 2. Identify the folder of the project
here::i_am('run_all.R')

# 3. Create folder for output
dir.create(
  here::here('output', 'mc'),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  here::here('output', 'forecasting'),
  recursive = TRUE,
  showWarnings = FALSE
)

# 4. Load the BAR function
Rcpp::sourceCpp(
  here::here('src', 'bar.cpp')
)

# 5. Run the replication scripts
source(here::here('R', 'simulations.R'))

# 6. Register the used version
writeLines(
  capture.output(sessionInfo()),
  here::here('output', 'sessionInfo.txt')
)