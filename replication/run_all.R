# ============================================================
# Replication code: Broken Adaptive Ridge
# Run after you opened "replication.Rproj"
# ============================================================

# 1. Load "here" package for directory
if (!requireNamespace('here', quietly = TRUE)) {
  install.packages('here', dependencies = TRUE)
}
library('here', character.only = TRUE)

# 2. Identify the folder of the project
here::i_am('run_all.R')

# 3. Create folders for output
dir.create(
  here::here('output', 'mc', 'csv'),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  here::here('output', 'mc', 'tables_tex'),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  here::here('output', 'mc', 'figures'),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  here::here('output', 'forecasting'),
  recursive = TRUE,
  showWarnings = FALSE
)

# 4. Run the Monte Carlo replication code
source(here::here('code', 'simulations.R'))

# 5. Run the Empirical replication code
source(here::here('code', 'empirical.R'))

# 6. Register the used version
writeLines(
  capture.output(sessionInfo()),
  here::here('output', 'sessionInfo.txt')
)