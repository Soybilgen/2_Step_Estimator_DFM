# Load FRED-MD Data Set
load(url("https://github.com/Soybilgen/2_Step_Estimator_DFM/raw/master/Example.RData"))

# Load necessary functions
source("https://raw.githubusercontent.com/Soybilgen/2_Step_Estimator_DFM/master/FactorExtraction.R")
source("https://raw.githubusercontent.com/Soybilgen/2_Step_Estimator_DFM/master/center.R")
source("https://raw.githubusercontent.com/Soybilgen/2_Step_Estimator_DFM/master/kalman_filter_diag.R")
source("https://raw.githubusercontent.com/Soybilgen/2_Step_Estimator_DFM/master/kalman_smoother_diag.R")
source("https://raw.githubusercontent.com/Soybilgen/2_Step_Estimator_DFM/master/kalman_update_diag.R")
source("https://raw.githubusercontent.com/Soybilgen/2_Step_Estimator_DFM/master/ricSW.R")
source("https://raw.githubusercontent.com/Soybilgen/2_Step_Estimator_DFM/master/smooth_update.R")

# Load necessary packages if you don't have
# install.packages("pracma") # Numeric Analysis Package
# install.packages("RSpectra") # Solvers for Large-Scale Eigenvalue and SVD Problems

# Determine the DFM specification
q <- 2 # dynamic factors
r <- 5 # static factors
p <- 4 # VAR lag

# Estimate the DFM
result <- FactorExtraction(x, q, r, p)

# Obtain dynamic factors
Dynamic_Factors <- result$F
