## These functions estimate a DFM with a 2-step estimator of Doz et. al. (2011).

These functions are based on the MATLAB functions of Giannone et. al. (2008).

### Simple Introduction

Let us define <img src="https://render.githubusercontent.com/render/math?math=x_{t}=(x_{1,t},x_{2,t},\dots,x_{n,t})', t_m=1,2,\dots,T_m"> as our standardized variables. Then, the factor model has the following representation:

<img src="https://render.githubusercontent.com/render/math?math=x_{t} = \Lambda f_{t}+\epsilon_{t}; \quad \epsilon_{t}\sim \mathbb{N}(0,\Sigma_{\epsilon_{t}}),">

where <img src="https://render.githubusercontent.com/render/math?math=\Lambda"> is an <img src="https://render.githubusercontent.com/render/math?math=n \times r"> matrix of factor loadings, <img src="https://render.githubusercontent.com/render/math?math=\epsilon_{t}"> is the idiosyncratic component, <img src="https://render.githubusercontent.com/render/math?math=f_{t_m}=(f_{1,t_m},f_{2,t_m},\dots,f_{r,t_m})'"> represents the unobserved common factors following a vector autoregression process (VAR) as follows:

<img src="https://render.githubusercontent.com/render/math?math=f_{t_m} = \sum_{i=1}^{p}A_i f_{t_m-i}+B\eta_{t}; \quad \eta_{t_m}\sim \mathbb{N}(0,I_{q}),">

where <img src="https://render.githubusercontent.com/render/math?math=B"> is an <img src="https://render.githubusercontent.com/render/math?math=r \times q"> matrix of full rank <img src="https://render.githubusercontent.com/render/math?math=q"> with <img src="https://render.githubusercontent.com/render/math?math=q \leqslant r">, <img src="https://render.githubusercontent.com/render/math?math=A_1,A_2,\dots,A_p"> are <img src="https://render.githubusercontent.com/render/math?math=r\times r"> matrices of autoregressive coefficients, and <img src="https://render.githubusercontent.com/render/math?math=\eta_{t}"> is the <img src="https://render.githubusercontent.com/render/math?math=q"> dimensional vector of common shocks, which follows a white-noise process.

Set of equations defined above is estimated using the two step estimation procedure of Doz et. al. (2011) as follows:
1. All variables are standardized in the first step.
2. First <img src="https://render.githubusercontent.com/render/math?math=r"> principal components are extracted from the balanced part of the data set where all observations are present and obtain the initial factor estimates and parameters using the function `ricSW(z, q, r, p)`. <img src="https://render.githubusercontent.com/render/math?math=z_{t}"> represents the balanced part of the data set.
3. Before moving to the Kalman smoother part, the ragged edge part of the data set are incorporated into the procedure by assigning an extremely large value to the variance of the idiosyncratic component where there is missing observations and replacing missing values in <img src="https://render.githubusercontent.com/render/math?math=x_{t}"> with arbitrary values.
4. Factors are re-estimated using one run of the Kalman filter function `kalman_filter_diag(y, A, C, Q, R, init_x, init_V, model)` and one run of the Kalman smoother function `kalman_smoother_diag(y, A, C, Q, R, init_x, init_V, model)` while incorporating the unbalanced part of the data set.

### Simple Example

Load FRED-MD Data Set
```
load(url("https://github.com/Soybilgen/2_Step_Estimator_DFM/raw/master/Example.RData"))
```
Load necessary functions
```
source("https://raw.githubusercontent.com/Soybilgen/2_Step_Estimator_DFM/master/FactorExtraction.R")
source("https://raw.githubusercontent.com/Soybilgen/2_Step_Estimator_DFM/master/center.R")
source("https://raw.githubusercontent.com/Soybilgen/2_Step_Estimator_DFM/master/kalman_filter_diag.R")
source("https://raw.githubusercontent.com/Soybilgen/2_Step_Estimator_DFM/master/kalman_smoother_diag.R")
source("https://raw.githubusercontent.com/Soybilgen/2_Step_Estimator_DFM/master/kalman_update_diag.R")
source("https://raw.githubusercontent.com/Soybilgen/2_Step_Estimator_DFM/master/ricSW.R")
source("https://raw.githubusercontent.com/Soybilgen/2_Step_Estimator_DFM/master/smooth_update.R")
```
Load necessary packages
```
install.packages("pracma") # Numeric Analysis Package
install.packages("RSpectra") # Solvers for Large-Scale Eigenvalue and SVD Problems
```
Determine the DFM specification
```
q <- 2 # dynamic factors
r <- 5 # static factors
p <- 4 # VAR lag
```
Estimate the DFM and obtain dynamic factors
```
result <- FactorExtraction(x, q, r, p)
Dynamic_Factors <- result$F
```
### References

Doz, C., Giannone, D., & Reichlin, L. (2011). A two-step estimator for large approximate dynamic factor models based on Kalman filtering. Journal of Econometrics, 164(1), 188-205.

Giannone, D., Reichlin, L., & Small, D. (2008). Nowcasting: The real-time informational content of macroeconomic data. Journal of Monetary Economics, 55(4), 665-676.
