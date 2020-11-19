# 2_Step_Estimator_DFM
These functions estimate a DFM with a 2-step estimator of Doz et. al. (2011).

Let us define <img src="https://render.githubusercontent.com/render/math?math=x_{t}=(x_{1,t},x_{2,t},\dots,x_{n,t})', t_m=1,2,\dots,T_m"> as our standardized variables. Then, the factor model has the following representation:

<img src="https://render.githubusercontent.com/render/math?math=x_{t} = \Lambda f_{t}+\epsilon_{t}; \quad \epsilon_{t}\sim \mathbb{N}(0,\Sigma_{\epsilon_{t}}),">

where <img src="https://render.githubusercontent.com/render/math?math=\Lambda"> is an <img src="https://render.githubusercontent.com/render/math?math=n \times r"> matrix of factor loadings, <img src="https://render.githubusercontent.com/render/math?math=\epsilon_{t}"> is the idiosyncratic component <img src="https://render.githubusercontent.com/render/math?math=f_{t_m}=(f_{1,t_m},f_{2,t_m},\dots,f_{r,t_m})'"> represents the unobserved common factors following a vector autoregression process (VAR) as follows:

<img src="https://render.githubusercontent.com/render/math?math=f_{t_m} = \sum_{i=1}^{p}A_i f_{t_m-i}+B\eta_{t}; \quad \eta_{t_m}\sim \mathbb{N}(0,I_{q}),">

where <img src="https://render.githubusercontent.com/render/math?math=B"> is an <img src="https://render.githubusercontent.com/render/math?math=r \times q"> matrix of full rank <img src="https://render.githubusercontent.com/render/math?math=q"> with <img src="https://render.githubusercontent.com/render/math?math=q \leqslant r">, <img src="https://render.githubusercontent.com/render/math?math=A_1,A_2,\dots,A_p"> are <img src="https://render.githubusercontent.com/render/math?math=r\times r"> matrices of autoregressive coefficients, and <img src="https://render.githubusercontent.com/render/math?math=\eta_{t}"> is the <img src="https://render.githubusercontent.com/render/math?math=q"> dimensional vector of common shocks, which follows a white-noise process.

Set of equations defined above can be estimated using the two step estimation procedure of Doz et. al. (2011) as follows:
