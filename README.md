# Bayesian Ensemble Trees for Causal Inference on Longitudinal Data (LongBet)

## About

This package implements the Bayesian Ensemble Trees for Causal Inference on 
Longitudinal Data for estimating time-varying conditional average treatment effect 
estimation; the manuscript will be available soon. 
This approach builds on the methodology behind Bayesian Causal Forests outlined 
in [Hahn et al.](https://projecteuclid.org/euclid.ba/1580461461) (2020) and 
incorporates several improvements to Bayesian Additive Regression Trees 
implemented by [He et al.](http://proceedings.mlr.press/v89/he19a.html) (2019).

This package is based on the source code of the [XBCF](https://github.com/socket778/XBCF) 
package.

## Installation
It can be installed from GitHub directly using devtools package in R. The CRAN version will be submitted soon.

```R
library(devtools)
install_github("google/longbet")
```

## Usage
```R
longBet(y, x, z, pcat, 
num_sweeps = 60, num_burnin = 20,
num_trees_pr = 50, num_trees_trt = 20,
mtry = 0L, n_min = 1L,
sig_knl = 1, lambda_knl = 5)
```

### Arguments
`y`: An n by t matrix of outcome variables.

`x`: n by p input matrix of covariates. (If the covariates matrix is different for the prognostic and treatment term, please use longBet_full).

`z`: An n by t matrix of treatment assignments.

`t`: time variable (post-treatment time for treatment term will be inferred based on input t and z).

`pcat`: The number of categorical inputs in matrix x.

`num_sweeps`: The total number of sweeps over data (default is 60).

`num_burnin`: The number of burn-in sweeps (default is 20).

`num_trees_pr`: The number of trees in the prognostic forest (default is 50).

`num_trees_trt`: The number of trees in the treatment forest (default is 20).

`mtry`: number of variables to be sampled as split candidate per tree.

`n_min`: The minimum node size. (default is 1)

`sig_knl`: variance parameter for squared exponential kernel (default is 1).

`lambda_knl`: lengthscale parameter for squared exponential kernel (default is 5).

### See Also
'longBet_full' for fitting with more hyperparameters.

'predict' will be available in the next update.

### Example
```R
require(longBet)

set.seed(1)
n <- 100
t1 <- 4
t0 <- 3

# generate dcovariates
x1 <- rnorm(n)
x2 <- sample(1:3,n,replace=TRUE,prob = c(0.4,0.3,0.3))
# TODO: memory bug occurs when there's no categorical variable 
x <- cbind(x1, x2)

# untreated outcome
mu <- outer(x1 * x2 , rnorm(t1, 5), '*')
# treatment effect
te <- outer(x1 + x2, rnorm(t1, 1), '*')

# generate observations
z <- rbinom(n,1,0.6)
y0 <- mu + 0.2 * sd(mu) * matrix(rnorm(n * t1), n, t1)
y1 <- y0 + te
y <- z * y1 + (1 - z) * y0

z_mat <- cbind(matrix(0, n, (t0 - 1)),  matrix(rep(z, t1 - t0 + 1), n, t1 - t0 + 1))

t_longbet <- proc.time()
longbet.fit <- longBet(y = y, x = x, z = z_mat, t = 1:t1, p_cat = 1,
num_trees_pr =  50, num_trees_trt = 50)
tau_hat_longBet <- apply(longbet.fit$tauhats.adjusted, c(1, 2), mean)
t_longbet <- proc.time() - t_longbet

print(paste0("longBet RMSE: ", sqrt(mean((as.vector(tau_longBet[, t0:t1]) - as.vector(te[,t0:t1]))^2))))
print(paste0("longBet runtime: ", round(as.list(t_longbet)$elapsed,2)," seconds"))

```