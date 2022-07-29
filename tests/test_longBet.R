# simple demonstration of longBet with default parameters
library(longBet)

#### 1. DATA GENERATION PROCESS
n <- 1000 # number of observations
t1 <- 20
t0 <- 1

# generate dcovariates
x1 <- rnorm(n)
x2 <- rnorm(n)
x3 <- rbinom(n,1,0.2)
x4 <- rbinom(n,1,0.7)
x5 <- factor(sample(1:3,n,replace=TRUE,prob = c(0.1,0.6,0.3)))
x5_reshape <- data.frame(x5 = x5)
x5_reshape <- model.matrix(~ . + 0, data = x5_reshape, 
contrasts.arg = lapply(x5_reshape, contrasts, contrasts=FALSE))

x <- as.matrix(cbind(x1,x2,x3,x4,x5_reshape))

# define average time-varying treatment effect
t <- 1:t1
ft <- log(2.5*t)*(t < 10) +  8*exp(-0.1*t)*(t >= 10)
# define heterogeneous treatment effects
tau <- 2 + 0.5*x[,2]*(2*x[,4]-1)
# time-varyiing heterogeneous treatment effect
tau_mat <- outer(tau, ft, "*")


# ## define prognostic function (RIC)
mu = function(x){
  lev = c(-0.5,0.75,0)
  result = 1 + x[,1]*(2*x[,4] - 2*(1-x[,4])) + lev[x5]
  return(result)
}

# compute propensity scores and treatment assignment
pi <- pnorm(-0.5 + mu(x) - x[,4] + 0.5*x[,2],0,3)
z <- rbinom(n,1,pi)
z_mat <- matrix(rep(z, t1), n, t1)

# generate outcome variable
# Ey = mu(x) + tau*z
Ey = tau_mat * z_mat
sig = 0.25*sd(Ey)
y = Ey + matrix(sig*rnorm(n*t1), n, t1)

# If you didn't know pi, you would estimate it here
pi_mat <- as.matrix(rep(pi, t1), n, t1)
pi_vec <- as.vector(pi_mat)

# Turn input into nt*1 vector and nt*p matrix
y_vec <- as.vector(y)
x_con <- c()
x_mod <- c()
t_vec <- c()
for (i in 1:t1){
  t_vec <- c(t_vec, rep(i, nrow(x)))
  x_con <- rbind(x_con, cbind(x, i))
  x_mod <- rbind(x_mod, cbind(x, max(0, i - t0 + 1)))
}
z_vec <- as.vector(z_mat)

ft_mat <- t(matrix(rep(ft, n), t1, n))
ft_vec <- as.vector(ft_mat)

#### 2. longBet
# run longBet
xbcf.fit <- longBet(y_vec, z_vec, x_con, x_mod, t_vec, pi_vec, 
n_trees_con =  0, n_trees_mod = 40,
alpha_mod = 0.95, beta_mod = 1.25,
pcat_con = 5,  pcat_mod = 5)

# get treatment individual-level estimates
tauhats <- getTaus(xbcf.fit)
betahats <- xbcf.fit$beta_draws
print(dim(betahats))

# main model parameters can be retrieved below
#print(xbcf.fit$model_params)

# compare results to inference
print(paste0("longBet RMSE with ft: ", sqrt(mean((tauhats*ft_vec - as.vector(tau_mat))^2))))
print(paste0("longBet pct error: ", mean( 100 * abs(tauhats * ft_vec - as.vector(tau_mat)) / as.vector(tau_mat) )))

print(paste0("longBet RMSE with estimated beta: ", sqrt(mean((tauhats*betahats - as.vector(tau_mat))^2))))
print(paste0("longBet pct error: ", mean( 100 * abs(tauhats * betahats - as.vector(tau_mat)) / as.vector(tau_mat) )))
