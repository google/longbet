# simple demonstration of longbet with default parameters
library(dbarts)
library(dplyr)
library(ggplot2)
library(tidyr)
library(longbet)
# DATA GENERATION PROCESS -------------------------------------------------


set.seed(1)
n <- 500
t1 <- 12
t0 <- 6

# generate dcovariates
x1 <- rnorm(n)
x2 <- rnorm(n)
x3 <- rnorm(n)
x4 <- rbinom(n,1,0.5)
x5 <- factor(sample(1:3,n,replace=TRUE,prob = c(0.4,0.3,0.3)))
x <- cbind(x1, x2, x3, x4, x5)

# define average time-varying treatment effect
post_t <- 1:(t1 - t0 + 1) 
beta_t <- dgamma(post_t, 2, 1)
# define heterogeneous treatment effects
tau <- 1 + 2 * x[,2] * x[,5]
# time-varyiing heterogeneous treatment effect
tau_mat <- 1 + 5 * outer(tau, beta_t , "*")


# ## define prognostic function (RIC)
mu = function(x){
  lev = c(2, -1, -4)
  result = -6 + lev[x[,5]] + 6 * abs(x[,3] - 1)
  return(result)
}

# generate auto regressive time trend
mu_t <- matrix(NA, n, t1)
mu_t[, 1] <- 1
alpha_t <- rnorm(t1, 0, 0.1) # generate pct change param
for (i in 2:t1){
  mu_t[, i] <- (1 + rnorm(n, alpha_t[i], 0.1)) * mu_t[, i-1]
}
mu_mat <- mu_t * matrix(rep(mu(x), t1), n, t1)

# compute propensity scores and treatment assignment
# pi <- pnorm(-0.5 + mu(x) - x[,4] + 0.5*x[,2],0,3)]
s <- sd(mu(x))
pi <- 0.8 * pnorm(1 * mu(x) / s - 0.5 * x[,1]) + 0.05 + runif(n) / 10
z <- rbinom(n,1,pi)
# z_mat <- matrix(rep(z, t1 - t0 + 1), n, t1 - t0 + 1)
z_mat <- matrix(rep(z, t1 - t0 + 1), n, t1 - t0 + 1)

# generate outcome variable
Ey = mu_mat
# Ey = matrix(0, nrow(mu_mat), ncol(mu_mat))
Ey[,t0:t1] <- Ey[, t0:t1] + tau_mat * z_mat
# Ey = tau_mat * z_mat
sig = 0.2*sd(Ey)
y = Ey + matrix(sig*rnorm(n*t1), n, t1)

# If you didn't know pi, you would estimate it here
pi_mat <- as.matrix(rep(pi, t1), n, t1)
pi_vec <- as.vector(pi_mat)



# Turn input into nt*1 vector and nt*p matrix
y_vec <- as.vector(y)
x_bart <- c()
t_vec <- c()
for (i in 1:t1){
  t_vec <- c(t_vec, rep(i, nrow(x)))
  x_bart <- rbind(x_bart, cbind(max(0, i - t0 + 1), x))
}
z_vec <- c(rep(0, n * (t0 - 1)), as.vector(z_mat))

expand_z_mat <- cbind(matrix(0, n, (t0 - 1)), z_mat)


# longbet -----------------------------------------------------------------
t_longbet <- proc.time()
longbet.fit <- longbet(y = y, x = x, z = expand_z_mat, t = 1:t1,
                       num_trees_pr =  50, num_trees_trt = 50 ,
                       num_sweeps = 100, num_burnin = 20,
                       pcat = ncol(x) - 3,  sig_knl = 1, lambda_knl = 1)
# TODO: lambda_knl is quite sensitve, need better understanding

# assume all unit get treated at t0 for test set to get CATE
z_test <- c(rep(0, t0 - 1), rep(1, t1 - t0 + 1)) %>% rep(times = n) %>%  matrix(nrow = t1, ncol = n) %>% t

longbet.pred <- predict.longbet(longbet.fit, x, z_test)
# longbet.pred <- predict.longbet(longbet.fit, x, 1:t1)
mu_hat_longbet <- apply(longbet.pred$muhats, c(1, 2), mean)
tau_hat_longbet <- apply(longbet.pred$tauhats, c(1, 2), mean)
tau_longbet <- tau_hat_longbet[,t0:t1]
t_longbet <- proc.time() - t_longbet

# results -----------------------------------------------------------------
# check ate
ate <- tau_mat %>% colMeans
ate_longbet <- tau_longbet %>% colMeans

print(paste0("longbet CATE RMSE: ", round( sqrt(mean((as.vector(tau_longbet) - as.vector(tau_mat))^2)), 2 ) ))
print(paste0("longbet ate rmse: ", round( sqrt(mean((ate_longbet - ate)^2)), 2)))
print(paste0("longbet runtime: ", round(as.list(t_longbet)$elapsed,2)," seconds"))


# # visualize ---------------------------------------------------------------
# ATE
ate_df <- data.frame(
  time = t0:t1,
  true = ate,
  longbet = ate_longbet
)

ate_df %>%
  gather("method", "ate", -time) %>%
  ggplot(aes(time, ate)) +
  geom_line(aes(color = method)) +
  ylab(labs(title = "Average Treatment Effect"))


# CATE
cate_df <- data.frame(
  true = as.vector(t(tau_mat)),
  lonbet = as.vector(t(tau_longbet)),
  time = rep(c(t0:t1), nrow(tau_mat)),
  id = as.vector(sapply(1:nrow(tau_longbet), rep, (t1 - t0 + 1)))
)

cate_plot <-  cate_df %>%
  gather("method", "cate", -time, -id) %>%
  ggplot() +
  geom_line(aes(time, cate, group = id, color = id)) +
  facet_wrap(~method)
plot(cate_plot)

# # muhat
# y0hat <- longbet.pred$muhats %>% apply(MARGIN = c(1, 2), mean) %>% data.frame
# y0_df <- data.frame(
#   true = as.vector(t(mu_mat)),
#   lonbet = as.vector(t(y0hat)),
#   time = rep(c(1:t1), nrow(mu_mat)),
#   id = as.vector(sapply(1:nrow(y0hat), rep, (t1)))
# )
# 
# y0_plot <-  y0_df %>%
#   gather("method", "cate", -time, -id) %>%
#   ggplot() +
#   geom_line(aes(time, cate, group = id, color = id)) +
#   facet_wrap(~method)
# plot(y0_plot)


# check convergence 
n_sweeps <- dim(longbet.pred$tauhats)[3]
rmse <- rep(NA, n_sweeps)
for(i in 1:n_sweeps){
  rmse[i] <- round( sqrt(mean((as.vector(longbet.pred$tauhats[,t0:t1,i]) - as.vector(tau_mat))^2)), 2)
}
rmse_df <- data.frame(rmse = rmse, sweeps = 1:n_sweeps)
rmse_trace <- rmse_df %>%
  ggplot(aes(sweeps, rmse)) + geom_point() + geom_line()
plot(rmse_trace)

# # sigma convergence
# sigma0 <- longbet.fit$sigma0_draws %>% as.vector
# sigma1 <- longbet.fit$sigma1_draws %>% as.vector
# sigma_df <- data.frame(iter = 1:length(sigma0), sigma0 = sigma0, sigma1 = sigma1)
# sigma_trace <- sigma_df %>% gather(key = "param", value = "value", -iter) %>%
#   ggplot(aes(iter, value, color = param)) + geom_point() + geom_line()
# plot(sigma_trace)
# 
# 
# param_df <- data.frame(iter = 1:longbet.fit$model_params$num_sweeps, 
#                        a = as.vector(longbet.fit$a_draws),
#                        b0 = as.vector(longbet.fit$b_draws[,1]),
#                        b1 = as.vector(longbet.fit$b_draws[,2]))
# param_trace <- param_df %>% gather(key = "param", value = "value", -iter) %>%
#   ggplot(aes(iter, value, color = param)) + geom_point() + geom_line()
# plot(param_trace)




