# simple demonstration of longbet with default parameters
library(dbarts)
library(dplyr)
library(ggplot2)
library(tidyr)
library(longbet)
# DATA GENERATION PROCESS -------------------------------------------------


set.seed(1)
n <- 2000
t0 <- 6 # treatment start time
t1 <- 12 # observed response period
t2 <- 18 # predict response period
alpha = 0.05

# generate dcovariates
x1 <- rnorm(n)
x2 <- rnorm(n)
x3 <- rnorm(n)
x4 <- rbinom(n,1,0.5)
x5 <- factor(sample(1:3,n,replace=TRUE,prob = c(0.4,0.3,0.3)))
x <- cbind(x1, x2, x3, x4, x5)

# define average time-varying treatment effect
post_t <- 1:(t2 - t0 + 1) 
beta_t <- dgamma(post_t, 2, 1)
# define heterogeneous treatment effects
tau <- 1 + 2 * abs(x[,2] * x[,5])
# time-varyiing heterogeneous treatment effect
tau_mat <- 2 + 2 * outer(tau, beta_t , "*")
# expand tau_mat
tau_mat <- cbind(matrix(0, nrow = n, ncol = t0 - 1), tau_mat)

# ## define prognostic function (RIC)
mu = function(x){
  lev = c(2, -1, -4)
  result = -6 + lev[x[,5]] + 6 * abs(x[,3] - 1)
  return(result)
}

# generate auto regressive time trend
mu_t <- matrix(NA, n, t2)
mu_t[, 1] <- 1
alpha_t <- rnorm(t2, 0, 0.1) # generate pct change param
for (i in 2:t2){
  mu_t[, i] <- (1 + rnorm(n, alpha_t[i], 0.1)) * mu_t[, i-1]
}
mu_mat <- mu_t * matrix(rep(mu(x), t2), n, t2)

# compute propensity scores and treatment assignment
# pi <- pnorm(-0.5 + mu(x) - x[,4] + 0.5*x[,2],0,3)]
s <- sd(mu(x))
pi <- 0.8 * pnorm(1 * mu(x) / s - 0.5 * x[,1]) + 0.05 + runif(n) / 10
z <- rbinom(n,1,pi)
z_mat <- matrix(rep(z, t2 - t0 + 1), n, t2 - t0 + 1)

# generate outcome variable
y0 <- mu_mat
y1 <- mu_mat + tau_mat
y <- z * y1 + (1 -z) * y0 +  0.2 * sd(mu_mat) * matrix(rnorm(n*t2), n, t2)

# If you didn't know pi, you would estimate it here
pi_mat <- as.matrix(rep(pi, t1), n, t1)
pi_vec <- as.vector(pi_mat)

# get training data
ytrain = y[, 1:t1]
ztrain <- cbind(matrix(0, n, (t0 - 1)),  matrix(rep(z, t1 - t0 + 1), n, t1 - t0 + 1))

# Turn input into nt*1 vector and nt*p matrix
y_vec <- as.vector(ytrain)
x_bart <- c()
t_vec <- c()
for (i in 1:t1){
  t_vec <- c(t_vec, rep(i, nrow(x)))
  x_bart <- rbind(x_bart, cbind(max(0, i - t0 + 1), x))
}
# z_vec <- c(rep(0, n * (t0 - 1)), as.vector(z_mat))
z_vec <- as.vector(ztrain)


# longbet -----------------------------------------------------------------
t_longbet <- proc.time()
longbet.fit <- longbet(y = ytrain, x = x, z = ztrain, t = 1:t1,
                       num_sweeps = 60,
                       num_trees_pr =  40, num_trees_trt = 40,
                       pcat = ncol(x) - 3, sig_knl =  1, lambda_knl = 1)

z_test <- c(rep(0, t0 - 1), rep(1, t2 - t0 + 1)) %>% rep(times = n) %>%  matrix(nrow = t2, ncol = n) %>% t
longbet.pred <- predict.longbet(longbet.fit, x, z_test, sigma = NULL, lambda = NULL)
# mu_hat_longbet <- apply(longbet.pred$muhats, c(1, 2), mean)
# tau_hat_longbet <- apply(longbet.pred$tauhats, c(1, 2), mean)
# tau_longbet <- tau_hat_longbet[,t0:t2]
t_longbet <- proc.time() - t_longbet
# 
# ate_longbet_fit <- apply(longbet.pred$tauhats, c(2, 3), mean)[t0:t2, ]
# ate <- tau_mat %>% colMeans
# ate_longbet <- ate_longbet_fit %>% rowMeans

longbet.ate <- get_ate(longbet.pred, alpha = 0.05)
longbet.att <- get_att(longbet.pred, z = ztrain, alpha = 0.05)
longbet.cate <- get_cate(longbet.pred, alpha = 0.05)

print(paste0("longbet CATE RMSE in-sample: ", sqrt(mean((as.vector(longbet.cate$cate[,t0:t1]) - as.vector(tau_mat[, t0:t1]))^2))))
print(paste0("longbet CATE RMSE extrapolate: ", sqrt(mean((as.vector(longbet.cate$cate[,(t1 + 1):t2]) - as.vector(tau_mat[, (t1 + 1):t2]))^2))))

ate <- tau_mat %>% colMeans
print(paste0("longbet ATE RMSE in-sample: ", round( sqrt(mean((longbet.ate$ate[t0:t1] - ate[t0:t1])^2)), 2)))
print(paste0("longbet ATE RMSE extrapolate: ", round( sqrt(mean((longbet.ate$ate[(t1 + 1):t2] - ate[(t1 + 1):t2])^2)), 2)))

print(paste0("longbet runtime: ", round(t_longbet[3],2)," seconds"))

# visualize ---------------------------------------------------------------
colors <- c("black", "#FFA500", "#00BFFF")
labels <- c("True", "LongBet", "BART")
names(colors) <- labels
# ATE
ate_df <- data.frame(
  time = t0:t2,
  true = ate[t0:t2],
  longbet = longbet.ate$ate[t0:t2],
  lower = longbet.ate$interval[1, t0:t2],
  upper = longbet.ate$interval[2, t0:t2]
)

ate_plot <- 
  ggplot(ate_df , aes(x = time, y = true)) +
  geom_line(aes(y = true, color = "True")) +
  geom_line(aes(y = longbet, color = "LongBet")) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = "LongBet"), alpha = 0.15, fill = colors[2]) +
  geom_vline(xintercept = t1, linetype = "dashed", color = "grey") +
  labs(x = "Time", y = "ATE", color = "Legend") +
  scale_color_manual(name = "Legend", values = colors, labels = labels)
print(ate_plot)
readline(prompt="Press [enter] to continue")

# CATE
cate_df <- data.frame(
  true = as.vector(t(tau_mat[, t0:t2])),
  lonbet = as.vector(t(longbet.cate$cate[, t0:t2])),
  time = rep(c(t0:t2), nrow(tau_mat)),
  id = as.vector(sapply(1:n, rep, (t2 - t0 + 1)))
)

cate_plot <- cate_df %>%
  gather("method", "cate", -time, -id) %>%
  ggplot() +
  geom_line(aes(time, cate, group = id, color = id)) +
  geom_vline(xintercept = t1, linetype = "dashed", color = "grey") +
  facet_wrap(~method, ncol = 1)
print(cate_plot)

# analyse tauhat ----------------------------------------------------------
# library(reshape2)
# beta_draws = as.data.frame(t(longbet.pred$beta_draws))
# beta_draws$trial = rownames(beta_draws)
# mdat = melt(beta_draws, id.vars="trial")
# mdat$time = as.numeric(gsub("V", "", mdat$variable))
# 
# beta_draws_plot = ggplot(mdat, aes(x=time, y=value, group=trial)) +
#   theme_bw() +
#   theme(panel.grid=element_blank()) +
#   geom_line(linewidth=0.5, alpha=0.5)
# plot(beta_draws_plot)
# 
# # per sweep
# tauhat <- longbet.pred$tauhats
# tau_draws = as.data.frame(matrix(tauhat[,100], n,t2))
# tau_draws$trial = rownames(tau_draws)
# mdat = melt(tau_draws, id.vars="trial")
# mdat$time = as.numeric(gsub("V", "", mdat$variable))
# 
# tau_draws_plot = ggplot(mdat, aes(x=time, y=value, group=trial)) +
#   theme_bw() +
#   theme(panel.grid=element_blank()) +
#   geom_line(linewidth=0.5, alpha=0.5) + 
#   labs(y = "Tauhat", title = "Tauhats at 100 sweep")
# plot(tau_draws_plot)
# 
# num_sweeps <- longbet.fit$model_params$num_sweeps 
# tau_draws <- matrix(NA, num_sweeps, t2)
# for (i in 1:num_sweeps){
#   temp <- matrix(tauhat[,i], n, t2)
#   tau_draws[i, ] <- temp[1,]
# }
# tau_draws = as.data.frame(tau_draws)
# tau_draws$trial = rownames(tau_draws)
# mdat = melt(tau_draws, id.vars="trial")
# mdat$time = as.numeric(gsub("V", "", mdat$variable))
# 
# tau_draws_plot = ggplot(mdat, aes(x=time, y=value, group=trial)) +
#   theme_bw() +
#   theme(panel.grid=element_blank()) +
#   geom_line(linewidth=0.5, alpha=0.5) + 
#   labs(y = "Tauhat", title = "Tauhats of 1 obs across sweeps")
# plot(tau_draws_plot)


