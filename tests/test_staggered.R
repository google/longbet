# simple demonstration of longbet with default parameters
library(dbarts)
library(dplyr)
library(ggplot2)
library(tidyr)
library(longBet)
# DATA GENERATION PROCESS -------------------------------------------------


set.seed(1)
n <- 2000
t0 <- 3 # treatment start time
t1 <- 12 # observed response period
alpha = 0.05

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
tau <- 1 + 2 * abs(x[,2] * x[,5])
# time-varyiing heterogeneous treatment effect
tau_mat <- 2 + 2 * outer(tau, beta_t , "*")
# expand tau_mat
tau_mat_full <- cbind(matrix(0, nrow = n, ncol = t0 - 1), tau_mat)

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
s <- sd(mu(x))
pi <- 0.8 * pnorm(1 * mu(x) / s - 0.5 * x[,1]) + 0.05 + runif(n) / 10
# z <- rbinom(n,1,pi)
# z_mat <- matrix(rep(z, t1 - t0 + 1), n, t1 - t0 + 1)
# higher propensity score means getting treated earlier
z_mat <- matrix(0, n, t0 - 1)
for (i in t0:t1){
  treated <- (z_mat[, i - 1] == 1)
  treatment <- rbinom(n, 1, pi^2) # adjust the probability
  z <- apply(cbind(treatment, treated), 1, max) # draw new treatment for time i
  z_mat <- cbind(z_mat, z)
}
colnames(z_mat) <- NULL
plot(1:t1, colMeans(z_mat), main = "Treated ratio over time", ylab = "Treated ratio")

# get post treatment time matrix
trt_time <- t(apply(z_mat, 1, cumsum))

get_te <- function(tau, trt){
  te <- rep(0, length(trt))
  te[trt > 0] = tau[trt[trt > 0]]
  return(te)
}
te <- t(mapply(get_te, data.frame(t(tau_mat)), data.frame(t(trt_time))))

# generate outcome variable
y0 <- mu_mat
y1 <- mu_mat + te
y <- z * y1 + (1 -z) * y0 +  0.2 * sd(mu_mat) * matrix(rnorm(n*t1), n, t1)

# # If you didn't know pi, you would estimate it here
# pi_mat <- as.matrix(rep(pi, t1), n, t1)
# pi_vec <- as.vector(pi_mat)

# get training data
ytrain <- y[, 1:t1]
ztrain <- z_mat
xtrain <- x

# longbet -----------------------------------------------------------------
t_longbet <- proc.time()
longbet.fit <- longbet(y = ytrain, x = xtrain, z = ztrain, t = 1:t1,
                       num_sweeps = 60,
                       num_trees_pr =  20, num_trees_trt = 20,
                       pcat = ncol(x) - 3)

longbet.pred <- predict.longBet(longbet.fit, x, 1:t1)


# mu_hat_longbet <- apply(longbet.pred$muhats, c(1, 2), mean)
# tau_hat_longbet <- apply(longbet.pred$tauhats, c(1, 2), mean)
# tau_longbet <- tau_hat_longbet[,t0:t1]
# t_longbet <- proc.time() - t_longbet
# 
# ate_longbet_fit <- apply(longbet.pred$tauhats, c(2, 3), mean)[t0:t1, ]
# ate <- tau_mat %>% colMeans
# ate_longbet <- ate_longbet_fit %>% rowMeans

longbet.ate <- get_ate(longbet.pred, alpha = 0.05)
longbet.att <- get_att(longbet.pred, z = ztrain, alpha = 0.05)
longbet.cate <- get_cate(longbet.pred, alpha = 0.05)

print(paste0("longbet CATE RMSE in-sample: ", sqrt(mean((as.vector(longbet.cate$cate[,t0:t1]) - as.vector(tau_mat[, t0:t1]))^2))))
print(paste0("longbet CATE RMSE extrapolate: ", sqrt(mean((as.vector(longbet.cate$cate[,(t1 + 1):t1]) - as.vector(tau_mat[, (t1 + 1):t1]))^2))))

ate <- tau_mat %>% colMeans
print(paste0("longbet ATE RMSE in-sample: ", round( sqrt(mean((longbet.ate$ate[t0:t1] - ate[t0:t1])^2)), 2)))
print(paste0("longbet ATE RMSE extrapolate: ", round( sqrt(mean((longbet.ate$ate[(t1 + 1):t1] - ate[(t1 + 1):t1])^2)), 2)))

print(paste0("longbet runtime: ", round(as.list(t_longbet)$elapsed,2)," seconds"))

# visualize ---------------------------------------------------------------
colors <- c("black", "#FFA500", "#00BFFF")
labels <- c("True", "LongBet", "BART")
names(colors) <- labels
# ATE
ate_df <- data.frame(
  time = t0:t1,
  true = ate[t0:t1],
  longbet = longbet.ate$ate[t0:t1],
  lower = longbet.ate$interval[1, t0:t1],
  upper = longbet.ate$interval[2, t0:t1]
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
  true = as.vector(t(tau_mat[, t0:t1])),
  lonbet = as.vector(t(longbet.cate$cate[, t0:t1])),
  time = rep(c(t0:t1), nrow(tau_mat)),
  id = as.vector(sapply(1:n, rep, (t1 - t0 + 1)))
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
# tau_draws = as.data.frame(matrix(tauhat[,100], n,t1))
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
# tau_draws <- matrix(NA, num_sweeps, t1)
# for (i in 1:num_sweeps){
#   temp <- matrix(tauhat[,i], n, t1)
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


