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
staggered = TRUE

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
gamma_i <- 1 + 2 * abs(x[,2] * x[,5])
# time-varyiing heterogeneous treatment effect
tau_mat <- 2 * outer(gamma_i, beta_t , "*")
# expand tau_mat
tau_mat_full <- cbind(matrix(0, nrow = n, ncol = t0 - 1), tau_mat)

# ## define prognostic function (RIC)
mu = function(x){
  result =  2 * abs(x[,5] + 0.5 * x[,3] - 1)
  return(result)
}
arima_order <- c(1, 0, 1) 
alpha_t <- arima.sim(model = list(order = arima_order, ar = 0.7, ma = -0.4), n = t2) + 1
nu_i <- mu(x)
mu_mat <- outer(nu_i, alpha_t, "*")

if (staggered){
  # compute propensity scores and treatment assignment
  s <- sd(mu(x))
  pi <- 0.2 * pnorm(0.5* mu(x) - 0.5 * x[,1])^2 + runif(n) / 10
  z_mat <- matrix(0, n, t0 - 1)
  for (i in t0:t1){
    treated <- (z_mat[, i - 1] == 1)
    treatment <- rbinom(n, 1, pi) # adjust the probability
    z <- apply(cbind(treatment, treated), 1, max) # draw new treatment for time i
    z_mat <- cbind(z_mat, z)
  }
  z_mat <- cbind(z_mat, matrix(rep(z_mat[,t1], t2 - t1), ncol = t2 -t1, byrow = FALSE))
  colnames(z_mat) <- NULL
} else {
  # compute propensity scores and treatment assignment
  s <- sd(mu(x))
  pi <- 0.8 * pnorm(1 * mu(x) / s - 0.5 * x[,1]) + 0.05 + runif(n) / 10
  z <- rbinom(n,1,pi)
  z_mat <- matrix(rep(z, t2 - t0 + 1), n, t2 - t0 + 1)
}

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
y <- z * y1 + (1 -z) * y0 +  0.2 * sd(mu_mat) * matrix(rnorm(n*t2), n, t2)

# get training data
ytrain = y[, 1:t1]
ztrain <- z_mat[,1:t1]

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
longbet.fit <- longbet(y = ytrain, x = x, x_trt = x, z = ztrain, t = 1:t1,
                       num_sweeps = 60,
                       num_trees_pr =  40, num_trees_trt = 40,
                       pcat = ncol(x) - 3, sig_knl =  1, lambda_knl = 1)

ztest <- z_mat
longbet.pred <- predict.longbet(longbet.fit, x, x, ztest, sigma = NULL, lambda = NULL)
t_longbet <- proc.time() - t_longbet

longbet.att <- get_att(longbet.pred, alpha = 0.05)
longbet.catt <- get_catt(longbet.pred, alpha = 0.05)

treated_in_sample <- ztest[, t0:t1] == 1
treated_extrapolate <- ztest[, (t1+1):t2] == 1
print(paste0("longbet CATE RMSE in-sample: ", sqrt(mean((as.vector(longbet.catt$catt[,t0:t1][treated_in_sample]) - as.vector(te[, t0:t1][treated_in_sample]))^2))))
print(paste0("longbet CATE RMSE extrapolate: ", sqrt(mean((as.vector(longbet.catt$catt[,(t1 + 1):t2][treated_extrapolate]) - as.vector(te[, (t1 + 1):t2][treated_extrapolate]))^2))))

# Align treatment effect 
treatment_period <- t2 - t0 + 1
align_te <- matrix(NA, nrow = n, ncol = treatment_period)
align_catt <- matrix(NA, nrow = n, ncol = treatment_period)
for (i in 1:n){
  if (sum(ztest[i,]) == 0) {next}
  post_t <- 1:sum(ztest[i,])
  align_te[i, post_t] = te[i, ztest[i,] == 1]
  align_catt[i, post_t] = longbet.catt$catt[i, ztest[i,]==1]
} 
att <- colMeans(align_te, na.rm = T)

print(paste0("longbet ATT RMSE in_sample: ", sqrt(mean(( att[1:(t1 - t0 + 1)] - longbet.att$att[1:(t1 - t0 + 1)])^2))))
print(paste0("longbet ATT RMSE extrapolate: ", sqrt(mean(( att[(t1 -t0 + 2): (t2 -t0 + 1)] -  longbet.att$att[(t1 -t0 + 2): (t2 -t0 + 1)])^2)) ))

print(paste0("longbet runtime: ", round(t_longbet[3],2)," seconds"))

# visualize ---------------------------------------------------------------
colors <- c("black", "#FFA500", "#00BFFF")
labels <- c("True", "LongBet", "BART")
names(colors) <- labels
# ATE
att_df <- data.frame(
  time = t0:t2,
  true = att,
  longbet = longbet.att$att,
  lower = longbet.att$interval[1,],
  upper = longbet.att$interval[2,]
)

att_plot <- 
  ggplot(att_df , aes(x = time, y = true)) +
  geom_line(aes(y = true, color = "True")) +
  geom_line(aes(y = longbet, color = "LongBet")) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = "LongBet"), alpha = 0.15, fill = colors[2]) +
  geom_vline(xintercept = t1, linetype = "dashed", color = "grey") +
  labs(x = "Time", y = "ATT", color = "Legend") +
  scale_color_manual(name = "Legend", values = colors, labels = labels)
print(att_plot)
readline(prompt="Press [enter] to continue")

# CATE
catt_df <- data.frame(
  true = as.vector(t(align_te)),
  lonbet = as.vector(t(align_catt)),
  time = rep(c(t0:t2), n),
  id = as.vector(sapply(1:n, rep, (t2 - t0 + 1)))
)


catt_plot <- catt_df %>%
  gather("method", "catt", -time, -id) %>%
  ggplot() +
  geom_line(aes(time, catt, group = id, color = id)) +
  geom_vline(xintercept = t1, linetype = "dashed", color = "grey") +
  facet_wrap(~method, ncol = 1)
print(catt_plot)

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


