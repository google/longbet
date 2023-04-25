# simple demonstration of longbet with default parameters
library(dbarts)
library(dplyr)
library(ggplot2)
library(tidyr)
library(longBet)
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
y = mu_mat + 0.2 * sd(mu_mat) * matrix(rnorm(n*t2), n, t2)
y[, t0:t2] = y[, t0:t2] + tau_mat * z_mat

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
                       num_sweeps = 100,
                       num_trees_pr =  20, num_trees_trt = 20,
                       pcat = ncol(x) - 3,  sig_knl = 1, lambda_knl = 1,
                       b_scaling = TRUE)
# TODO: lambda_knl is quite sensitve, need better understanding
sigma_knl = mean( sqrt( apply(longbet.fit$beta_draws[t0:t1,], 2, var) ))
lambda_knl = 4 #(t1 - t0 + 1) / 2

longbet.pred <- predict.longBet(longbet.fit, x, 1:t2, sigma = sigma_knl, lambda = lambda_knl)
mu_hat_longbet <- apply(longbet.pred$muhats.adjusted, c(1, 2), mean)
tau_hat_longbet <- apply(longbet.pred$tauhats.adjusted, c(1, 2), mean)
tau_longbet <- tau_hat_longbet[,t0:t2]
t_longbet <- proc.time() - t_longbet

cat("          beta draws: ", round(rowMeans(longbet.fit$beta_draws),3), "\n")
cat("predicted beta draws: ", round(rowMeans(longbet.pred$beta_draws),3), "\n")

ate_longbet_fit <- apply(longbet.pred$tauhats.adjusted, c(2, 3), mean)[t0:t2, ]
ate <- tau_mat %>% colMeans
ate_longbet <- ate_longbet_fit %>% rowMeans

print(paste0("longbet CATE RMSE: ", sqrt(mean((as.vector(tau_longbet) - as.vector(tau_mat))^2))))
print(paste0("longbet ate rmse: ", round( sqrt(mean((ate_longbet - ate)^2)), 2)))
print(paste0("longbet runtime: ", round(as.list(t_longbet)$elapsed,2)," seconds"))

mu_hat_longbet_fit <- apply(longbet.fit$muhats.adjusted, c(1, 2), mean)
tau_hat_longbet_fit <- apply(longbet.fit$tauhats.adjusted, c(1, 2), mean)
tau_longbet_fit <- tau_hat_longbet_fit[,t0:t1]
print(paste0("longbet CATE in-sample: ", sqrt(mean((as.vector(tau_longbet_fit) - as.vector(tau_mat[, 1:(t1 - t0 + 1)]))^2))))

# # # bart --------------------------------------------------------------------
# xtr <- cbind(z_vec, x_bart)
# xte <- cbind(1 - z_vec, x_bart)
# ytr <- y_vec
# 
# t_bart = proc.time()
# ce_bart <- list()
# bartps<-bart(x.train = xtr, y.train = ytr, x.test = xte, ndpost = 200)
# ppd_test<-t(apply(bartps$yhat.test,1,function(x) rnorm(n=length(x),mean=x,sd=bartps$sigma)))
# ppd_test_mean<-apply(ppd_test,2,mean)
# 
# ## individual causal effects ##
# ce_bart$ite<-rep(NA,length(ytr))
# ce_bart$ite[which(xtr[,1]==1)] <- ytr[which(xtr[,1]==1)] - ppd_test_mean[which(xtr[,1]==1)]
# ce_bart$ite[which(xtr[,1]==0)] <- ppd_test_mean[which(xtr[,1]==0)] - ytr[which(xtr[,1]==0)]
# tau_bart <- matrix(ce_bart$ite, n, t1)[,t0:t1]
# 
# ce_bart$itu<-apply(ppd_test,2,quantile,probs=0.975)
# ce_bart$itl<-apply(ppd_test,2,quantile,probs=0.025)
# 
# ## average causal effects ##
# ce_bart$ate <- mean(ce_bart$ite)
# t_bart = proc.time() - t_bart



# # results -----------------------------------------------------------------
# # check ate
# ate_longbet_fit <- apply(longbet.fit$tauhats.adjusted, c(2, 3), mean)[t0:t1, ]
# 
# ate <- tau_mat %>% colMeans
# ate_bart <- tau_bart %>% colMeans
# ate_longbet <- ate_longbet_fit %>% rowMeans
# 
# print(paste0("longbet CATE RMSE: ", round( sqrt( mean(  as.vector(tau_longbet - tau_mat)^2 ) ), 2 ) ))
# print(paste0("longbet ate rmse: ", round( sqrt(mean((ate_longbet - ate)^2)), 2)))
# print(paste0("longbet runtime: ", round(as.list(t_longbet)$elapsed,2)," seconds"))
# 
# print(paste0("bart CATE RMSE: ", round( sqrt( mean( (tau_bart - tau_mat)^2 ) ), 2)))
# print(paste0("BART ate rmse: ", round( sqrt(mean((ate_bart - ate)^2)), 2)))
# print(paste0("bart runtime: ", round(as.list(t_bart)$elapsed,2)," seconds"))


# visualize ---------------------------------------------------------------
colors <- c("black", "#FFA500", "#00BFFF")
labels <- c("True", "LongBet", "BART")
names(colors) <- labels
# ATE
ate_df <- data.frame(
  time = t0:t2,
  true = ate,
  longbet = ate_longbet,
  longbet_up = apply(ate_longbet_fit, 1, quantile, probs = 1 - alpha / 2),
  longbet_low = apply(ate_longbet_fit, 1, quantile, probs = alpha / 2),
  longbet_beta = rowMeans(longbet.pred$beta_draws)[t0:t2]
)

ate_plot <- 
  ggplot(ate_df , aes(x = time, y = true)) +
  geom_line(aes(y = true, color = "True")) +
  geom_line(aes(y = longbet, color = "LongBet")) +
  geom_ribbon(aes(ymin = longbet_low, ymax = longbet_up, fill = "LongBet"), alpha = 0.15, fill = colors[2]) +
  labs(x = "Time", y = "ATE", color = "Legend") +
  scale_color_manual(name = "Legend", values = colors, labels = labels)
print(ate_plot)
# readline(prompt="Press [enter] to continue")

# CATE
cate_df <- data.frame(
  true = as.vector(t(tau_mat)),
  lonbet = as.vector(t(tau_longbet)),
  time = rep(c(t0:t2), nrow(tau_mat)),
  id = as.vector(sapply(1:nrow(tau_longbet), rep, (t2 - t0 + 1)))
)

cate_plot <- cate_df %>%
  gather("method", "cate", -time, -id) %>%
  ggplot() +
  geom_line(aes(time, cate, group = id, color = id)) +
  facet_wrap(~method)
print(cate_plot)

# readline(prompt="Press [enter] to continue")

# CATE error
cate_error <- data.frame(
  lonbet = as.vector(t(tau_longbet - tau_mat)),
  time = rep(c(t0:t2), nrow(tau_mat)),
  id = as.vector(sapply(1:nrow(tau_longbet), rep, (t2 - t0 + 1)))
)

cate_error_plot <- cate_error %>%
  gather("method", "cate", -time, -id) %>%
  ggplot() +
  geom_line(aes(time, cate, group = id, color = id)) +
  labs(x = "Time", y = "CATE Error") +
  facet_wrap(~method)
print(cate_error_plot)

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


