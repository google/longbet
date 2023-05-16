# simple demonstration of longbet with default parameters
library(dbarts)
library(dplyr)
library(ggplot2)
library(tidyr)
library(longBet)
# DATA GENERATION PROCESS -------------------------------------------------


set.seed(1)
n <- 2000
t1 <- 12
t0 <- 6
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
  # x_bart <- rbind(x_bart, cbind(i, z * max(0, i - t0 + 1), x))
  x_bart <- rbind(x_bart, cbind(max(0, i - t0 + 1) , x))
}
z_vec <- c(rep(0, n * (t0 - 1)), as.vector(z_mat))

expand_z_mat <- cbind(matrix(0, n, (t0 - 1)), z_mat)


# longbet -----------------------------------------------------------------
t_longbet <- proc.time()
longbet.fit <- longbet(y = y, x = x, z = expand_z_mat, t = 1:t1,
                       num_sweeps = 100,
                       num_trees_pr =  20, num_trees_trt = 20 ,
                       pcat = ncol(x) - 3,  sig_knl = 1, lambda_knl = 1)
# TODO: lambda_knl is quite sensitve, need better understanding

longbet.pred <- predict.longBet(longbet.fit, x, 1:t1)

longbet.ate <- get_ate(longbet.pred, alpha = 0.05)
longbet.att <- get_att(longbet.pred, z = expand_z_mat, alpha = 0.05)
longbet.cate <- get_cate(longbet.pred, alpha = 0.05)
t_longbet <- proc.time() - t_longbet


tau_longbet <- longbet.cate$cate[,t0:t1]
print(paste0("longbet CATE RMSE: ", sqrt(mean((as.vector(tau_longbet) - as.vector(tau_mat))^2))))
print(paste0("longbet CATT RMSE: ", sqrt(mean((as.vector(tau_longbet[z == 1, ]) - as.vector(tau_mat[z==1,]))^2))))
print(paste0("longbet runtime: ", round(as.list(t_longbet)$elapsed,2)," seconds"))

# # bart --------------------------------------------------------------------
xtr <- cbind(z_vec, x_bart)
xte <- cbind(1 - z_vec, x_bart)
ytr <- y_vec

# t_bart = proc.time()
# ce_bart <- list()
# bartps<-bart(x.train = xtr, y.train = ytr, x.test = xte, ndpost = 200)
# ppd_test<-t(apply(bartps$yhat.test,1,function(x) rnorm(n=length(x),mean=x,sd=bartps$sigma)))
# # ppd_test <- bartps$yhat.test
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
# t_bart = proc.time() - t_bart

t_bart = proc.time()

xtrt <- cbind(rep(1, nrow(x_bart)), x_bart)
xctrl <- cbind(rep(0, nrow(x_bart)), x_bart)
bartps <- bart(x.train = xtr, y.train = ytr, x.test = rbind(xtrt, xctrl), ndpost = 200)
y1_test_mean <- apply(bartps$yhat.test[, 1:nrow(xtrt)], 2, mean)
y0_test_mean <- apply(bartps$yhat.test[, (nrow(xtrt) + 1):(nrow(xtrt) + nrow(xctrl))], 2, mean)
tau_bart <- matrix(y1_test_mean - y0_test_mean, n, t1)[, t0:t1]

t_bart = proc.time() - t_bart
# results -----------------------------------------------------------------
# check ate

ate <- tau_mat %>% colMeans
ate_bart <- tau_bart %>% colMeans
ate_longbet <- tau_longbet %>% colMeans

print(paste0("longbet CATE RMSE: ", round( sqrt( mean(  as.vector(tau_longbet - tau_mat)^2 ) ), 2 ) ))
print(paste0("longbet ate rmse: ", round( sqrt(mean((ate_longbet - ate)^2)), 2)))
print(paste0("longbet runtime: ", round(as.list(t_longbet)$elapsed,2)," seconds"))

print(paste0("bart CATE RMSE: ", round( sqrt( mean( as.vector(tau_bart - tau_mat)^2 ) ), 2)))
print(paste0("BART ate rmse: ", round( sqrt(mean((ate_bart - ate)^2)), 2)))
print(paste0("bart runtime: ", round(as.list(t_bart)$elapsed,2)," seconds"))


colMeans(tau_longbet) - ate_longbet
colMeans(tau_bart) - ate_bart
dim(tau_bart)
# visualize ---------------------------------------------------------------
colors <- c("black", "#FFA500", "#00BFFF")
labels <- c("True", "LongBet", "BART")
names(colors) <- labels
# ATE
ate_longbet_fit <- apply(longbet.pred$tauhats, c(2, 3), mean)[t0:t1, ]
ate_df <- data.frame(
  time = t0:t1,
  true = ate,
  bart = ate_bart,
  longbet = ate_longbet,
  longbet_up = apply(ate_longbet_fit, 1, quantile, probs = 1 - alpha / 2),
  longbet_low = apply(ate_longbet_fit, 1, quantile, probs = alpha / 2)
)

ate_plot <- 
  ggplot(ate_df , aes(x = time, y = true)) +
  geom_line(aes(y = true, color = "True")) +
  geom_line(aes(y = longbet, color = "LongBet")) +
  geom_line(aes(y = bart, color = "BART")) +
  geom_ribbon(aes(ymin = longbet_low, ymax = longbet_up, fill = "LongBet"), alpha = 0.15, fill = colors[2]) +
  # geom_ribbon(aes(ymin = longbet_low, ymax = longbet_up, fill = "BART"), alpha = 0.15, fill = colors[3]) +
  labs(x = "Time", y = "ATE", color = "Legend") +
  scale_color_manual(name = "Legend", values = c("black", "#FFA500", "#00BFFF"), labels = c("True", "LongBet", "BART"))
print(ate_plot)

readline(prompt="Press [enter] to continue")

# CATE
cate_df <- data.frame(
  true = as.vector(t(tau_mat)),
  lonbet = as.vector(t(tau_longbet)),
  bart = as.vector(t(tau_bart)),
  time = rep(c(t0:t1), nrow(tau_mat)),
  id = as.vector(sapply(1:nrow(tau_longbet), rep, (t1 - t0 + 1)))
)

cate_plot <- cate_df %>%
  gather("method", "cate", -time, -id) %>%
  ggplot() +
  geom_line(aes(time, cate, group = id, color = id)) +
  facet_wrap(~method, ncol = 2)
plot(cate_plot)

readline(prompt="Press [enter] to continue")
