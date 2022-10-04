# simple demonstration of longbet with default parameters
library(longBet)
library(dbarts)
library(dplyr)
library(ggplot2)
library(tidyr)

#### 1. DATA GENERATION PROCESS
# n <- 200 # number of observations
# t1 <- 20
# t0 <- 11

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

t_longbet <- proc.time()
longbet.fit <- longBet(y = y, x = x, z = expand_z_mat, t = 1:t1,
num_trees_pr =  50, num_trees_trt = 50 ,
pcat = ncol(x) - 3,  sig_knl = 1, lambda_knl = 1)
# TODO: lambda_knl is quite sensitve, need better understanding

mu_hat_longbet <- apply(longbet.fit$muhats.adjusted, c(1, 2), mean)
tau_hat_longbet <- apply(longbet.fit$tauhats.adjusted, c(1, 2), mean)
tau_longbet <- tau_hat_longbet[,t0:t1]
t_longbet <- proc.time() - t_longbet



cat("beta draws: ", rowMeans(longbet.fit$beta_draws), "\n")
print(paste0("longbet RMSE: ", sqrt(mean((as.vector(tau_longbet) - as.vector(tau_mat))^2))))
print(paste0("longbet runtime: ", round(as.list(t_longbet)$elapsed,2)," seconds"))

att <- colMeans(tau_mat[z == 1, ])
rep_z <- matrix(rep(z, t1), n, t1)
ytilde <- y - matrix(rowMeans(y)) %*% rep(1,t1)- rep(1,n) %*% t(matrix(colMeans(y))) + mean(y)
recenter <- mean((colMeans(ytilde[matrix(rep_z[,1])==1,]) - colMeans(ytilde[matrix(rep_z[,1])==0,]))[1:t0-1])
tau_hat_freq <- colMeans(ytilde[matrix(rep_z[,1])==1,]) - colMeans(ytilde[matrix(rep_z[,1])==0,])-recenter
print(paste0("Frequentist att rmse: ", sqrt( mean((tau_hat_freq[t0:t1] - att)^2) ) ))
att_longbet <- colMeans(tau_longbet[z == 1, ])
print(paste0("longbet att rmse: ", sqrt(mean( (att_longbet - att)^2 ))))


## create covariates for bart ##
xtr <- cbind(z_vec, x_bart)
xte <- cbind(1 - z_vec, x_bart)
ytr <- y_vec

t_bart = proc.time()
ce_bart <- list()
bartps<-bart(x.train = xtr, y.train = ytr, x.test = xte)
ppd_test<-t(apply(bartps$yhat.test,1,function(x) rnorm(n=length(x),mean=x,sd=bartps$sigma)))
ppd_test_mean<-apply(ppd_test,2,mean)

## individual causal effects ##
ce_bart$ite<-rep(NA,length(ytr))
ce_bart$ite[which(xtr[,1]==1)] <- ytr[which(xtr[,1]==1)] - ppd_test_mean[which(xtr[,1]==1)]
ce_bart$ite[which(xtr[,1]==0)] <- ppd_test_mean[which(xtr[,1]==0)] - ytr[which(xtr[,1]==0)]
tau_bart <- matrix(ce_bart$ite, n, t1)[,t0:t1]

ce_bart$itu<-apply(ppd_test,2,quantile,probs=0.975)
ce_bart$itl<-apply(ppd_test,2,quantile,probs=0.025)

## average causal effects ##
ce_bart$ate <- mean(ce_bart$ite)

t_bart = proc.time() - t_bart


# compare results to inference
print(paste0("longbet RMSE: ", round( sqrt(mean((as.vector(tau_longbet) - as.vector(tau_mat))^2)), 2 ) ))
print(paste0("longbet runtime: ", round(as.list(t_longbet)$elapsed,2)," seconds"))


print(paste0("bart RMSE: ", round( sqrt(mean((tau_bart - tau_mat)^2)), 2)))
print(paste0("bart runtime: ", round(as.list(t_bart)$elapsed,2)," seconds"))

# check att
att <- tau_mat[z == 1,] %>% colMeans
att_bart <- tau_bart[z == 1,] %>% colMeans
att_longbet <- tau_longbet[z == 1,] %>% colMeans
print(paste0("longbet att rmse: ", round( sqrt(mean((att_longbet - att)^2)), 2)))
print(paste0("BART att rmse: ", round( sqrt(mean((att_bart - att)^2)), 2)))


# visualize ---------------------------------------------------------------

att_df <- data.frame(
  time = t0:t1,
  true = att,
  bart = att_bart,
  longbet = att_longbet
)

att_df %>% 
  gather("method", "att", -time) %>%
  ggplot(aes(time, att)) + 
  geom_line(aes(color = method)) + 
  ylab(labs(title = "Average Treatment on Treated"))
