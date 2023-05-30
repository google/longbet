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
                       pcat = ncol(xtrain) - 3)

longbet.pred <- predict.longBet(longbet.fit, x, ztrain)
longbet.ate <- get_ate(longbet.pred, alpha = 0.05)
longbet.cate <- get_cate(longbet.pred, alpha = 0.05)

#TODO estimate CATT
# post_trt_time <- t(apply(ztrain, 1, get_post_trt_time, t = 1:t1))
treatment_period <- t1 - t0 + 1
# Align treatment effect 
align_te <- matrix(NA, nrow = n, ncol = treatment_period)
align_catt <- matrix(NA, nrow = n, ncol = treatment_period)

for (i in 1:n){
  if (sum(ztrain[i,]) == 0) {next}
  post_t <- 1:sum(ztrain[i,])
  align_te[i, post_t] = te[i, ztrain[i,] == 1]
  align_catt[i, post_t] = longbet.cate$cate[i, ztrain[i,] == 1]
} 

att <- colMeans(align_te, na.rm = T)
att_hat <- colMeans(align_catt, na.rm = T)

print(paste0("longBet ATT RMSE: ", sqrt(mean(( att - att_hat)^2))))
print(paste0("longBet CATT RMSE: ", sqrt(mean((align_te - align_catt)^2, na.rm = T))))
print("longbet CATT RMSE by time: ")
print(sqrt(colMeans( (align_te - align_catt)^2, na.rm = T)))
print(paste0("longBet runtime: ", round(as.list(t_longbet)$elapsed,2)," seconds"))


# Visualize ---------------------------------------------------------------
# ATE
att_df <- data.frame(
  time = t0:t1,
  true = att,
  longbet = att_hat
)

plot_att <- att_df %>%
  gather("method", "att", -time) %>%
  ggplot(aes(time, att)) +
  geom_line(aes(color = method)) +
  ylab(labs(title = "Average Treatment Effect on Treated"))
print(plot_att)

# CATE
catt_df <- data.frame(
  true = as.vector(t(align_te)),
  lonbet = as.vector(t(align_catt)),
  time = rep(c(t0:t1), nrow(tau_mat)),
  id = as.vector(sapply(1:nrow(align_te), rep, (t1 - t0 + 1)))
)

plot_catt <- catt_df %>%
  filter(id < 100) %>%
  gather("method", "cate", -time, -id) %>%
  ggplot() +
  geom_line(aes(time, cate, group = id, color = id)) +
  facet_wrap(~method)
print(plot_catt)


