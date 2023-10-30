# simple demonstration of longbet with default parameters
library(dbarts)
library(dplyr)
library(ggplot2)
library(tidyr)
library(longbet)
library(forecast)
# DATA GENERATION PROCESS -------------------------------------------------


set.seed(1)
n <- 2000
t0 <- 6 # treatment start time
t1 <- 12 # observed response period
alpha = 0.05

# generate covariates
x1 <- rnorm(n)
x2 <- rnorm(n)
x3 <- rnorm(n)
x4 <- rbinom(n,1,0.5)
x5 <- factor(sample(1:3,n,replace=TRUE,prob = c(1/3, 1/3,1/3)))
x <- cbind(x1, x2, x3, x4, x5)

# define average time-varying treatment effect
post_t <- 1:(t1 - t0 + 1) 
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
alpha_t <- arima.sim(model = list(order = arima_order, ar = 0.7, ma = -0.4), n = t1) + 1
nu_i <- mu(x)
mu_mat <- outer(nu_i, alpha_t, "*")

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
colnames(z_mat) <- NULL
treated <- apply(z_mat, 1, function(x) any(x != 0))
treatment_group <- apply(z_mat, 1, function(x) t1 - t0 - sum(x) + 2)
treatment_group[treatment_group == max(treatment_group)] <- 0 # Control group

plot(1:t1, colMeans(z_mat), main = "Treated ratio over time", ylab = "Treated ratio")

# get post treatment time matrix
trt_time <- t(apply(z_mat, 1, cumsum))


early_bird_effect <- c(1.2, 1.2, 1, 1, 0.9, 0.9, 0.8, 0.8)
get_treatment_effect <- function(ground_truth, trt_time){
  trt_effect <- rep(0, length(trt_time))
  eb_effect <- early_bird_effect[sum(trt_time == 0) - t0 + 2]
  trt_effect[trt_time > 0] = ground_truth[trt_time[trt_time > 0]] * eb_effect
  return(trt_effect)
}
te <- t(mapply(get_treatment_effect, data.frame(t(tau_mat)), data.frame(t(trt_time))))

# generate outcome variable
y0 <- mu_mat
y1 <- mu_mat + te
y <- z * y1 + (1 -z) * y0 +  0.2 * sd(mu_mat) * matrix(rnorm(n*t1), n, t1)


# Visualize observed data -------------------------------------------------
# TODO: Visualize observed data and treatment effect
treated_obs <- data.frame(
  time = 1:t1,
  avg = colMeans(y[treated,]),
  upper = apply(y[treated,], 2, quantile, prob = 0.975),
  lower = apply(y[treated,], 2, quantile, prob = 0.025),
  group = "Treated"
)

control_obs <- data.frame(
  time = 1:t1,
  avg = colMeans(y[!treated,]),
  upper = apply(y[!treated,], 2, quantile, prob = 0.975),
  lower = apply(y[!treated,], 2, quantile, prob = 0.025),
  group = "Control"
)

gg_obs <- rbind(treated_obs, control_obs) %>% ggplot(aes(time, avg, group = group, colour = group, fill = group)) +
  geom_line(aes(color = group)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .3, colour = NA) + 
  facet_wrap(~group, nrow = 1)
plot(gg_obs)

# Visualize treatment effect by group
trt_effect_df <- data.frame()
for (group in 1:max(treatment_group)){
  temp_y <- te[treatment_group == group, t0:t1]
  temp_df <- data.frame(
    time = t0:t1,
    avg = colMeans(temp_y),
    upper = apply(temp_y, 2, quantile, prob = 0.975),
    lower = apply(temp_y, 2, quantile, prob = 0.025),
    group = group
  )
  trt_effect_df <- rbind(trt_effect_df, temp_df)
}

gg_trt_effect <- trt_effect_df %>% ggplot(aes(time, avg, group = group, colour = group, fill = group)) +
  geom_pointrange(aes(ymin = lower, ymax = upper)) + 
  facet_wrap(~group, nrow = 4, dir = "v")
plot(gg_trt_effect)




# Training ----------------------------------------------------------------
# TODO: use estimated propensity score in training data
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

longbet.pred <- predict.longbet(longbet.fit, x, ztrain)
longbet.att <- get_att(longbet.pred, alpha = 0.05)
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

print(paste0("longbet ATT RMSE: ", sqrt(mean(( att - att_hat)^2))))
print(paste0("longbet CATT RMSE: ", sqrt(mean((align_te - align_catt)^2, na.rm = T))))
print("longbet CATT RMSE by time: ")
print(sqrt(colMeans( (align_te - align_catt)^2, na.rm = T)))
print(paste0("longbet runtime: ", round(as.list(t_longbet)$elapsed,2)," seconds"))


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

# muhat
y0hat <- longbet.pred$muhats %>% apply(MARGIN = c(1, 2), mean) %>% data.frame
y0_df <- data.frame(
  true = as.vector(t(mu_mat)),
  lonbet = as.vector(t(y0hat)),
  time = rep(c(1:t1), nrow(mu_mat)),
  id = as.vector(sapply(1:nrow(y0hat), rep, (t1)))
)

y0_plot <-  y0_df %>%
  gather("method", "y0", -time, -id) %>%
  ggplot() +
  geom_line(aes(time, y0, group = id, color = id)) +
  facet_wrap(~method)
plot(y0_plot)


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


