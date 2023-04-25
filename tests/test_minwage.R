# compare longbet to did method on the minimum wage dataset.
# estimate the effect of minimum wage on teen employment
# did paper: https://reader.elsevier.com/reader/sd/pii/S0304407620303948?token=967C9A6B23A76EC6C271F99CC54C0CAC8829014C025C3AB9F600C1B5F39F0A7EB6CD905A583911978C1025A798AD05D6&originRegion=us-east-1&originCreation=20230425183217
# required package: did

library(longBet)
library(did)
library(tidyr)
library(dplyr)
data(mpdta)

# to start with, only compare the control group and treated group in 2004
data <- mpdta[mpdta$first.treat < 2006, ]
data$treat[data$first.treat > data$year] = 0
ttrain <- unique(data$year)

xtrain <- data[, c("countyreal", "lpop")] %>% 
  group_by(countyreal) %>%
  summarise(lpop = mean(lpop))

ytrain <- data[, c("year", "lemp", "countyreal")] %>%
  spread(key = "year", value = "lemp")

ztrain <- data[, c("year", "treat", "countyreal")] %>%
  spread(key = "year", value = "treat")

# check x, y, z id are correct
all(xtrain$countyreal == ytrain$countyreal)
all(xtrain$countyreal == ztrain$countyreal)
countyreal <- xtrain$countyreal
xtrain$countyreal <- NULL
ytrain$countyreal <- NULL
ztrain$countyreal <- NULL

xtrain <- as.matrix(xtrain)
ytrain <- as.matrix(ytrain)
ztrain <- as.matrix(ztrain)

t_longbet <- proc.time()
longbet.fit <- longbet(y = ytrain, x = xtrain, z = ztrain, t = 1:5,
                       num_sweeps = 100,
                       num_trees_pr =  20, num_trees_trt = 20,
                       pcat = 0,  b_scaling = TRUE)
# TODO: lambda_knl is quite sensitve, need better understanding
# sigma_knl = mean( sqrt( apply(longbet.fit$beta_draws[t0:t1,], 2, var) ))
# lambda_knl = 4 #(t1 - t0 + 1) / 2
# 
# longbet.pred <- predict.longBet(longbet.fit, xtrain, ttrain, sigma = sigma_knl, lambda = lambda_knl)
# mu_hat_longbet <- apply(longbet.pred$muhats.adjusted, c(1, 2), mean)
# tau_hat_longbet <- apply(longbet.pred$tauhats.adjusted, c(1, 2), mean)
# tau_longbet <- tau_hat_longbet[,t0:t2]
# t_longbet <- proc.time() - t_longbet
