# simple demonstration of longbet with default parameters
require(dplyr)
require(ggplot2)
require(tidyr)
require(longbet)
require(panelView)
require(stringr)
require(MetBrewer)

set.seed(100)
n <- 2000      # number of observations
t0 <- 6        # earliest treatment adoption time
t1 <- 12       # total time period
alpha <- 0.05
pr_types <- c("linear", "non-linear")
trt_types <- c("homogeneous", "heterogeneous")
pcat <- 2    # number of categorical variable

source('dgp.R')

att.results <- data.frame(
  pr = character(),
  trt = character(),
  trees = double(),
  RMSE = double(),
  Bias = double(),
  Coverage = double(),
  Cover0 = double(),
  I.L = double()
)

catt.results <- data.frame(
  pr = character(),
  trt = character(),
  trees = double(),
  RMSE = double(),
  Bias = double(),
  Coverage = double(),
  Cover0 = double(),
  I.L = double()
)

for (pr in pr_types){
  for (trt in trt_types){
    # Data generating
    data <- dgp(n, t0, t1, pr_type = pr, trt_type = trt)
    xtrain <- data$x
    ytrain <- data$y
    ztrain <- data$z
    
    # align treatment effect
    align_tau <- matrix(NA, nrow = n, ncol = t1 - t0 + 1)
    for (i in 1:n){
      if (sum(ztrain[i,]) == 0) {next}
      align_tau[i, 1:sum(ztrain[i,])] = data$tau[i, ztrain[i,] == 1]
    } 
    att <- colMeans(align_tau, na.rm = T)
    
    for (num_trees in c(20, 30, 40, 60)){
      longbet.fit <- longbet(y = ytrain, x = xtrain, x_trt = xtrain, z = ztrain, t = 1:t1,
                             num_sweeps = 500,
                             num_trees_pr =  num_trees, num_trees_trt = num_trees,
                             pcat = ncol(xtrain) - 3)
      
      longbet.pred <- predict.longbet(longbet.fit, xtrain, xtrain, ztrain)
      # align catt
      num_sweeps <- dim(longbet.pred$tauhats)[3]
      longbet.catt.sweeps <- array(NA, dim = c(n, t1 - t0 + 1, num_sweeps))
      for (i in 1:n){
        if (sum(ztrain[i,]) == 0) {next}
        longbet.catt.sweeps[i, 1:sum(ztrain[i,]), ] = longbet.pred$tauhats[i, ztrain[i,] == 1, ]
      } 
      
    
      longbet.att <- longbet.catt.sweeps %>%
        apply(c(2, 3), mean, na.rm = T) %>% t() %>%
        data.frame() %>%
        gather("t", "CATT") %>%
        mutate(t =  as.double(str_replace_all(t, c("X" = ""))) - 1) %>%
        group_by(t) %>%
        summarise(
          estimate = mean(CATT),
          conf.low = quantile(CATT, prob = alpha / 2),
          conf.high = quantile(CATT, prob = 1 - alpha / 2),
          method = "LongBet"
        )
      
      longbet.catt <- apply(longbet.catt.sweeps, c(1, 2), mean, na.rm = T)
      longbet.catt.low <- apply(longbet.catt.sweeps, c(1, 2), quantile, prob = alpha /2 , na.rm = T)
      longbet.catt.high <- apply(longbet.catt.sweeps, c(1, 2), quantile, prob = 1 - alpha /2 , na.rm = T)
      
      att.results[nrow(att.results) + 1,] <- c(pr, trt, num_trees, att.metric(att, longbet.att))
      catt.results[nrow(catt.results) + 1, ] <- c(pr, trt, num_trees, catt.metric(align_tau, longbet.catt, longbet.catt.low, longbet.catt.high))
    
    }
  }
}

write.csv(att.results, 'trees_sim.csv', row.names = F)

att.results <- read.csv('trees_sim.csv')
att.results %>%
  select(-c(Bias, Cover0, I.L)) %>%
  gather("metric", "value", -pr, -trt, -trees) %>% 
  mutate_at(c('trees', 'value'), as.numeric) %>%
  mutate(across(where(is.numeric), ~ round(., 3))) %>%
  ggplot(aes(x = trees, y = value, group = interaction(pr, trt), color = interaction(pr, trt))) +
  geom_line() + 
  facet_wrap(~metric, scales = "free_y", nrow = 1)



