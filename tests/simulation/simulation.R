
# Description -------------------------------------------------------------
# This script conduct simulation studies to evaluate the performance of Longbet
# and its baseline methods on panel data with staggered adoption

# Required libraries ------------------------------------------------------
require(longbet)
require(dplyr)
require(tidyr)
require(stringr)
require(forecast) # for data generating process
require(did)

source('dgp.R') # script for data generating process

# Set up ------------------------------------------------------------------
set.seed(100)
mc <- 10 #1000 # monte carlo iterations
n <- 2000      # number of observations
t0 <- 6        # earliest treatment adoption time
t1 <- 12       # total time period
alpha <- 0.05
pr_types <- c("linear", "non-linear")
trt_types <- c("homogeneous", "heterogeneous")
pcat <- 2    # number of categorical variable

if (file.exists('att.csv') & file.exists('catt.csv')){
  att.results <- read.csv('att.csv')
  catt.results <- read.csv('catt.csv')
} else {
  att.results <- data.frame(
    iter = double(),
    pr = character(),
    trt = character(),
    method = character(),
    RMSE = double(),
    Bias = double(),
    Coverage = double(),
    I.L = double(),
    Time = double()
  )
  
  catt.results <- data.frame(
    iter = double(),
    pr = character(),
    trt = character(),
    method = character(),
    RMSE = double(),
    Bias = double(),
    Coverage = double(),
    I.L = double()
  )
}

# simulation ---------------------------------------------------------------
for (pr in pr_types){
  for (trt in trt_types){
    for (iter in 1:mc){
      if (any((att.results$pr == pr) & (att.results$trt == trt))){
        if ( iter <= max(att.results$iter[(att.results$pr == pr) & (att.results$trt == trt)]) ){
          # skip this iteration
          next 
        }
      }
      
      # Data generating
      data <- dgp(n, t0, t1, pr_type = pr, trt_type = trt)
      xtrain <- data$x
      ytrain <- data$y
      ztrain <- data$z
      
      # panel view
      first.treat <- t1 - rowSums(ztrain) + 1
      first.treat[first.treat == t1 + 1] <- 0
      panel.data <- data.frame(
        ytrain = as.vector( t(ytrain) ),
        ztrain = as.vector( t(ztrain) ),
        first.treat = as.vector( sapply(first.treat, rep, t1)),
        id = as.vector( sapply(1:n, rep, t1)),
        time = rep(1:t1, n),
        X1 = as.vector( sapply(xtrain[,1], rep, t1)),
        X2 = as.vector( sapply(xtrain[,2], rep, t1)),
        X3 = as.vector( sapply(xtrain[,3], rep, t1)),
        X4 = as.vector( sapply(xtrain[,4], rep, t1)),
        X5 = as.vector( sapply(xtrain[,5], rep, t1))
      )
      
      # align treatment effect
      align_tau <- matrix(NA, nrow = n, ncol = t1 - t0 + 1)
      for (i in 1:n){
        if (sum(ztrain[i,]) == 0) {next}
        align_tau[i, 1:sum(ztrain[i,])] = data$tau[i, ztrain[i,] == 1]
      } 
      att <- colMeans(align_tau, na.rm = T)
      
      
      # Longbet -----------------------------------------------------------------
      longbet.time <- proc.time()
      longbet.fit <- longbet(y = ytrain, x = xtrain, z = ztrain, t = 1:t1,
                             num_sweeps = 100, num_trees_pr =  20, num_trees_trt = 20,
                             pcat = pcat)

      longbet.pred <- predict.longbet(longbet.fit, xtrain, ztrain)
  
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
      longbet.time <- proc.time() - longbet.time
      
      att.results[nrow(att.results) + 1,] <- c(iter, pr, trt, 'LongBet100', att.metric(att, longbet.att), as.numeric(longbet.time[3]))
      catt.results[nrow(catt.results) + 1, ] <- c(iter, pr, trt, 'LongBet100', catt.metric(align_tau, longbet.catt, longbet.catt.low, longbet.catt.high))
      
      
      # Longbet 500 sweeps-----------------------------------------------------------------
      longbet.time <- proc.time()
      longbet.fit <- longbet(y = ytrain, x = xtrain, z = ztrain, t = 1:t1,
                             num_sweeps = 500, num_trees_pr =  20, num_trees_trt = 20,
                             pcat = pcat)
      
      longbet.pred <- predict.longbet(longbet.fit, xtrain, ztrain)
      
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
      longbet.time <- proc.time() - longbet.time
      
      att.results[nrow(att.results) + 1,] <- c(iter, pr, trt, 'LongBet500', att.metric(att, longbet.att), as.numeric(longbet.time[3]))
      catt.results[nrow(catt.results) + 1, ] <- c(iter, pr, trt, 'LongBet500', catt.metric(align_tau, longbet.catt, longbet.catt.low, longbet.catt.high))
      
      
      # Baseline: DiD with multiple periods --------------------------------------------------------------
      did.time <- proc.time()
      did.out <- att_gt(yname = "ytrain",
                        gname = "first.treat",
                        idname = "id",
                        tname = "time",
                        xformla = ~ X1 + X3 + X2 + X4 + X5,
                        data = panel.data,
                        est_method = "dr",
                        control_group = "notyettreated"
      )
      DiD <- aggte(did.out, type = "dynamic", na.rm = TRUE) %>% 
        tidy() %>% 
        rename(t = event.time) %>% 
        filter(t >= 0 & t < 8) %>% 
        select(t, estimate, conf.low, conf.high) %>% 
        mutate(method = "DiD")
      did.time <- proc.time() - did.time
      att.results[nrow(att.results) + 1,] <- c(iter, pr, trt, 'DiD', att.metric(att, DiD), as.numeric(did.time[3]))
      
      # Baseline: DiD Non-linear ------------------------------------------------
      did_nl.time <- proc.time()
      if (pr == "linear"){
        did_nl.out <- att_gt(yname = "ytrain",
                             gname = "first.treat",
                             idname = "id",
                             tname = "time",
                             xformla = ~ X1 + X3 + X2 + X4 + X5,
                             data = panel.data,
                             est_method = "dr",
                             control_group = "notyettreated"
        )
      } else {
        did_nl.out <- att_gt(yname = "ytrain",
                             gname = "first.treat",
                             idname = "id",
                             tname = "time",
                             xformla = ~ X1 * X3 + X2 + X4 + X5,
                             data = panel.data,
                             est_method = "dr",
                             control_group = "notyettreated")
      }
      
      DiD_nl <- aggte(did_nl.out, type = "dynamic", na.rm = TRUE) %>% 
        tidy() %>% 
        rename(t = event.time) %>% 
        filter(t >= 0 & t < 8) %>% 
        select(t, estimate, conf.low, conf.high) %>% 
        mutate(method = "Non-linear DiD")
      did_nl.time = proc.time() - did_nl.time
      att.results[nrow(att.results) + 1,] <- c(iter, pr, trt, 'Non-linear DiD', att.metric(att, DiD_nl), as.numeric(did_nl.time[3]))
      
      write.csv(att.results, file = 'att.csv', row.names= F)
      write.csv(catt.results, file = 'catt.csv', row.names= F)
      
    }
  }
}

att.results <- read.csv('att.csv')
print(att.results)



