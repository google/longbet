
# Description -------------------------------------------------------------
# This script conduct simulation studies to evaluate the performance of Longbet
# and its baseline methods on panel data with staggered adoption

# Required libraries ------------------------------------------------------
require(longBet)
require(dplyr)
require(forecast) # for data generating process

source('dgp.R') # script for data generating process

# Set up ------------------------------------------------------------------
set.seed(100)
mc <- 1 #1000 # monte carlo iterations
n <- 5000      # number of observations
t0 <- 6        # earliest treatment adoption time
t1 <- 12       # total time period
pr_types <- c("linear", "non-linear")
trt_types <- c("homogeneous", "heterogeneous")
pcat <- 2    # number of categorical variable


# simulation ---------------------------------------------------------------
for (pr in pr_types){
  for (trt in trt_types){
    # check if simulation results exist
    filename <- paste(pr, " ", trt, ".csv", sep = "")
    if (file.exists(filename)) {
      results <- read.csv(filename)
      max_iter <- max(results$iter)
    } else { 
      max_iter <- 0 
      results <- data.frame(
        iter = integer(),
        method = character(),
        ATT.RMSE = double(),
        ATT.Bias = double(),
        ATT.Coverage = double(),
        CATT.RMSE = double(),
        CATT.Bias = double(),
        CATT.Coverage = double(),
        Time = double()
      )
    }
    
    for (iter in 1:mc){
      if (iter <= max_iter){
        # skip this iteration
        break 
      }
      
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
      
      
      # Longbet -----------------------------------------------------------------
      t_longbet <- proc.time()
      longbet.fit <- longbet(y = ytrain, x = xtrain, z = ztrain, t = 1:t1,
                             num_sweeps = 60, num_trees_pr =  20, num_trees_trt = 20,
                             pcat = pcat)

      longbet.pred <- predict.longBet(longbet.fit, xtrain, ztrain)
      
      # align catt
      num_sweeps <- dim(longbet.pred$tauhats)[3]
      align_catt <- array(NA, dim = c(n, t1 - t0 + 1, num_sweeps))
      for (i in 1:n){
        if (sum(ztrain[i,]) == 0) {next}
        # align_catt[i, 1:sum(ztrain[i,])] = longbet.cate$cate[i, ztrain[i,] == 1]
        align_catt[i, 1:sum(ztrain[i,]), ] = longbet.pred$tauhats[i, ztrain[i,] == 1, ]
      } 
      longbet.att.sweeps <- apply(align_catt, c(2, 3), mean, na.rm = T)
      longbet.att <- longbet.att.sweeps %>% rowMeans(na.rm = T)
      longbet.att.lower <- longbet.att.sweeps %>% apply(1, quantile, prob = alpha / 2, na.rm = T)
      longbet.att.upper <- longbet.att.sweeps %>% apply(1, quantile, prob = 1 - alpha / 2, na.rm = T)
      
      longbet.catt <- apply(align_catt, c(1, 2), mean, na.rm = T)
      longbet.catt.lower <- apply(align_catt, c(1, 2), quantile, prob = alpha /2 , na.rm = T)
      longbet.catt.upper <- apply(align_catt, c(1, 2), quantile, prob = 1 - alpha /2 , na.rm = T)
      t_longbet <- proc.time() - t_longbet
      
      longbet.results <- list()
      longbet.results$iter <- iter
      longbet.results$method <- 'LongBet'
      longbet.results$ATT.RMSE <-  sqrt(mean( (att - longbet.att)^2 ) )
      longbet.results$ATT.Bias <- mean(abs(att - longbet.att))
      longbet.results$ATT.Coverage <- mean( (att >= longbet.att.lower) & (att <= longbet.att.upper) )

      longbet.results$CATT.RMSE <- sqrt(mean((align_tau - longbet.catt)^2, na.rm = T))
      longbet.results$CATT.Bias <- mean( abs( align_tau - longbet.catt ), na.rm = T)
      longbet.results$CATT.Coverage <- mean( (align_tau >= longbet.catt.lower) & (align_tau <= longbet.catt.upper), na.rm = T )
      results <- rbind(results, data.frame(longbet.results))
      
      write.csv(results, file = filename, row.names= F)
      
    }
  }
}




