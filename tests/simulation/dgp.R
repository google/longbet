# Required packages
require(forecast)

# Data generating function ----------------------------------------------------- 
# Input: 
# n: number of observations
# t0: treatment start time
# t1: total time period
# pr_type: prognostic effect types ("parallel" or "non-parallel")
# trt_type: treatment effect types ("homogeneous" or "heterogeneous")

# Output: 
# x: covariate matrix
# y: response panel
# z: treatment panel
# y0: potential outcomes 
# tau: treatment effect

dgp <- function(n, t0 = 6, t1 = 12, t2 = NULL, 
                pr_type = "non-parallel", trt_type = "heterogeneous",
                staggered_effect = NULL){
  # t2: Extended period for prediction
  
  # screening
  if (t0 > t1){
    stop("Treatment start time (t0) can not be greater than t1")
  }
  if(!(pr_type %in% c("parallel", "non-parallel"))){
    stop("pr_type should be either parallel or non-parallel")
  }
  if(!(trt_type %in% c("homogeneous", "heterogeneous"))){
    stop("trt_type should be either homogeneous or heterogeneous")
  }
  
  if (is.null(t2)) {
    t2 <- t1
  }
  
  results <- list()
  # generate covariates
  x <- matrix(data = NA, nrow = n, ncol = 5)
  x[,1] <- rnorm(n)
  x[,2] <- rnorm(n)
  x[,3] <- rnorm(n)
  x[,4] <- rbinom(n,1,0.5)
  x[,5] <- factor(sample(1:3,n,replace=TRUE,prob = c(1/3, 1/3, 1/3)))
  results$x <- x
  
  # generate prognostic effect
  f_t <- arima.sim(model = list(order = c(1, 0, 1), ar = 0.7, ma = -0.4), n = t2) + 1
  g <- c(2, -1)
  if(pr_type == "parallel"){
    gamma_x <-  g[x[, 4] + 1] + 1 * x[, 1] * abs(x[, 3] - 1)
    y0 <- t( matrix(rep(f_t, n), t2, n) ) + matrix(rep(gamma_x, t2), n, t2)
  } else if (pr_type == "non-parallel"){
    gamma_x <-  g[x[, 4] + 1] + 1 * x[, 1] * abs(x[, 3] - 1)
    y0 <- outer(gamma_x, f_t, "*")
  }
  results$y0 <- y0
  
  # generate treatment effect 
  s <- 1:(t2 - t0 + 1) 
  h_s <- s * exp(-s)
  if(trt_type == "homogeneous"){
    nu_x <- rep(2, n)
  } else if (trt_type == "heterogeneous"){
    nu_x <- 2 + x[, 2] * x[, 5]
  }
  tau <- outer(nu_x, h_s , "*")
  
  # compute propensity scores and generate treatment assignment
  pi <- 0.2 * pnorm(0.5* gamma_x - 0.5 * x[,1])^2 + runif(n) / 10
  z <- matrix(0, n, t2)
  for (i in t0:t1){ # only assign treatment in observed period
    treated <- (z[, i - 1] == 1)
    treatment <- rbinom(n, 1, pi) # draw new treatment for time i
    z[, i] <- apply(cbind(treated, treatment), 1, max) 
  }
  if (t2 > t1){
    z[,(t1 + 1): t2] <- matrix(rep(z[,t1], t2 - t1), n, t2 - t1)
  }
  results$z <- z
  
  # get post treatment time matrix
  post_t <- t(apply(z, 1, cumsum))
  
  # add staggered effect
  if (is.null(staggered_effect)){
    staggered_effect <- rep(1, t2 - t0 + 1)
  } else {
    # check staggered effect size
    if (length(staggered_effect) < t1 - t0 + 1){
      print("Staggered effect should have length equal or greater than t1 - t0")
      stop()
    }
  }
  get_treatment_effect <- function(ground_truth, trt_time){
    trt_effect <- rep(0, length(trt_time))
    eb_effect <- staggered_effect[sum(trt_time == 0) - t0 + 2]
    trt_effect[trt_time > 0] = ground_truth[trt_time[trt_time > 0]] * eb_effect
    return(trt_effect)
  }
  tau <- t(mapply(get_treatment_effect, data.frame(t(tau)), data.frame(t(post_t))))
  results$tau <- tau
  
  # generate outcome variable
  results$y <-  y0 + tau + matrix(rnorm(n*t2, mean = 0, sd = 0.2), n, t2)
  
  
  return(results)
}

getCohort <- function(z){
  # z: n by t treatment indicator at time t
  # get cohorts that adopt the treatment at the same time
  # 0 for control group
  z.rowSums <- rowSums(z)
  cohort <- max(z.rowSums) - z.rowSums + 1
  cohort[z.rowSums == 0] <- 0
  return(cohort)
}

# Calculate metrics -------------------------------------------------------

getAttHat <- function(tauhat){
  # tauhat: tauhats from longbet prediction, n by t by sweeps
  # return: return estimated att across time
  att.hat <- apply(tauhat, c(2, 3), mean) %>% rowMeans
  return(att.hat)
}

getAttCI <- function(tauhat, alpha = 0.05){
  # tauhat: tauhats from longbet prediction, n by t by sweeps
  # return: return credible intervals across 
  att.hats <- apply(tauhat, c(2, 3), mean)
  lower <- apply(att.hats, 1, quantile, probs = alpha/2)
  upper <- apply(att.hats, 1, quantile, probs = 1 - alpha/2)
  return(list(lower = lower, upper = upper))
}

getCohortAtt <- function(tau, cohort){
  # return dataframe of ATT for each group at each time
  tau.df <- data.frame(tau)
  tau.df$cohort <- cohort
  df <- tau.df %>% 
    group_by(cohort) %>%
    summarize(across(starts_with("X"), mean)) %>%
    gather(key = "Time", value = "ATT", -cohort) %>%
    mutate(Time = as.numeric(gsub("X", "", Time))) %>% 
    filter(cohort != 0) %>%
    filter(Time >= cohort)
  return(df)
}

getCohortAttBackup <- function(tau, cohort){
  # return dataframe of ATT for each group at each time
  tau.df <- data.frame(tau)
  tau.df$cohort <- cohort
  df <- tau.df %>% 
    group_by(cohort) %>%
    summarize(across(starts_with("X"), mean))
  return(df)
}

getCohortAttHat <- function(tauhat, cohort){
  # tauhat: tauhats from longbet prediction, n by t by sweeps
  # cohort: a vector of length n indicating group for each unit
  # return: dataframe of ATT for each group at each time
  
  unique.cohort <- sort(unique(cohort))
  df <- data.frame(cohort = unique.cohort, matrix(NA, nrow = length(unique.cohort), ncol = dim(tauhat)[2]))
  for (i in 1:length(unique.cohort)){
    df[i,2:ncol(df)] <- getAttHat(tauhat[cohort == unique.cohort[i],,])
  }
  mutate_df <-  df %>%
    gather(key = "Time", value = "True", -cohort) %>%
    mutate(Time = as.numeric(gsub("X", "", Time))) %>% 
    filter(cohort != 0) %>%
    filter(Time >= cohort)
  return(mutate_df)
}

getCohortAttHatBackup <- function(tauhat, cohort){
  # tauhat: tauhats from longbet prediction, n by t by sweeps
  # cohort: a vector of length n indicating group for each unit
  # return: dataframe of ATT for each group at each time
  
  unique.cohort <- sort(unique(cohort))
  df <- data.frame(cohort = unique.cohort, matrix(NA, nrow = length(unique.cohort), ncol = dim(tauhat)[2]))
  for (i in 1:length(unique.cohort)){
    df[i,2:ncol(df)] <- getAttHat(tauhat[cohort == unique.cohort[i],,])
  }
  return(df)
}

getCohortAttLongBet<- function(tauhat, cohort, alpha = 0.05){
  # tauhat: tauhats from longbet prediction, n by t by sweeps
  # cohort: a vector of length n indicating group for each unit
  # return: list of dataframe of credible intervals for each group at each time
  
  unique.cohort <- sort(unique(cohort))
  estimate <- data.frame(cohort = unique.cohort, matrix(NA, nrow = length(unique.cohort), ncol = dim(tauhat)[2]))
  lower <- data.frame(cohort = unique.cohort, matrix(NA, nrow = length(unique.cohort), ncol = dim(tauhat)[2]))
  upper <- data.frame(cohort = unique.cohort, matrix(NA, nrow = length(unique.cohort), ncol = dim(tauhat)[2]))
  for (i in 1:length(unique.cohort)){
    estimate[i,2:ncol(estimate)] <- getAttHat(tauhat[cohort == unique.cohort[i],,])
    
    CI <- getAttCI(tauhat[cohort == unique.cohort[i],,], alpha)
    lower[i,2:ncol(lower)] <- CI$lower
    upper[i,2:ncol(upper)] <- CI$upper
  }
  df <-  estimate %>%
    gather(key = "Time", value = "ATT", -cohort) %>%
    mutate(Time = as.numeric(gsub("X", "", Time))) %>% 
    left_join(lower %>%
                gather(key = "Time", value = "conf.low", -cohort) %>%
                mutate(Time = as.numeric(gsub("X", "", Time))),
              by = c("cohort", "Time")) %>%
    left_join(upper %>%
                gather(key = "Time", value = "conf.high", -cohort) %>%
                mutate(Time = as.numeric(gsub("X", "", Time))),
              by = c("cohort", "Time")) %>%
    filter(cohort != 0) %>%
    filter(Time >= cohort) 
  return(df)
}


getCohortAttDiD <- function(did, cohort.att, alpha, t0 = 1){
  # cohort.att: take reference for which cohort is of interest
  
  # adjust did cohort and time
  did$group <- did$group - t0 + 1
  did$t <- did$t - t0 + 1
  did.df <- data.frame(cohort = did$group, Time = did$t, ATT = did$att, SE = did$se)
  
  # merge
  cohort.att$ATT <- NULL
  df <- merge(cohort.att, did.df, by = c("cohort", "Time"), how = "left")
  
  # get 95% simultaneous confidence band
  df$conf.low <- df$ATT +  df$SE * qnorm(alpha / 2 / nrow(df))
  df$conf.high <- df$ATT +  df$SE * qnorm(1 - alpha / 2 / nrow(df))
  df$SE <- NULL
  return(df)
}


getCohortAttCIBackup <- function(tauhat, cohort, alpha = 0.05){
  # tauhat: tauhats from longbet prediction, n by t by sweeps
  # cohort: a vector of length n indicating group for each unit
  # return: list of dataframe of credible intervals for each group at each time
  
  unique.cohort <- sort(unique(cohort))
  lower <- data.frame(cohort = unique.cohort, matrix(NA, nrow = length(unique.cohort), ncol = dim(tauhat)[2]))
  upper <- data.frame(cohort = unique.cohort, matrix(NA, nrow = length(unique.cohort), ncol = dim(tauhat)[2]))
  for (i in 1:length(unique.cohort)){
    CI <- getAttCI(tauhat[cohort == unique.cohort[i],,], alpha)
    lower[i,2:ncol(lower)] <- CI$lower
    upper[i,2:ncol(upper)] <- CI$upper
  }
  return(list(lower = lower, upper = upper))
}

getCohortMetrics <- function(att, att.hat, att.ci = NULL){
  # att: ATT for each cohort, the first column is cohort
  # att.hat: Estimated att for each cohort
  # att.ci: A list of credible intervals of estimated att for each cohort
  # return: A table of RMSE and Coverage (if provided) on estimated ATT for each cohort
  
  # check if cohort in att match the one in att.hat
  if (names(att)[1] != 'cohort'){
    print("The first column in att should be cohort")
    stop()
  }
  if (names(att.hat)[1] != 'cohort'){
    print("The first column in att.hat should be cohort")
    stop()
  }
  if (any(att$cohort != att.hat$cohort)){
    print("The cohort column in att.hat should match with att")
    stop()
  }
  if (ncol(att) != ncol(att.hat)){
    print("The dimension of att.hat does not match with att table")
    stop()
  }
   
  df <- data.frame(cohort = att$cohort, RMSE = NA, coverage = NA)
  for (i in 1:nrow(df)){
    df$RMSE[i] <- sqrt( mean( as.numeric( (att[i,-1] - att.hat[i, -1])^2 )) ) # get rmse
  }
  
  if (is.null(att.ci)){
    df$coverage <- NULL
  } else {
    if (any(names(att.ci) != c('lower', 'upper'))){
      print("att.ci should have two tables named 'lower' and 'upper'.")
      stop()
    }
    for (i in 1:nrow(df)){
      df$coverage[i] <- mean( (att[i, 2:ncol(att)] >= att.ci$lower[i, 2:ncol(att)]) & (att[i, 2:ncol(att)] <= att.ci$upper[i, 2:ncol(att)]) )
    }
  }
  return(df)
}

att.metric <- function(att, estimate){
  # check estimate has the right form
  if (!is.data.frame(estimate)){
    stop("Estimate needs to be a data frame")
  }
  if (!all(c("t", "estimate", "conf.low", "conf.high", "method")%in% colnames(estimate))){
    stop("Estimate not having the right columns")
  }
  if (any (estimate$t != 0:(length(att) - 1))){
    stop(paste0("Estimate t should be ", toString( 0:(length(att) - 1) ) ))
  }
  metrics <- c()
  metrics['RMSE'] <- sqrt(mean( (att - estimate$estimate)^2 ))
  metrics['Bias'] <- mean(abs(att - estimate$estimate))
  metrics['Coverage'] <- mean( (att >= estimate$conf.low) & (att <= estimate$conf.high))
  metrics['Cover0'] <-  mean( (0 >= estimate$conf.low) & (0 <= estimate$conf.high))
  metrics['I.L'] <- mean(estimate$conf.high - estimate$conf.low)
  return(metrics)
}

cohort.att.metric <- function(att, estimate){
  # check estimate has the right form
  if (!is.data.frame(estimate)){
    stop("Estimate needs to be a data frame")
  }
  if (!all(c("cohort", "Time", "ATT", "conf.low", "conf.high")%in% colnames(estimate))){
    stop(cat("Estimate table must have the following columns: ", c("cohort", "Time", "ATT", "conf.low", "conf.high")))
  }
  # merge and check
  df <- merge(att, estimate, by = c("cohort", "Time"), suffixes = c(".true", ".est"), how = "left")
  
  if (any(is.null(df[, c("ATT.est", "conf.low", "conf.high")]))){
    print("estimate table is not complete")
    stop()
  }
  
  metrics <- c()
  metrics['RMSE'] <- sqrt(mean( (df$ATT.true - df$ATT.est)^2 ))
  metrics['Bias'] <- mean(abs(df$ATT.true - df$ATT.est))
  metrics['Coverage'] <- mean( (df$ATT.true >= df$conf.low) & (df$ATT.true <= df$conf.high))
  metrics['Cover0'] <-  mean( (0 >= df$conf.low) & (0 <= df$conf.high))
  metrics['I.L'] <- mean(df$conf.high - df$conf.low)
  return(metrics)
}

catt.metric.aligned <- function(align_catt, estimate, conf.low, conf.high){
  # check dimensions
  if (any(dim(align_catt) != dim(estimate))){
    stop("Dimension of estimate does not match aligned catt")
  }
  if (any(dim(align_catt) != dim(conf.low))){
    stop("Dimension of conf.low does not match aligned catt")
  }
  if (any(dim(align_catt) != dim(conf.high))){
    stop("Dimension of conf.high does not match aligned catt")
  }
  metrics <- c()
  metrics['RMSE'] <- sqrt(mean((align_catt - estimate)^2, na.rm = T))
  metrics['Bias'] <- mean( abs( align_catt - estimate ), na.rm = T)
  metrics['Coverage'] <- mean( (align_catt >= conf.low) & (align_tau <= conf.high), na.rm = T )
  metrics['Cover0'] <-  mean( (0 >= conf.low) & (0 <= conf.high), na.rm = T)
  metrics['I.L'] <- mean(conf.high - conf.low, na.rm = T)
  return(metrics)
}

catt.metric.longbet <- function(tau, tauhats, z, alpha = 0.05){
  # z: treatment metric, only measure treatment effect when treatment is assigned
  
  # check dimension
  if (any(dim(tau) != dim(tauhats)[1:2])){
    print("Dimensions of tau and tauhats do not match")
    stop()
  }
  if(any(dim(tau) != dim(z))){
    print(" Dimensions of tau and z do not match")
    stop()
  }
  
  treated <- z==1
  estimate <- apply(tauhats, c(1, 2), mean)
  conf.low <- apply(tauhats, c(1, 2), quantile, probs = alpha/2)
  conf.high <- apply(tauhats, c(1, 2), quantile, probs = 1 - alpha/2)
  
  metrics <- c()
  metrics['RMSE'] <- sqrt(mean((tau[treated] - estimate[treated])^2, na.rm = T))
  metrics['Bias'] <- mean( abs( tau[treated] - estimate[treated] ), na.rm = T)
  metrics['Coverage'] <- mean( (tau[treated] >= conf.low[treated]) & (tau[treated] <= conf.high[treated]), na.rm = T )
  metrics['Cover0'] <-  mean( (0 >= conf.low[treated]) & (0 <= conf.high[treated]), na.rm = T)
  metrics['I.L'] <- mean(conf.high[treated] - conf.low[treated], na.rm = T)
  return(metrics)
}

