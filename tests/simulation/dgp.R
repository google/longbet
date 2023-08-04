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

dgp <- function(n, t0 = 6, t1 = 12, pr_type = "non-parallel", trt_type = "heterogeneous"){
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
  f_t <- arima.sim(model = list(order = c(1, 0, 1), ar = 0.7, ma = -0.4), n = t1) + 1
  g <- c(2, -1)
  if(pr_type == "parallel"){
    gamma_x <- 1 + g[x[, 4] + 1] +  5 * x[, 1] + 2 * x[, 3]
    y0 <- t( matrix(rep(f_t, n), t1, n) ) + matrix(rep(gamma_x, t1), n, t1)
  } else if (pr_type == "non-parallel"){
    gamma_x <-  1 + g[x[, 4] + 1] + 6 * x[, 1] * abs(x[, 3] - 1)
    y0 <- outer(gamma_x, f_t, "*")
  }
  results$y0 <- y0
  
  # generate treatment effect 
  s <- 1:(t1 - t0 + 1) 
  h_s <- s * exp(-s)
  trt_type = "homogeneous"
  if(trt_type == "homogeneous"){
    nu_x <- rep(2, n)
  } else if (trt_type == "heterogeneous"){
    nu_x <- 2 + 5 * x[, 2] * x[, 5]
  }
  tau <- outer(nu_x, h_s , "*")
  
  # compute propensity scores and generate treatment assignment
  pi <- 0.2 * pnorm(0.5* gamma_x - 0.5 * x[,1])^2 + runif(n) / 10
  z <- matrix(0, n, t1)
  for (i in t0:t1){
    treated <- (z[, i - 1] == 1)
    treatment <- rbinom(n, 1, pi) # draw new treatment for time i
    z[, i] <- apply(cbind(treated, treatment), 1, max) 
  }
  results$z <- z
  
  # get post treatment time matrix
  post_t <- t(apply(z, 1, cumsum))
  
  get_tau <- function(trt, post_t){
    tau <- rep(0, length(post_t))
    tau[post_t > 0] = trt[post_t[post_t > 0]]
    return(tau)
  }
  tau <- t(mapply(get_tau, data.frame(t(tau)), data.frame(t(post_t))))
  results$tau <- tau
  
  # generate outcome variable
  results$y <-  y0 + tau + matrix(rnorm(n*t1, mean = 0, sd = 0.2), n, t1)
  
  return(results)
}


# Calculate metrics -------------------------------------------------------
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

catt.metric <- function(align_catt, estimate, conf.low, conf.high){
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

