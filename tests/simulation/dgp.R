# Required packages
require(forecast)

# Data generating function ----------------------------------------------------- 
# Input: 
# n: number of observations
# t0: treatment start time
# t1: total time period
# pr_type: prognostic effect types ("linear" or "non-linear")
# trt_type: treatment effect types ("homogeneous" or "heterogeneous")

# Output: 
# x: covariate matrix
# y: response panel
# z: treatment panel
# y0: potential outcomes 
# tau: treatment effect

dgp <- function(n, t0 = 6, t1 = 12, pr_type = "non-linear", trt_type = "heterogeneous"){
  # screening
  if (t0 > t1){
    stop("Treatment start time (t0) can not be greater than t1")
  }
  if(!(pr_type %in% c("linear", "non-linear"))){
    stop("pr_type should be either linear or non-linear")
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
  if(pr_type == "linear"){
    gamma_x <- 1 + g[x[, 4] + 1] +  2 * x[, 3]
  } else if (pr_type == "non-linear"){
    gamma_x <- -6 + g[x[, 4] + 1] + 6 * x[, 1] * abs(x[, 3] - 1)
  }
  y0 <- outer(gamma_x, f_t, "*")
  results$y0 <- y0
  
  # generate treatment effect 
  s <- 1:(t1 - t0 + 1) 
  h_s <- s * exp(-s)
  if(trt_type == "homogeneous"){
    nu_x <- rep(1, n)
  } else if (trt_type == "heterogeneous"){
    nu_x <- 1 + 5 * x[, 2] * x[, 5]
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
  results$y <-  y0 + tau + matrix(rnorm(n*t1, mean = 0, sd = 0.5), n, t1)
  
  return(results)
}
