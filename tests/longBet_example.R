require(longBet)

set.seed(1)
n <- 100
t1 <- 4
t0 <- 3

# generate dcovariates
x1 <- rnorm(n)
x2 <- sample(1:3,n,replace=TRUE, prob = c(0.4,0.3,0.3))
# TODO: memory bug occurs when there's no categorical variable 
x <- cbind(x1, x2)

# untreated outcome
mu <- outer(x1 * x2 , rnorm(t1, 5), '*')
# treatment effect
te <- outer(x1 + x2, rnorm(t1, 1), '*')

# generate observations
z <- rbinom(n,1,0.6)
y0 <- mu + 0.2 * sd(mu) * matrix(rnorm(n * t1), n, t1)
y1 <- y0 + te
y <- z * y1 + (1 - z) * y0

z_mat <- cbind(matrix(0, n, (t0 - 1)),  matrix(rep(z, t1 - t0 + 1), n, t1 - t0 + 1))

t_longbet <- proc.time()
longbet.fit <- longbet(y = y, x = x, z = z_mat, t = 1:t1, pcat = 1,
num_trees_pr =  50, num_trees_trt = 50)
longbet.pred <- predict.longBet(longbet.fit, x, z_mat)
longbet.ate <- get_ate(longbet.pred, alpha = 0.05)
longbet.cate <- get_cate(longbet.pred, alpha = 0.05)
t_longbet <- proc.time() - t_longbet

ate <- colMeans(te)
print(paste0("longBet CATE RMSE: ", sqrt(mean((as.vector(longbet.cate$cate[, t0:t1]) - as.vector(te[,t0:t1]))^2))))
print(paste0("longBet ATE RMSE: ", sqrt(mean((as.vector(longbet.ate$ate[t0:t1]) - as.vector(ate[t0:t1]))^2))))
print(paste0("longBet runtime: ", round(as.list(t_longbet)$elapsed,2)," seconds"))
