# compare longbet to did method on the minimum wage dataset.
# estimate the effect of minimum wage on teen employment
# did paper: https://reader.elsevier.com/reader/sd/pii/S0304407620303948?token=967C9A6B23A76EC6C271F99CC54C0CAC8829014C025C3AB9F600C1B5F39F0A7EB6CD905A583911978C1025A798AD05D6&originRegion=us-east-1&originCreation=20230425183217
# required package: did https://github.com/bcallaway11/did

library(longBet)
library(did)
library(tidyr)
library(dplyr)
data(mpdta)

alpha <- 0.05
gp_year <- c(2004, 2006, 2007)
longbet_att <- data.frame(
  group = numeric(),
  t = numeric(),
  att = numeric(),
  upper = numeric(),
  lower = numeric()
)
for (year in gp_year){
  data <- mpdta[(mpdta$first.treat == 0) | (mpdta$first.treat == year), ]
  data$treat[data$first.treat > data$year] = 0
  
  t0 <- max(data$first.treat) - 2003 + 1 
  t1 <- length(unique(data$year))
  ttrain <- unique(data$year)
  ttreated <- ttrain[t0:t1]
  
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
  longbet.fit <- longbet(y = ytrain, x = xtrain, z = ztrain, t = ttrain,
                         num_sweeps = 100,
                         num_trees_pr =  20, num_trees_trt = 20,
                         pcat = 0,  sig_knl = 1, lambda_knl = 1)
  # TODO: lambda_knl is quite sensitve, need better understanding
  if (t1 > t0){
  sigma_knl = mean( sqrt( apply(longbet.fit$beta_draws[t0:t1,], 2, var) ))
  } else {
    sigma_knl = 1
  }
  lambda_knl = 1
  
  longbet.pred <- predict.longBet(longbet.fit, xtrain, ttrain, sigma = sigma_knl, lambda = lambda_knl)
  t_longbet <- proc.time() - t_longbet
  
  treated <- ztrain[,ncol(ztrain)]
  att_longbet_fit <- apply(longbet.pred$tauhats.adjusted[treated,,], c(2, 3), mean)[t0:t1, ]
  if (t1 > t0){
    att_df <- data.frame(
      group = rep(year, length(ttreated)),
      t = ttreated,
      att = att_longbet_fit %>% rowMeans,
      upper = apply(att_longbet_fit, 1, quantile, probs = 1 - alpha / 2),
      lower = apply(att_longbet_fit, 1, quantile, probs = alpha / 2)
    )
  } else {
    att_df <- data.frame(
      group = rep(year, length(ttreated)),
      t = ttreated,
      att = att_longbet_fit %>% mean,
      upper = quantile(att_longbet_fit, probs = 1 - alpha / 2),
      lower = quantile(att_longbet_fit, probs = alpha / 2)
    )
  }
  longbet_att <- rbind(longbet_att, att_df)
}

# did ---------------------------------------------------------------------

out <- att_gt(yname = "lemp",
              gname = "first.treat",
              idname = "countyreal",
              tname = "year",
              xformla = ~1,
              data = mpdta,
              est_method = "reg"
)
ggdid(out, ylim = c(-.25,.1))

did_att <- data.frame(
  group = out$group,
  t = out$t,
  att = out$att,
  se = out$se
)
k <- 1.96
did_att$upper <- did_att$att + k * did_att$se
did_att$lower <- did_att$att - k * did_att$se
did_att$se <- NULL

# visualization -----------------------------------------------------------

require(ggplot2)
# Plot ggdid with longbet results
longbet_att$method <- rep("LongBet", nrow(longbet_att))
did_att$method <- rep("DiD", nrow(did_att))

att_plot <- rbind(longbet_att, did_att) %>%
  ggplot(aes(x = t, y = att)) +
  geom_point(aes(color = method), position = position_dodge(.3)) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = method), width = .2,  position = position_dodge(.3)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey")+
  facet_wrap(~group, ncol = 1) +
  labs(y = "ATT")
plot(att_plot)
