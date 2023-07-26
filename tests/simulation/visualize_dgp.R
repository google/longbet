# simple demonstration of longbet with default parameters
require(dplyr)
require(ggplot2)
require(tidyr)
require(longbet)
require(panelView)
require(stringr)
require(MetBrewer)
# baseline approach
require(did) # https://github.com/bcallaway11/did
require(didimputation)
require(fixest)

set.seed(1)
n <- 2000
t0 <- 6 # treatment start time
t1 <- 12 # observed response period
alpha <- 0.05
coefs_total <- c()

source('dgp.R')

par(mfrow = c(2, 2))
for (pr_type in c("linear", "non-linear")){
  for (trt_type in c("homogeneous", "heterogeneous")){
    data <- dgp(n, t0, t1, pr_type = pr_type, trt_type = trt_type)
  
    # get training data
    ytrain <- data$y
    ztrain <- data$z
    xtrain <- data$x
    
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
    att.df <- data.frame(
      t = 0:(t1 - t0),
      estimate = att,
      conf.low = att,
      conf.high = att,
      method = "True",
      pr = pr_type,
      trt = trt_type
    )
  
   # panelview(ytrain ~ ztrain, data = panel.data, index = c("id","time"), xlab = "Time", ylab = "Unit", axis.lab.gap = 5, display.all = T)
  
    # longbet -----------------------------------------------------------------
    longbet.time <- proc.time()
    longbet.fit <- longbet(y = ytrain, x = xtrain, z = ztrain, t = 1:t1,
                           num_sweeps = 60,
                           num_trees_pr =  20, num_trees_trt = 20,
                           pcat = ncol(xtrain) - 3)
    
    longbet.pred <- predict.longbet(longbet.fit, xtrain, ztrain)
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
    longbet.time <- proc.time() - longbet.time
    
    longbet.results <- list()
    longbet.results$method <- 'LongBet'
    longbet.results$ATT.RMSE <-  sqrt(mean( (att - longbet.att)^2 ) )
    longbet.results$ATT.Bias <- mean(abs(att - longbet.att))
    longbet.results$ATT.Coverage <- mean( (att >= longbet.att.lower) & (att <= longbet.att.upper) )
    
    longbet.results$CATT.RMSE <- sqrt(mean((align_tau - longbet.catt)^2, na.rm = T))
    longbet.results$CATT.Bias <- mean( abs( align_tau - longbet.catt ), na.rm = T)
    longbet.results$CATT.Coverage <- mean( (align_tau >= longbet.catt.lower) & (align_tau <= longbet.catt.upper), na.rm = T )
    longbet.results$Time <- longbet.time[1]
    print(data.frame(longbet.results))
    
    longbet.results <- data.frame(
      t = 0:(t1 - t0),
      estimate = longbet.att,
      conf.low = longbet.att.lower,
      conf.high = longbet.att.upper,
      method = "LongBet",
      pr = pr_type,
      trt = trt_type
    )
    
    
    # Baseline: DiD with multiple periods --------------------------------------------------------------
    did.time <- proc.time()
    did.out <- att_gt(yname = "ytrain",
                      gname = "first.treat",
                      idname = "id",
                      tname = "time",
                      xformla = ~ X1 + X2 + X3 + X4 + X5,
                      data = panel.data,
                      est_method = "dr",
                      control_group = "notyettreated"
    )
    did.es <- aggte(did.out, type = "dynamic")
    did.att <- did.es$att.egt[did.es$egt >= 0]
    did.att.lower <- (did.es$att.egt + qnorm(alpha / 2) * did.es$se.egt)[did.es$egt >= 0]
    did.att.upper <- (did.es$att.egt + qnorm(1 - alpha / 2) * did.es$se.egt)[did.es$egt >= 0]
    # ggdid(es)
    did.time <- proc.time() - did.times
    
    DiD <- aggte(did.out, type = "dynamic", na.rm = TRUE) %>% 
      tidy() %>% 
      rename(t = event.time) %>% 
      filter(t >= 0) %>% 
      select(t, estimate, conf.low, conf.high) %>% 
      mutate(method = "DiD", pr = pr_type, trt = trt_type)
    
    
    # Baseline: DiD Non-linear ------------------------------------------------
    did_nl.time <- proc.time()
    if (pr_type == "linear"){
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
      mutate(method = "Non-linear DiD", pr = pr_type, trt = trt_type)
    did_nl.time = proc.time() - did_nl.time
    
    # results -----------------------------------------------------------------
    # coefs <- bind_rows(twfe, stacked, CSA, SA, coef_imp, asyn, fect, PM) 
    coefs <- bind_rows(att.df, longbet.results, DiD, DiD_nl)
    coefs_total <- rbind(coefs_total, coefs)
  }
}

plot <- 
  coefs_total %>% 
  transform(pr = factor(pr,levels=c("linear", "non-linear")), 
            trt = factor(trt, levels = c("homogeneous", "heterogeneous") ) ) %>% 
  arrange(pr, trt) %>%
  # filter(method %in% c("True", "LongBet", "DiD")) %>%
  ggplot(aes(x = t, y = estimate, color = method)) + 
  geom_point(aes(x = t, y = estimate), position = position_dodge2(width = 0.8), size = 1) +
  geom_linerange(aes(x = t, ymin = conf.low, ymax = conf.high), position = position_dodge2(width = 0.8), size = 0.75) +
  geom_line(data = subset(coefs_total, method %in% "True"), aes(x = t, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = .25, alpha = 0.75) + 
  geom_vline(xintercept = -0.5, linetype = "dashed", size = .25) +
  facet_wrap(vars(pr, trt), scales = "free_y") + 
  scale_color_manual(name="Estimation Method", values= met.brewer("Cross", 8, "discrete")) +
  theme(legend.position= 'bottom') +
  labs(title = 'Event Time Estimates', y="ATT", x = "Relative Time") + 
  guides(col = guide_legend(nrow = 3)) 
print(plot)


