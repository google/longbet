
# Description -------------------------------------------------------------
# This script conduct simulation studies to evaluate the performance of Longbet
# and its baseline methods on panel data with staggered adoption
# Added staggered effect and evaluate by cohort and time.

# Required libraries ------------------------------------------------------
require(longbet)
require(dplyr)
require(ggplot2)
require(tidyr)
require(panelView)
require(stringr)
require(MetBrewer)
# baseline approach
require(did) # https://github.com/bcallaway11/did
require(fixest)
require(didimputation)

source('dgp.R') # script for data generating process

# Set up ------------------------------------------------------------------
set.seed(100)
mc <- 100 #1000 # monte carlo iterations
n <- 2000      # number of observations
t0 <- 6        # earliest treatment adoption time
t1 <- 12       # total time period
alpha <- 0.05
pr_types <- c("parallel", "non-parallel")
trt_types <- c("homogeneous", "heterogeneous")
staggered_effect <- seq(1.3, 0.8, length.out = t1 - t0 + 1)
pcat <- 2    # number of categorical variable

if (file.exists('cohort_att.csv') & file.exists('cohort_catt.csv')){
  att.results <- read.csv('cohort_att.csv')
  catt.results <- read.csv('cohort_catt.csv')
} else {
  att.results <- data.frame(
    iter = double(),
    pr = character(),
    trt = character(),
    method = character(),
    RMSE = double(),
    Bias = double(),
    Coverage = double(),
    Cover0 = double(),
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
    Cover0 = double(),
    I.L = double()
  )
}

# simulation ---------------------------------------------------------------
for (pr in pr_types){
  for (trt in trt_types){
    for (iter in 1:mc){
      if (any((att.results$pr == pr) & (att.results$trt == trt))){
        if ( iter %in% (att.results$iter[(att.results$pr == pr) & (att.results$trt == trt)]) ){
          # skip this iteration
          next 
        }
      }
      
      # Data generating
      data <- dgp(n, t0, t1, pr_type = pr, trt_type = trt, staggered_effect = staggered_effect)
      xtrain <- data$x
      ytrain <- data$y
      ztrain <- data$z
      
      cohort <- getCohort(ztrain)
      xmod <- cbind(cohort, xtrain)
      # get att by cohort
      cohort.att <- getCohortAtt(data$tau[,t0:t1], cohort)
      
      
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
      
      
      # Longbet -----------------------------------------------------------------
      longbet.time <- proc.time()
      longbet.fit <- longbet(y = ytrain, x = xtrain, x_trt = xmod, z = ztrain, t = 1:t1,
                             num_sweeps = 100, num_trees_pr =  50, num_trees_trt = 50,
                             pcat = pcat, pcat_trt = pcat)
      
      longbet.pred <- predict.longbet(longbet.fit, xtrain, xmod, ztrain)
      longbet.att.cohort <- getCohortAttLongBet(longbet.pred$tauhats[, t0:t1, ], cohort, alpha = alpha / nrow(cohort.att))
      longbet.time <- proc.time() - longbet.time
      
      att.results[nrow(att.results) + 1,] <- c(iter, pr, trt, 'LongBet', cohort.att.metric(cohort.att,longbet.att.cohort), as.numeric(longbet.time[3]))
      catt.results[nrow(catt.results) + 1, ] <- c(iter, pr, trt, 'LongBet', catt.metric.longbet(data$tau, longbet.pred$tauhats, data$z, alpha))
      
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
      did.time <- proc.time() - did.time
      DiD <- getCohortAttDiD(did.out, cohort.att, alpha, t0)
      att.results[nrow(att.results) + 1,] <- c(iter, pr, trt, 'DiD', cohort.att.metric(cohort.att, DiD), as.numeric(did.time[3]))
      
      # Baseline: DiD Non-linear ------------------------------------------------
      did_nl.time <- proc.time()
      did_nl.out <- att_gt(yname = "ytrain",
                           gname = "first.treat",
                           idname = "id",
                           tname = "time",
                           xformla = ~ X1 * X3 + X2 + X4 + X5,
                           data = panel.data,
                           est_method = "dr",
                           control_group = "notyettreated")
      
      
      DiD_nl <- aggte(did_nl.out, type = "dynamic", na.rm = TRUE) %>% 
        tidy() %>% 
        rename(t = event.time) %>% 
        filter(t >= 0 & t < 8) %>% 
        select(t, estimate, conf.low, conf.high) %>% 
        mutate(method = "Non-linear DiD")
      did_nl.time = proc.time() - did_nl.time
      
      DiD_nl <- getCohortAttDiD(did_nl.out, cohort.att, alpha, t0)
      att.results[nrow(att.results) + 1,] <- c(iter, pr, trt, 'Non-linear DiD', cohort.att.metric(cohort.att, DiD_nl), as.numeric(did_nl.time[3]))
      
      # Baseline: DiD Imputation ------------------------------------------------
      imp.tm <- proc.time()
      imp.data <- panel.data
      imp.data$cohort <- imp.data$first.treat - t0 + 1
      # get weight for each cohort and time
      wtr <- rep(NA, nrow(cohort.att))
      for(i in 1:nrow(cohort.att)){
        c <- cohort.att$cohort[i]
        t <- cohort.att$Time[i] + t0 - 1
        imp.data[, paste("cohort", c, "_time", t, sep = "")] <- (imp.data$time == t) & (imp.data$cohort == c)
        wtr[i] <- paste("cohort", c, "_time", t, sep = "")
      }
      did_imp.fit <- did_imputation(data = imp.data, yname = "ytrain", gname = "first.treat",
                                    tname = "time", idname = "id",  wtr = wtr,
                                    first_stage = ~ 0 | id + time) 
      
      
      did_imp <- did_imp.fit %>%
        separate(term, into = c("cohort", "Time"), sep = "_") %>% # get cohort and time
        mutate(cohort = gsub("cohort", "", cohort) %>% as.numeric,
               Time = (gsub("time", "", Time) %>% as.numeric) - t0 + 1, 
               ATT = estimate,
               conf.low = estimate + std.error * qnorm(alpha / 2 / nrow(did_imp.fit)), 
               conf.high = estimate + std.error * qnorm(1 - alpha / 2 / nrow(did_imp.fit))) %>% # get 95% simult band
        select(-lhs, -estimate, -std.error)
      
      imp.tm <- proc.time() - imp.tm
      att.results[nrow(att.results) + 1,] <- c(iter, pr, trt, 'DiD Imputation', cohort.att.metric(cohort.att, did_imp), as.numeric(imp.tm[3]))
      
    }
  }
}

att.results <- read.csv('cohort_att.csv')
att.summary <- att.results %>% 
  group_by(pr, trt, method) %>%
  summarise(RMSE = mean(RMSE),
            Coverage = mean(Coverage),
            Cover0 = mean(Cover0),
            Time = mean(Time)) 

catt.results <- read.csv('cohort_catt.csv')
catt.summary <- catt.results %>% 
  group_by(pr, trt, method) %>%
  summarise(RMSE = mean(RMSE),
            Coverage = mean(Coverage),
            Cover0 = mean(Cover0)) 

summary <- merge(att.summary, catt.summary, by = c('pr', 'trt', 'method'), 
                 suffixes = c(".ATT",".CATT"), all.x = T) %>%
  arrange(
    factor(pr, levels = c("parallel", 'non-parallel')),
    factor(trt, levels = c("homogeneous", "heterogeneous")),
    factor(method, levels = c("LongBet", "DiD", "Non-linear DiD", "DiD Imputation"))) %>%
  relocate(Time, .after = last_col()) %>%
  mutate(across(where(is.numeric), ~ round(., 3)))
print(summary)
write.csv(summary, 'cohort_results.csv', row.names = F)


