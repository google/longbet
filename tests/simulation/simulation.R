
# Description -------------------------------------------------------------
# This script conduct simulation studies to evaluate the performance of Longbet
# and its baseline methods on panel data with staggered adoption

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
      
      # Longbet 500 sweeps-----------------------------------------------------------------
      longbet.time <- proc.time()
      longbet.fit <- longbet(y = ytrain, x = xtrain, x_trt = xtrain, z = ztrain, t = 1:t1,
                             num_sweeps = 500, num_trees_pr =  20, num_trees_trt = 20,
                             pcat = pcat)
      
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
      longbet.time <- proc.time() - longbet.time
      
      att.results[nrow(att.results) + 1,] <- c(iter, pr, trt, 'LongBet', att.metric(att, longbet.att), as.numeric(longbet.time[3]))
      catt.results[nrow(catt.results) + 1, ] <- c(iter, pr, trt, 'LongBet', catt.metric(align_tau, longbet.catt, longbet.catt.low, longbet.catt.high))
      
      
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
      att.results[nrow(att.results) + 1,] <- c(iter, pr, trt, 'Non-linear DiD', att.metric(att, DiD_nl), as.numeric(did_nl.time[3]))
      
      # Baseline: DiD Imputation ------------------------------------------------
      tm <- proc.time()
      did_imp.fit <- did_imputation(data = panel.data, yname = "ytrain", gname = "first.treat",
                                    tname = "time", idname = "id", 
                                    first_stage = ~ 0 | id + time,
                                    horizon=TRUE) 
      did_imp <- did_imp.fit %>% 
        select(t = term, estimate, std.error) %>%
        mutate(
          conf.low = estimate - 1.96 * std.error,
          conf.high = estimate + 1.96 * std.error,
          t = as.numeric(t)
        ) %>%
        mutate(method = "DiD Imputation") %>% 
        select(c(t, estimate, conf.low, conf.high, method)) %>% 
        filter(t >= 0 & t < 10)
      tm <- proc.time() - tm
      att.results[nrow(att.results) + 1,] <- c(iter, pr, trt, 'DiD Imputation', att.metric(att, did_imp), as.numeric(tm[3]))
      
      
      # Baseline: DiD Imputation with Covariates --------------------------------
      imp.data <- data.frame(
        outcome = as.vector( t(ytrain) ),
        treat_status = as.vector( t(ztrain) ),
        uniq_id = as.vector( sapply(1:n, rep, t1)),
        panel_time = rep(1:t1, n),
        move_date = as.vector( sapply(first.treat, rep, t1)),
        cov1 = as.vector( sapply(xtrain[,1], rep, t1)),
        cov2 = as.vector( sapply(xtrain[,2], rep, t1)),
        cov3 = as.vector( sapply(xtrain[,3], rep, t1)),
        cov4 = as.vector( sapply(xtrain[,4], rep, t1)),
        cov5 = as.vector( sapply(xtrain[,5], rep, t1))
      )
      imp.data$event_time <- sapply(imp.data$panel_time - imp.data$move_date + 1, function(x) max(x, 0))
      imp.data$event_time[imp.data$move_date == 0] <- 0
      
      imp.data$X1 <- cut(imp.data$cov1, breaks = 20, labels = 1:20)
      imp.data$X2 <- cut(imp.data$cov2, breaks = 20, labels = 1:20)
      imp.data$X3 <- cut(imp.data$cov3, breaks = 20, labels = 1:20)
      
      
      estimate_event_time_effects_weighted<-function(control_obs,treated_obs,weights_vec, include_cov){
        control_ids<-unique(control_obs$uniq_id)
        treated_ids<-unique(treated_obs$uniq_id)
        
        status_changers<-intersect(control_ids, treated_ids)
        treated_obs%>%filter(uniq_id %in% status_changers)->imputed_group
        
        #Add more interactions of panel_time with covariates here as needed
        if (include_cov){
          fixed_effect_object<-fixest::feols(outcome ~ -1 | uniq_id + panel_time^cov4^cov5^X1^X2^X3, data=control_obs, weights = weights_vec)
        } else{
          fixed_effect_object<-fixest::feols(outcome ~ -1 | uniq_id + panel_time, data=control_obs, weights = weights_vec)
        }
        predicted_outcome<-predict(fixed_effect_object,newdata=imputed_group)
        imputed_group$counterfactual_outcome<-predicted_outcome
        imputed_group$tau_est<-imputed_group$outcome-imputed_group$counterfactual_outcome
        
        imputed_group%>%
          group_by(event_time)%>%
          summarise(period_effect=mean(tau_est,na.rm = TRUE))->event_time_estimates
        
        return(list(event=event_time_estimates,individual=imputed_group))
      }
      
      #Main function to do bootstrap inference for effects
      boot_BJS_effects<-function(dataset,period_number,boot_number, alpha, include_cov){
        parentids<-unique(dataset$uniq_id)
        individual_frame<-data.frame()
        event_frame<-data.frame()
        control_obs<-dataset[dataset$treat_status==0,]
        treated_obs<-dataset[dataset$treat_status ==1,]
        control_obs%>%arrange(uniq_id)->control_obs
        control_obs%>%select(uniq_id)%>%table()->obs_counts
        for(i in (1:boot_number)){
          weight_vector<-rep(rexp(nrow(obs_counts)), obs_counts)
          BJS_results<-estimate_event_time_effects_weighted(control_obs,treated_obs,weight_vector, include_cov)
          event_study_estimates<-BJS_results$event
          event_study_estimates$boot_id<-i
          event_study_estimates<-dplyr::bind_rows(event_study_estimates,
                                                  data.frame(event_time=9999,period_effect=sum(event_study_estimates$period_effect)))
          individual_estimates<-BJS_results$individual
          individual_estimates$boot_id<-i
          event_frame<-dplyr::bind_rows(event_frame,event_study_estimates)
          individual_frame<-dplyr::bind_rows(individual_frame,individual_estimates)
        }
        
        event_frame %>%
          group_by(event_time) %>% 
          summarise(conf.low=quantile(period_effect,alpha / 2,na.rm = TRUE),
                    estimate = mean(period_effect, na.rm = TRUE),
                    conf.high=quantile(period_effect, 1 - alpha / 2,na.rm = TRUE)) %>%
          rename(t = event_time) %>%
          mutate(method = "DiD Imputation Cov")  %>%
          mutate(t = t - 1) %>%
          filter(t >= 0 & t < 999) ->event_effects
        
        # individual_frame%>%group_by(move_date,event_time)%>%summarise(low=quantile(tau_est,alpha / 2,na.rm = TRUE),
        #                                                               high=quantile(tau_est,1-alpha/2,na.rm = TRUE))->cohort_per_period_effects
        # 
        # individual_frame%>%group_by(move_date)%>%summarise(low=quantile(total_effect,alpha/2,na.rm = TRUE),
        #                                                    high=quantile(total_effect,1-alpha/2,na.rm = TRUE))->cohort_total_effects
        # 
        # results_list=list(
        #   event_times=event_effects,
        #   cohort_period=cohort_per_period_effects,
        #   cohort_total=cohort_total_effects
        # )
        return (event_effects)}
      
      # tm <- proc.time()
      # did_imp_check <- boot_BJS_effects(imp.data, t1, 500, alpha, include_cov = FALSE)
      # tm <- proc.time() - tm
      # did_imp_check %>% mutate(method = "DiD Imputation Check")
      # att.results[nrow(att.results) + 1,] <- c(iter, pr, trt, 'DiD Imputation Check', att.metric(att, did_imp_check), as.numeric(tm[3]))
      
      tm <- proc.time()
      did_imp_cov <- boot_BJS_effects(imp.data, t1, 500, alpha, include_cov = TRUE)
      tm <- proc.time() - tm
      att.results[nrow(att.results) + 1,] <- c(iter, pr, trt, 'DiD Imputation Cov', att.metric(att, did_imp_cov), as.numeric(tm[3]))
      
      
      
      write.csv(att.results, file = 'att.csv', row.names= F)
      write.csv(catt.results, file = 'catt.csv', row.names= F)
    }
  }
}

att.results <- read.csv('att.csv')
att.summary <- att.results %>% 
  group_by(pr, trt, method) %>%
  summarise(RMSE = mean(RMSE),
            Coverage = mean(Coverage),
            Cover0 = mean(Cover0),
            Time = mean(Time)) 

catt.results <- read.csv('catt.csv')
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
    factor(method, levels = c("LongBet", "DiD", "Non-linear DiD", "DiD Imputation", "DiD Imputation Cov"))) %>%
  relocate(Time, .after = last_col()) %>%
  mutate(across(where(is.numeric), ~ round(., 3)))
print(summary)
write.csv(summary, 'results.csv', row.names = F)


