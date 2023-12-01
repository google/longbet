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

# DATA GENERATION PROCESS -------------------------------------------------
set.seed(3)
n <- 2000
t0 <- 7 # treatment start time
t1 <- 12 # observed response period
alpha <- 0.05

source('dgp.R')
pr_type = "non-parallel"
trt_type = "heterogeneous"
staggered_effect <- seq(1.3, 0.8, length.out = t1 - t0 + 1)
# staggered_effect <- rep(1, length.out = t1 - t0 + 1)
data <- dgp(n, t0, t1, pr_type = pr_type, trt_type = trt_type, 
            staggered_effect = staggered_effect )

# get training data
ytrain <- data$y
ztrain <- data$z
xtrain <- data$x

cohort <- getCohort(ztrain)
xmod <- cbind(xtrain, cohort)
# xmod <- xtrain

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

att.results <- data.frame(
  method = character(),
  RMSE = double(),
  Bias = double(),
  Coverage = double(),
  Cover0 = double(),
  I.L = double(),
  Time = double()
)

catt.results <- data.frame(
  method = character(),
  RMSE = double(),
  Bias = double(),
  Coverage = double(),
  Cover0 = double(),
  I.L = double()
)


# longbet -----------------------------------------------------------------
longbet.time <- proc.time()
longbet.fit <- longbet(y = ytrain, x = xtrain, x_trt = xmod, z = ztrain, t = 1:t1,
                       num_sweeps = 100,
                       num_trees_pr =  50, num_trees_trt = 50,
                       pcat = ncol(xtrain) - 3, pcat_trt = ncol(xmod) - 3)

longbet.pred <- predict.longbet(longbet.fit, xtrain, xmod, ztrain)
longbet.att.cohort <- getCohortAttLongBet(longbet.pred$tauhats[, t0:t1, ], cohort, alpha = alpha / nrow(cohort.att))
longbet.time <- longbet.time - proc.time()

att.results[nrow(att.results) + 1,] <- c('LongBet', cohort.att.metric(cohort.att,longbet.att.cohort), as.numeric(longbet.time[3]))
catt.results[nrow(catt.results) + 1, ] <- c('LongBet', catt.metric.longbet(data$tau, longbet.pred$tauhats, data$z, alpha))


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
did.time <- did.time - proc.time()

DiD <- getCohortAttDiD(did.out, cohort.att, alpha, t0)

att.results[nrow(att.results) + 1,] <- c('DiD', cohort.att.metric(cohort.att, DiD), as.numeric(did.time[3]))

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
did_nl.time <- did_nl.time - proc.time()

DiD_nl <- getCohortAttDiD(did_nl.out, cohort.att, alpha, t0)
att.results[nrow(att.results) + 1,] <- c('Non-linear DiD', cohort.att.metric(cohort.att, DiD_nl), as.numeric(did_nl.time[3]))


# Baseline: twfe ----------------------------------------------------------
# panel.data$time_to_treatment <- panel.data$time - panel.data$first.treat
# panel.data$time_to_treatment[panel.data$first.treat==0] <- -1000
# # How can I add continuous pre-treatment covariates X1, X2, X3 in the fixed effect model?
# twfe <- panel.data %>% 
#   do(broom::tidy(feols(ytrain ~ + i(time_to_treatment, ref = c(-1, -1000)) | id + time + X4 + X5, 
#                        data = .), conf.int = TRUE)) %>% 
#   mutate(t =  as.double(str_replace_all(term, c("time_to_treatment::" = "", ":treated" = "")))) %>% 
#   filter(t > -12 & t < 8) %>%
#   select(t, estimate, conf.low, conf.high) %>% 
#   # add in data for year -1
#   bind_rows(tibble(t = -1, estimate = 0, 
#                    conf.low = 0, conf.high = 0
#   )) %>% 
#   mutate(method = "TWFE")

# Baseline: DiD imputation ----------------------------------------------------------
# TODO: Figure out how to fit did imputation to get grouped effect
tm <- proc.time()
panel.data$cohort <- panel.data$first.treat - t0 + 1
wtr <- rep(NA, nrow(cohort.att))
for(i in 1:nrow(cohort.att)){
  c <- cohort.att$cohort[i]
  t <- cohort.att$Time[i] + t0 - 1
  panel.data[, paste("cohort", c, "_time", t, sep = "")] <- (panel.data$time == t) & (panel.data$cohort == c)
  wtr[i] <- paste("cohort", c, "_time", t, sep = "")
}

did_imp.fit <- did_imputation(data = panel.data, yname = "ytrain", gname = "first.treat",
                              tname = "time", idname = "id", wtr = wtr,
                              first_stage = ~ 0 | id + time,
                              horizon=T)

did_imp <- did_imp.fit %>%
  separate(term, into = c("cohort", "Time"), sep = "_") %>% # get cohort and time
  mutate(cohort = gsub("cohort", "", cohort) %>% as.numeric,
         Time = (gsub("time", "", Time) %>% as.numeric) - t0 + 1, 
         ATT = estimate,
         conf.low = estimate + std.error * qnorm(alpha / 2 / nrow(did_imp.fit)), 
         conf.high = estimate + std.error * qnorm(1 - alpha / 2 / nrow(did_imp.fit))) %>% # get 95% simult band
  select(-lhs, -estimate, -std.error)
  
tm <- proc.time() - tm
att.results[nrow(att.results) + 1,] <- c('DiD Imputation', cohort.att.metric(cohort.att, did_imp), as.numeric(tm[3]))


# # Baseline: DiD Imputation with Covariates --------------------------------
# imp.data <- data.frame(
#   outcome = as.vector( t(ytrain) ),
#   treat_status = as.vector( t(ztrain) ),
#   uniq_id = as.vector( sapply(1:n, rep, t1)),
#   panel_time = rep(1:t1, n),
#   move_date = as.vector( sapply(first.treat, rep, t1)),
#   cov1 = as.vector( sapply(xtrain[,1], rep, t1)),
#   cov2 = as.vector( sapply(xtrain[,2], rep, t1)),
#   cov3 = as.vector( sapply(xtrain[,3], rep, t1)),
#   cov4 = as.vector( sapply(xtrain[,4], rep, t1)),
#   cov5 = as.vector( sapply(xtrain[,5], rep, t1))
# )
# imp.data$event_time <- sapply(imp.data$panel_time - imp.data$move_date + 1, function(x) max(x, 0))
# imp.data$event_time[imp.data$move_date == 0] <- 0
# 
# imp.data$X1 <- cut(imp.data$cov1, breaks = 20, labels = 1:20)
# imp.data$X2 <- cut(imp.data$cov2, breaks = 20, labels = 1:20)
# imp.data$X3 <- cut(imp.data$cov3, breaks = 20, labels = 1:20)
# 
# 
# estimate_event_time_effects_weighted<-function(control_obs,treated_obs,weights_vec, include_cov){
#   control_ids<-unique(control_obs$uniq_id)
#   treated_ids<-unique(treated_obs$uniq_id)
#   
#   status_changers<-intersect(control_ids, treated_ids)
#   treated_obs%>%filter(uniq_id %in% status_changers)->imputed_group
#   
#   #Add more interactions of panel_time with covariates here as needed
#   if (include_cov){
#     fixed_effect_object<-fixest::feols(outcome ~ -1 | uniq_id + panel_time^cov4^cov5^X1^X2^X3, data=control_obs, weights = weights_vec)
#   } else{
#     fixed_effect_object<-fixest::feols(outcome ~ -1 | uniq_id + panel_time, data=control_obs, weights = weights_vec)
#   }
#   predicted_outcome<-predict(fixed_effect_object,newdata=imputed_group)
#   imputed_group$counterfactual_outcome<-predicted_outcome
#   imputed_group$tau_est<-imputed_group$outcome-imputed_group$counterfactual_outcome
#   
#   imputed_group%>%
#     group_by(event_time)%>%
#     summarise(period_effect=mean(tau_est,na.rm = TRUE))->event_time_estimates
#   
#   return(list(event=event_time_estimates,individual=imputed_group))
# }
# 
# #Main function to do bootstrap inference for effects
# boot_BJS_effects<-function(dataset,period_number,boot_number, alpha, include_cov){
#   parentids<-unique(dataset$uniq_id)
#   individual_frame<-data.frame()
#   event_frame<-data.frame()
#   control_obs<-dataset[dataset$treat_status==0,]
#   treated_obs<-dataset[dataset$treat_status ==1,]
#   control_obs%>%arrange(uniq_id)->control_obs
#   control_obs%>%select(uniq_id)%>%table()->obs_counts
#   for(i in (1:boot_number)){
#     weight_vector<-rep(rexp(nrow(obs_counts)), obs_counts)
#     BJS_results<-estimate_event_time_effects_weighted(control_obs,treated_obs,weight_vector, include_cov)
#     event_study_estimates<-BJS_results$event
#     event_study_estimates$boot_id<-i
#     event_study_estimates<-dplyr::bind_rows(event_study_estimates,
#                                             data.frame(event_time=9999,period_effect=sum(event_study_estimates$period_effect)))
#     individual_estimates<-BJS_results$individual
#     individual_estimates$boot_id<-i
#     event_frame<-dplyr::bind_rows(event_frame,event_study_estimates)
#     individual_frame<-dplyr::bind_rows(individual_frame,individual_estimates)
#   }
#   
#   event_frame %>%
#     group_by(event_time) %>% 
#     summarise(conf.low=quantile(period_effect,alpha / 2,na.rm = TRUE),
#               estimate = mean(period_effect, na.rm = TRUE),
#               conf.high=quantile(period_effect, 1 - alpha / 2,na.rm = TRUE)) %>%
#     rename(t = event_time) %>%
#     mutate(method = "DiD Imputation Cov")  %>%
#     mutate(t = t - 1) %>%
#     filter(t >= 0 & t < 999) ->event_effects
#   
#   # individual_frame%>%group_by(move_date,event_time)%>%summarise(low=quantile(tau_est,alpha / 2,na.rm = TRUE),
#   #                                                               high=quantile(tau_est,1-alpha/2,na.rm = TRUE))->cohort_per_period_effects
#   # 
#   # individual_frame%>%group_by(move_date)%>%summarise(low=quantile(total_effect,alpha/2,na.rm = TRUE),
#   #                                                    high=quantile(total_effect,1-alpha/2,na.rm = TRUE))->cohort_total_effects
#   # 
#   # results_list=list(
#   #   event_times=event_effects,
#   #   cohort_period=cohort_per_period_effects,
#   #   cohort_total=cohort_total_effects
#   # )
#   return (event_effects)}
# 
# # tm <- proc.time()
# # did_imp_check <- boot_BJS_effects(imp.data, t1, 500, alpha, include_cov = FALSE)
# # tm <- proc.time() - tm
# # did_imp_check %>% mutate(method = "DiD Imputation Check")
# # att.results[nrow(att.results) + 1,] <- c(iter, pr, trt, 'DiD Imputation Check', att.metric(att, did_imp_check), as.numeric(tm[3]))
# 
# tm <- proc.time()
# did_imp_cov <- boot_BJS_effects(imp.data, t1, 500, alpha, include_cov = TRUE)
# tm <- proc.time() - tm
# att.results[nrow(att.results) + 1,] <- c(iter, pr, trt, 'DiD Imputation Cov', att.metric(att, did_imp_cov), as.numeric(tm[3]))
# 
# 
# 
# write.csv(att.results, file = 'att.csv', row.names= F)
# write.csv(catt.results, file = 'catt.csv', row.names= F)
# results -----------------------------------------------------------------
coefs <- cohort.att %>% 
  mutate(conf.low = ATT, conf.high = ATT, method = "True")  %>%
  bind_rows(longbet.att.cohort %>% mutate(method = "LongBet")) %>%
  bind_rows(DiD %>% mutate(method = "DiD")) %>%
  bind_rows(DiD_nl %>% mutate(method = "Non-linear DiD"))

plot <- coefs %>% 
  ggplot(aes(x = Time, y = ATT, color = method)) + 
  geom_point(aes(x = Time, y = ATT), position = position_dodge2(width = 0.8), size = 1) +
  geom_linerange(aes(x = Time, ymin = conf.low, ymax = conf.high), position = position_dodge2(width = 0.8), size = 0.75) +
  # geom_line(data = att.df, aes(x = t, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = .25, alpha = 0.75) + 
  geom_vline(xintercept = -0.5, linetype = "dashed", size = .25) +
  scale_color_manual(name="Estimation Method", values= met.brewer("Cross", 8, "discrete")) +
  theme(legend.position= 'bottom') +
  labs(title = 'Event Time Estimates', y="ATT", x = "Relative Time") + 
  guides(col = guide_legend(nrow = 3)) + 
  facet_wrap(~cohort)
print(plot)


print(att.results)
print(catt.results)
