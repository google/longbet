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
t0 <- 6 # treatment start time
t1 <- 12 # observed response period
alpha <- 0.05

source('dgp.R')
pr_type = "non-parallel"
trt_type = "heterogeneous"
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
  method = rep("True", t1 - t0 + 1)
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
longbet.fit <- longbet(y = ytrain, x = xtrain, z = ztrain, t = 1:t1,
                       num_sweeps = 100,
                       num_trees_pr =  20, num_trees_trt = 20,
                       pcat = ncol(xtrain) - 3)

longbet.pred <- predict.longbet(longbet.fit, xtrain, ztrain)
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

att.results[nrow(att.results) + 1,] <- c('LongBet', att.metric(att, longbet.att), as.numeric(longbet.time[3]))
catt.results[nrow(catt.results) + 1, ] <- c('LongBet', catt.metric(align_tau, longbet.catt, longbet.catt.low, longbet.catt.high))


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
att.results[nrow(att.results) + 1,] <- c('DiD', att.metric(att, DiD), as.numeric(did.time[3]))

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
  mutate(method = "Non-linear DiD")
did_nl.time = proc.time() - did_nl.time
att.results[nrow(att.results) + 1,] <- c('Non-linear DiD', att.metric(att, DiD_nl), as.numeric(did_nl.time[3]))


# Baseline: twfe ----------------------------------------------------------
panel.data$time_to_treatment <- panel.data$time - panel.data$first.treat
panel.data$time_to_treatment[panel.data$first.treat==0] <- -1000
# How can I add continuous pre-treatment covariates X1, X2, X3 in the fixed effect model?
twfe <- panel.data %>% 
  do(broom::tidy(feols(ytrain ~ + i(time_to_treatment, ref = c(-1, -1000)) | id + time + X4 + X5, 
                       data = .), conf.int = TRUE)) %>% 
  mutate(t =  as.double(str_replace_all(term, c("time_to_treatment::" = "", ":treated" = "")))) %>% 
  filter(t > -12 & t < 8) %>%
  select(t, estimate, conf.low, conf.high) %>% 
  # add in data for year -1
  bind_rows(tibble(t = -1, estimate = 0, 
                   conf.low = 0, conf.high = 0
  )) %>% 
  mutate(method = "TWFE")


# results -----------------------------------------------------------------
# coefs <- bind_rows(longbet.att, DiD, twfe)
coefs <- bind_rows(att.df, longbet.att, DiD, DiD_nl)

plot <- coefs %>% 
  ggplot(aes(x = t, y = estimate, color = method)) + 
  geom_point(aes(x = t, y = estimate), position = position_dodge2(width = 0.8), size = 1) +
  geom_linerange(aes(x = t, ymin = conf.low, ymax = conf.high), position = position_dodge2(width = 0.8), size = 0.75) +
  # geom_line(data = att.df, aes(x = t, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = .25, alpha = 0.75) + 
  geom_vline(xintercept = -0.5, linetype = "dashed", size = .25) +
  scale_color_manual(name="Estimation Method", values= met.brewer("Cross", 8, "discrete")) +
  theme(legend.position= 'bottom') +
  labs(title = 'Event Time Estimates', y="ATT", x = "Relative Time") + 
  guides(col = guide_legend(nrow = 3)) 
print(plot)


print(att.results)
print(catt.results)