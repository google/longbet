# simple demonstration of longbet with default parameters
require(dplyr)
require(ggplot2)
require(tidyr)
require(longbet)
require(panelView)
require(stringr)
require(MetBrewer)

# DATA GENERATION PROCESS -------------------------------------------------
set.seed(3)
n <- 2000
t0 <- 7 # treatment start time
t1 <- 12 # observed response period
t2 <- 16 # extended treatment effect
alpha <- 0.05

source('dgp.R')
pr_type = "non-parallel"
trt_type = "homogeneous"
staggered_effect <- c(1.2, 1.2, 1, 0.9, 0.8, 0.8)
data <- dgp(n, t0, t1, t2, pr_type = pr_type, trt_type = trt_type, staggered_effect = staggered_effect)

# get training data
ytrain <- data$y[, 1:t1]
ztrain <- data$z[, 1:t1]
xtrain <- data$x

ypred <- data$y[, 1:t2]
zpred <- data$z[, 1:t2]

cohort <- getCohort(zpred) # get treatment group
cohort.att <- getCohortAtt(data$tau[,t0:t2], cohort)
xmod <- cbind(xtrain, cohort) # covariates for treatment arm

# longbet -----------------------------------------------------------------
longbet.time <- proc.time()
longbet.fit <- longbet(y = ytrain, x = xtrain, x_trt = xmod, z = ztrain, t = 1:t1,
                       num_sweeps = 100,
                       num_trees_pr =  40, num_trees_trt = 40,
                       pcat = ncol(xtrain) - 3)

longbet.ext <- predict.longbet(longbet.fit, xtrain, xmod, zpred)
longbet.time <- proc.time() - longbet.time

# get ATT by cohort
longbet.att.hat <- getCohortAttHat(longbet.ext$tauhats[,t0:t2,], cohort)
longbet.att.ci <- getCohortAttCI(longbet.ext$tauhats[,t0:t2,], cohort, alpha)

# calculate metrics for each group
longbet.metrics <- getCohortMetrics(cohort.att, longbet.att.hat, longbet.att.ci)
print(round(longbet.metrics, 3))

# results -----------------------------------------------------------------
colors <- c("black", "#FFA500", "#00BFFF")
labels <- c("True", "LongBet", "BART")
# names(colors) <- labels
# Visualize results by cohort

cohort.att.combined <- cohort.att %>%
  gather(key = "Time", value = "True", -cohort) %>%
  mutate(Time = as.numeric(gsub("X", "", Time)) + t0 - 1) %>%
  left_join(longbet.att.hat %>%
           gather(key = "Time", value = "LongBet", -cohort) %>%
           mutate(Time = as.numeric(gsub("X", "", Time)) + t0 - 1),
           by = c("cohort", "Time")) %>%
  left_join(longbet.att.ci$lower %>%
              gather(key = "Time", value = "Lower", -cohort) %>%
              mutate(Time = as.numeric(gsub("X", "", Time)) + t0 - 1),
            by = c("cohort", "Time")) %>%
  left_join(longbet.att.ci$upper %>%
            gather(key = "Time", value = "Upper", -cohort) %>%
            mutate(Time = as.numeric(gsub("X", "", Time)) + t0 - 1),
          by = c("cohort", "Time"))

cohort.att.plot <- cohort.att.combined %>%
  filter(cohort > 0) %>%
  filter(Time - t0 + 1 >= cohort) %>%
  ggplot(aes(x = Time)) +
  geom_line(aes(y = True, color = "True")) +  # Map color to legend
  geom_line(aes(y = LongBet, color = "LongBet")) +  # Map color to legend
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = "LongBet"), alpha = 0.15, fill = colors[2]) +
  geom_vline(xintercept = t1, linetype = "dashed", color = "grey") +
  labs(x = "Time", y = "ATT per cohort", color = "Legend") +
  scale_color_manual(
    values = c("True" = colors[1], "LongBet" = colors[2])
  ) +
  facet_wrap(~cohort)

print(cohort.att.plot)
# readline(prompt="Press [enter] to continue")
