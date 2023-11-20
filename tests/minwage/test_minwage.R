# compare longbet to did method on the minimum wage dataset.
# estimate the effect of minimum wage on teen employment
# did paper: https://reader.elsevier.com/reader/sd/pii/S0304407620303948?token=967C9A6B23A76EC6C271F99CC54C0CAC8829014C025C3AB9F600C1B5F39F0A7EB6CD905A583911978C1025A798AD05D6&originRegion=us-east-1&originCreation=20230425183217
# required package: did https://github.com/bcallaway11/did

library(longbet)
library(XBART)
library(did)
library(tidyr)
library(dplyr)
data(mpdta)

income <- read.csv('income.csv') # https://www.kaggle.com/datasets/thedevastator/2013-irs-us-income-data-by-zip-code
mpdta <- merge(mpdta, income, by = "countyreal", all.x = T)
# fill missing value
mpdta$income[is.na(mpdta$income)] <- mean(mpdta$income, na.rm = T)


alpha <- 0.05
gp_year <- c(2004, 2006, 2007)
longbet_att <- data.frame(
  group = numeric(),
  t = numeric(),
  att = numeric(),
  upper = numeric(),
  lower = numeric()
)

# longbet staggered adoption ----------------------------------------------
data <- mpdta %>%  spread(key = "year", value = "lemp")

xtrain <- as.matrix(data[, c("lpop", "income")])
# xtrain <- cbind(xtrain, as.numeric(data$first.treat))
ytrain <- as.matrix(data[, c("2003", "2004", "2005", "2006", "2007")])
get_z <- function(first.treat){
  if (first.treat == 0) {
    return(rep(0, 5))
  } else {
      return(first.treat <= c(2003:2007))}
}
ztrain <- sapply(data$first.treat, get_z) %>% t 

# get propensity score
yclass <- factor(data$first.treat, labels = c(0, 1, 2, 3)) %>% as.numeric()
yclass <- yclass - 1
fit.ps <- XBART.multinomial(y = matrix(yclass), num_class = 4, X = xtrain,
                          # num_trees = 20, num_sweeps = 100, burnin = 20,
                          p_categorical = 0)
ps.hat <- predict.XBARTmultinomial(fit.ps, xtrain)
xtrain <- cbind(ps.hat$prob[,2:4], xtrain)

# add group in treatment covariates
x_trt <- cbind(xtrain, data$first.treat)

longbet.fit <- longbet(y = ytrain, x = xtrain, x_trt = x_trt, z = ztrain, t = 1:ncol(ztrain),
                       num_sweeps = 200, num_burnin = 80, 
                       num_trees_pr =  20, num_trees_trt = 10,
                       pcat = 0, lambda_knl = 1,
                       split_time_ps = T, split_time_trt = T)

longbet.pred <- predict.longbet(longbet.fit, xtrain, x_trt, ztrain)
longbet.att <- get_att(longbet.pred, alpha = 0.05)
longbet.catt <- get_catt(longbet.pred, alpha = 0.05)

# reshape tauhats to get att credible interval per group
n <- dim(longbet.pred$tauhats)[1]
t <- dim(longbet.pred$tauhats)[2]
num_sweeps <- dim(longbet.pred$tauhats)[3]

tauhats <- longbet.pred$tauhats %>% 
  aperm(perm = c(1,3,2)) %>% 
  matrix(nrow = n * num_sweeps, ncol = t, byrow = F) %>%
  data.frame
colnames(tauhats) <- c(2003:2007)
tauhats$group <- rep(data$first.treat, num_sweeps)
tauhats$sweeps <- sapply(1:num_sweeps, rep, times = n) %>% as.vector

staggered_att <- 
  tauhats %>% 
  filter(group != 0) %>%
  gather(key = "t", value = "catt", -group, -sweeps) %>%
  group_by(group, t, sweeps) %>%
  summarise(catt = mean(catt)) %>%
  ungroup() %>%
  group_by(group, t) %>%
  summarise(att = mean(catt), upper = quantile(catt, 0.975), lower = quantile(catt, 0.025)) %>%
  filter(t >= group)
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
# longbet_att$method <- rep("LongBet", nrow(longbet_att))
staggered_att$method <- rep("Staggered", nrow = (staggered_att))
did_att$method <- rep("DiD", nrow(did_att))


att_plot <- rbind(did_att, staggered_att) %>% #longbet_att
  ggplot(aes(x = t, y = att)) +
  geom_point(aes(color = method), position = position_dodge(.3)) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = method), width = .2,  position = position_dodge(.3)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  facet_wrap(~group, ncol = 1) +
  labs(y = "ATT")
plot(att_plot)

# Plot estimated muhat and tauhat of each group, compare to groun --------
muhat <- longbet.pred$muhats %>% apply(MARGIN = c(1, 2), mean) %>% data.frame
tauhat <- longbet.pred$tauhats %>% apply(MARGIN = c(1, 2), mean) %>% data.frame
yhat <- muhat + tauhat
colnames(muhat) <- c(2003:2007)
colnames(yhat) <- c(2003:2007)
muhat$group <-data$first.treat %>% as.factor
yhat$group <- data$first.treat %>% as.factor

mu_df <- muhat %>% 
  gather(key = "year", value = "lemp", -group) %>%
  group_by(group, year) %>% summarise(lemp = mean(lemp))
mu_df$year <- as.numeric(mu_df$year)

yhat_df <- yhat %>% 
  gather(key = "year", value = "lemp", -group) %>%
  group_by(group, year) %>% summarise(lemp = mean(lemp))
yhat_df$year <- as.numeric(mu_df$year)

ground_truth <- data %>%
  select('first.treat', '2003', '2004', '2005', '2006', '2007') %>%
  gather(key = "year", value = "lemp", -first.treat) %>%
  group_by(first.treat, year) %>%
  summarise(lemp = mean(lemp))
ground_truth$group <- as.factor(ground_truth$first.treat)
ground_truth$year <- as.numeric(ground_truth$year)
ground_truth$first.treat <- NULL

ground_truth$label <- "true"
mu_df$label <- "y0hat"
yhat_df$label <- "y1hat"

yhat_plot <- rbind(ground_truth, mu_df, yhat_df) %>% 
  ggplot(aes(x = year, y = lemp, color = label)) +
  geom_point() +  geom_line() + 
  facet_wrap(~group)
plot(yhat_plot)

# check counfoundings by varifying impute accuracy on y0_hat --------------
y0_hat <- longbet.pred$muhats %>% apply(MARGIN = c(1, 2), mean) %>% data.frame
for (gp in 0:3){
  gp_z <- ztrain[yclass == gp, ]
  rmse <- mean( (ytrain[yclass == gp, ][gp_z == 0] - y0_hat[yclass == gp, ][gp_z == 0])^2) %>% sqrt
  print(paste0("Group ", gp, " RMSE for y0 = ", round(rmse, 3), sep =""))
}

# check observed confoundings explain the observed outcomes ---------------
check_df <- data.frame(x_trt)
check_df$group <- data$treat %>% as.factor()
check_df$lemp <- ytrain[,1]

check_lpop <- check_df %>% 
  ggplot(aes(x = lpop, y = lemp, group = group, color = group)) +
  geom_point()

check_income <- check_df %>% 
  ggplot(aes(x = income, y = lemp, group = group, color = group)) +
  geom_point()

combined_check <- gridExtra::grid.arrange(check_lpop, check_income, ncol = 2)  # You can specify the layout
print(combined_check)


# check att per sweep -----------------------------------------------------
att_sweeps <- 
  tauhats %>% 
  filter(group != 0) %>%
  gather(key = "t", value = "catt", -group, -sweeps) %>%
  group_by(group, t, sweeps) %>%
  summarise(catt = mean(catt)) %>%
  ungroup() %>%
  filter(t >= group)

sweeps_plt <- att_sweeps %>%
  ggplot(aes(x = sweeps, y = catt)) +
  geom_point() + geom_line() +
  facet_wrap(~group + t, nrow = 2)
plot(sweeps_plt)

# conditional treatment effect --------------------------------------------
# feed treatment effect in a CART tree
# for simplicity feed the the treatment effect in 2007
library(rpart)
library(rpart.plot)
tree.x <- data.frame(
  lpop = data$lpop,
  income = data$income,
  tauhat = longbet.catt$catt[,5]
)
fit.tree <- rpart(tauhat ~ ., data=tree.x, cp=0.008)
rpart.plot(fit.tree)

library(partykit)
tree_party <- as.party(fit.tree) #coercing the rpart.object into partykit
tree_fit <- fitted(tree_party) #a two-column data.frame with the fitted node numbers and the observed responses on the training data.
subgroups <- predict(tree_party, newdata = tree.x, type = "node")

tauhats$subgroups <- rep(subgroups, num_sweeps)

catt <- 
  tauhats %>% 
  filter(group != 0) %>%
  gather(key = "t", value = "catt", -group, -sweeps, -subgroups) %>%
  group_by(group, t, sweeps, subgroups) %>%
  summarise(catt = mean(catt)) %>%
  ungroup() %>%
  group_by(group, t, subgroups) %>%
  summarise(att = mean(catt), upper = quantile(catt, 0.975), lower = quantile(catt, 0.025)) %>%
  filter(t >= group)

catt_plot <- catt %>%
  ggplot(aes(x = t, y = att)) +
  geom_point(position = position_dodge(.3)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .2,  position = position_dodge(.3)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  facet_wrap(~group + subgroups, nrow = 3) +
  labs(y = "CATT")
plot(catt_plot)
