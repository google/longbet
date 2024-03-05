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

x <- mpdta %>%
  select(c("countyreal","lpop", "first.treat", "income")) %>%
  unique()

# check overlap assumption
x %>% mutate(first.treat = as.factor(first.treat)) %>%
  ggplot() + 
  geom_violin(aes(first.treat, lpop, group = first.treat))


x %>% mutate(first.treat = as.factor(first.treat)) %>%
  ggplot() + 
  geom_violin(aes(first.treat, income, group = first.treat))

