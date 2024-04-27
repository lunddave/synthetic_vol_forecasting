
library(highfrequency)
library(xts)

RVSPY <- as.xts(SPYRM$RV5, order.by = SPYRM$DT)

day_of <- ifelse( index(RVSPY) == '2016-11-07' | index(RVSPY) == '2016-11-08' , 1, 0)

day_after <- ifelse(index(RVSPY) == '2016-11-09' | index(RVSPY) == '2016-11-10'  , 1, 0)

names(RVSPY) <- c('ret')

RVSPY$day_of <- day_of
RVSPY$day_after <- day_after


head(RVSPY)

x <- HARmodel(data = RVSPY , periods = c(1,5,22), RVest = c("rCov"),
              type = "HAR", h = 1, transform = NULL, inputType = "RM"
              ,externalRegressor = RVSPY$day_of
              )

class(x)
x
summary(x)
plot(x)
predict(x)


#Dan's work from week beginning April 15th, 2024

library(highfrequency)
library(tidyverse)

# Example 1: HAR
# Forecasting daily Realized volatility for the
# S&P 500 using the basic HARmodel: HAR
library(xts)
RVSPY <- as.xts(SPYRM$RV5, order.by = SPYRM$DT)
x <- HARmodel(data = RVSPY,
              periods = c(1,5,22),
              RVest = c("rCov"),
              type = "HAR",
              h = 1,
              transform = NULL,
              inputType = "RM")
class(x)

x
summary(x)
plot(x)
predict(x)

foo = as.data.frame(RVSPY)
head(foo)

RVSPY <- as.xts(SPYRM$RV5, order.by = SPYRM$DT) #WHEN? "2014-01-02" "2019-12-31"

foo = data.frame(
  Date = SPYRM$DT,
  RV5 = SPYRM$RV5
)
head(foo)
foo %>% arrange(desc(RV5))


## https://www.forbes.com/advisor/investing/fed-funds-rate-history/

Fed_rate_cuts =
  c(#"2019-10-31", "2019-09-19", "2019-08-01",
    "2018-12-20", "2018-09-27", "2018-06-14",
    "2018-03-22", "2017-12-14", "2017-06-15",
    "2017-03-16", "2016-12-15", "2015-12-17")

foo %>% filter(Date %in% Fed_rate_cuts)

#December 20, 2018	+25	2.25% to 2.50%
#  Sept. 27, 2018	+25	2.0% to 2.25%
#  Jun. 14, 2018	+25	1.75% to 2.0%
#  March 22, 2018	+25	1.50% to 1.75%
#  Dec. 14, 2017	+25	1.25% to 1.50%
#  June 15, 2017	+25	1.00% to 1.25%
#  March 16, 2017	+25	0.75% to 1.00%
#  Dec. 15, 2016	+25	0.5% to 0.75%
#  Dec. 17, 2015	+25	0.25% to 0.50%


library(RcppRoll)
bar = foo %>%
  mutate(RV5_5 = roll_sum(RV5, 5, align = "right", fill = NA)/5) %>%
  mutate(RV5_22 = roll_sum(RV5, 22, align = "right", fill = NA)/22)
baz = bar[complete.cases(bar), ]
head(baz)

qux = data.frame(
  Y = baz$RV5[-nrow(baz)],
  baz[-1, ]) %>%
  mutate(rate_cut = ifelse(Date %in% Fed_rate_cuts,1,0))
qux[qux$rate_cut == 1, "rate_cut"] = 1:9
qux$rate_cut = as.factor(qux$rate_cut)

m1 = lm(Y ~ RV5 + RV5_5 + RV5_22 + rate_cut, data = qux)
summary(m1)

rate_cut_9 = qux[qux$Date == "2018-12-20", ]
rate_cut_9

newdat = data.frame(RV5 = 0.0004136769,
                    RV5_5 = 0.0002954981,
                    RV5_22 = 0.0001593874,
                    rate_cut = 0)
newdat$rate_cut = as.factor(newdat$rate_cut)

rate_cut_9_pred = predict(m1, newdata = newdat)

rate_cut_9$Y - rate_cut_9_pred

alpha_mean = mean(coef(m1)[5:12])
rate_cut_9$Y - (rate_cut_9_pred + mean(coef(m1)[5:12]))


#MY RIFF Dan's work from week beginning April 15th, 2024

# Example 1: HAR
# Forecasting daily Realized volatility for the
# S&P 500 using the basic HARmodel: HAR
library(xts)
RVSPY <- as.xts(SPYRM$RV5, order.by = SPYRM$DT)
x <- HARmodel(data = RVSPY,
              periods = c(1,5,22),
              RVest = c("rCov"),
              type = "HAR",
              h = 1,
              transform = NULL,
              inputType = "RM")
class(x)

x
summary(x)
plot(x)
predict(x)

foo = as.data.frame(RVSPY)
head(foo)

RVSPY <- as.xts(SPYRM$RV5, order.by = SPYRM$DT) #WHEN? "2014-01-02" "2019-12-31"

foo = data.frame(
  Date = SPYRM$DT,
  RV5 = SPYRM$RV5
)
head(foo)
foo %>% arrange(desc(RV5))


## https://www.forbes.com/advisor/investing/fed-funds-rate-history/

Fed_rate_cuts =
  c("2014-11-05", "2016-06-23", "2016-11-09", "2018-11-07")

foo %>% filter(Date %in% Fed_rate_cuts)



library(RcppRoll)
bar = foo %>%
  mutate(RV5_5 = roll_sum(RV5, 5, align = "right", fill = NA)/5) %>%
  mutate(RV5_22 = roll_sum(RV5, 22, align = "right", fill = NA)/22)
baz = bar[complete.cases(bar), ]
head(baz)

qux = data.frame(
  Y = baz$RV5[-nrow(baz)],
  baz[-1, ]) %>%
  mutate(rate_cut = ifelse(Date %in% Fed_rate_cuts,1,0))
qux[qux$rate_cut == 1, "rate_cut"] = 1:4
qux$rate_cut = as.factor(qux$rate_cut)

wo = lm(Y ~ RV5 + RV5_5 + RV5_22, data = qux)
summary(wo)

m1 = lm(Y ~ RV5 + RV5_5 + RV5_22 + rate_cut, data = qux)
summary(m1)

rate_cut_9 = qux[qux$Date == "2018-12-20", ]
rate_cut_9

newdat = data.frame(RV5 = 0.0004136769,
                    RV5_5 = 0.0002954981,
                    RV5_22 = 0.0001593874,
                    rate_cut = 0)
newdat$rate_cut = as.factor(newdat$rate_cut)

rate_cut_9_pred = predict(m1, newdata = newdat)

rate_cut_9$Y - rate_cut_9_pred

alpha_mean = mean(coef(m1)[5:12])
rate_cut_9$Y - (rate_cut_9_pred + mean(coef(m1)[5:12]))

#MY SECOND RIFF Dan's work from week beginning April 15th, 2024

# Example 1: HAR
# Forecasting daily Realized volatility for the
# S&P 500 using the basic HARmodel: HAR
library(xts)
RVSPY <- as.xts(SPYRM$RV5, order.by = SPYRM$DT)
x <- HARmodel(data = RVSPY,
              periods = c(1,5,22),
              RVest = c("rCov"),
              type = "HAR",
              h = 1,
              transform = NULL,
              inputType = "RM")
class(x)

x
summary(x)
plot(x)
predict(x)

foo = as.data.frame(RVSPY)
head(foo)

RVSPY <- as.xts(SPYRM$RV5, order.by = SPYRM$DT) #WHEN? "2014-01-02" "2019-12-31"

foo = data.frame(
  Date = SPYRM$DT,
  RV5 = SPYRM$RV5
)
head(foo)
foo %>% arrange(desc(RV5))


## https://www.forbes.com/advisor/investing/fed-funds-rate-history/

Fed_rate_cuts =
  c("2014-11-27"
    , "2015-07-14"
    , "2014-11-05"
    , "2016-06-23"
    , "2016-11-09"
    , "2018-11-07")


foo %>% filter(Date %in% Fed_rate_cuts)



library(RcppRoll)
bar = foo %>%
  mutate(RV5_5 = roll_sum(RV5, 5, align = "right", fill = NA)/5) %>%
  mutate(RV5_22 = roll_sum(RV5, 22, align = "right", fill = NA)/22)
baz = bar[complete.cases(bar), ]
head(baz)

qux = data.frame(
  Y = baz$RV5[-nrow(baz)],
  baz[-1, ]) %>%
  mutate(rate_cut = ifelse(Date %in% Fed_rate_cuts,1,0))
qux[qux$rate_cut == 1, "rate_cut"] = 1:5
qux$rate_cut = as.factor(qux$rate_cut)

wo = lm(Y ~ RV5 + RV5_5 + RV5_22, data = qux)
summary(wo)

m1 = lm(Y ~ RV5 + RV5_5 + RV5_22 + rate_cut, data = qux)
summary(m1)

rate_cut_9 = qux[qux$Date == "2018-12-20", ]
rate_cut_9

newdat = data.frame(RV5 = 0.0004136769,
                    RV5_5 = 0.0002954981,
                    RV5_22 = 0.0001593874,
                    rate_cut = 0)
newdat$rate_cut = as.factor(newdat$rate_cut)

rate_cut_9_pred = predict(m1, newdata = newdat)

rate_cut_9$Y - rate_cut_9_pred

alpha_mean = mean(coef(m1)[5:12])
rate_cut_9$Y - (rate_cut_9_pred + mean(coef(m1)[5:12]))

