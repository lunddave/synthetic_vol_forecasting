library(highfrequency)
library(tidyverse)
library(xts)

QL <- function(pred, vol) {pred/vol - log(pred/vol) - 1}

RVSPY <- as.xts(SPYRM$RV1, order.by = SPYRM$DT)

foo = data.frame(
  Date = SPYRM$DT,
  RV1 = SPYRM$RV1
)

head(foo)

## https://www.forbes.com/advisor/investing/fed-funds-rate-history/

Fed_rate_cuts =
  c(#"2018-12-20",
    "2018-09-27", "2018-06-14",
    "2018-03-22", "2017-12-14", "2017-06-15",
    "2017-03-16", "2016-12-15", "2015-12-17")

foo %>% filter(Date %in% Fed_rate_cuts)

library(RcppRoll)

bar = foo %>%
  mutate(tomorrow_RV1 = lead(RV1)) %>%
  mutate(RV1 = RV1) %>%
  mutate(RV5 = roll_mean(RV1, 5, align = "right", fill = NA)) %>%
  mutate(RV5_22 = roll_mean(RV1, 22, align = "right", fill = NA))

baz = bar[complete.cases(bar), ]
head(baz)

qux = data.frame(baz %>% mutate(rate_cut = ifelse(Date %in% Fed_rate_cuts,1,0)))

qux$log_tomorrow_RV1<- log(qux$tomorrow_RV1)

qux[qux$rate_cut == 1, "rate_cut"] = Fed_rate_cuts

qux$rate_cut = as.factor(qux$rate_cut)

m1 = lm(tomorrow_RV1 ~ RV1 + RV5 + RV5_22 + rate_cut, data = qux)
summary(m1)
plot(m1)

forecast_period = qux[qux$Date == as.Date("2018-12-20") - 1, ]

newdat <- forecast_period[,c(2,4,5,6)]

newdat$rate_cut = as.factor(newdat$rate_cut)
newdat

pred = predict(m1, newdata = newdat)

no_events <- length(Fed_rate_cuts)
no_coef <- length(coef(m1))

FE_mean <- mean(coef(m1)[(no_coef-no_events+1):no_coef])

adjusted_pred <- pred + FE_mean

QL_loss_adjusted <- QL(adjusted_pred, forecast_period$tomorrow_RV1)
QL_loss_unadjusted <- QL(pred, forecast_period$tomorrow_RV1)
QL_loss_unadjusted > QL_loss_adjusted

#log

log_lin = lm(log_tomorrow_RV1 ~ RV1 + RV5 + RV5_22 + rate_cut, data = qux)
summary(log_lin)
plot(log_lin)

forecast_period = qux[qux$Date == as.Date("2018-12-20") - 1, ]

pred = predict(log_lin, newdata = newdat)

FE_mean <- mean(coef(log_lin)[(no_coef-no_events+1):no_coef])

adjusted_pred <- exp(pred + FE_mean)

QL_loss_adjusted <- QL(adjusted_pred, forecast_period$tomorrow_RV1)
QL_loss_unadjusted <- QL(exp(pred), forecast_period$tomorrow_RV1)
QL_loss_unadjusted > QL_loss_adjusted
