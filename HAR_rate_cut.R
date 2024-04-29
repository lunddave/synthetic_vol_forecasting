library(highfrequency)
library(tidyverse)

QL <- function(adj, unadj) {adj/unadj - log(adj/unadj) - 1}



RVSPY <- as.xts(SPYRM$RV1, order.by = SPYRM$DT)

foo = data.frame(
  Date = SPYRM$DT,
  RV1 = SPYRM$RV1
)

head(foo)

foo %>% arrange(desc(RV1))


## https://www.forbes.com/advisor/investing/fed-funds-rate-history/

Fed_rate_cuts =
  c("2018-12-20", "2018-09-27", "2018-06-14",
    "2018-03-22", "2017-12-14", "2017-06-15",
    "2017-03-16", "2016-12-15", "2015-12-17")

Fed_rate_cuts =
  c(#"2018-12-20",
    "2018-09-27", "2018-06-14",
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
  mutate(Y = lead(RV1)) %>%
  mutate(RV1 = RV1) %>%
  mutate(RV5 = roll_mean(RV1, 5, align = "right", fill = NA)) %>%
  mutate(RV5_22 = roll_mean(RV1, 22, align = "right", fill = NA))

baz = bar[complete.cases(bar), ]
head(baz)

qux = data.frame(baz %>% mutate(rate_cut = ifelse(Date %in% Fed_rate_cuts,1,0)))

qux$log_Y <- log(qux$Y)

qux[qux$rate_cut == 1, "rate_cut"] = 1:8

qux$rate_cut = as.factor(qux$rate_cut)

m1 = lm(Y ~ RV1 + RV5 + RV5_22 + rate_cut + 0, data = qux)
summary(m1)
plot(m1)

log_lin = lm(log_Y ~ RV1 + RV5 + RV5_22 + rate_cut + 0, data = qux)
summary(log_lin)
plot(log_lin)

rate_cut_9 = qux[qux$Date == as.Date("2018-12-20") - 1, ]
rate_cut_9

newdat = data.frame(RV5 = 0.0004136769,
                    RV5_5 = 0.0002954981,
                    RV5_22 = 0.0001593874,
                    rate_cut = 0)

newdat <- rate_cut_9[,c(2,4,5,6)]

newdat$rate_cut = as.factor(newdat$rate_cut)

rate_cut_9_pred = predict(m1, newdata = newdat)

rate_cut_9$Y - rate_cut_9_pred

alpha_mean = mean(coef(m1)[5:12]) #use arithmetic mean

resid_adjusted <- rate_cut_9$Y - (rate_cut_9_pred + alpha_mean)
resid_adjusted_SE <- resid_adjusted**2

resid_unadjusted <- rate_cut_9$Y - rate_cut_9_pred
resid_unadjusted_SE <- resid_unadjusted**2

resid_unadjusted_SE > resid_adjusted_SE
resid_unadjusted_SE - resid_adjusted_SE


QL_loss_adjusted <- QL(rate_cut_9_pred + mean(coef(m1)[5:12]), rate_cut_9$Y)
QL_loss_unadjusted <- QL(rate_cut_9_pred, rate_cut_9$Y)
QL_loss_unadjusted > QL_loss_adjusted

# We exponentiate predictions of log values
QL_loss_adjusted <- QL(exp(rate_cut_9_pred + mean(coef(m1)[5:12])), exp(rate_cut_9$Y))
QL_loss_unadjusted <- QL(exp(rate_cut_9_pred), exp(rate_cut_9$Y))
QL_loss_unadjusted > QL_loss_adjusted

