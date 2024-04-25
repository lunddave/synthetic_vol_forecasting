
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
