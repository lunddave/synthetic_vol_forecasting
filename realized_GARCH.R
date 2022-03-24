

library(rugarch)
require(xts)

# https://www.r-bloggers.com/2014/01/the-realized-garch-model/

data(spyreal)

head(spyreal)
tail(spyreal)

#Let's plot the realized vol column
plot(spyreal$SPY_RK)
lines(spyreal$SPY_OC, col = 'green')

spec <- ugarchspec(mean.model = list(armaOrder = c(0, 0), 
                                     include.mean = FALSE), 
                                      variance.model = 
                     list(model = 'realGARCH', garchOrder = c(1, 1)))
# fit <- ugarchfit(spec, out.sample = 25, spyreal[, 1] * 100, 
#                  solver = 'hybrid', 
#                  realizedVol = spyreal[,2] * 100)

fit <- ugarchfit(spec, spyreal[, 1] * 100, 
                 solver = 'hybrid', 
                 realizedVol = spyreal[,2] * 100)

spec
fit

tail(spyreal)
tail(sigma(fit))
tail(fitted(fit))

# https://stats.stackexchange.com/questions/58339/rugarch-understand-forecasting
# https://quant.stackexchange.com/questions/7932/forecasting-using-rugarch-package

fc1 <- ugarchforecast(fit, data = as.matrix(1), 
                     n.ahead = 5, 
                     n.roll = 10, 
                     #out.sample = 2
                     )

fc1 <- ugarchforecast(fit, data = as.matrix(1), 
                      n.ahead = 5
)

# What is out.sample?
# Optional. If a specification object is supplied, indicates how many data points 
# to keep for out of sample testing.

sigma(fc1)
fitted(fc1)

# I most likely will not need much of this sophisticated functionality.
# I will need only..

fit <- ugarchfit(spec, spyreal[, 1] * 100, 
                 solver = 'hybrid', 
                 realizedVol = spyreal[,2] * 100)
fc2 <- ugarchforecast(fit, data = as.matrix(10,12,13), 
                      n.ahead = 5)
fc2
fc2@forecast$realizedFor
plot.ts(fc2@forecast$realizedFor)
fc2@forecast$sigmaFor
fc2@forecast$sigmaDF
sigma(fc2)
fitted(fc2)
