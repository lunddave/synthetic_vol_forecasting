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

## Now using data gleaned from Wharton (WRDS) in December 2023
wharton <- read.csv('/home/david/Desktop/synthetic_vol_forecasting/implied_vol_datasets/uw5rmndubuhtpvbt.csv')
names(wharton)
subset <- wharton[(wharton$days == 30) & wharton$cp_flag == 'C',]
subset$date <- as.Date(subset$date)
subset$impl_volatility <- as.numeric(subset$impl_volatility) * 100

IYG <- read.csv('/home/david/Desktop/synthetic_vol_forecasting/implied_vol_datasets/IYG.csv')
IYG$log_ret <- c(NA,diff(log(IYG$Adj.Close)) * 100)

merged <- merge(subset, IYG, by.x = 'date', by.y = 'Date')

merged$exog_var <- rnorm(nrow(merged))

subset_xts <- xts(merged, order.by = merged$date)

subset_xts <- subset_xts[,c('impl_volatility','log_ret','exog_var')]

storage.mode(subset_xts) <- "numeric"

#Let's plot the realized vol column
plot(subset_xts$impl_volatility)
lines(subset_xts$log_ret, col = 'green')

spec <- ugarchspec(mean.model = list(armaOrder = c(0, 0),
                                     include.mean = FALSE),
                   variance.model = list(model = 'realGARCH'
                                         ,garchOrder = c(1, 1)
                                         ,external.regressors = subset_xts$exog_var

                                         ))

fit <- ugarchfit(spec
                 ,subset_xts$log_ret
                 ,solver = 'hybrid'
                 ,realizedVol = subset_xts$impl_volatility * 100
                 )


fit

tail(sigma(fit))

ni = newsimpact(fit, z = seq(-2, 2, length.out = 100))
plot(ni$zx, (ni$zy), ylab = ni$yexpr, xlab = ni$xexpr, type = 'l', main = 'News Impact realGARCH')
abline(v = 0)
abline(h = 0)
grid()

forc1 <- ugarchforecast(fit
                        , n.ahead = 5
                        , n.roll = 0
                        , externalforecasts = list(vregfor = (5)))
#https://stats.stackexchange.com/questions/566840/how-do-i-forecast-with-external-regressors-in-the-rugarch-package-in-r-the-regr
forc1
