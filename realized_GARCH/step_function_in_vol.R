require(garchx)

set.seed(2)

# Rig the data the way we want it
n <- 200
step_function_series_low_vol <- rnorm(n/2, sd = .1)
step_function_series_high_vol <- rnorm(n/2, sd = .2)
series <- c(step_function_series_low_vol,step_function_series_high_vol)
plot.ts(series)

indicator_vec_low_vol <- rep(0,n/2)
indicator_vec_high_vol <- rep(1,n/2)
indicator_series <- c(indicator_vec_low_vol, indicator_vec_high_vol)

#Fit GARCH(1,1) models and plot the fitted sigma values
par(mfrow = c(1,2))

mod <- garchx(series, xreg = indicator_series)
mod
plot.ts(mod$fitted, ylab = 'Sigma^2'
        , main = 'An Indicator Variable at the Changepoint\n Allows Quick Reaction')
abline(v = n/2, col = 'green')

mod_wo_indicator <- garchx(series)
mod_wo_indicator
plot.ts(mod_wo_indicator$fitted, ylab = 'Sigma^2'
        ,main = 'No Indicator Variable Results\n in Delayed Reaction\n and Erratic Estimation')
abline(v = n/2, col = 'red')
