options(digits = 7)

Tee <- 96
n <- 6
shock_time_vec <- rep(Tee/2, n)
Y <- matrix(rep(rnorm(Tee), n), ncol = n)


X <- list()

#we make X into a list
for (i in 1:n){
  X[[i]] <- cbind(rnorm(Tee)
             ,rnorm(Tee)
             ,rnorm(Tee))
}

temp <- SynthVolForecast(Y
                         ,X
                         ,shock_time_vec
                         ,rep(1, n)
                         ,garch_order = c(1,1)
                         ,plots = TRUE)
temp$convex_combination
temp$forecast

## Example

library(quantmod)
library(garchx)
library(lmtest)

# time series of interest
getSymbols('COP',src='yahoo',from = '2020-01-08' , to = '2020-03-10')

# covariates
getSymbols('CL',src='yahoo',from = '2020-01-07' , to = '2020-03-09')
getSymbols('SPY',src='yahoo',from = '2020-01-07' , to = '2020-03-09')
getSymbols('DX-Y.NYB',src='yahoo',from = '2020-01-07' , to = '2020-03-09')
getSymbols('^VIX',src='yahoo',from = '2020-01-07' , to = '2020-03-09')
getSymbols('^IRX',src='yahoo',from = '2020-01-07' , to = '2020-03-09')

Y <- cbind(COP$COP.Close
           ,VIX$VIX.Close)

stock_list <- c("CL", "SPY", 'DX-Y.NYB', '^VIX', '^IRX')
start_date <- '2020-01-07'
end_date <- '2020-03-09'
master_df <- NULL

for (idx in seq(length(stock_list))){
  stock_index = stock_list[idx]
  getSymbols(stock_index, verbose = TRUE, src = "yahoo",
             from=start_date,to=end_date)
  temp_df = as.data.frame(get(stock_index))
  temp_df$Date = row.names(temp_df)
  temp_df$Index = stock_index
  row.names(temp_df) = NULL
  colnames(temp_df) = c("Open", "High", "Low", "Close",
                        "Volume", "Adjusted", "Date", "Index")
  temp_df = temp_df[c("Date", "Index", "Open", "High",
                      "Low", "Close", "Volume", "Adjusted")]
  master_df = rbind(master_df, temp_df)
}

# Build an indicator variable with a 1 at only T*+1
post_shock_indicator <- c(rep(0, nrow( na.omit(diff(log(CL$CL.Adjusted)))) - 1), 1)

# Throw the external regressors into a matrix
X <- as.matrix(
  cbind(
    na.omit(diff(log(CL$CL.Adjusted))),
    na.omit(diff(log(SPY$SPY.Adjusted))),
    diff(   log(   na.omit( as.data.frame(get('DX-Y.NYB'))$"DX-Y.NYB.Adjusted")  )),
    diff(log(   na.omit( as.data.frame(get('VIX'))$"VIX.Adjusted")  )),
    diff(log(   na.omit(  as.data.frame(get('IRX'))$"IRX.Adjusted")  )),
    post_shock_indicator
  )
)

mymod <- garchx( as.numeric(na.omit(diff(log(COP$COP.Adjusted)))) , order = c(1,1),
                 xreg = X, control = list(eval.max = 10000, iter.max = 15000, rel.tol = 1e-6))
mymod
coeftest(mymod)
BIC(mymod)
predict(mymod, n.ahead = 1, newxreg = matrix(c(-.01,0), nrow = 1))




temp <- SynthVolForecast <- function(series_matrix
                                     ,covariates_series_list
                                     ,shock_time_vec
                                     ,shock_length_vec
                                     ,geometric_sets
                                     ,p_dbw
                                     ,covariate_vector_list #tk
                                     ,days_before_shocktime_vec #tk I may want to remove this
                                     ,garch_order
                                     ,plots = TRUE
)

  temp$forecast
