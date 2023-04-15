# Purpose herein is verifying that a GARCH model can be fit
# with external regressors such that at least one regressor
# is a vector with only 1 non-zero value

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

mymod <- garchx( as.numeric(na.omit(diff(log(COP$COP.Adjusted)))) , order = c(2,1),
                  xreg = X, control = list(eval.max = 10000, iter.max = 15000, rel.tol = 1e-6))
mymod
coeftest(mymod)
BIC(mymod)
predict(mymod, n.ahead = 1, newxreg = matrix(c(-.01,0), nrow = 1))

mymod_differenced <- garchx(  as.numeric(diff(diff(log(COP$COP.Adjusted)))), order = c(0,1), 
                  xreg = X, control = list(eval.max = 10000, iter.max = 6250, rel.tol = 1e-6))
mymod_differenced
coeftest(mymod_differenced)
AIC(mymod_differenced)


# Let's look at other two shock dates
# 
# 2008-09-08
# 2014-11-27

# loads COP and West Texas Intermediate
getSymbols('COP',src='yahoo',from = '2008-08-08' , to = '2008-09-09')
getSymbols('CL',src='yahoo',from = '2008-08-07' , to = '2008-09-08')

# Build an indicator variable with a 1 at only T*+1
post_shock_indicator <- c(rep(0, nrow(diff(log(CL$CL.Adjusted))) - 1), 1)
# Throw the external regressors into a matrix
X <- as.matrix(cbind(diff(log(CL$CL.Adjusted)), post_shock_indicator))

mymod <- garchx(  as.numeric(diff(log(COP$COP.Adjusted))), order = c(0,1),
                  xreg = X, control = list(eval.max = 10000, iter.max = 15000, rel.tol = 1e-6))
mymod
coeftest(mymod)
AIC(mymod)
predict(mymod, n.ahead = 1, newxreg = matrix(c(-.01,0), nrow = 1))

mymod_differenced <- garchx(  as.numeric(diff(diff(log(COP$COP.Adjusted)))), order = c(1,1), 
                              xreg = X, control = list(eval.max = 10000, iter.max = 6250, rel.tol = 1e-6))
mymod_differenced
coeftest(mymod_differenced)
AIC(mymod_differenced)


####

# loads COP and West Texas Intermediate
getSymbols('COP',src='yahoo',from = '2014-10-09' , to = '2014-12-07')
getSymbols('CL',src='yahoo',from = '2014-10-09' , to = '2014-12-07')
getSymbols('IRX',src='yahoo',from = '2014-10-09' , to = '2014-12-07')

# Build an indicator variable with a 1 at only T*+1
post_shock_indicator <- c(rep(0, nrow(diff(log(CL$CL.Adjusted))) - 6), rep(1,6))
# Throw the external regressors into a matrix
X <- as.matrix( cbind( diff(log(CL$CL.Open)), diff(log(IRX$IRX.Open)),  post_shock_indicator))

mymod <- garchx(  as.numeric(diff(log(COP$COP.Adjusted))), order = c(1,1),
                  xreg = X, control = list(eval.max = 10000, iter.max = 25000, rel.tol = 1e-6))
mymod
coeftest(mymod)
BIC(mymod)

predict(mymod, n.ahead = 1, newxreg = matrix(c(-.01,1), nrow = 1))










# time series of interest
getSymbols('COP',src='yahoo',from = '2014-10-09' , to = '2014-12-05')

# covariates
getSymbols('CL',src='yahoo',from = '2014-10-09' , to = '2014-12-05')
getSymbols('SPY',src='yahoo',from = '2014-10-09' , to = '2014-12-05')
getSymbols('DX-Y.NYB',src='yahoo',from = '2014-10-09', to = '2014-12-05')
getSymbols('^VIX',src='yahoo',from = '2014-10-09', to = '2014-12-05')
getSymbols('^IRX',src='yahoo',from = '2014-10-09', to = '2014-12-05')

# Build an indicator variable with a 1 at only T*+1

# Let k denote the number of days the shock last IN addition to Nov 28 2014
k <- 4
post_shock_indicator <- c(rep(0, nrow( na.omit(diff(log(CL$CL.Adjusted)))) - (k + 1)), rep(1,k+1))
# Throw the external regressors into a matrix
X <- as.matrix(
  cbind(
    #na.omit(diff(log(CL$CL.Adjusted))),
    #na.omit(diff(log(SPY$SPY.Adjusted))),
    #diff( log(   na.omit( as.data.frame(get('DX-Y.NYB'))$"DX-Y.NYB.Adjusted")  )),
    #diff(log(   na.omit( as.data.frame(get('VIX'))$"VIX.Adjusted")  )),
    #diff(log(   na.omit(  as.data.frame(get('IRX'))$"IRX.Adjusted")  )),
    post_shock_indicator
  )
)

mymod <- garchx(  as.numeric(na.omit(diff(log(COP$COP.Adjusted)))) , order = c(1,0),
                  xreg = X, control = list(eval.max = 10000, iter.max = 15000, rel.tol = 1e-6))
mymod
coeftest(mymod)
BIC(mymod)

predict(mymod, n.ahead = 1, newxreg = matrix(c(0), nrow = 1)) #predict as if shock effect is gone
predict(mymod, n.ahead = 1, newxreg = matrix(c(1), nrow = 1)) #predict as if shock effect remains



# March 2008

# time series of interest
getSymbols('COP',src='yahoo',from = '2008-01-02' , to = '2008-03-21')

# covariates
getSymbols('CL',src='yahoo',from = '2008-01-02' , to = '2008-03-21')
getSymbols('SPY',src='yahoo',from = '2008-01-02' , to = '2008-03-21')
getSymbols('DX-Y.NYB',src='yahoo',from = '2008-01-02', to = '2008-03-21')
getSymbols('^VIX',src='yahoo',from = '2008-01-02', to = '2008-03-21')
getSymbols('^IRX',src='yahoo',from = '2008-01-02', to = '2008-03-21')

# Build an indicator variable with a 1 at only T*+1

# Let k denote the number of days the shock last IN addition to Nov 28 2014
k <- 3
post_shock_indicator <- c(rep(0, nrow( na.omit(diff(log(CL$CL.Adjusted)))) - (k + 1)), rep(1,k+1))
# Throw the external regressors into a matrix
X <- as.matrix(
  cbind(
    #na.omit(diff(log(CL$CL.Adjusted))),
    #na.omit(diff(log(SPY$SPY.Adjusted))),
    diff( log(   na.omit( as.data.frame(get('DX-Y.NYB'))$"DX-Y.NYB.Adjusted")  )),
    diff(log(   na.omit( as.data.frame(get('VIX'))$"VIX.Adjusted")  )),
    #diff(log(   na.omit(  as.data.frame(get('IRX'))$"IRX.Adjusted")  )),
    post_shock_indicator
  )
)

mymod <- garchx(  as.numeric(na.omit(diff(log(COP$COP.Adjusted)))) , order = c(1,0),
                  xreg = X, control = list(eval.max = 10000, iter.max = 15000, rel.tol = 1e-6))
mymod
coeftest(mymod)
BIC(mymod)

predict(mymod, n.ahead = 1, newxreg = matrix(c(0), nrow = 1)) #predict as if shock effect is gone
predict(mymod, n.ahead = 1, newxreg = matrix(c(1), nrow = 1)) #predict as if shock effect remains

