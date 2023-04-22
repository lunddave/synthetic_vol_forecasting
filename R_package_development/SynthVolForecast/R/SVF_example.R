options(digits = 7)

### BEGIN EXAMPLE WITH SIMULATED DATA

Tee <- 192 #Specify length of time series
n <- 4 #Specify number of donors
shock_time_vec <- rep(Tee/2, n+1)

Y <- list()

#we populate the covariate list
for (i in 1:(n+1)){
  Y[[i]] <- rnorm(Tee)
}

X <- list()

#we populate the covariate list
for (i in 1:(n+1)){
  X[[i]] <- cbind(rnorm(Tee)
             ,rnorm(Tee)
             ,rnorm(Tee))
}

system.time(

  {

    temp <- SynthVolForecast(Y
                             ,X
                             ,shock_time_vec
                             ,rep(1, n+1)
                             ,garch_order = c(1,1)
                             ,plots = TRUE)
  }

  )

temp$convex_combination
temp$forecast

### END EXAMPLE WITH SIMULATED DATA

### BEGIN EXAMPLE WITH ConocoPhillips DATA

library(quantmod)
library(garchx)
library(lmtest)
library(bizdays)
library(lubridate)

## BEGIN USER DATA INPUTS##
log_ret_covariates <- c("COP"
                        ,"CL=F"
                        ,'^VIX'
                        ,'^IRX'
                        ,'SPY'
                        ,'XOM'
                        ,'CVX'
                        ,'SHEL'
                        ,'TTE'
                        ,'BP') #first should be TSUS
level_covariates <- c('^VIX')

shock_dates <- c("2020-03-06"
                 ,"2014-11-26"
                 , "2008-09-25"
                 , "2008-09-12"
                 , "2008-09-05"
                 , "2008-03-14")

k <- 1
## END USER DATA INPUTS##

nyse <- timeDate::holidayNYSE(2000:year(Sys.Date()) +1)
create.calendar(name='NYSE', holidays=nyse, weekdays=c('saturday', 'sunday'))
shock_dates_as_dates <- as.Date(shock_dates)
start_dates <- offset(shock_dates_as_dates, -1*252, "NYSE")
k_periods_after_shock <- offset(shock_dates_as_dates, k, "NYSE")

shock_dates
shock_dates_as_dates
start_dates
k_periods_after_shock

market_data_list <- vector("list", length(shock_dates))
names(market_data_list) <- shock_dates

for (i in 1:length(shock_dates)){

  to_add <- lapply(log_ret_covariates, function(sym) {
    dailyReturn(na.omit(getSymbols(sym
                                   ,from=start_dates[i]
                                   ,to=k_periods_after_shock[i]+10 #tk +10
                                   ,auto.assign=FALSE))
                                   ,type='log')})

  to_add_2 <- lapply(level_covariates, function(sym) {
    na.omit(getSymbols(sym
                                   ,from=start_dates[i]
                                   ,to=k_periods_after_shock[i]+10 #tk +10
                                   ,auto.assign=FALSE))[,4]})

  to_add <- c(to_add, to_add_2)

  market_data_list[[i]] <- do.call(merge, to_add)

}

##################################

#now build Y
Y <- list()
for (i in 1:length(start_dates)){
  Y[[i]] <- market_data_list[[i]][,1]
}

#Now build X
X <- list()
for (i in 1:length(start_dates)){
  X[[i]] <- market_data_list[[i]][,-1]
}

n <- length(Y) - 1

#Now run the algorithm
temp <- SynthVolForecast(Y
                         ,X
                         ,shock_time_vec = shock_dates
                         ,rep(k, n+1)
                         ,1:length(X)
                         ,covariate_indices = c(length(X))
                         ,garch_order = c(1,1,1)
                         ,plots = TRUE)



temp$convex_combination
temp$forecast

## GARCH on
mod <- garchx(Y[[1]][1:253])
fitted(mod)
plot.ts(fitted(mod))
