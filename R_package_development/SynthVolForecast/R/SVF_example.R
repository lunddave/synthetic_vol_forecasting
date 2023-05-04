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
                             ,rep(2, n+1)
                             ,garch_order = c(1,1)
                             ,plots = TRUE)
  }

  )

temp$linear_combinations
temp$predictions

### END EXAMPLE WITH SIMULATED DATA

### BEGIN EXAMPLE WITH ConocoPhillips DATA

library(quantmod)
library(bizdays)
library(lubridate)

## BEGIN USER DATA INPUTS##
log_ret_covariates <- c("COP" #first should be time series under study
                        ,"CL=F"
                        ,'^VIX'
                        ,'^IRX'
                        ,'SPY'
                        ,'XOM'
                        ,'CVX'
                        ,'SHEL'
                        ,'TTE'
                        ,'BP'
                        )

level_covariates <- c('^VIX','^IRX')

volume_covariates <- c("COP" #first should be time series under study
                        ,'SPY'
                        ,'XOM'
                        ,'CVX'
                        ,'SHEL'
                        ,'TTE'
                        ,'BP'
)

shock_dates <- c("2020-03-06"
                 ,"2014-11-26"
                 , "2008-09-25"
                 , "2008-09-12"
                 , "2008-09-05"
                 , "2008-03-14"
                 )

k <- 1
## END USER DATA INPUTS##

nyse <- timeDate::holidayNYSE(2000:year(Sys.Date())+1)
create.calendar(name='NYSE', holidays=nyse, weekdays=c('saturday', 'sunday'))
shock_dates_as_dates <- as.Date(shock_dates)
start_dates <- offset(shock_dates_as_dates, round(-.2*252), "NYSE")
k_periods_after_shock <- offset(shock_dates_as_dates, k, "NYSE")

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
                       ,auto.assign=FALSE))[,6]})

  to_add_3 <- lapply(volume_covariates, function(sym) {
    na.omit(getSymbols(sym
                       ,from=start_dates[i]
                       ,to=k_periods_after_shock[i]+10 #tk +10
                       ,auto.assign=FALSE))[,6]})

  to_add <- c(to_add, to_add_2, to_add_3)
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

n <- length(start_dates) - 1

#Now run the algorithm
temp <- SynthVolForecast(Y
                         ,X
                         ,shock_time_vec = shock_dates
                         ,rep(k, n+1)
                         ,dwb_indices = NULL
                         #,covariate_indices = c((length(X)-1):length(X))
                         ,garch_order = c(1,1,1)
                         ,plots = TRUE)

temp$linear_combinations
temp$predictions

#Now run the algorithm
temp <- SynthPrediction(Y
                         ,X
                         ,shock_time_vec = shock_dates
                         ,rep(k, n+1)
                         ,dwb_indices = NULL
                         ,covariate_indices = (length(X)-1):length(X)
                         ,plots = TRUE
                         ,display_ground_truth_choice = TRUE)


## GARCH on
mod <- garchx(Y[[1]][1:253])
plot.ts(Y[[1]][1:253])
fitted(mod)
plot.ts(fitted(mod))

COP
