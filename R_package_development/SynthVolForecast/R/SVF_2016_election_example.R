source('SynthVolForecast_functions.R',
       echo = FALSE,
       verbose = FALSE)

options(digits = 7, scipen = 7)

### BEGIN 2016 election example

library(quantmod)
library(bizdays)
library(lubridate)

## BEGIN USER DATA INPUTS##
ground_truth <- c(0.000712, 0.000976)
ground_truth <- c(0.000712)

log_ret_covariates <- c("IYG" #first should be time series under study
                        ,"GBP=X"
                        ,"6B=F"
                        ,"CL=F"
                        ,"^VIX"
                        ,"^IRX"
                        ,"^FVX"
                        ,"^TNX"
                       # ,"^TYX"
                        )

level_covariates <- c(#'^VIX'
                      #,"GBP=X"
                      #,'^IRX'
                      )

volume_covariates <- c()

shock_dates <- c("2016-11-08",
                   "2016-06-23"
                 , "2014-11-04"
                 , "2012-11-06"
                 , "2010-11-02"
                 , "2008-11-04"
                 , "2006-11-07"
                 #, "2004-11-02"
                 #, "2002-11-05"
                 #, "2000-11-07"
                 )

k <- 1
## END USER DATA INPUTS##

nyse <- timeDate::holidayNYSE(2000:year(Sys.Date())+1)
create.calendar(name='NYSE', holidays=nyse, weekdays=c('saturday', 'sunday'))

shock_dates_as_dates <- as.Date(shock_dates)

start_dates <- offset(shock_dates_as_dates, round(-5*252), "NYSE")

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

  merged_data <- do.call(merge, to_add)
  complete_cases_merged_data <- merged_data[complete.cases(merged_data),]
  market_data_list[[i]] <- complete_cases_merged_data

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

png('SVF_2016.png')

#Now run the algorithm
temp <- SynthVolForecast(Y
                         ,X
                         ,shock_time_vec = shock_dates
                         ,rep(k, n+1)
                         ,dwb_indices = NULL
                         #,covariate_indices = length(X)
                         ,garch_order = c(1,0,1)
                         ,plots = TRUE
                         ,ground_truth_vec = ground_truth)

dev.off()

#Now run the algorithm
#temp <- SynthPrediction(Y
#                         ,X
#                         ,shock_time_vec = shock_dates
#                         ,rep(k, n+1)
#                         ,dwb_indices = NULL
#                         ,covariate_indices = length(X)
#                         ,plots = TRUE
#                        ,display_ground_truth_choice = TRUE)



