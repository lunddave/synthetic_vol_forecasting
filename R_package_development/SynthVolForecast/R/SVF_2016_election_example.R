source('SynthVolForecast_functions.R',
       echo = FALSE,
       verbose = FALSE)

options(digits = 7, scipen = 7)

### BEGIN 2016 election example
packs <- c('quantmod'
          , 'bizdays'
          , 'lubridate'
          )

suppressPackageStartupMessages(lapply(packs, require, character.only = TRUE))

## BEGIN USER DATA INPUTS##
ground_truth <- c(0.000712, 0.000976)
ground_truth <- c(0.000712)
#ground_truth <- c(0.000712, 0.000976, .0006)


k <- 1

TSUS <- 'IYG'

log_ret_covariates <- c(#"GBP=X",
                        "6B=F"
                        ,"CL=F"
                        ,"^VIX"
                        ,"^IRX"
                        ,"^FVX"
                        ,"^TNX"
                        ,"^TYX"
                        )

level_covariates <- c('^VIX'
                      #,"GBP=X"
                      #,'^IRX'
                      )

volume_covariates <- c('IYG')

shock_dates <- c("2016-11-08",
                 "2016-06-23"
               , "2014-11-04"
               , "2012-11-06"
               , "2010-11-02"
               , "2008-11-04"
               , "2006-11-07"
               , "2004-11-02"
               , "2002-11-05"
               #, "2000-11-07"
               )

## END USER DATA INPUTS##

nyse <- timeDate::holidayNYSE(2000:year(Sys.Date())+1)
create.calendar(name='NYSE', holidays=nyse, weekdays=c('saturday', 'sunday'))

shock_dates_as_dates <- as.Date(shock_dates)

start_dates <- offset(shock_dates_as_dates, round(-1.8*252), "NYSE")

k_periods_after_shock <- offset(shock_dates_as_dates, k, "NYSE")

market_data_list <- vector("list", length(shock_dates))
names(market_data_list) <- shock_dates

for (i in 1:length(shock_dates)){

  data_TSUS <- lapply(TSUS, function(sym) {
    dailyReturn(na.omit(getSymbols(sym
                                   ,from=start_dates[i]
                                   ,to=k_periods_after_shock[i]+20 #tk +10
                                   ,auto.assign=FALSE))[,6]
                                   ,type='log')})

  data_log_ret_covariates <- lapply(log_ret_covariates, function(sym) {
    dailyReturn(na.omit(getSymbols(sym
                                   ,from=start_dates[i]
                                   ,to=k_periods_after_shock[i]+20 #tk +10
                                   ,auto.assign=FALSE))[,6]
                                   ,type='log')})

  data_level_covariates <- lapply(level_covariates, function(sym) {
    na.omit(getSymbols(sym
                       ,from=start_dates[i]
                       ,to=k_periods_after_shock[i]+20 #tk +10
                       ,auto.assign=FALSE))[,6]})

  data_volume_covariates <- lapply(volume_covariates, function(sym) {
    dailyReturn(na.omit(getSymbols(sym
                       ,from=start_dates[i]
                       ,to=k_periods_after_shock[i]+20 #tk +10
                       ,auto.assign=FALSE))[,6])})

  data_absolute_return_covariates <- lapply(data_log_ret_covariates, abs)

  to_add <- c(data_TSUS
            , data_log_ret_covariates
            , data_level_covariates
            , data_volume_covariates
            , data_absolute_return_covariates
            )

  merged_data <- do.call(merge, to_add)
  market_data_list[[i]] <- merged_data

}

##################################

#now build Y
Y <- list()
for (i in 1:length(start_dates)){
  Y_i <- market_data_list[[i]][,1]
  Y_i_drop_NA <- Y_i[complete.cases(Y_i)]
  #print('Here is the type of object we are working with:')
  #print(class(Y_i_drop_NA))
  #print('Here are the rownames')
  #print(index(Y_i_drop_NA))
  if (shock_dates[i] %in% index(Y_i_drop_NA)){
    print('The shock date is in the series.')
  }
  else{print('Shock date NOT in series.')}
  Y[[i]] <- Y_i_drop_NA
}


#Now build X
X <- list()
for (i in 1:length(start_dates)){
  X[[i]] <- market_data_list[[i]][,-1]
}

n <- length(start_dates) - 1

time_date <- gsub(" ", "", format(Sys.time(), "%a%b%d%X%Y"), fixed = TRUE)
png_save_name <- paste("/home/david/Desktop/synthetic_vol_forecasting/R_package_development/SynthVolForecast/R/real_data_output_plots/savetime_"
                       ,time_date
                       ,'_'
                       ,TSUS
                       ,'_'
                       ,paste(log_ret_covariates,collapse='-')
                       ,'_'
                       ,paste(level_covariates,collapse='-')
                       ,'_'
                       ,paste(volume_covariates,collapse='-')
                       ,'_'
                       ,paste(data_absolute_return_covariates,collapse='-')
                       ,'_'
                       ,paste(shock_dates,collapse='-')
                       ,".png"
                       ,sep="")

png(png_save_name)

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
