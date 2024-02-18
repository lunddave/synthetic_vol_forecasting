#library(SynthVolForecast)

options(digits = 7, scipen = 7)

sysname <- Sys.info()["sysname"]

if(sysname == "Darwin") {
  setwd("~/Desktop/PhD/synthetic_vol_forecasting/") # example on mac machine
} else if(sysname == "Linux") {
  setwd('~/Desktop/synthetic_vol_forecasting/') # example on linux machine
} 

### BEGIN 2016 election example
packs <- c('quantmod'
           ,'bizdays'
           ,'lubridate'
           ,'rlist'
           ,'Rsolnp'
           ,'garchx'
           ,'lmtest'
           ,'RColorBrewer'
           ,'forecast'
)

suppressPackageStartupMessages(lapply(packs, require, character.only = TRUE))

## BEGIN USER DATA INPUTS##
ground_truth <- c(0.000712, 0.000976)
ground_truth <- c(0.000712)

ground_truth <- c(0.000712, 0.000976, .0006, .0001, .00001)

ground_truth <- c(0.000712, 0.000976, .0006)

ground_truth <- c(0.000712, 0.000976)
ground_truth <- c(0.000712)

k <- 1

TSUS <- "SPY"

log_ret_covariates <- c("^VIX"
                        ,"^IRX"
                        ,"^FVX"
                        ,"^TNX"
                        ,"^TYX"
                        ,"CL=F"
                        #,"DX-Y.NYB"
                      )

level_covariates <- c('^VIX'
                      #,"GBP=X"
                      #,'^IRX'
)

volume_covariates <- c()

FRED_covariates <- c('AAA', 'BAA')

shock_dates <- list('Oct 7' = "2023-10-09"
                 ,'Rus-Ukr' = '2022-02-24'
                 ,'2021' = "2021-05-10"
                 ,'2018' = "2018-04-12"
                 ,'2014' = "2014-07-07"
                 , 'Lebanon' = '2006-07-12'
                 #, '2000' ="2000-10-01"
)

shock_dates <- c(shock_dates[1], list.reverse(shock_dates[2:length(shock_dates)]))

## END USER DATA INPUTS##

nyse <- timeDate::holidayNYSE(2000:year(Sys.Date())+1)
create.calendar(name='NYSE', holidays=nyse, weekdays=c('saturday', 'sunday'))

shock_dates_as_dates <- as.Date(as.Date(unlist(shock_dates)))

start_dates <- offset(shock_dates_as_dates, round(-3*252), "NYSE")

k_periods_after_shock <- offset(shock_dates_as_dates, k, "NYSE")

market_data_list <- vector("list", length(shock_dates))
names(market_data_list) <- names(shock_dates)

# Now we loop through shock dates

for (i in 1:length(shock_dates)){

  print(shock_dates_as_dates[i])

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

  data_absolute_return_covariates <- c()

  if (length(FRED_covariates) > 0){

    #Get FRED data (requires subtracting one series from another)
    data_FRED_covariates <- lapply(FRED_covariates, function(sym) {
      na.omit(getSymbols(sym
                         ,src = 'FRED'
                         ,from=start_dates[i]
                         ,to=k_periods_after_shock[i]+20 #tk +10
                         ,auto.assign=FALSE))[,1]})

    # Now we are going manually and carefully add the spread between BAA and AAA credit

    diff_log_FRED_spread <- diff(log(data_FRED_covariates[[2]] - data_FRED_covariates[[1]]))

    month_before_shock <- month(shock_dates_as_dates[i]) - 1
    year_of_shock <- year(shock_dates_as_dates[i]) #This is not necessarily robust to a shock in January

    first_of_month_before_shock <- paste(year_of_shock
                                         ,'-'
                                         ,month_before_shock
                                         ,'-'
                                         ,'01'
                                         ,sep = '')

    first_of_month_before_shock <- as.Date(first_of_month_before_shock)

    diff_log_FRED_spread_latest_before_shock <- diff_log_FRED_spread[first_of_month_before_shock]

    # date <- shock_dates_as_dates[i]
    # data <- as.numeric(diff_log_FRED_spread_latest_before_shock)
    # row_to_add <- xts(data, order.by = date)

    #Drop NA
    diff_log_FRED_spread_latest_before_shock <- na.omit(diff_log_FRED_spread_latest_before_shock)

    } #end lapply for FRED covariates

    to_add <- c(data_TSUS
                , data_log_ret_covariates
                , data_level_covariates
                , data_volume_covariates
                , data_absolute_return_covariates
    )

    merged_data <- do.call(merge, to_add)

    #We add column names to the data so we can analyze it more easily when we print it

    if (length(log_ret_covariates) > 0)
    {log_ret_covariates_colname <- paste(log_ret_covariates, "_log_ret", sep="")}
    else
    {log_ret_covariates_colname <- c()}

    if (length(level_covariates) > 0)
    {level_covariates_colname <- paste(level_covariates, "_raw", sep="")}
    else
    {level_covariates_colname <- c()}

    if (length(volume_covariates) > 0)
    {volume_covariates_colname <- paste(volume_covariates, "_volume", sep="")}
    else
    {volume_covariates_colname <- c()}

    if (length(data_absolute_return_covariates) > 0)
    {abs_log_ret_covariates_colname <- paste(log_ret_covariates, "_abs_log_ret", sep="")}
    else
    {abs_log_ret_covariates_colname <- c()}

    colnames(merged_data) <- c(TSUS
                               , log_ret_covariates_colname
                               , level_covariates_colname
                               , volume_covariates_colname
                               , abs_log_ret_covariates_colname
    )

    if (length(FRED_covariates) > 0){

      dates <- index(merged_data)
      times <- length(dates)
      data <- rep(as.numeric(diff_log_FRED_spread_latest_before_shock), times)
      Debt_Risk_Spread <- xts(data, order.by = dates)

      #Merge FRED data
      merged_data <- merge(merged_data, Debt_Risk_Spread)

    } #end conditional for FRED_covariates

    market_data_list[[i]] <- merged_data
  }

##################################

#now build Y
Y <- list()
for (i in 1:length(start_dates)){
  Y_i <- market_data_list[[i]][,1]
  Y_i_drop_NA <- Y_i[complete.cases(Y_i)]
  print('Here is the type of object we are working with:')
  print(class(Y_i_drop_NA))
  # print('Here are the rownames')
  # print(index(Y_i_drop_NA))
  print('Here is the shock date')
  print(shock_dates_as_dates[i])

  if (shock_dates[i] %in% index(Y_i_drop_NA)){
    print('The shock date is in the series.')
  }
  else{
    print(paste('Shock date ', i, ' NOT in series ',i,".", sep = ''))
  }

  Y[[i]] <- Y_i_drop_NA
} #end for loop for building Y

#Now build X
X <- list()
for (i in 1:length(start_dates)){
  X[[i]] <- market_data_list[[i]][,-1]
}

n <- length(start_dates) - 1

time_date <- gsub(" ", "", format(Sys.time(), "%a%b%d%X%Y"), fixed = TRUE)
png_save_name <- paste("real_data_output_plots/savetime_"
                       ,time_date
                       ,'_'
                       ,TSUS
                       ,'_'
                       ,paste(log_ret_covariates,collapse='-')
                       ,'_'
                       ,paste(level_covariates,collapse='-')
                       ,'_'
                       # ,paste(volume_covariates,collapse='-')
                       # ,'_'
                       # ,paste(data_absolute_return_covariates,collapse='-')
                       # ,'_'
                       ,paste(shock_dates,collapse='-')
                       ,".png"
                       ,sep="")

png(png_save_name,width = 800, height = 600)

#Now run the algorithm
temp <- SynthVolForecast(Y
                         ,X
                         ,shock_time_vec = unlist(shock_dates)
                         ,rep(k, n+1)
                         ,dbw_scale = TRUE
                         ,dbw_center = TRUE
                         ,dbw_indices = NULL
                         # ,dbw_princ_comp_input = 3
                         #,covariate_indices = length(X)
                         ,garch_order = c(1,1,1)
                         ,plots = TRUE
                         ,shock_time_labels = names(shock_dates)
                         ,ground_truth_vec = ground_truth
                         )

dev.off()
