options(digits = 7, scipen = 7)

sysname <- Sys.info()["sysname"]

if(sysname == "Darwin") {
  setwd("~/Desktop/PhD/synthetic_vol_forecasting/") # example on mac machine
} else if(sysname == "Linux") {
  setwd('~/Desktop/synthetic_vol_forecasting/synthVolForecast.R') # example on linux machine
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
           ,'highfrequency'
           ,'tidyverse'
           ,'xts'
           ,'dplyr'
           ,'RcppRoll'
)

suppressPackageStartupMessages(lapply(packs, require, character.only = TRUE))

## BEGIN USER DATA INPUTS##

k <-1

log_ret_covariates <- c(
                       # "CL=F"
                        "^VIX"
                        ,"^IRX"
                        # ,"^FVX"
                        # ,"^TNX"
                        # ,"^TYX"
                        ,'^XAU'
                        #,"DX-Y.NYB"
                      )

level_covariates <- c('^VIX'
                      #,'^IRX'
)

volume_covariates <- c('SPY')

FRED_covariates <- c('AAA', 'BAA')

shock_dates <- list('Powell Hikes' = "2018-12-18"
                ,'Dec. 2015' = "2015-12-15"
                 ,'Dec. 2016' = "2016-12-13"
                 ,'Mar. 2017' = "2017-03-14"
                 ,'Jun. 2017' = "2017-06-13"
                 , 'Dec. 2017' ="2017-12-12"
                 ,'Mar. 2018' = "2018-03-20"
                 ,'Jun. 2018' = "2018-06-12"
                 ,'Sept. 2018' =  "2018-09-25"
)

# shock_dates <- list('Powell Hikes' = "2018-12-19"
#                     ,'Dec. 2015' = "2015-12-16"
#                     ,'Dec. 2016' = "2016-12-14"
#                     ,'Mar. 2017' = "2017-03-15"
#                     ,'Jun. 2017' = "2017-06-14"
#                     , 'Dec. 2017' ="2017-12-13"
#                     ,'Mar. 2018' = "2018-03-21"
#                     ,'Jun. 2018' = "2018-06-13"
#                     ,'Sept. 2018' =  "2018-09-26"
# )


#shock_dates <- c(shock_dates[1], list.reverse(shock_dates[2:length(shock_dates)]))

## END USER DATA INPUTS##

nyse <- timeDate::holidayNYSE(2000:year(Sys.Date())+1)
create.calendar(name='NYSE', holidays=nyse, weekdays=c('saturday', 'sunday'))

shock_dates_as_dates <- as.Date(as.Date(unlist(shock_dates)))

start_dates <- offset(shock_dates_as_dates, round(-.1*252), "NYSE")

k_periods_after_shock <- offset(shock_dates_as_dates, k, "NYSE")

market_data_list <- vector("list", length(shock_dates))
names(market_data_list) <- names(shock_dates)

# Now we loop through shock dates

for (i in 1:length(shock_dates)){

  print(shock_dates_as_dates[i])

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

    to_add <- c(data_log_ret_covariates
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

    colnames(merged_data) <- c(log_ret_covariates_colname
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

#Build Y
RVSPY <- as.xts(SPYRM$RV1, order.by = SPYRM$DT)

RVSPY_subset <- data.frame(
  Date = SPYRM$DT
  ,RV1 = SPYRM$RV1
  #,close = SPYRM$CLOSE
)

RVSPY_subset_2 <- RVSPY_subset %>%
  mutate(tomorrow_RV1 = lead(RV1))

RVSPY_subset_3 = RVSPY_subset_2 %>%
  mutate(RV5 = roll_mean(RV1, 5, align = "right", fill = NA)) %>%
  mutate(RV22 = roll_mean(RV1, 22, align = "right", fill = NA))

RVSPY_complete <- RVSPY_subset_3[complete.cases(RVSPY_subset_3), ]
head(RVSPY_complete)

#RVSPY_complete <- data.frame(RVSPY_complete %>% mutate(donor = ifelse(Date %in% shock_dates,1,0)))

RVSPY_complete$log_tomorrow_RV1 <- log(RVSPY_complete$tomorrow_RV1)
RVSPY_complete$difflog <- c(NA,diff(log(RVSPY_complete$tomorrow_RV1)))

#RVSPY_complete[RVSPY_complete$rate_move == 1, "donor"] <- shock_dates_as_dates

#RVSPY_complete$donor = as.factor(RVSPY_complete$donor)

RVSPY_complete$RV1_neg <- ifelse(RVSPY_complete$difflog > 0, 0, RVSPY_complete$RV1)
RVSPY_complete$RV5_neg <- ifelse(RVSPY_complete$difflog > 0, 0, RVSPY_complete$RV5)
RVSPY_complete$RV22_neg <- ifelse(RVSPY_complete$difflog > 0, 0, RVSPY_complete$RV22)

rownames(RVSPY_complete) <- RVSPY_complete$Date

head(RVSPY_complete)

# RVSPY_final <- xts(RVSPY_complete, order.by = RVSPY_complete$Date)

RVSPY_final <- RVSPY_complete[,c('tomorrow_RV1'
                                 , 'RV1'
                                 , 'RV5'
                                 , 'RV22'
                                 , 'RV1_neg'
                                 ,'RV5_neg'
                                 ,'RV22_neg')]

head(RVSPY_final)

#Date index
#First column is outcome
#Last column is donor indicator
#Columns in between are regressors


#Now build X
X <- list()
for (i in 1:length(start_dates)){
  X[[i]] <- market_data_list[[i]][,-1]
}

n <- length(start_dates) - 1

time_date <- gsub(" ", "", gsub(':', '', format(Sys.time(), "%a%b%d%X%Y")), fixed = TRUE)
png_save_name <- paste("real_data_output_plots/savetime_"
                       ,time_date
                       ,'_'
                       #,TSUS
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

#png(png_save_name,width = 800, height = 600)

#Now run the algorithm
output <- HAR(RVSPY_final
            ,X
            ,shock_time_vec = unlist(shock_dates)
            ,shock_length_vec
            ,k=1 #does it make sense for a k-step ahead?
            ,dbw_scale = TRUE
            ,dbw_center = TRUE
            ,dbw_indices = 1:ncol(X[[1]])
            ,dbw_Y_lookback = c(0)
            ,dbw_princ_comp_input = NULL
            ,covariate_indices = NULL
            ,geometric_sets = NULL #tk
            ,plots = TRUE #tk
            ,shock_time_labels = names(shock_dates)
            ,ground_truth_vec = NULL
            ,Y_lookback_indices_input = list(seq(1,15,1))
            ,X_lookback_indices_input = rep(list(c(10)),length(1:ncol(X[[1]])))
)

output$predictions
output$linear_combinations

#dev.off()
