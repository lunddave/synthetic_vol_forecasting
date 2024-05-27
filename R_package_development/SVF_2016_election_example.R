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
)

suppressPackageStartupMessages(lapply(packs, require, character.only = TRUE))

## BEGIN USER DATA INPUTS##
ground_truth <- c(0.000712, 0.000976)
ground_truth <- c(0.000712)
#ground_truth <- c(0.000712, 0.000976, .0006)

k <- length(ground_truth)

TSUS <- 'IYG'

log_ret_covariates <- c(#"GBP=X",
                         #"6B=F",
                        "CL=F"
                        ,"^VIX"
                        ,"^IRX"
                        ,"^FVX"
                        ,"^TNX"
                        ,"^TYX"
                        #,"DX-Y.NYB"
                      )

level_covariates <- c('^VIX'
                      #,"GBP=X"
                      #,'^IRX'
)

volume_covariates <- c()

FRED_covariates <- c('AAA', 'BAA')

shock_dates_outside_loop <- list('2016 Election' = "2016-11-08"
                                 ,'Brexit' = "2016-06-22"
                               # ,'2014 Midterm' = "2014-11-04"
                               ,'2012 Election' = "2012-11-06"
                               # , '2010 Midterm' ="2010-11-02"
                               ,'2008 Election' = "2008-11-04"
                               # , '2006 Midterm' ="2006-11-07"
                               ,'2004 Election' = "2004-11-02"
                               #,'2002 Midterm' =  "2002-11-05"
                               #,'2000 Election' = "2000-11-07"
              )

shock_dates_outside_loop <- c(shock_dates_outside_loop[1]
                              , list.reverse(shock_dates_outside_loop[2:length(shock_dates_outside_loop)]))
## END USER DATA INPUTS##

#We now loop through donors 2...n+1 as part of LOO procedure

number_of_covariates <- length(log_ret_covariates) +
                        length(level_covariates) +
                        length(volume_covariates) +
                        2 #The +2 is for the TSUS and the FRED covariates'

list_from_looping <- list()

#NOTE: to run without an Leave-one-out procedures, set...
# z <- u <- 0

for (u in 0:number_of_covariates){

  for (z in c(0,2:length(shock_dates_outside_loop))){

    if (z > 0){
      shock_dates <- shock_dates_outside_loop[-z]
      } else {shock_dates <- shock_dates_outside_loop}

    nyse <- timeDate::holidayNYSE(2000:year(Sys.Date())+1)
    create.calendar(name='NYSE', holidays=nyse, weekdays=c('saturday', 'sunday'))

    shock_dates_as_dates <- as.Date(as.Date(unlist(shock_dates)))

    start_dates <- offset(shock_dates_as_dates, round(-3.5*252), "NYSE")

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

          print('Number of rows in merged_data')
          print(ncol(merged_data))

          print('Value of number_of_covariates')
          print(number_of_covariates)

          if (u > 0){
            merged_data_uth_covariate_dropped <- merged_data[,-u]
          }
          else {
            merged_data_uth_covariate_dropped <- merged_data
          }

          covariates_col_names <- colnames(merged_data_uth_covariate_dropped)

          covariate_string <- paste(covariates_col_names,collapse="_")

        } #end conditional for FRED_covariates

        market_data_list[[i]] <- merged_data_uth_covariate_dropped
    }

    ##################################

    #now build Y
    Y <- list()
    for (i in 1:length(start_dates)){
      Y_i <- market_data_list[[i]][,1]
      Y_i_drop_NA <- Y_i[complete.cases(Y_i)]

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

    time_date <- gsub(" ", "", gsub(':', '', format(Sys.time(), "%a%b%d%X%Y")), fixed = TRUE)
    png_save_name <- paste("real_data_output_plots/"
                           ,time_date
                           ,'_'
                           ,TSUS
                           ,'_'
                           ,covariate_string
                            ,'_'
                           ,paste(shock_dates,collapse='-')
                           #,".png"
                           ,sep="")

    png_save_name <- gsub('.', '', png_save_name, fixed = TRUE)

    png_save_name <- paste(png_save_name, '.png', sep = "")

    png(png_save_name,width = 800, height = 600)

    #Now run the algorithm
    temp <- SynthVolForecast(Y
                             ,X
                             ,shock_time_vec = unlist(shock_dates)
                             ,rep(k, n+1)
                             ,dbw_scale = TRUE
                             ,dbw_center = TRUE
                             ,dbw_indices = NULL
                             #,covariate_indices = length(X)
                             ,garch_order = c(1,1,0)
                             ,plots = TRUE
                             ,shock_time_labels = names(shock_dates)
                             ,ground_truth_vec = ground_truth
                             ,Y_lookback_indices = list(seq(1,30,1))
                             ,X_lookback_indices = rep(list(c(1)),ncol(X[[1]]))
                             )

    list_from_looping <- append(list_from_looping, list(temp))

    dev.off()

  }#end loop for dropping donors

}#end loop for dropping covariates

prediction_matrix <- sapply(list_from_looping,"[[",2)
prediction_matrix <- matrix(as.numeric(t(prediction_matrix)), ncol = 3)

median_unadj_pred <- round(median(prediction_matrix[,1]),5)
indices_we_want <- round(prediction_matrix[,1],5) == median_unadj_pred
#we exclude indices where a non-2016 election is the TSUS

apply(prediction_matrix[indices_we_want,], 2, mean)

lin_combs <- sapply(list_from_looping,"[[",1)
df_lin_combs <- do.call(rbind.data.frame, lin_combs)
names(df_lin_combs) <- c('d1','d2','d3')
apply(df_lin_combs, 1, sum)
lin_combs_averages <- apply(df_lin_combs[indices_we_want,], 2, mean)
lin_combs_averages #these are no meaningful averages, since we drop donors
#tk but what about calculating averages for sets with same donors?



#And what's the loss on these averaged forecasts?
forecast_averages <- apply(prediction_matrix[indices_we_want,], 2, mean)

losses_from_averaging <- QL_loss_function(forecast_averages
                                , ground_truth)
losses_from_averaging

loss_matrix <- sapply(list_from_looping,"[[",3)
loss_matrix <- matrix(as.numeric(t(loss_matrix)), ncol = 3)
loss_matrix <- as.data.frame(loss_matrix)
apply(loss_matrix[indices_we_want,], 2, mean)

nrow_loss_matrix <- nrow(loss_matrix)

loss_matrix$loss_from_avg <- rep(losses_from_averaging[2], nrow_loss_matrix)

names(loss_matrix) <- c('unadj'
                        , 'adj'
                        , 'arith_mean'
                        ,'LOSS_from_averaging_all_forecasts')


loss_matrix$CH_wins <- ifelse(loss_matrix_with_forecast_averages$unadj > loss_matrix_with_forecast_averages$adj
                              , TRUE, FALSE)
loss_matrix$AM_wins <- ifelse(loss_matrix_with_forecast_averages$unadj > loss_matrix_with_forecast_averages$arith_mean
                              , TRUE, FALSE)
loss_matrix$average_adjusted_wins <- ifelse(loss_matrix_with_forecast_averages$unadj > loss_matrix_with_forecast_averages$arith_mean
                              , TRUE, FALSE)

win_df <- t(apply(loss_matrix, 1, function(x) x == min(x)))

colSums(win_df)



#Things we want from the loop:
# 1) each convex comb
# 2) unadj, adj, and arithmetic_mean
# 3) loss of unadj, adj, and arithmetic_mean
# 4) the FE estimates? this is not easy to justify


