set.seed(1986)

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
           ,'MCS'

)

suppressPackageStartupMessages(lapply(packs, require, character.only = TRUE))

## BEGIN USER DATA INPUTS##
ground_truth <- c(0.000712, 0.000976)
ground_truth <- c(0.000712)
#ground_truth <- c(0.000712, 0.000976, .0006)

k <- length(ground_truth)

TSUS <- 'IYG'

log_ret_covariates <- c(#"GBP=X",
                        #"6B=F"
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

#FRED_covariates <- c('AAA', 'BAA')

FRED_covariates <- c()

shock_dates_outside_loop <- list('2016 Election' = c("2016-11-08", list(1, c(1,1,0)))
                                 ,'Brexit' = c("2016-06-23", list(1, c(1,1,0)))
                                ,'Brexit Poll Released' = c("2016-06-13", list(1, c(1,1,0)))
                                #   ,'Boris Johnson' = "2016-02-20"
                                ,'Brexit Announced' = c("2016-02-19", list(1, c(1,1,0)))
                                # ,'2014 Midterm' = "2014-11-04"
                                ,'2012 Election' = c("2012-11-06", list(1, c(1,1,0)))
                                # , '2010 Midterm' ="2010-11-02"
                                ,'2008 Election' = c("2008-11-04", list(1, c(1,1,0)))
                                # , '2006 Midterm' ="2006-11-07"
                                ,'2004 Election' = c("2004-11-02", list(1, c(1,1,0)))
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
                        1 #The +2 is for the TSUS and the FRED covariates'

list_from_looping <- list()

#NOTE: to run without an Leave-one-out procedures, set...
# z <- u <- 0

for (u in c(-1,0,2:number_of_covariates)){

  for (z in c(0,2:length(shock_dates_outside_loop))){

    if (z > 0){
      shock_dates <- shock_dates_outside_loop[-z]
      dropped_donor <- shock_dates_outside_loop[z]
    } else {
        shock_dates <- shock_dates_outside_loop
        dropped_donor <- 'None'
        }

    nyse <- timeDate::holidayNYSE(2000:year(Sys.Date())+1)
    create.calendar(name='NYSE', holidays=nyse, weekdays=c('saturday', 'sunday'))

    shock_dates_as_dates <- as.Date(sapply(shock_dates,"[[",1))

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
                                       ,to=k_periods_after_shock[i]+40 #tk +10
                                       ,auto.assign=FALSE))[,6]
                    ,type='log')})

      data_log_ret_covariates <- lapply(log_ret_covariates, function(sym) {
        dailyReturn(na.omit(getSymbols(sym
                                       ,from=start_dates[i]
                                       ,to=k_periods_after_shock[i]+40 #tk +10
                                       ,auto.assign=FALSE))[,6]
                    ,type='log')})

      data_level_covariates <- lapply(level_covariates, function(sym) {
        na.omit(getSymbols(sym
                           ,from=start_dates[i]
                           ,to=k_periods_after_shock[i]+40 #tk +10
                           ,auto.assign=FALSE))[,6]})

      data_volume_covariates <- lapply(volume_covariates, function(sym) {
        dailyReturn(na.omit(getSymbols(sym
                                       ,from=start_dates[i]
                                       ,to=k_periods_after_shock[i]+40 #tk +10
                                       ,auto.assign=FALSE))[,6])})

      data_absolute_return_covariates <- lapply(data_log_ret_covariates, abs)

      data_absolute_return_covariates <- c()

      if (length(FRED_covariates) > 0){

        #Get FRED data (requires subtracting one series from another)
        data_FRED_covariates <- lapply(FRED_covariates, function(sym) {
          na.omit(getSymbols(sym
                             ,src = 'FRED'
                             ,from=start_dates[i]
                             ,to=k_periods_after_shock[i]+40 #tk +10
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
        {log_ret_covariates_colname <- paste("Log_Return_of_", log_ret_covariates, sep="")}
        else
        {log_ret_covariates_colname <- c()}

        if (length(level_covariates) > 0)
        {level_covariates_colname <- paste("Raw", level_covariates, sep="")}
        else
        {level_covariates_colname <- c()}

        if (length(volume_covariates) > 0)
        {volume_covariates_colname <- paste("Volume", volume_covariates, sep="")}
        else
        {volume_covariates_colname <- c()}

        if (length(data_absolute_return_covariates) > 0)
        {abs_log_ret_covariates_colname <- paste("Absolute_Log_Return_of", log_ret_covariates, sep="")}
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

        print('Number of rows in merged_data')
        print(ncol(merged_data))

        print('Value of number_of_covariates')
        print(number_of_covariates)

        if (u > 0){
          merged_data_uth_covariate_dropped <- merged_data[,-u]
          covariate_string <- names(merged_data[,u])
          Y_lookback_indices_u_loop <- list(seq(1,3,1))
        }
        else if (u == 0){
          merged_data_uth_covariate_dropped <- merged_data
          covariates_col_names <- colnames(merged_data_uth_covariate_dropped)
          covariate_string <- 'None'
          Y_lookback_indices_u_loop <- list(seq(1,3,1))
        }
        else{ #This case is for the time series under study not being used in the DBW
          merged_data_uth_covariate_dropped <- merged_data
          covariates_col_names <- colnames(merged_data_uth_covariate_dropped)[-1]
          covariate_string <- TSUS
          Y_lookback_indices_u_loop <- NULL
        }

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

      if (shock_dates_as_dates[i] %in% index(Y_i_drop_NA)){
        print('The shock date is in the series.')
      }
      else{
        print(paste(shock_dates_as_dates[i], ' is NOT in series Y ',i,".", sep = ''))
      }

      Y[[i]] <- Y_i_drop_NA
    } #end for loop for building Y

    print('We have now finished building Y.')

    #Now build X
    X <- list()
    for (i in 1:length(start_dates)){
      X[[i]] <- market_data_list[[i]][,-1]
    }

    n <- length(start_dates) - 1

    time_date <- gsub(" ", "", gsub(':', '', format(Sys.time(), "%a%b%d%X%Y")), fixed = TRUE)

    print('We have now accessed the system time from the computer.')

    #We have two naming conventions:

    #1 What donors and covariates are INCLUDED?
    png_save_name <- paste("real_data_output_plots/"
                           ,time_date
                           ,'_'
                           ,TSUS
                           ,'_'
                           ,covariate_string
                            ,'_'
                           ,dropped_donor
                           #,".png"
                           ,sep="")

    png_save_name <- gsub('.', '', png_save_name, fixed = TRUE)

    png_save_name <- paste(png_save_name, '.png', sep = "")

    print('We have created a file name for the png to output.')

    png(png_save_name,width = 800, height = 600)

    print('We have used the png function in R.')



    #Now run the algorithm
    function_output <- SynthVolForecast(Y
                             ,X
                             ,shock_time_vec = sapply(shock_dates,"[[",1)
                             ,shock_length_vec = sapply(shock_dates,"[[",2)
                             ,dbw_scale = TRUE
                             ,dbw_center = TRUE
                             ,dbw_indices = NULL
                             #,covariate_indices = length(X)
                             ,garch_order = lapply(shock_dates,"[[",3)
                             ,plots = TRUE
                             ,shock_time_labels = names(shock_dates)
                             ,ground_truth_vec = ground_truth
                             ,Y_lookback_indices = Y_lookback_indices_u_loop
                             ,X_lookback_indices = rep(list(c(1)),ncol(X[[1]]))
                             )

    function_output[[4]] <- covariate_string
    function_output[[5]] <- dropped_donor

    list_from_looping <- append(list_from_looping, list(function_output))

    dev.off()

  }#end loop for dropping donors

}#end loop for dropping covariates

library(xtable)

covs <- sapply(list_from_looping,"[[",4)

covs <- gsub("_", " ", covs, fixed = TRUE)
covs <- gsub(" .", " ", covs, fixed = TRUE)
covs

dons <- sapply(sapply(list_from_looping,"[[",5),'[[',1)



prediction_matrix <- sapply(list_from_looping,"[[",2)
prediction_matrix <- matrix(as.numeric(t(prediction_matrix)), ncol = 3)

indices_we_want <- 1:nrow(prediction_matrix)

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

forecast_medians <- apply(prediction_matrix[indices_we_want,], 2, median)

losses_from_averaging <- QL_loss_function(forecast_averages
                                , ground_truth)
losses_from_averaging

losses_from_medians <- QL_loss_function(forecast_medians
                                          , ground_truth)
losses_from_medians

loss_matrix <- sapply(list_from_looping,"[[",3)
loss_matrix <- matrix(as.numeric(t(loss_matrix)), ncol = 3)
loss_matrix <- as.data.frame(loss_matrix)
apply(loss_matrix[indices_we_want,], 2, mean)

nrow_loss_matrix <- nrow(loss_matrix)

loss_matrix$loss_from_avg <- rep(losses_from_averaging[2], nrow_loss_matrix)

loss_matrix$loss_from_median <- rep(losses_from_medians[2], nrow_loss_matrix)

names(loss_matrix) <- c('unadj'
                        ,'QLLoss(AdjustedForecast)'
                        ,'arith_mean'
                        ,'QLLoss(AverageAdjustedForecast)'
                        ,'QLLoss(MedianAdjustedForecast)')

loss_matrix$dropped_covariate <- covs
loss_matrix$dropped_donor <- dons

loss_matrix <- as.data.frame(lapply(loss_matrix, unlist))

print(loss_matrix[order(loss_matrix[,2]),], row.names = FALSE)

print(xtable(loss_matrix[order(loss_matrix[,2]),c(2,6,7)],digits=4)
      , include.rownames = FALSE
      , size="\\fontsize{8pt}{10pt}\\selectfont")

lm1 <- lm(as.numeric(loss_matrix[,2]) ~ 0 + as.factor(loss_matrix$dropped_covariate) + as.factor(loss_matrix$dropped_donor))
summary(lm1)

lm1 <- lm(as.numeric(loss_matrix[,2]) ~ 0 + as.factor(loss_matrix$dropped_donor))
summary(lm1)

lm1 <- lm(as.numeric(loss_matrix[,2]) ~ 0 + as.factor(loss_matrix$dropped_covariate) )
summary(lm1)

plot(lm1)

win_df <- t(apply(loss_matrix, 1, function(x) x == min(x)))

colSums(win_df)


print(xtable(colSums(win_df),digits=4)
      , include.rownames = FALSE
      , size="\\fontsize{9pt}{10pt}\\selectfont")


#MCS
loss_only <- as.data.frame(t(as.matrix(loss_matrix[,2], ncol = 1)))
loss_only[2,] <- rgamma(50,1,1)
cols <- as.vector(paste(dons, covs, sep = ''))

colnames(loss_only) <- cols

loss_xts <- xts(x=loss_only, order.by=as.Date(c("2016-11-08","2016-11-09")))

dim(loss_xts)
head(loss_xts)

dimnames(loss_xts)

class(dimnames(loss_xts)[[2]])

dimnames(Loss)
class(dimnames(Loss)[[2]])

MCS <- MCSprocedure(Loss=loss_xts,alpha=0.1,B=5000,statistic='Tmax',cl=NULL)

MCS
