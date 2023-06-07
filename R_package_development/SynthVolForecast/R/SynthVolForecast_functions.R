## Required packages
packs <- c('Rsolnp'
          , 'garchx'
          , 'lmtest'
          , 'forecast'
          )
suppressPackageStartupMessages(lapply(packs, require, character.only = TRUE))


### START QL_loss_function
QL_loss_function <- function(pred, gt){pred/gt - log(pred/gt) - 1}
### END QL_loss_function


### START dbw
dbw <- function(X
                ,dbw_indices
                ,shock_time_vec
                ,scale = FALSE
                ,sum_to_1 = 1
                ,bounded_below_by = 0
                ,bounded_above_by = 1
                ,normchoice = c('l1', 'l2')[2]
                ,penalty_normchoice = c('l1', 'l2')[1]
                ,penalty_lambda = 0
) { # https://github.com/DEck13/synthetic_prediction/blob/master/prevalence_testing/numerical_studies/COP.R
  # X is a list of covariates for the time series
  # X[[1]] should be the covariate of the time series to predict
  # X[[k]] for k = 2,...,n+1 are covariates for donors

  # T^* is a vector of shock-effects time points
  # shock effect point must be > 2

  # number of time series for pool
  n <- length(X) - 1

  # COVARIATE FOR TIME SERIES UNDER STUDY AT shock_time_vec
  X1 <- X[[1]][shock_time_vec[1], dbw_indices, drop = FALSE] # we get only 1 row

  #We notify user if p > n, i.e. if linear system is overdetermined
  p <- length(dbw_indices)

  # LOOP for grab shock_time_vec covariate vector for each donor
  X0 <- c()
  for (i in 1:n) {
    X0[[i]] <- X[[i+1]][shock_time_vec[i+1], dbw_indices
                          , drop = FALSE] #get 1 row from each donor
  }

  #################################
  #begin if statement
  if (scale == TRUE) {
    dat <- rbind(X1, do.call('rbind', X0)) # do.call is for cluster computing?
    dat <- apply(dat, 2, function(x) scale(x, center = TRUE, scale = TRUE))
    X1 <- dat[1, dbw_indices
              , drop = FALSE]
    X0 <- c()
    for (i in 1:n) {
      X0[[i]] <- dat[i+1, dbw_indices, drop = FALSE] #we are repopulating X0[[i]] with scaled+centered data
    } #end loop
  } #end if statement
  #################################

  # objective function
  weightedX0 <- function(W) {
    # W is a vector of weight of the same length of X0
    n <- length(W)
    p <- ncol(X1)
    XW <- matrix(0, nrow = 1, ncol = p)
    for (i in 1:n) {
      XW <- XW + W[i] * X0[[i]]
    } #end of loop

    #normchoice
    if (normchoice == 'l1') {
      norm <- as.numeric(norm(matrix(X1 - XW), type = "1"))
    }
    else {
      norm <- as.numeric(crossprod(matrix(X1 - XW)))
    }

    #now add penalty
    if (penalty_normchoice == 'l1' & penalty_lambda > 0) {
      norm <- norm + penalty_lambda * norm(as.matrix(W), type = "1")
    }
    else if (penalty_normchoice == 'l2' & penalty_lambda > 0) {
      norm <- norm + penalty_lambda * as.numeric(crossprod(matrix(W)))
    }
    else {norm <- norm}

    return(norm)
  } #end objective function

  # optimization and return statement

  # I have added features
  # 1) The option to remove the sum-to-1 constraint
  # 2) The option to change the lower bound to -1 or NA
  # 3) option to change the upper bound to NA'
  # 4) option to choose l1 or l2 norm as distance function

  #Thus I need if statements to implement these...

  # conditional for sum to 1
  if (is.na(sum_to_1) == FALSE) {eq_constraint <- function(W) sum(W) - 1 }
  else{eq_constraint = NULL}

  # conditional for bounding below
  if (is.na(bounded_below_by) == FALSE)
  {
    lower_bound = rep(bounded_below_by, n)
  }
  else if (is.na(bounded_below_by) == TRUE)  {
    lower_bound = NULL
  }

  #conditional for bounding above
  if (is.na(bounded_above_by) == FALSE)
  {
    upper_bound = rep(1, n)
  }
  else if (is.na(bounded_above_by) == TRUE)  {
    upper_bound = NULL
  }

  object_to_return <- solnp(par = rep(1/n, n),
                            fun = weightedX0,
                            eqfun = eq_constraint,
                            eqB = 0,
                            LB = lower_bound, UB = upper_bound,
                            control = list(trace = 0
                                           , 1.0e-8
                                           , tol = 1e-19
                                           , outer.iter = 100000
                                           , inner.iter = 100000))
  return(object_to_return$pars)

} #END dbw function
### END dbw

### START GARCH plot_maker_garch
plot_maker_garch <- function(fitted_vol
                      ,shock_times_for_barplot
                      ,shock_time_vec #mk
                      ,shock_length_vec
                      ,unadjusted_pred
                      ,w_hat
                      ,omega_star_hat
                      ,adjusted_pred
                      ,arithmetic_mean_based_pred
                      ,ground_truth_vec){

  if (is.character(shock_times_for_barplot) == FALSE){
    shock_times_for_barplot <- 1:length(shock_times_for_barplot)
  }

  par(mfrow = c(1,2))

  #PLOT ON THE LEFT:
  #Plot donor weights
  barplot(w_hat
          ,  main = 'Donor Pool Weights'
          , names.arg = shock_times_for_barplot[-1]
          , cex.names=.8
          , las=2)

  #Plot target series and prediction

  thing_to_get_max_of <- c(as.numeric(fitted_vol)
                        , unadjusted_pred
                        , adjusted_pred
                        , ground_truth_vec
                        , arithmetic_mean_based_pred
                        )

  max_for_y_lim <- max(thing_to_get_max_of)

  #PLOT ON THE RIGHT:
  plot.ts(fitted_vol[1:shock_time_vec[1]], #mk
       main = 'Two Forecasts Following T*', #mk can improve this title
       ylab = '',
       xlab = "Trading Days",
       xlim = c(0, shock_time_vec[1] + 5), #mk
       ylim = c(min(0, fitted_vol),  max_for_y_lim))

  title(ylab = expression(sigma^2), line = 2.05, cex.lab = 1.99) # Add y-axis text

  # Here is the color scheme we will use
  colors_for_adjusted_pred <- c('red'
                              , "green"
                              , "purple"
                              , 'blue'
                              )

  # Let's add the plain old GARCH prediction
  points(y = unadjusted_pred
         ,x = (shock_time_vec[1]+1):(shock_time_vec[1]+shock_length_vec[1])
         ,col = colors_for_adjusted_pred[1]
         ,cex = 1.1
         ,pch = 15)

  # Now plot the adjusted predictions
  points(y = adjusted_pred
         ,x = (shock_time_vec[1]+1):(shock_time_vec[1]+shock_length_vec[1])
         ,col = colors_for_adjusted_pred[2]
         ,cex = 1.1
         ,pch = 23)

  # Now plot the arithmetic mean-based predictions
  points(y = arithmetic_mean_based_pred
         ,x = (shock_time_vec[1]+1):(shock_time_vec[1]+shock_length_vec[1])
         ,col = colors_for_adjusted_pred[3]
         ,cex = 1.1
         ,pch = 23)

  # Now plot Ground Truth
  points(y = ground_truth_vec
         ,x = (shock_time_vec[1]+1):(shock_time_vec[1]+shock_length_vec[1])
         ,col = colors_for_adjusted_pred[4]
         ,cex = 1.1
         ,pch = 23)

  labels_for_legend <- c('GARCH (unadjusted)'
                        , 'Adjusted'
                        , 'Arithmetic Mean'
                        , 'Ground Truth'
                        )

  legend(x = "topleft",  # Coordinates (x also accepts keywords) #mk
         legend = labels_for_legend,
         1:length(labels_for_legend), # Vector with the name of each group
         colors_for_adjusted_pred,   # Creates boxes in the legend with the specified colors
         title = 'Prediction Method',      # Legend title,
         cex = .9)

}
### END plot_maker_garch


### START plot_maker_synthprediction
plot_maker_synthprediction <- function(Y
                       ,shock_times_for_barplot
                       ,shock_time_vec #mk
                       ,shock_length_vec
                       ,unadjusted_pred
                       ,w_hat
                       ,omega_star_hat
                       ,adjusted_pred
                       ,display_ground_truth = FALSE){

  if (is.character(shock_times_for_barplot) == FALSE){
    shock_times_for_barplot <- 1:length(shock_times_for_barplot)
  }

  #First print donor series
  par(mfrow = c(round(sqrt(n)),ceiling(sqrt(n))))

  for (i in 2:(n+1)){
    plot.ts(Y[[i]][1:shock_time_vec[i]]
            ,xlab = 'Trading Days'
            ,ylab = 'Differenced Logarithm'
            ,main = paste('Donor ', i,': ', shock_times_for_barplot[i], sep = '')
            ,xlim = c(0, shock_time_vec[i] + 5)
            ,ylim = c(min(Y[[i]]),  max(Y[[i]]))
            )

    if (display_ground_truth == TRUE){

      lines(y = Y[[i]][shock_time_vec[i]:(shock_time_vec[i] + shock_length_vec[i])]
            ,x = shock_time_vec[i]:(shock_time_vec[i] + shock_length_vec[i])
            ,col = 'purple'
            ,cex = 1.1
            ,lty = 3)

      points(y = Y[[i]][(shock_time_vec[i]+1):(shock_time_vec[i] + shock_length_vec[i])]
             ,x = (shock_time_vec[i]+1):(shock_time_vec[i] + shock_length_vec[i])
             ,col = 'red'
             ,cex = 1.1
             ,pch = 24)

    }
  }

  #Now print time series under study
  par(mfrow = c(1,2))

  #PLOT ON THE LEFT:
  #Plot donor weights
  barplot(w_hat
          ,  main = 'Donor Pool Weights'
          , names.arg = shock_times_for_barplot[-1]
          , cex.names=.8)

  #Plot target series and prediction

  thing_to_get_max_of <- c(as.numeric(Y[[1]]), unadjusted_pred, adjusted_pred)

  max_for_y_lim <- max(thing_to_get_max_of)

  #PLOT ON THE RIGHT:
  plot.ts(Y[[1]][1:shock_time_vec[1]], #mk
          main = 'Two Forecasts Following T*', #mk can improve this title
          ylab = '',
          xlab = "Trading Days",
          xlim = c(0, shock_time_vec[1] + 5), #mk
          ylim = c(min(0, Y[[1]]),  max_for_y_lim))

  title(ylab = 'Differenced Logarithm', line = 2.05, cex.lab = 1.99) # Add y-axis text

  # Here is the color scheme we will use
  colors_for_adjusted_pred <- c('red', "green",'purple')

  # Let's add the plain old GARCH prediction
  points(y = unadjusted_pred
         ,x = (shock_time_vec[1]+1):(shock_time_vec[1]+shock_length_vec[1])
         ,col = colors_for_adjusted_pred[1]
         ,cex = .9
         ,pch = 15)

  # Now plot the adjusted predictions
  points(y = adjusted_pred
         ,x = (shock_time_vec[1]+1):(shock_time_vec[1]+shock_length_vec[1])
         ,col = colors_for_adjusted_pred[2]
         ,cex = 1.1
         ,pch = 23)

  if (display_ground_truth == TRUE){

    lines(y = Y[[1]][shock_time_vec[1]:(shock_time_vec[1] + shock_length_vec[1])]
          ,x = shock_time_vec[1]:(shock_time_vec[1] + shock_length_vec[1])
          ,col = colors_for_adjusted_pred[3]
          ,cex = 1.1
          ,lty = 3)

    points(y = Y[[1]][(shock_time_vec[1]+1):(shock_time_vec[1] + shock_length_vec[1])]
          ,x = (shock_time_vec[1]+1):(shock_time_vec[1] + shock_length_vec[1])
          ,col = colors_for_adjusted_pred[3]
          ,cex = 1.1
          ,pch = 24)

  }

  labels_for_legend <- c('ARIMA (unadjusted)', 'Adjusted Prediction', 'Actual')

  legend(x = "topleft",  # Coordinates (x also accepts keywords) #mk
         legend = labels_for_legend,
         1:length(labels_for_legend), # Vector with the name of each group
         colors_for_adjusted_pred,   # Creates boxes in the legend with the specified colors
         title = 'Prediction Method',      # Legend title,
         cex = .9)

}
### END plot_maker_synthprediction


### START SynthVolForecast
SynthVolForecast <- function(Y_series_list
                             ,covariates_series_list
                             ,shock_time_vec
                             ,shock_length_vec
                             ,dwb_indices = NULL
                             ,covariate_indices = NULL
                             ,geometric_sets = NULL #tk
                             ,days_before_shocktime_vec = NULL #tk I may want to remove this
                             ,garch_order = NULL
                             ,plots = TRUE
                             ,ground_truth_vec
){
  ### BEGIN Doc string
  #tk
  ### END Doc string

  ### BEGIN Populate defaults
  n <- length(Y_series_list) - 1

  if (is.null(garch_order) == TRUE) {
    garch_order <- c(1,1,1)
  }

  if (is.null(dwb_indices) == TRUE) {
    dwb_indices <- 1:ncol(X[[1]])
  }

  ### END Populate defaults




  ## BEGIN Check that inputs are all comformable/acceptable
  n <- length(Y_series_list) - 1 #tk
  ## END Check that inputs are all comformable/acceptable

  integer_shock_time_vec <- c() #mk

  ## BEGIN Check whether shock_time_vec is int/date

  for (i in 1:(n+1)){

    if (is.character(shock_time_vec[i]) == TRUE){
      integer_shock_time_vec[i] <- which(index(Y[[i]]) == shock_time_vec[i]) #mk
    }
    else{ integer_shock_time_vec <- shock_time_vec}

  }

  ## END Check whether shock_time_vec is int/date

  ## BEGIN calculate weight vector
  w_hat <- dbw(X, #tk
               dwb_indices,
               integer_shock_time_vec,
               scale = TRUE,
               sum_to_1 = TRUE, #tk
               bounded_below_by = 0, #tk
               bounded_above_by = 1, #tk
               # normchoice = normchoice, #tk
               # penalty_normchoice = penalty_normchoice,
               # penalty_lambda = penalty_lambda
               )
  ## END calculate weight vector

  ## BEGIN estimate fixed effects in donors
  omega_star_hat_vec <- c()

  for (i in 2:(n+1)){

    # Make indicator variable w/ a 1 at only T*+1, T*+2,...,T*+shock_length_vec[i]
    vec_of_zeros <- rep(0, integer_shock_time_vec[i])
    vec_of_ones <- rep(1, shock_length_vec[i])
    post_shock_indicator <- c(vec_of_zeros, vec_of_ones)
    last_shock_point <- integer_shock_time_vec[i] + shock_length_vec[i]

    #subset X_i
    if (is.null(covariate_indices) == TRUE) {
      X_i_penultimate <- cbind(Y_series_list[[i]][1:last_shock_point]
                               , post_shock_indicator)
      X_i_final <- X_i_penultimate[,2]
    }
    else {
      X_i_subset <- X[[i]][1:last_shock_point,covariate_indices]
      X_i_with_indicator <- cbind(X_i_subset, post_shock_indicator)
      X_i_final <- X_i_with_indicator
    }

    print(paste('Now fitting a GARCH model to donor series number ', i,'.', sep = ''))
    fitted_garch <- garchx(Y_series_list[[i]][1:last_shock_point] #tk
                   , order = garch_order
                   , xreg = X_i_final
                   , backcast.values = NULL
                   , control = list(eval.max = 10000
                   , iter.max = 15000
                   , rel.tol = 1e-6))

    coef_test <- coeftest(fitted_garch)
    extracted_fixed_effect <- coef_test[dim(coeftest(fitted_garch))[1], 1]
    omega_star_hat_vec <- c(omega_star_hat_vec, extracted_fixed_effect)

  } ## END loop for computing fixed effects

  ## END estimate fixed effects in donors

  ## BEGIN compute linear combination of fixed effects
  omega_star_hat <- w_hat %*% omega_star_hat_vec
  ## END compute linear combination of fixed effects


  ## BEGIN fit GARCH to target series

  if (is.null(covariate_indices) == TRUE){

    fitted_garch <- garchx(Y_series_list[[1]][1:integer_shock_time_vec[1]]
                           , order = garch_order
                           , xreg = NULL 
                           , backcast.values = NULL
                           , control = list(eval.max = 10000
                                            , iter.max = 15000
                                            , rel.tol = 1e-6))

    print('Outputting the fitted GARCH for time series under study.')
    print(fitted_garch)

    unadjusted_pred <- predict(fitted_garch, n.ahead = shock_length_vec[1])
  }
  else{
    ## BEGIN fit GARCH to target series
    fitted_garch <- garchx(Y_series_list[[1]][1:integer_shock_time_vec[1]]
                           , order = garch_order
                           , xreg = X[[1]][1:integer_shock_time_vec[1],covariate_indices]
                           , control = list(eval.max = 10000
                                            , iter.max = 15000
                                            , rel.tol = 1e-6))

    print('Outputting the fitted GARCH for the time series under study.')
    print(fitted_garch)

    #Note: for forecasting, we use last-observed X value
    X_to_use_in_forecast <- X[[1]][integer_shock_time_vec[1],covariate_indices]

    X_replicated_for_forecast_length <- matrix(rep(X_to_use_in_forecast, k)
                                               , nrow = shock_length_vec[1]
                                               , byrow = TRUE)


    forecast_period <- (integer_shock_time_vec[1]+1):(integer_shock_time_vec[1]+shock_length_vec[1])
    mat_X_for_forecast <- cbind(Y_series_list[[1]][forecast_period]
                           , X_replicated_for_forecast_length)

    unadjusted_pred <- predict(fitted_garch
                               , n.ahead = shock_length_vec[1]
                               , newxreg = mat_X_for_forecast[,-1])
  }

  adjusted_pred <- unadjusted_pred + rep(omega_star_hat, k)

  arithmetic_mean_based_pred <- rep(mean(omega_star_hat_vec), k) + unadjusted_pred

  QL_loss_unadjusted_pred <- sum(QL_loss_function(unadjusted_pred, ground_truth))

  QL_loss_adjusted_pred <- sum(QL_loss_function(adjusted_pred, ground_truth))

  list_of_linear_combinations <- list(w_hat)
  list_of_forecasts <- list(unadjusted_pred, adjusted_pred)
  names(list_of_forecasts) <- c('unadjusted_pred', 'adjusted_pred')

  output_list <- list(list_of_linear_combinations
                      , list_of_forecasts)

  names(output_list) <- c('linear_combinations', 'predictions')

  ## tk OUTPUT
  cat('--------------------------------------------------------------\n',
      '-------------------SynthVolForecast Results-------------------','\n',
      '--------------------------------------------------------------\n',
      'Donors:', n, '\n',  '\n',
      'Shock times:', shock_time_vec, '\n', '\n',
      'Lengths of shock times:', shock_length_vec, '\n', '\n',
      'Convex combination:',w_hat,'\n', '\n',
      'Shock estimates:', omega_star_hat_vec, '\n', '\n',
      'Aggregate estimated shock effect:', omega_star_hat, '\n', '\n',
      'Unadjusted Forecast:', unadjusted_pred,'\n', '\n',
      'Adjusted Forecast:', adjusted_pred,'\n', '\n',
      'Arithmetic-Mean-Based Forecast:',arithmetic_mean_based_pred,'\n','\n',
      'Ground Truth (estimated by realized volatility):', ground_truth_vec,'\n', '\n',
      'QL Loss of unadjusted:', QL_loss_unadjusted_pred,'\n', '\n',
      'QL Loss of adjusted:', QL_loss_adjusted_pred,'\n', '\n'
  )

  ## PLOTS

  if (plots == TRUE){
    cat('\n User has opted to produce plots.','\n')
    plot_maker_garch(fitted(fitted_garch)
               ,shock_time_vec
               ,integer_shock_time_vec
               ,shock_length_vec
               ,unadjusted_pred
               ,w_hat
               ,omega_star_hat
               ,adjusted_pred
               ,arithmetic_mean_based_pred
               ,ground_truth_vec)

  }

  return(output_list)

} ### END SynthVolForecast

### START SynthPrediction
SynthPrediction <- function(Y_series_list
                             ,covariates_series_list
                             ,shock_time_vec
                             ,shock_length_vec
                             ,dwb_indices = NULL
                             ,covariate_indices = NULL
                             ,geometric_sets = NULL #tk
                             ,days_before_shocktime_vec = NULL #tk I may want to remove this
                             ,arima_order = NULL
                             ,user_ic_choice = c('aicc','aic','bic')[1]
                             ,plots = TRUE
                             ,display_ground_truth_choice = FALSE
){
  ### BEGIN Doc string
  #tk
  ### END Doc string

  ### BEGIN Populate defaults
  n <- length(Y_series_list) - 1

  if (is.null(arima_order) == TRUE) {
    arima_order <- c(1,1,1)
  }

  if (is.null(dwb_indices) == TRUE) {
    dwb_indices <- 1:ncol(X[[1]])
  }

  ### END Populate defaults

  ## BEGIN Check that inputs are all comformable/acceptable
  n <- length(Y_series_list) - 1 #tk
  ## END Check that inputs are all comformable/acceptable

  integer_shock_time_vec <- c() #mk

  ## BEGIN Check whether shock_time_vec is int/date

  for (i in 1:(n+1)){

    if (is.character(shock_time_vec[i]) == TRUE){
      integer_shock_time_vec[i] <- which(index(Y[[i]]) == shock_time_vec[i]) #mk
    }
    else{ integer_shock_time_vec <- shock_time_vec}

  }

  ## END Check whether shock_time_vec is int/date

  ## BEGIN estimate fixed effects in donors
  omega_star_hat_vec <- c()

  order_of_arima <- list()

  for (i in 2:(n+1)){

    # Make indicator variable w/ a 1 at only T*+1, T*+2,...,T*+shock_length_vec[i]
    vec_of_zeros <- rep(0, integer_shock_time_vec[i])
    vec_of_ones <- rep(1, shock_length_vec[i])
    post_shock_indicator <- c(vec_of_zeros, vec_of_ones)
    last_shock_point <- integer_shock_time_vec[i] + shock_length_vec[i]

    #subset X_i
    if (is.null(covariate_indices) == TRUE) {
      X_i_penultimate <- cbind(Y_series_list[[i]][1:last_shock_point]
                               , post_shock_indicator)
      X_i_final <- X_i_penultimate[,2]
    }
    else {
      X_i_subset <- X[[i]][1:last_shock_point,covariate_indices]
      X_i_with_indicator <- cbind(X_i_subset, post_shock_indicator)
      X_i_final <- X_i_with_indicator
    }

    print('Now fitting the donor ARIMA models')

    arima <- auto.arima(Y_series_list[[i]][1:last_shock_point]
                        ,xreg=X_i_final
                        ,ic = user_ic_choice)

    print(arima)

    order_of_arima[[i]] <- arima$arma #tk

    coef_test <- coeftest(arima)
    extracted_fixed_effect <- coef_test[nrow(coef_test),1]
    omega_star_hat_vec <- c(omega_star_hat_vec, extracted_fixed_effect)

  } ## END loop for computing fixed effects

  ## END estimate fixed effects in donors

  ## BEGIN compute linear combination of fixed effects
  w_hat <- dbw(X, #tk
               dwb_indices,
               integer_shock_time_vec,
               scale = TRUE,
               sum_to_1 = TRUE, #tk
               bounded_below_by = 0, #tk
               bounded_above_by = 1, #tk
               # normchoice = normchoice, #tk
               # penalty_normchoice = penalty_normchoice,
               # penalty_lambda = penalty_lambda
  )

  omega_star_hat <- as.numeric(w_hat %*% omega_star_hat_vec)
  ## END compute linear combination of fixed effects

  ## BEGIN fit GARCH to target series

  if (is.null(covariate_indices) == TRUE){

    arima <- auto.arima(Y_series_list[[1]][1:integer_shock_time_vec[1]]
                        ,xreg = NULL
                        ,ic = user_ic_choice)

    unadjusted_pred <- predict(arima, n.ahead = shock_length_vec[1])

  }
  else{
    ## BEGIN fit GARCH to target series

    X_lagged <- lag.xts(X[[1]][1:integer_shock_time_vec[1],covariate_indices])

    arima <- auto.arima(Y_series_list[[1]][1:integer_shock_time_vec[1]]
                        ,xreg = X_lagged
                        ,ic = user_ic_choice)

    print(arima)

    #Note: for forecasting, we use last-observed X value
    X_to_use_in_forecast <- X[[1]][integer_shock_time_vec[1],covariate_indices]

    X_replicated_for_forecast_length <- matrix(rep(X_to_use_in_forecast, k)
                                               , nrow = shock_length_vec[1]
                                               , byrow = TRUE)

    forecast_period <- (integer_shock_time_vec[1]+1):(integer_shock_time_vec[1]+shock_length_vec[1])
    mat_X_for_forecast <- cbind(Y_series_list[[1]][forecast_period]
                                , X_replicated_for_forecast_length)

    unadjusted_pred <- predict(arima
                               , n.ahead = shock_length_vec[1]
                               , newxreg = mat_X_for_forecast[,-1])
  }


  ##We take care of housekeeping
  #tk
  order_of_arima[[1]] <- arima$arma

  print('now we print dataframe with orders...')
  order_matrix <- matrix(unlist(order_of_arima), byrow = TRUE, nrow = length(order_of_arima))

  print(order_matrix)

  if( length(unique(order_matrix[,2])) > 1 ) {
    message <- paste('NOT all series are I(', order_matrix[,2], ')')
    warning(message)
  }

  ##

  adjusted_pred <- unadjusted_pred$pred + omega_star_hat

  list_of_linear_combinations <- list(w_hat)
  list_of_forcasts <- list(unadjusted_pred, adjusted_pred)
  names(list_of_forecasts) <- c('unadjusted_pred', 'adjusted_pred')

  output_list <- list(list_of_linear_combinations
                      , list_of_forecasts)

  names(output_list) <- c('linear_combinations', 'predictions')

  ## tk OUTPUT
  cat('SynthVolForecast Details','\n',
      '-------------------------------------------------------------\n',
      'Donors:', n, '\n',
      'Shock times:', shock_time_vec, '\n',
      'Lengths of shock times:', shock_length_vec, '\n',
      'Convex combination',w_hat,'\n',
      'Shock estimates provided by donors:', omega_star_hat_vec, '\n',
      'Aggregate estimated shock effect:', omega_star_hat, '\n',
      'Actual change in stock price at T* + 1:', Y_series_list[[1]][integer_shock_time_vec[1]+1],'\n',
      'Adjusted forecasted change in stock price at T* + 1:', unadjusted_pred$pred,'\n',
      'MSE unadjusted:', (Y_series_list[[1]][integer_shock_time_vec[1]+1]-unadjusted_pred$pred)**2,'\n',
      'Adjusted forecasted change in stock price at T* + 1:', adjusted_pred,'\n',
      'MSE adjusted:', (Y_series_list[[1]][integer_shock_time_vec[1]+1]-adjusted_pred)**2,'\n'
  )

  ## PLOTS

  if (plots == TRUE){
    cat('User has opted to produce plots.','\n')
    plot_maker_synthprediction(Y_series_list
               ,shock_time_vec
               ,integer_shock_time_vec
               ,shock_length_vec
               ,unadjusted_pred$pred
               ,w_hat
               ,omega_star_hat
               ,adjusted_pred
               ,display_ground_truth = display_ground_truth_choice)

  }

  return(output_list)

} ### END SynthPrediction

