####################### BEGIN Auxiliary functions #######################

### START QL_loss_function
QL_loss_function <- function(pred, gt){gt/pred - log(gt/pred) - 1}
### END QL_loss_function

#We specify some functions for transforming y

#mean_square_y will be used in garch models
mean_square_y <- function(y){return((y-mean(y))**2)}

#identity function will be used in synthetic prediction models
id <- function(y){return(y)}

### START dbw
dbw <- function(X
                ,dbw_indices
                ,shock_time_vec
                ,scale = FALSE
                ,center = FALSE
                ,sum_to_1 = 1
                ,bounded_below_by = 0
                ,bounded_above_by = 1
                ,princ_comp_count = NULL
                ,normchoice = c('l1', 'l2')[2]
                ,penalty_normchoice = c('l1', 'l2')[1]
                ,penalty_lambda = 0
                ,Y = NULL
                ,Y_lookback_indices = list(seq(1,1,1))
                ,X_lookback_indices = rep(list(c(1)),length(dbw_indices))
                ,inputted_transformation
) { # https://github.com/DEck13/synthetic_prediction/blob/master/prevalence_testing/numerical_studies/COP.R
  # X is a list of covariates for the time series
  # X[[1]] should be the covariate of the time series to predict
  # X[[k]] for k = 2,...,n+1 are covariates for donors

  # number of time series for pool
  n <- length(X) - 1

  normchoice_number <- unlist(strsplit(normchoice, split = ""))[2]

  #We notify user if p > n, i.e. if linear system is overdetermined
  p <- length(dbw_indices)
  print('Here is the number of covariates used in dbw:')
  print(p)
  if (p > n){cat('p > n, i.e. system is overdetermined from an unconstrained point-of-view.')}

  ## We now perform the complicated task of grabbing a specified
  ## 'lookback' length for each i=1,2,...,n+1 and each covariate.

  col_returner <- function(df){return(df[,dbw_indices])}
  print('col_returner succeeded')

  X_subset1 <- lapply(X, col_returner)
  print('col_returner with lapply succeeded')

  #Task: for each entry in Y, make it the first column of X's corresponding entry
  if (is.null(Y_lookback_indices) == FALSE){

    print('User has provided Y_lookback_indices, so we include them.')
    X_lookback_indices <- c(Y_lookback_indices, X_lookback_indices)

    X_Y_combiner <- function(y,x) {

          print('We print the transformation and its class')
          print(inputted_transformation)
          print(class(inputted_transformation))

          transformed_series <- inputted_transformation(y) #tk

          return(merge(transformed_series,x, all = FALSE))

    } #end X_Y_combiner

    print('We are about to combine X_subset1 and Y...')
    print(class(Y))
    print(length(Y))
    print(class(X_subset1))
    print(length(X_subset1))

    combined_X <- mapply(X_Y_combiner
                         , y = Y
                         , x = X_subset1
                         , SIMPLIFY = FALSE)

  }
  else{combined_X <- X_subset1}

  row_returner <- function(df, stv){
    print(paste('Shock occurs at ', stv, sep = ''))
    print(paste('Row count of the series is ', nrow(df), sep = ''))
    return(df[1:(stv),])
  }

  X_subset2 <- mapply(row_returner, df = combined_X, stv = shock_time_vec, SIMPLIFY=FALSE)

  cov_extractor <- function(X_df){

        #Get the row count of X_df
        len <- nrow(X_df)

        #Function that maps a vector of indices to a padded vector
        # of length len, where the vector is TRUE at the indices and FALSE otherwise
          padded_vector_maker <- function(x)
            {
            vec <- rep(FALSE,len)
            vec[x] <- TRUE
            vec_reversed <- rev(vec)
            return(vec_reversed)
            }

          covariate_vals_in_list <- lapply(X_lookback_indices, padded_vector_maker)
          boolmat <- as.matrix(do.call(data.frame, covariate_vals_in_list))

          return(as.matrix(X_df)[boolmat])

  }

  X_subset <- lapply(X_subset2, cov_extractor)

  # Now bind the TSUS covariates to the donor covariates
  dat <- do.call('rbind', X_subset)
  print('Pre-scaling')
  print(dat)

  if (scale == TRUE) {cat('User has chosen to scale covariates.')}
  if (center == TRUE) {cat('User has chosen to center covariates.')}

    dat <- apply(dat, 2, function(x) scale(x, center = center, scale = scale))
    print('Post-scaling (nothing will happen if center and scale are set to FALSE).')
    print(dat)

    ## We output details of SVD
    dat.svd <- svd(dat)
    sing_vals <- dat.svd$d / sum(dat.svd$d)
    print('Singular value percentages for the donor pool X data:')
    print(paste(round(100 * sing_vals,2), "%", sep = ""))

    #Now project in direction of first princ_comp_count principal components
    if (is.null(princ_comp_count) == FALSE){
      print(paste('We are using ', princ_comp_count, ' principal components.', sep = ''))
      print(dat.svd$v)
      dat <- dat %*% dat.svd$v[,1:princ_comp_count]
    }

    X1 <- dat[1, , drop = FALSE]
    X0 <- split(dat[-1,],seq(nrow(dat[-1,]))) #tk what is split doing?

    print('We inspect X1.')
    print(X1)
    print('We inspect X0.')
    print(X0)

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
    norm_output <- as.numeric(norm(matrix(X1 - XW), type = normchoice_number))

    #now add penalty
    if (penalty_normchoice == 'l1' & penalty_lambda > 0) {
      norm_output <- norm_output + penalty_lambda * norm(as.matrix(W), type = "1")
    }
    else if (penalty_normchoice == 'l2' & penalty_lambda > 0) {
      norm_output <- norm_output + penalty_lambda * as.numeric(crossprod(matrix(W)))
    }
    else {norm_output <- norm_output}

    return(norm_output)
  } #end objective function

  # optimization and return statement

  # I have added features
  # 1) The option to remove the sum-to-1 constraint
  # 2) The option to change the lower bound to -1 or NA
  # 3) option to change the upper bound to NA'
  # 4) option to choose l1 or l2 norm as distance function

  #Thus I need if statements to implement these.

  # conditional for sum to 1
  if (is.na(sum_to_1) == FALSE) {eq_constraint <- function(W) sum(W) - 1}
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
    upper_bound = rep(bounded_above_by, n)
  }
  else if (is.na(bounded_above_by) == TRUE)  {
    upper_bound = NULL
  }

  object_to_return <- Rsolnp::solnp(par = rep(1/n, n),
                            fun = weightedX0,
                            eqfun = eq_constraint,
                            eqB = 0,
                            LB = lower_bound, UB = upper_bound,
                            control = list(trace = 1
                                           , 1.0e-12
                                           , tol = 1e-27
                                           , outer.iter = 1000000000
                                           , inner.iter = 10000000))

  #We print the loss from the optimization
  loss <- round(norm(X1 - object_to_return$pars %*% dat[-1,], type = normchoice_number),3)

  print(paste('The loss of distanced-based weighting is ', loss, ',',
              ' which is ',
              100*round(loss/norm(X1, type=normchoice_number),3),
              "% of the norm of the vector we are trying to approximate.", sep = ""))

  if (object_to_return$convergence == 0){convergence <- 'convergence'}
  else {convergence <- 'failed_convergence'}

  pair_to_return <- list(object_to_return$pars, convergence)

  names(pair_to_return) <- c('opt_params', 'convergence')

  return(pair_to_return)

} #END dbw function
### END dbw

### START GARCH plot_maker_garch
plot_maker_garch <- function(fitted_vol
                            ,shock_time_labels = NULL
                            ,shock_time_vec #mk
                            ,shock_length_vec
                            ,unadjusted_pred
                            ,w_hat
                            ,omega_star_hat
                            ,omega_star_hat_vec
                            ,omega_star_std_err_hat_vec
                            ,adjusted_pred
                            ,arithmetic_mean_based_pred
                            ,ground_truth_vec){

  if (is.character(shock_time_labels) == FALSE | is.null(shock_time_labels) == TRUE){
    shock_time_labels <- 1:length(shock_time_vec)
  }

  par(mfrow = c(1,3), mar=c(15,9,4,2))

  barplot_colors <- RColorBrewer::brewer.pal(length(w_hat),'Set3')

  #PLOT ON THE LEFT:
print('We plot the weights.')
  # Plot donor weights
  barplot(w_hat
          , main = 'Donor Pool Weights'
          , names.arg = shock_time_labels[-1]
          , cex.names=1.3
          , cex.main=1.5
          , las=2
          , col = barplot_colors
          )

  #PLOT IN THE MIDDLE
  print('We plot the FE estimates.')
  #Plot FE estimates
    bp <- barplot(omega_star_hat_vec
          , main = 'Donor-Pool-Supplied\n FE Estimates\nand Standard Errors Estimates'
          , names.arg = shock_time_labels[-1]
          , cex.names=1.4
          , cex.main=1.5
          , las=2
          , col = barplot_colors
          , ylim = c(0, 1.4 * max(omega_star_hat_vec)) )

    # Add the labels with some offset to be above the bar
    print('We print the std errors')
    print(omega_star_std_err_hat_vec)
    #omega_star_std_err_hat_vec <- ifelse(omega_star_std_err_hat_vec, nan)

    #https://stackoverflow.com/questions/65057352/how-to-add-labels-above-the-bar-of-barplot-graphics
    text(x = bp,
                ,y = omega_star_hat_vec + .00029
                ,cex = 1.3
                ,labels = round(omega_star_std_err_hat_vec, 5)
                , srt= 90)

  title(ylab = expression(sigma^2), line = 3.05, cex.lab = 1.99) # Add y-axis text

  #Plot target series and prediction

  thing_to_get_max_of <- c(as.numeric(fitted_vol)
                        , unadjusted_pred
                        , adjusted_pred
                        , ground_truth_vec
                        , arithmetic_mean_based_pred
                        )

  max_for_y_lim <- max(thing_to_get_max_of)

  x_ax_first_point_of_shock <- index(fitted_vol)[shock_time_vec[1]-1] + 1 #do I use this?
  x_ax_end_point <- index(fitted_vol)[shock_time_vec[1]-1] + length(adjusted_pred)

  #PLOT ON THE RIGHT:
  print('We plot the fitted volatility series.')
  plot(y = fitted_vol[1:shock_time_vec[1]], #mk
          x = index(fitted_vol)[1:shock_time_vec[1]],
       main = 'Post-Shock Volatility Forecast', #mk can improve this title
       cex.main=1.5,
       ylab = '',
       type="l",
       xlab = '',
       xlim =  as.Date(c(index(fitted_vol)[1], x_ax_end_point)),
       ylim = c(min(0, fitted_vol),  max_for_y_lim))

  title(ylab = expression(sigma^2), line = 2.05, cex.lab = 1.99) # Add y-axis text

  # Here is the color scheme we will use
  # https://colorbrewer2.org/?type=diverging&scheme=RdYlBu&n=4
  colors_for_adjusted_pred <- c('#d7191c','#fdae61','#abd9e9','#2c7bb6')

  # Let's add the plain old GARCH prediction
  points(y = unadjusted_pred
         ,x = x_ax_first_point_of_shock:x_ax_end_point
         ,col = colors_for_adjusted_pred[1]
         ,cex = 3.5
         ,pch = 15)

  # Now plot the adjusted predictions
  points(y = adjusted_pred
         ,x = x_ax_first_point_of_shock:x_ax_end_point
         ,col = colors_for_adjusted_pred[2]
         ,cex = 3.5
         ,pch = 18)

  # Now plot the arithmetic mean-based predictions
  points(y = arithmetic_mean_based_pred
         ,x = x_ax_first_point_of_shock:x_ax_end_point
         ,col = colors_for_adjusted_pred[3]
         ,cex = 3.5
         ,pch = 19)

  # Now plot Ground Truth tk
  if (is.null(ground_truth_vec) == FALSE)
    {
    points(y = ground_truth_vec
           ,x = x_ax_first_point_of_shock:x_ax_end_point
           ,col = colors_for_adjusted_pred[4]
           ,cex = 3.5
           ,pch = 17)
  }

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

  #par(mfrow = c(1,1), mar=c(15,6,4,2))

}
### END plot_maker_garch

### START plot_maker_synthprediction
plot_maker_synthprediction <- function(Y
                                     ,shock_time_labels = NULL
                                     ,shock_time_vec #mk
                                     ,shock_length_vec
                                     ,unadjusted_pred
                                     ,w_hat
                                     ,omega_star_hat
                                     ,omega_star_hat_vec
                                     ,adjusted_pred
                                     ,display_ground_truth = FALSE){


  if (is.character(shock_time_labels) == FALSE | is.null(shock_time_labels) == TRUE){
    shock_time_labels <- 1:length(shock_time_vec)
  }

  n <- length(Y) - 1

  #First print donor series
  par(mfrow = c(round(sqrt(n)),ceiling(sqrt(n))))

  for (i in 2:(n+1)){
    plot.ts(Y[[i]][1:shock_time_vec[i]]
            ,xlab = ' '
            ,ylab = 'Log Return'
            ,main = paste('Donor ', i,': ', shock_time_labels[i], sep = '')
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
             # ,col = 'red'
             ,cex = 1.1
             ,pch = 17)

    }
  }

  #Now print time series under study
  par(mfrow = c(1,3), mar=c(15,6,4,2))

  barplot_colors <- RColorBrewer::brewer.pal(length(w_hat),'Set3')

  #PLOT ON THE LEFT:
  #Plot donor weights
  barplot(w_hat
          , main = 'Donor Pool Weights'
          , names.arg = shock_time_labels[-1]
          , cex.names=.95
          , las=2
          , col = barplot_colors)

  #PLOT IN THE MIDDLE

  #Plot FE estimates
  barplot(omega_star_hat_vec
          , main = 'Donor-Pool-Supplied \n FE Estimates'
          , names.arg = shock_time_labels[-1]
          , cex.names=.95
          , las=2
          , col = barplot_colors)

  #Plot target series and prediction

  thing_to_get_max_of <- c(as.numeric(Y[[1]]), unadjusted_pred, adjusted_pred)

  max_for_y_lim <- max(thing_to_get_max_of)

  #PLOT ON THE RIGHT:
  plot.ts(Y[[1]][1:shock_time_vec[1]], #mk
          main = 'Post-shock Forecasts',
          ylab = '',
          xlab = ' ',
          xlim = c(0, shock_time_vec[1] + 5), #mk
          ylim = c(min(0, Y[[1]]),  max_for_y_lim))

  title(ylab = 'Log Return', line = 2.05, cex.lab = 1.99) # Add y-axis text

  # Here is the color scheme we will use
  #https://colorbrewer2.org/?type=diverging&scheme=RdYlBu&n=4
  colors_for_adjusted_pred <- c('#d7191c','#fdae61','#abd9e9')

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
         ,cex = 2
         ,pch = 19)

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

### START plot_maker_HAR
plot_maker_HAR <- function(Y
                           ,shock_time_labels = NULL
                           ,shock_time_vec #mk
                           ,shock_length_vec = 1 #tk
                           ,unadjusted_pred
                           ,w_hat
                           ,omega_star_hat
                           ,omega_star_hat_vec
                           #,omega_star_hat_vec_std_err
                           ,adjusted_pred
                           ,display_ground_truth = FALSE){

  n <- length(shock_time_vec) - 1

  if (is.character(shock_time_labels) == FALSE | is.null(shock_time_labels) == TRUE){
    shock_time_labels <- 1:(n+1)
  }

  par(mfrow = c(round(sqrt(n)),ceiling(sqrt(n))))

  for (i in 2:(n+1)){

    plot.ts(Y[[i]][1:shock_time_vec[i]]
            ,xlab = ' '
            ,ylab = 'Realized Measure of Volatility'
            ,main = paste('Donor ', i,': ', shock_time_labels[i], sep = '')
            ,xlim = c(0, shock_time_vec[i] + 5)
            ,ylim = c(min(Y[[i]]),  max(Y[[i]]))
    )

    if (display_ground_truth == TRUE){

      lines(y = Y[[i]][shock_time_vec[i]:(shock_time_vec[i] + shock_length_vec[i])]
            ,x = shock_time_vec[i]:(shock_time_vec[i] + shock_length_vec[i])
            ,col = 'purple'
            ,cex = 1.1
            ,lty = 3)

      #tk points not working?

      points(y = Y[[i]][(shock_time_vec[i]+1):(shock_time_vec[i] + shock_length_vec[i])]
             ,x = (shock_time_vec[i]+1):(shock_time_vec[i] + shock_length_vec[i])
             # ,col = 'red'
             ,cex = 1.1
             ,pch = 17)

    }
  }

  #Now print time series under study
  par(mfrow = c(1,3), mar=c(15,6,4,2))

  barplot_colors <- RColorBrewer::brewer.pal(length(w_hat),'Set3')

  #PLOT ON THE LEFT:
  #Plot donor weights
  barplot(w_hat
          , main = 'Donor Pool Weights'
          , names.arg = shock_time_labels[-1]
          , cex.names=.95
          , las=2
          , col = barplot_colors)

  #PLOT IN THE MIDDLE

  #Plot FE estimates
  barplot(omega_star_hat_vec
          , main = 'Donor-Pool-Supplied \n FE Estimates'
          , names.arg = shock_time_labels[-1]
          , cex.names=.95
          , las=2
          , col = barplot_colors)

  #Plot target series and prediction

  Y_to_plot <- Y[[1]][(shock_time_vec[1]-45):(shock_time_vec[1]-1)] #tk 30?

  x_ax_first_point_of_shock <- index(Y[[1]])[(shock_time_vec[1]-45)] #tk rename
  x_ax_end_point <- index(Y[[1]])[shock_time_vec[1]] + 5 #tk why 5

  thing_to_get_max_of <- c(as.numeric(Y[[1]]) #tk should we shorten time frame?
                           , unadjusted_pred
                           , adjusted_pred)

  max_for_y_lim <- max(thing_to_get_max_of)

  #PLOT ON THE RIGHT:
  plot(y = Y_to_plot, #mk
       x = index(Y_to_plot),
       main = 'Post-Shock Volatility Forecast', #mk can improve this title
       cex.main=1.5,
       ylab = '',
       type="l",
       xlab = '',
       xlim =  as.Date(c(x_ax_first_point_of_shock, x_ax_end_point)),
       ylim = c(min(0, Y_to_plot),  max_for_y_lim)
       )

  title(ylab = 'Realized Measure of Volatility', line = 2.05, cex.lab = 1.99) # Add y-axis text

  # Here is the color scheme we will use
  #https://colorbrewer2.org/?type=diverging&scheme=RdYlBu&n=4
  colors_for_adjusted_pred <- c('#d7191c','#fdae61','#abd9e9')

  # Let's add the HAR prediction
  points(y = unadjusted_pred
         ,x =  index(Y[[1]])[shock_time_vec[1]]
         ,col = colors_for_adjusted_pred[1]
         ,cex = 2
         ,pch = 19)

  # Now plot the adjusted predictions
  points(y = adjusted_pred
         ,x =  index(Y[[1]])[shock_time_vec[1]]
         ,col = colors_for_adjusted_pred[2]
         ,cex = 2
         ,pch = 19)

  if (display_ground_truth == TRUE){ #tk what is this doing? printing the actual, if it exists

    Y_not_yet_plotted <- Y[[1]][shock_time_vec[1]:(shock_time_vec[1]+1)] #tk need to fix
    connecting_the_two <- Y[[1]][(shock_time_vec[1]-1):shock_time_vec[1]]


    lines(y = connecting_the_two
          ,x = index(connecting_the_two)
          ,col = 'black'
          ,cex = 1.1
          ,lty = 3
          ,lwd = 3)

    points(y = Y_not_yet_plotted[1:length(adjusted_pred)]
           ,x = index(Y_not_yet_plotted)[1:length(adjusted_pred)]
           ,col = colors_for_adjusted_pred[3]
           ,cex = 2
           ,pch = 19)
  }

  labels_for_legend <- c('HAR (unadjusted)', 'Adjusted HAR Prediction', 'Actual')

  legend(x = "topleft",  # Coordinates (x also accepts keywords) #mk
         legend = labels_for_legend,
         1:length(labels_for_legend), # Vector with the name of each group
         colors_for_adjusted_pred,   # Creates boxes in the legend with the specified colors
         title = 'Prediction Method',      # Legend title,
         cex = 1.3)

}
### END plot_maker_HAR

####################### END Auxiliary functions #######################

### START SynthPrediction
SynthPrediction <- function(Y_series_list
                            ,covariates_series_list
                            ,shock_time_vec
                            ,shock_length_vec
                            ,k = 1
                            ,dbw_scale = TRUE
                            ,dbw_center = TRUE
                            ,dbw_indices = NULL
                            ,princ_comp_input = min(length(shock_time_vec), ncol(covariates_series_list[[1]]))
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
  print('The indices used for distance-based-weight are ')
  print(dbw_indices)
  if (is.null(dbw_indices) == TRUE) {
    dbw_indices <- 1:ncol(covariates_series_list[[1]]) #tk
  }
  ### END Populate defaults
  ## BEGIN Check that inputs are all comformable/acceptable
  n <- length(Y_series_list) - 1 #tk
  ## END Check that inputs are all comformable/acceptable
  integer_shock_time_vec <- c() #mk
  integer_shock_time_vec_for_convex_hull_based_optimization <- c() #mk
  ## BEGIN Check whether shock_time_vec is int/date
  for (i in 1:(n+1)){
    if (is.character(shock_time_vec[i]) == TRUE){
      integer_shock_time_vec[i] <- which(index(Y[[i]]) == shock_time_vec[i]) #mk
      integer_shock_time_vec_for_convex_hull_based_optimization[i] <- which(index(covariates_series_list[[i]]) == shock_time_vec[i]) #mk
    }
    else{
      integer_shock_time_vec[i] <- shock_time_vec[i]
      integer_shock_time_vec_for_convex_hull_based_optimization[i] <- shock_time_vec[i]
    }
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
      X_i_subset <- covariates_series_list[[i]][1:last_shock_point,covariate_indices]
      X_i_with_indicator <- cbind(X_i_subset, post_shock_indicator)
      X_i_final <- X_i_with_indicator
    }
    print('Now fitting the donor ARIMA models')
    arima <- forecast::auto.arima(Y_series_list[[i]][1:last_shock_point]
                                  ,xreg=X_i_final
                                  ,ic = user_ic_choice)
    print(arima)
    order_of_arima[[i]] <- arima$arma #tk
    coef_test <- lmtest::coeftest(arima)
    extracted_fixed_effect <- coef_test[nrow(coef_test),1]
    omega_star_hat_vec <- c(omega_star_hat_vec, extracted_fixed_effect)
  } ## END loop for computing fixed effects
  ## END estimate fixed effects in donors
  ## BEGIN compute linear combination of fixed effects
  dbw_output <- dbw(covariates_series_list, #tk
                    dbw_indices,
                    integer_shock_time_vec,
                    scale = TRUE,
                    center = TRUE,
                    sum_to_1 = TRUE, #tk
                    bounded_below_by = 0, #tk
                    bounded_above_by = 1, #tk
                    # normchoice = normchoice, #tk
                    # penalty_normchoice = penalty_normchoice,
                    # penalty_lambda = penalty_lambda
                    Y = Y_series_list,
                    inputted_transformation = id
  )
  w_hat <- dbw_output[[1]]
  omega_star_hat <- as.numeric(w_hat %*% omega_star_hat_vec)
  ## END compute linear combination of fixed effects
  ## BEGIN fit GARCH to target series
  if (is.null(covariate_indices) == TRUE){
    arima <- forecast::auto.arima(Y_series_list[[1]][1:(integer_shock_time_vec[1])]
                                  ,xreg = NULL
                                  ,ic = user_ic_choice)
    unadjusted_pred <- predict(arima, n.ahead = shock_length_vec[1])
  }
  else{
    ## BEGIN fit GARCH to target series
    X_lagged <- lag.xts(covariates_series_list[[1]][1:integer_shock_time_vec[1],covariate_indices])
    arima <- forecast::auto.arima(Y_series_list[[1]][1:integer_shock_time_vec[1]]
                                  ,xreg = X_lagged
                                  ,ic = user_ic_choice)
    print(arima)
    #Note: for forecasting, we use last-observed X value
    X_to_use_in_forecast <- covariates_series_list[[1]][integer_shock_time_vec[1],covariate_indices]
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
  list_of_forecasts <- list(unadjusted_pred, adjusted_pred)
  names(list_of_forecasts) <- c('unadjusted_pred', 'adjusted_pred')
  output_list <- list(list_of_linear_combinations
                      , list_of_forecasts)
  names(output_list) <- c('linear_combinations', 'predictions')
  ## tk OUTPUT
  cat('SynthPrediction Details','\n',
      '-------------------------------------------------------------\n',
      'Donors:', n, '\n',
      'Shock times:', shock_time_vec, '\n',
      'Lengths of shock times:', shock_length_vec, '\n',
      'Optimization Success:', dbw_output[[2]], '\n', '\n',
      'Convex combination',w_hat,'\n',
      'Shock estimates provided by donors:', omega_star_hat_vec, '\n',
      'Aggregate estimated shock effect:', omega_star_hat, '\n',
      'Actual change in stock price at T* + 1:', Y_series_list[[1]][integer_shock_time_vec[1]+1],'\n',
      'Unadjusted forecasted change in stock price at T*+1:', unadjusted_pred$pred,'\n',
      'MSE unadjusted:', (as.numeric(Y_series_list[[1]][integer_shock_time_vec[1]+1])-unadjusted_pred$pred)**2,'\n',
      'Adjusted forecasted change in stock price at T*+1:', adjusted_pred,'\n',
      'MSE adjusted:', (as.numeric(Y_series_list[[1]][integer_shock_time_vec[1]+1])-adjusted_pred)**2,'\n'
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
                               ,omega_star_hat_vec
                               ,adjusted_pred
                               ,display_ground_truth = display_ground_truth_choice
    )
  }
  return(output_list)
} ### END SynthPrediction

### START SynthVolForecast
SynthVolForecast <- function(Y_series_list
                             ,covariates_series_list
                             ,shock_time_vec
                             ,shock_length_vec
                             ,k=1
                             ,dbw_scale = TRUE
                             ,dbw_center = TRUE
                             ,dbw_indices = NULL
                             ,dbw_Y_lookback = c(0)
                             ,dbw_princ_comp_input = NULL
                             ,covariate_indices = NULL
                             ,geometric_sets = NULL #tk
                             ,days_before_shocktime_vec = NULL #tk I may want to remove this
                             ,garch_order = NULL
                             ,common_series_assumption = FALSE
                             ,plots = TRUE
                             ,shock_time_labels = NULL
                             ,ground_truth_vec = NULL
                             ,Y_lookback_indices_input = list(seq(1,30,1))
                             ,X_lookback_indices_input = rep(list(c(1)),length(dbw_indices))
){
  ### BEGIN Doc string
  #tk
  ### END Doc string

  ### BEGIN Populate defaults
  n <- length(Y_series_list) - 1

  if (is.null(garch_order) == TRUE) {

    garch_order <- list(length = n+1)

    for (i in 1:(n+1)){garch_order[[i]] <- c(1,1,0)}
}

  if (is.null(dbw_indices) == TRUE) {dbw_indices <- 1:ncol(covariates_series_list[[1]])}

  ### END Populate defaults

  ## BEGIN Check that inputs are all comformable/acceptable
  n <- length(Y_series_list) - 1 #tk
  ## END Check that inputs are all comformable/acceptable

  integer_shock_time_vec <- c() #mk
  integer_shock_time_vec_for_convex_hull_based_optimization <- c() #mk

  ## BEGIN Check whether shock_time_vec is int/date

  for (i in 1:(n+1)){

    if (is.character(shock_time_vec[i]) == TRUE){
      integer_shock_time_vec[i] <- which(index(Y[[i]]) == shock_time_vec[i]) #mk
      integer_shock_time_vec_for_convex_hull_based_optimization[i] <- which(index(covariates_series_list[[i]]) == shock_time_vec[i]) #mk
    }
    else{
      integer_shock_time_vec[i] <- shock_time_vec[i]
      integer_shock_time_vec_for_convex_hull_based_optimization[i] <- shock_time_vec[i]
    }

  }

  ## END Check whether shock_time_vec is int/date

  ## BEGIN calculate weight vector
  dbw_output <- dbw(covariates_series_list, #tk
                    dbw_indices,
                    integer_shock_time_vec_for_convex_hull_based_optimization,
                    scale = dbw_scale,
                    center = dbw_center,
                    sum_to_1 = TRUE, #tk
                    princ_comp_count = dbw_princ_comp_input,
                    bounded_below_by = 0, #tk
                    bounded_above_by = 1, #tk
                    # normchoice = normchoice, #tk
                    # penalty_normchoice = penalty_normchoice,
                    # penalty_lambda = penalty_lambda
                    Y = Y_series_list,
                    Y_lookback_indices = Y_lookback_indices_input,
                    X_lookback_indices = X_lookback_indices_input,
                    inputted_transformation = mean_square_y
  )

  w_hat <- dbw_output[[1]]

  ## END calculate weight vector

  ## BEGIN estimate fixed effects in donors
  omega_star_hat_vec <- c()
  omega_star_std_err_hat_vec <- c()

  #tk
  if (common_series_assumption == TRUE){
    print('tk TODO')

    #step 1: create dummy vector with n+1 shocks
    #NOTA BENE: n different fixed effects, or
    # 1 fixed effect estimated at n shocks?
    vec_of_zeros <- rep(0, integer_shock_time_vec[i])
    vec_of_ones <- rep(1, shock_length_vec[i])
    vec_of_final_zeros <- rep(0, 20) #tk
    post_shock_indicator <- c(vec_of_zeros, vec_of_ones)
    last_shock_point <- integer_shock_time_vec[i] + shock_length_vec[i] + 20 #tk

    #step 2: fit model

  }

  else{

    for (i in 2:(n+1)){

      # Make indicator variable w/ a 1 at only T*+1, T*+2,...,T*+shock_length_vec[i]
      vec_of_zeros <- rep(0, integer_shock_time_vec[i])
      vec_of_ones <- rep(1, shock_length_vec[i])
      vec_of_final_zeros <- rep(0, 20)
      post_shock_indicator <- c(vec_of_zeros, vec_of_ones, vec_of_final_zeros)
      last_shock_point <- integer_shock_time_vec[i] + shock_length_vec[i] + 20 #tk

      #subset X_i
      if (is.null(covariate_indices) == TRUE) {
        X_i_penultimate <- cbind(Y_series_list[[i]][1:last_shock_point] #tk
                                 , post_shock_indicator)
        X_i_final <- X_i_penultimate[,2]
      }
      else {
        X_i_subset <- covariates_series_list[[i]][1:last_shock_point,covariate_indices]
        X_i_with_indicator <- cbind(X_i_subset, post_shock_indicator)
        X_i_final <- X_i_with_indicator
      }

      print('We print the tail of the covariate df we use in GARCH model:')
      print(tail(X_i_final))

      fitted_garch <- garchx::garchx(Y_series_list[[i]][1:last_shock_point] #tk
                                     , order = garch_order[[i]]
                                     , xreg = X_i_final
                                     , backcast.values = NULL
                                     , control = list(eval.max = 100000
                                                      , iter.max = 1500000
                                                      , rel.tol = 1e-8))

      cat('\n===============================================================\n')
      print(paste('Outputting GARCH estimates for donor series number ', i,'.', sep = ''))
      print(fitted_garch)
      print(paste('Outputting AIC for donor series number ', i,'.', sep = ''))
      print(AIC(fitted_garch))
      cat('\n===============================================================\n')

      coef_test <- lmtest::coeftest(fitted_garch)
      extracted_fixed_effect <- coef_test[dim(lmtest::coeftest(fitted_garch))[1], 1]
      extracted_fixed_effect_std_err <- coef_test[dim(lmtest::coeftest(fitted_garch))[1], 2]

      omega_star_hat_vec <- c(omega_star_hat_vec, extracted_fixed_effect)
      omega_star_std_err_hat_vec <- c(omega_star_std_err_hat_vec, extracted_fixed_effect_std_err)

    } ## END loop for computing fixed effects

  }

  ## END estimate fixed effects in donors

  ## BEGIN compute linear combination of fixed effects
  omega_star_hat <- w_hat %*% omega_star_hat_vec
  ## END compute linear combination of fixed effects

  ## BEGIN fit GARCH to target series

  if (is.null(covariate_indices) == TRUE){

    fitted_garch <- garchx::garchx(Y_series_list[[1]][1:(integer_shock_time_vec[1])]
                                   , order = garch_order[[1]]
                                   , xreg = NULL
                                   , backcast.values = NULL
                                   , control = list(eval.max = 100000
                                                    , iter.max = 1500000
                                                    , rel.tol = 1e-8))

    cat('\n===============================================================\n')
    print('Outputting the fitted GARCH for time series under study.')
    print(fitted_garch)
    print('Outputting AIC for time series under study.')
    print(AIC(fitted_garch))
    cat('\n===============================================================\n')

    unadjusted_pred <- predict(fitted_garch, n.ahead = shock_length_vec[1])
  }
  else{
    ## BEGIN fit GARCH to target series
    fitted_garch <- garchx::garchx(Y_series_list[[1]][1:(integer_shock_time_vec[1])]
                                   , order = garch_order
                                   , xreg = covariates_series_list[[1]][1:(integer_shock_time_vec[1]),covariate_indices]
                                   , backcast.values = NULL
                                   , control = list(eval.max = 100000
                                                    , iter.max = 1500000
                                                    , rel.tol = 1e-8))

    cat('\n===============================================================\n')
    print('Outputting the fitted GARCH for the time series under study.')
    print(fitted_garch)
    cat('\n===============================================================\n')

    #Note: for forecasting, we use last-observed X value
    X_to_use_in_forecast <- covariates_series_list[[1]][integer_shock_time_vec[1],covariate_indices]

    X_replicated_for_forecast_length <- matrix(rep(X_to_use_in_forecast, k)
                                               , nrow = shock_length_vec[1]
                                               , byrow = TRUE)

    forecast_period <- (integer_shock_time_vec[1]):(integer_shock_time_vec[1]+shock_length_vec[1])
    mat_X_for_forecast <- cbind(Y_series_list[[1]][forecast_period]
                                , X_replicated_for_forecast_length)

    unadjusted_pred <- predict(fitted_garch
                               , n.ahead = shock_length_vec[1]
                               , newxreg = mat_X_for_forecast[,-1])
  }

  print('Now we get the adjusted predictions.')
  adjusted_pred <- unadjusted_pred + rep(omega_star_hat, k)

  arithmetic_mean_based_pred <- rep(mean(omega_star_hat_vec), k) + unadjusted_pred

  if (is.null(ground_truth_vec) == TRUE){
    QL_loss_unadjusted_pred <- NA
    QL_loss_adjusted_pred <- NA
  }
  else {
    QL_loss_unadjusted_pred <- sum(QL_loss_function(unadjusted_pred
                                                    , ground_truth_vec))
    QL_loss_adjusted_pred <- sum(QL_loss_function(adjusted_pred
                                                  , ground_truth_vec))
    QL_loss_arithmetic_mean_based_pred <- sum(QL_loss_function(arithmetic_mean_based_pred
                                                   , ground_truth_vec))
  }

  list_of_linear_combinations <- list(w_hat)

  list_of_forecasts <- list(unadjusted_pred
                            , adjusted_pred
                            , arithmetic_mean_based_pred)

  list_of_losses <- list(QL_loss_unadjusted_pred
                         , QL_loss_adjusted_pred
                         , QL_loss_arithmetic_mean_based_pred)

  names(list_of_forecasts) <- c('unadjusted_pred'
                                , 'adjusted_pred'
                                , 'arithmetic_mean_based_pred')

  names(list_of_losses) <- c('unadjusted_pred'
                                , 'adjusted_pred'
                                , 'arithmetic_mean_based_pred')

  output_list <- list(list_of_linear_combinations
                      , list_of_forecasts
                      , list_of_losses)

  names(output_list) <- c('linear_combinations'
                          , 'predictions'
                          , 'loss')

  ## tk OUTPUT
  cat('--------------------------------------------------------------\n',
      '-------------------SynthVolForecast Results-------------------','\n',
      '--------------------------------------------------------------\n',
      'Donors:', n, '\n',  '\n',
      'Shock times:', shock_time_vec, '\n', '\n',
      'Lengths of shock times:', shock_length_vec, '\n', '\n',
      'Optimization Success:', dbw_output[[2]], '\n', '\n',
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
                     ,shock_time_labels
                     ,integer_shock_time_vec
                     ,shock_length_vec
                     ,unadjusted_pred
                     ,w_hat
                     ,omega_star_hat
                     ,omega_star_hat_vec
                     #,omega_star_std_err_hat_vec
                     ,adjusted_pred
                     ,arithmetic_mean_based_pred
                     ,ground_truth_vec)
  }

  return(output_list)

} ### END SynthVolForecast

### START HAR
HAR <- function(Y
                ,covariates_series_list
                ,shock_time_vec
                ,shock_length_vec
                ,k=1 #does it make sense for a k-step ahead?
                ,dbw_scale = TRUE
                ,dbw_center = TRUE
                ,dbw_indices = 1:ncol(covariates_series_list[[1]])
                ,dbw_Y_lookback = c(0)
                ,dbw_princ_comp_input = NULL
                ,covariate_indices = NULL
                ,geometric_sets = NULL #tk
                ,plots = TRUE
                ,shock_time_labels = NULL
                ,ground_truth_vec = NULL
                ,Y_lookback_indices_input = list(seq(1,3,1))
                ,X_lookback_indices_input = rep(list(c(1)),length(dbw_indices))
){
  ### BEGIN Doc string

  # Input:
  #   Y_series_list - a list of df, where in each df, realized measure is first col
  # ,covariates_series_list
  # ,shock_time_vec
  # ,shock_length_vec
  # ,k=1
  # ,dbw_scale = TRUE
  # ,dbw_center = TRUE
  # ,dbw_indices = NULL
  # ,dbw_Y_lookback = c(0)
  # ,dbw_princ_comp_input = NULL
  # ,covariate_indices = NULL
  # ,geometric_sets = NULL #tk
  # ,days_before_shocktime_vec = NULL #tk I may want to remove this
  # ,garch_order = NULL
  # ,plots = TRUE
  # ,shock_time_labels = NULL
  # ,ground_truth_vec = NULL
  # ,Y_lookback_indices_input = list(seq(1,3,1))
  # ,X_lookback_indices_input = rep(list(c(1)),length(dbw_indices))

  ### END Doc string

  #Big thing to consider: what needs to change when Y is a just a df?

  n <- length(X) - 1 #tk

  integer_shock_time_vec <- c() #mk
  integer_shock_time_vec_for_convex_hull_based_optimization <- c() #mk

  omega_star_hat_vec <- c()

  if (is.data.frame(Y) == TRUE){

    print('Y is a dataframe')

    ## BEGIN Check whether shock_time_vec is int/date

    for (i in 1:length(shock_time_vec)){

      print(paste('The ', i, 'th shock time is ', shock_time_vec[i], sep = ''))

      if (is.character(shock_time_vec[i]) == TRUE){
        integer_shock_time_vec[i] <- which(row.names(Y) == shock_time_vec[i]) #mk
        integer_shock_time_vec_for_convex_hull_based_optimization[i] <- which(index(covariates_series_list[[i]]) == shock_time_vec[i]) #mk
      }
      else{
        integer_shock_time_vec[i] <- shock_time_vec[i]
        integer_shock_time_vec_for_convex_hull_based_optimization[i] <- shock_time_vec[i]
      }

    }

    ## END Check whether shock_time_vec is int/date

    shock_dates_as_dates <- as.Date(as.Date(unlist(shock_dates)))

    shock_dates_as_dates_without_TSUS <- shock_dates_as_dates[-1]

    Y_with_donor_col <- data.frame(Y %>% mutate(donor = ifelse(row.names(RVSPY_final) %in% shock_dates_as_dates_without_TSUS,1,0)))

    print('Here are the rows with 1 in the indicator variable:')
    print(Y_with_donor_col[Y_with_donor_col$donor == 1,])

    Y_with_donor_col[Y_with_donor_col$donor == 1, "donor"] <- shock_dates_as_dates_without_TSUS

    Y_with_donor_col$donor = as.factor(Y_with_donor_col$donor)

    training_period <- Y_with_donor_col[row.names(Y_with_donor_col) <= as.Date(shock_dates_as_dates[1]) , ] #tk

    print('We inspect the last few rows of the training period:')
    print(tail(training_period))

    m1 <- lm(training_period[,1] ~. , data = training_period[,-c(1)])

    cat('We inspect the fitted linear model.')
    print(summary(m1))

    forecast_period <- Y_with_donor_col[row.names(Y_with_donor_col) == as.Date(shock_dates_as_dates[1]), ] #tk

    outcome <- forecast_period[,1]

    print('Inspect the forecast period:')
    print(head(forecast_period))

    forecast_period_w_outcome_dropped <- forecast_period[,-c(1)] #drop the outcome var

    print('Make the donor variable a factor:')
    forecast_period_w_outcome_dropped$donor <- as.factor(forecast_period_w_outcome_dropped$donor)

    unadjusted_pred <- predict(m1, newdata = forecast_period_w_outcome_dropped)

    no_events <- length(shock_dates_as_dates) - 1
    no_coef <- length(coef(m1))

    omega_star_hat_vec <- c(omega_star_hat_vec, coef(m1)[(no_coef-no_events+1):no_coef])

  }
  else{

    print('Y is not a dataframe.')

    for (i in 1:length(shock_time_vec)){

      print(paste('The ', i, 'th shock time is ', shock_time_vec[i], sep = ''))

      if (is.character(shock_time_vec[i]) == TRUE){
        integer_shock_time_vec[i] <- which(row.names(Y) == shock_time_vec[i]) #mk
        integer_shock_time_vec_for_convex_hull_based_optimization[i] <- which(index(covariates_series_list[[i]]) == shock_time_vec[i]) #mk
      }
      else{
        integer_shock_time_vec[i] <- shock_time_vec[i]
        integer_shock_time_vec_for_convex_hull_based_optimization[i] <- shock_time_vec[i]
      }

    } #end loop

    for (i in 2:(n+1)){

      # Make indicator variable w/ a 1 at only T*+1, T*+2,...,T*+shock_length_vec[i]
      vec_of_zeros <- rep(0, integer_shock_time_vec[i])
      vec_of_ones <- rep(1, shock_length_vec[i])
      post_shock_indicator <- c(vec_of_zeros, vec_of_ones)
      last_shock_point <- integer_shock_time_vec[i] + shock_length_vec[i]

      #subset X_i
      if (is.null(covariate_indices) == TRUE) {
        X_i_penultimate <- cbind(Y[[i]][1:last_shock_point] #tk
                                 , post_shock_indicator)
        X_i_final <- X_i_penultimate[,2]
      }
      else {
        X_i_subset <- covariates_series_list[[i]][1:last_shock_point,covariate_indices]
        X_i_with_indicator <- cbind(X_i_subset, post_shock_indicator)
        X_i_final <- X_i_with_indicator
      }

      print('We print the tail of the covariate df we use in the model:')
      print(tail(X_i_final))

      omega_star_hat_vec <- c(omega_star_hat_vec, HAR_lm)

    } #end loop?

  }

  print('We print the integers of the shock times.')
  print(integer_shock_time_vec)

  print('We print the booleans for the shock times.')
  print(integer_shock_time_vec_for_convex_hull_based_optimization) #tk seems wrong

  print('Here is the Y data we will use in dbw:')
  make_xts <- xts(Y[,1], order.by = as.Date(row.names(Y))) #tk

  colnames(make_xts) <- c('User_chosen_Y')

  Y_list <- rep(list(make_xts), n+1) #tk as.Date here is inelegant

  ## BEGIN calculate weight vector
  dbw_output <- dbw(covariates_series_list, #tk
                    dbw_indices,
                    integer_shock_time_vec_for_convex_hull_based_optimization,
                    scale = dbw_scale,
                    center = dbw_center,
                    sum_to_1 = TRUE, #tk
                    princ_comp_count = dbw_princ_comp_input,
                    bounded_below_by = 0, #tk
                    bounded_above_by = 1, #tk
                    # normchoice = normchoice, #tk
                    # penalty_normchoice = penalty_normchoice,
                    # penalty_lambda = penalty_lambda
                    Y = Y_list,
                    Y_lookback_indices = Y_lookback_indices_input,
                    X_lookback_indices = X_lookback_indices_input,
                    inputted_transformation = id
  )
  ## END calculate weight vector

  ## BEGIN compute linear combination of fixed effects
  w_hat <- dbw_output[[1]]
  omega_star_hat <- w_hat %*% omega_star_hat_vec
  ## END compute linear combination of fixed effects

  ## Forecasts and loss

  arithmetic_mean_based_pred <- mean(omega_star_hat_vec) + unadjusted_pred

  adjusted_pred <- unadjusted_pred + omega_star_hat

  if (min(Y[[1]]) < 0){
    QL_loss_arithmetic_mean <- QL_loss_function(exp(arithmetic_mean_based_pred), exp(outcome))
    QL_loss_adjusted <- QL_loss_function(exp(adjusted_pred), exp(outcome)) #tk
    QL_loss_unadjusted <- QL_loss_function(exp(unadjusted_pred), exp(outcome)) #tk
  }
  else{
    QL_loss_arithmetic_mean <- QL_loss_function(arithmetic_mean_based_pred, outcome)
    QL_loss_adjusted <- QL_loss_function(adjusted_pred, outcome) #tk
    QL_loss_unadjusted <- QL_loss_function(unadjusted_pred, outcome) #tk
  }

  list_of_linear_combinations <- list(w_hat)

  list_of_forecasts <- list(unadjusted_pred
                            , adjusted_pred
                            , arithmetic_mean_based_pred)

  list_of_losses <- list(QL_loss_unadjusted
                         , QL_loss_adjusted
                         , QL_loss_arithmetic_mean)

  names(list_of_forecasts) <- c('unadjusted_pred'
                                , 'adjusted_pred'
                                , 'arithmetic_mean_based_pred')

  output_list <- list(list_of_linear_combinations
                      , list_of_forecasts
                      , list_of_losses)

  names(output_list) <- c('linear_combinations', 'predictions', 'loss')

  ## tk OUTPUT
  cat('--------------------------------------------------------------\n',
      '-------------------HAR Results-------------------','\n',
      '--------------------------------------------------------------\n',
      'Donors:', n, '\n',  '\n',
      'Shock times:', shock_time_vec, '\n', '\n',
      'Lengths of shock times:', '#tk FIX', '\n', '\n',
      'Optimization Success:', dbw_output[[2]], '\n', '\n',
      'Convex combination:',w_hat,'\n', '\n',
      'Shock estimates:', omega_star_hat_vec, '\n', '\n',
      'Aggregate estimated shock effect:', omega_star_hat, '\n', '\n',
      'Unadjusted Forecast:', unadjusted_pred,'\n', '\n',
      'Adjusted Forecast:', adjusted_pred,'\n', '\n',
      'Arithmetic-Mean-Based Forecast:',arithmetic_mean_based_pred,'\n','\n',
      'Ground Truth (estimated by realized volatility):', outcome,'\n', '\n',
      'QL Loss of Arithmetic Mean', QL_loss_arithmetic_mean, '\n', '\n',
      'QL Loss of unadjusted:', QL_loss_unadjusted,'\n', '\n',
      'QL Loss of adjusted:', QL_loss_adjusted,'\n', '\n'
  )

  ## PLOTS

  #tk make plot_maker_HAR
  if (plots == TRUE){
    cat('\n User has opted to produce plots.','\n')
    plot_maker_HAR(Y_list
                   ,shock_time_labels = shock_time_labels
                   ,shock_time_vec = integer_shock_time_vec
                   ,shock_length_vec = rep(1, n+1) #tk
                   ,unadjusted_pred
                   ,w_hat
                   ,omega_star_hat
                   ,omega_star_hat_vec
                   #,omega_star_hat_vec_std_err
                   ,adjusted_pred
                   ,display_ground_truth = TRUE) #tk what is this set to true?
  }



  return(output_list)

} ### END HAR
