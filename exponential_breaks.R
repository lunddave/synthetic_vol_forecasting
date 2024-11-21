library(quantmod)
library(garchx)
library(lmtest)
library(extraDistr)
library(gnorm)
library(Rsolnp)
library(RColorBrewer)
library(DescTools)
library(LaplacesDemon)
library(mvtnorm)
library(latex2exp)
library(MASS)
library(tikzDevice)
library(lmtest)

options(scipen = 7)

####################### BEGIN Auxiliary functions #######################

#Decay maker
decay_maker <- function(t
                        ,eta = 4
                        ,psi = .1
                        ,T_i_star = 38 ){
  # DOC STRING HERE
  return(eta*(1 - exp(-psi*(t-(T_i_star )))))
}

#series maker
series_maker <- function(ser_len
                         ,sd
                         ,alpha
                         ,eta = 5
                         ,psi = .02
                         ,T_i_star = 38
                         ){

  # DOC STRING HERE

  preshock <- alpha + rnorm(T_i_star, sd = sd)
  
  sequence_to_use <- seq(T_i_star+1,ser_len)

  shock_and_onward <- alpha + mapply(decay_maker
                                     , sequence_to_use,
                                     MoreArgs = list(
                                       eta
                                       , psi
                                       , T_i_star)) + rnorm(n = length(sequence_to_use)
                                                            , sd = sd)

  return(c(preshock,shock_and_onward))

}

#A distance-based weighting method function adapted from code by Jilei Lin (PhD Student at GWU).
# It returns a vector determined by the user's choice of distance-based-weighting method.

dbw <- function(X,
                Tstar,
                scale = FALSE,
                sum_to_1 = 1,
                bounded_below_by = 0,
                bounded_above_by = 1,
                princ_comp_count = 2,
                normchoice = c('l1', 'l2')[2],
                penalty_normchoice = c('l1', 'l2')[1],
                penalty_lambda = 0
                ) { # https://github.com/DEck13/synthetic_prediction/blob/master/prevalence_testing/numerical_studies/COP.R
  # X is a list of covariates for the time series
  # X[[1]] should be the covariate of the time series to predict
  # X[[p]] for p = 2,...,n+1 are covariates for donors

  # T^* is a vector of shock-effects time points
  # shock effect point must be > 2

  # number of time series for pool
  n <- length(X) - 1

  # COVARIATE FOR TIME SERIES UNDER STUDY AT TSTAR
  X1 <- X[[1]][Tstar[1], , drop = FALSE] # we get only 1 row

  # LOOP for grab TSTAR covariate vector for each donor
  X0 <- c()
  for (i in 1:n) {
    X0[[i]] <- X[[i + 1]][Tstar[i + 1], , drop = FALSE] #get 1 row from each donor
  }

        if (scale == TRUE) { #begin if statement
          dat <- rbind(X1, do.call('rbind', X0)) # do.call is for cluster computing?
          dat <- apply(dat, 2, function(x) scale(x, center = TRUE, scale = TRUE))
          X1 <- dat[1, , drop = FALSE]
          X0 <- c()
          for (i in 1:n) {
            X0[[i]] <- dat[i + 1, , drop = FALSE] #we are repopulating X0[[i]] with scaled+centered data
          } #end loop
        } #end if statement

  # objective function
  weightedX0 <- function(W) {
    # W is a vector of weight of the same length of X0
    n <- length(W)
    p <- ncol(X1)
    XW <- matrix(0, nrow = 1, ncol = p)
    for (i in 1:n) {
      XW <- XW + W[i] * X0[[i]]
    } #end of loop

    if (normchoice == 'l1') {
      norm_output <- as.numeric(norm(matrix(X1 - XW), type = "1"))
    }
    else {
      norm_output <- as.numeric(crossprod(matrix(X1 - XW)))
    }

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
            upper_bound = rep(bounded_above_by, n)
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
                                           , tol = 1e-9
                                           , outer.iter = 10000
                                           , inner.iter = 10000))

  #We calculate the loss from the optimization
  if (normchoice == 'l1') {
    loss <- round(norm(X1 - object_to_return$pars %*% dat[-1,], type = '1'),3)
  }
  else {
    loss <- round(norm(X1 - object_to_return$pars %*% dat[-1,], type = '2'),3)
  }

  #And finally, get norm of X1
  if (normchoice == 'l1') {
    X1_norm <- as.numeric(norm(X1, type = "1"))
  }
  else {
    X1_norm <- round(norm(X1, type = '2'),3)
  }

  output <- list()

  output[[1]] <- object_to_return$pars

  output[[2]] <- loss

  output[[3]] <- X1_norm

  return(output)

} #END dbw function

####################### END Auxiliary functions #######################

#####################################################################
########### A function that produces exactly one simulation case for synthetic volatility
#####################################################################

exp_break_maker <- function(n
                            ,p
                            ,H
                            ,alpha
                            ,eta
                            ,a
                            ,b
                            ,shock_sd
                            ,mu_delta
){

  ## Doc String
  # donor_pool_size
  # ,p
  # ,H
  # ,eta
  # ,a
  # ,b
  # ,replication_number
  # ,optimization_norm
  # ,mu_eps_star
  # ,vol_shock_sd
  # ,M21_M22_vol_mu_delta

  #Simulate series lengths
  Tee <- rdunif(n+1, a, b)
  
  shock_time_vec <- a + rep(5, n+1)
  
  print('Here is the shock time vec:')
  print(shock_time_vec)
    
  ############ Simulate all n+1 series   ############

  #Create null lists for dependent variable and independent variable output
  Y <- vector(mode = "list", length = n+1)
  X <- vector(mode = "list", length = n+1)

  # This vector will be used to store psi estimates
  eta_estimate <- c()
  psi_estimate <- c()

  #Create covariate MVN mean and sigma parameters
  vector_M21_M22_vol_mu_delta <- rep(1, p)
  matrix_M21_M22_vol_sd_delta <- matrix(diag(1**2,p), ncol=p)

  #Create M21 level and M21 vol cross donor random effects vectors
  M21_vol_cross_donor_random_effect <- mu_delta * seq(1,p,1) / (p*(p+1)/2)

  # For each of n+1 series...
  for (i in 1:(n+1)){

    #Now add the design matrix (covariates) to the list X
    covariates <- matrix(NA, nrow = shock_time_vec[i], ncol = p)
    # print('Here is the NA matrix we just generated:')
    # print(covariates)
    covariate_vec <- rmvnorm(n=1, mean=vector_M21_M22_vol_mu_delta, sigma=matrix_M21_M22_vol_sd_delta)

    covariates[shock_time_vec[i], ] <- covariate_vec
    X[[i]] <- covariates

    print(paste('now we simulate a series of length', Tee[i]))
    
    print('Here is our generated psi value:')
    
    # print(M21_vol_cross_donor_random_effect)
    # 
    # print(as.matrix(covariate_vec))
    # 
    # print(dim(t(as.matrix(M21_vol_cross_donor_random_effect))))
    # 
    # print(dim(as.matrix(covariate_vec)))
    
    psi_generated <- as.numeric(t(as.matrix(M21_vol_cross_donor_random_effect)) %*% t(as.matrix(covariate_vec)))
    
    simulated_series <- series_maker(ser_len = Tee[i]
                                      ,sd = shock_sd
                                      ,alpha = alpha
                                      ,eta = eta
                                      ,psi = psi_generated
                                      ,T_i_star = shock_time_vec[i]
                                     )
    
    print('We generated the series.')
    
    plot.ts(simulated_series)
    
    # fit the model 
    start_values <- c(eta = eta
                      , psi = psi_generated
                      ) 
    
    decay_indices <- seq(1, length(simulated_series) - shock_time_vec[i])
    
    post_shock_indices <- seq(shock_time_vec[i] + 1, length(simulated_series))
    
    pre_shock_alpha <- mean(simulated_series[1:shock_time_vec[i]] )
    
    fit <- nls(simulated_series[post_shock_indices] - pre_shock_alpha ~ eta *(1 - exp(-psi * decay_indices)),
               start = start_values,
               algorithm = "port",
               control = nls.control(maxiter = 1000000, tol = 1e-09, minFactor = 1/(2**24),
                                     printEval = FALSE, warnOnly = FALSE, scaleOffset = 0,
                                     nDcentral = FALSE))
    
    coefficients <- coeftest(fit)
    
    eta_hat <- coefficients[nrow(coefficients)-1,1]
    psi_hat <- coefficients[nrow(coefficients),1]
    
    # alpha_hat <- mean(simulated_series[1:shock_time_vec[i]])
    # eta_hat <- simulated_series[Tee[i]]
    # epsilon_hat <- simulated_series - alpha_hat
    # epsilon_hat_post_shock <- epsilon_hat[shock_time_vec[i]+1:Tee[i]]
    # divisor_vec <- seq(1, Tee[i])
    # psi_hat <- (1/(Tee[i]-shock_time_vec[i]))*sum(-log(1 - epsilon_hat_post_shock/eta_hat))/divisor_vec

    eta_estimate <- c(eta_estimate,psi_hat)
    psi_estimate <- c(psi_estimate,psi_hat)

  } #end loop for n+1 series
  
  print('We are all done generating data.')

  #Now make xreg into a dataframe
  eta_psi_matrix <- cbind(eta_estimate, psi_estimate)
  
  #NEXT, we get the vectors w for all sensible methods
  w <- list() #initialize
  dbw_loss <- c()
  X1_norm <- c()

  matrix_of_specs <- matrix(c(rep(1,6),
                              rep(NA,6),

                              rep(c(0,-1,NA),4),

                              rep(c(1,NA), 6)),
                              byrow = FALSE, nrow = 12)

  #We drop the 4th row because it's functionally no different from the first OR have lower bound > upper bound
  matrix_of_specs <- matrix_of_specs[-4,]

  for (i in 1:nrow(matrix_of_specs)){

  dbw_output <- dbw(X,
                    shock_time_vec[i], #tk
                 scale = TRUE,
                 sum_to_1 = matrix_of_specs[i,1],
                 bounded_below_by = matrix_of_specs[i,2],
                 bounded_above_by = matrix_of_specs[i,3],
                 #normchoice = normchoice,
                 #penalty_normchoice = penalty_normchoice,
                 #penalty_lambda = penalty_lambda
                 )

  w[[i]] <- dbw_output[[1]]
  # dbw_loss[i] <- dbw_output[[2]]
  # X1_norm[i] <- dbw_output[[3]]

  }

  # Now we place these linear combinations into a matrix
  w_mat <- matrix(unlist(w), nrow = nrow(matrix_of_specs), byrow = TRUE)
  
  #add a row to the matrix corresponding to the arithmetic mean
  col_count <- dim(w_mat)[,1]
  w_mat <- rbind(w_mat, rep(1/col_count, col_count))
  
  #add a row corresponding to the unadjusted forecast
  w_mat <- rbind(w_mat, rep(0, col_count))

  #Second, we calculate omega_star_hat, which is the dot product of w and the estimated shock effects
  weighted_estimates <- as.numeric(w_mat %*% eta_psi_matrix)
  
  TSUS_prediction <- mean(simulated_series[1:shock_time_vec[1]]) + 
    
  decay_preds <- decay_maker(seq(1, Tee[1] - shock_time_vec[1])
                          ,eta = weighted_estimates[1,1]
                          ,psi = weighted_estimates[1,2]
                          ,T_i_star = 0 )
  
  pred_matrix <- rbind(TSUS_prediction + decay_preds
                       , rep(TSUS_prediction,length(decay_preds)))
  
  loss_matrix <- cbind(sum((pred_matrix[1,1:5]-pred_matrix[2,1:5])**2)
                       ,sum((pred_matrix[1,1:10]-pred_matrix[2,1:10])**2)
                       ,sum((pred_matrix[1,1:50]-pred_matrix[2,1:50])**2))
  
  return(list(pred_matrix,loss_matrix))
}


exp_break_maker(n = 10
              ,p = 3
              ,H = 1
              ,alpha = 100
              ,eta = -2
              ,a = 200
              ,b = 800
              #,optimization_norm
              ,shock_sd = .001
              ,mu_delta = .01
            )
