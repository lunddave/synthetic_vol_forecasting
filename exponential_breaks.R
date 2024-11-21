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
series_maker <- function(length
                         ,sd
                         ,alpha
                         ,eta = 5
                         ,psi = .02
                         ,T_i_star = 38
                         ){

  # DOC STRING HERE

  preshock <- alpha + rnorm(T_i_star, sd = sd)

  shock_and_onward <- alpha + mapply(decay_maker, seq(T_i_star+1,length,1),
                                     MoreArgs = list(
                                       eta
                                       , psi
                                       , T_i_star)) +
    rnorm(length-T_i_star, sd = sd)

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

exp_break_maker <- function(donor_pool_size
                            ,p
                            ,H
                            ,eta
                            ,a
                            ,b
                            ,replication_number
                            ,optimization_norm
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

  ############ Simulate all n+1 series   ############

  #Create null lists for dependent variable and independent variable output
  Y <- vector(mode = "list", length = n+1)
  X <- vector(mode = "list", length = n+1)

  # This vector will be used to store psi estimates
  psi_estimate <- c()

  #Create covariate MVN mean and sigma parameters
  vector_M21_M22_vol_mu_delta <- rep(mu_x, p)
  matrix_M21_M22_vol_sd_delta <- matrix(diag(sigma_x**2,p), ncol=p)

  #Create M21 level and M21 vol cross donor random effects vectors
  M21_vol_cross_donor_random_effect <- mu_delta * seq(1,p,1) / (p*(p+1)/2)

  # For each of n+1 series...
  for (i in 1:(n+1)){

    # #Epsilon vector for the VAR
    # innovations_matrix_entries <- rnorm( (shock_time_vec[i]) * p, sd = sigma_x)
    # sim_VAR_innovations <- matrix(innovations_matrix_entries, ncol = p, byrow = T)
    #
    # #Note: we need only covariate information up through (shock_time_vec[i] ), since
    # #(a) we assume the covariate to be a lag1 r.v. and (b) we model the shock as a function of
    # #what occurred at (shock_time_vec[i])
    #
    # VAR_process <- VAR.sim(B = simVAR_params,
    #                        lag = 1,
    #                        include = "none",
    #                        n = shock_time_vec[i], # we do not need more data points than this
    #                        innov = sim_VAR_innovations) + 1

    #   #As of March 22nd, 2023, a VAR(p) for the covariates has been deprecated.  Using MVN now.

    #Now add the design matrix (covariates) to the list X
    covariates <- matrix(NA, nrow = shock_time_vec[i], ncol = p)
    covariate_vec <- rmvnorm(n=1, mean=vector_M21_M22_vol_mu_delta, sigma=matrix_M21_M22_vol_sd_delta)
    covariates[shock_time_vec[i], ] <- covariate_vec
    X[[i]] <- covariates


    simulated_series <- series_maker(Tee[i]
                                        ,shock_sd
                                        ,100
                                        ,eta = -15
                                        ,psi = M21_vol_cross_donor_random_effect %*% covariates
                                        ,58)
    
    alpha_hat <- mean(simulated_series[1:shock_time_vec[i]])
    
    eta_hat <- simulated_series[Tee[i]]
    
    epsilon_hat <- simulated_series - alpha_hat
    epsilon_hat_post_shock <- epsilon_hat[shock_time_vec[i]+1:Tee[i]]
    
    divisor_vec <- seq(1, Tee[i])
    
    psi_hat <- (1/(Tee[i]-shock_time_vec[i]))*sum(-log(1 - epsilon_hat_post_shock/eta_hat))/divisor_vec
    
    psi_estimate <- c(psi_estimate,psi_hat)

  } #end loop for n+1 series

  #Now make xreg into a dataframe
  psi_estimate <- data.frame(matrix(psi_estimate, nrow = n+1, byrow = TRUE))
  
  # fit the model 
  start_values <- c(a=4, b=2)
  fit <- nls(y ~ a * exp(1 - b * x),
             start = start_values,
             algorithm = "port",
             control = nls.control(maxiter = 1000))
  summary(fit)

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
                 T_star,
                 scale = TRUE,
                 sum_to_1 = matrix_of_specs[i,1],
                 bounded_below_by = matrix_of_specs[i,2],
                 bounded_above_by = matrix_of_specs[i,3],
                 normchoice = normchoice,
                 penalty_normchoice = penalty_normchoice,
                 penalty_lambda = penalty_lambda)

  w[[i]] <- dbw_output[[1]]
  dbw_loss[i] <- dbw_output[[2]]
  X1_norm[i] <- dbw_output[[3]]

  }

  # Now we place these linear combinations into a matrix
  w_mat <- matrix(unlist(w), nrow = nrow(matrix_of_specs), byrow = TRUE)

  #Second, we calculate omega_star_hat, which is the dot product of w and the estimated shock effects
  omega_star_hat_vec <- as.numeric(w_mat %*% psi_estimate)

  #Third, we add the lin_reg_pred to omega_star_hat_vec
  omega_star_hat_vec <- c(omega_star_hat_vec, lin_reg_pred)

}
