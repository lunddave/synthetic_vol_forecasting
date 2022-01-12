# David Lundquist
# Simulations for Synthetic Prediction GARCH

library(quantmod)
library(garchx)
library(lmtest)
library(extraDistr)
library(gnorm)
library(tsDyn)
library(Rsolnp)

options(scipen = 6)

## Doc String

# synth_vol_sim: function that simulates (n+1)*(p+1) time series: 
# a response series and p covariate series for each of the n donors 
# and for the time series under study, as well.  The series
# must experience an exogenous shock at exactly one discrete time point
# in the series.
# --Input:
#   --n
#   --p
#   --model for the shock (M1, M2)
#   --sigma of the shock
#   --shock time vector (if specified by user. Otherwise, simulated.)
#   -- a,b, the parameters used to construct the discrete uniform from which we get series lengths

synth_vol_sim <- function(n, p, arch_param, garch_param, 
                          level_model, 
                          vol_model,
                          sigma_GARCH_innov, sigma_x, 
                          shock_time_vec, 
                          level_shock_length,
                          vol_shock_length,
                          a, b, 
                          mu_eps_star, M2_mu_eps_star, sigma_eps_star,
                          mu_omega_star, M2_mu_omega_star, vol_shock_multiplier,
                          vol_shock_sd,
                          omega_shape, omega_rate,
                          level_GED_alpha,
                          level_GED_beta,
                          ...){
  
  #Simulate series lengths
  Tee <- rdunif(n+1, a, b)
  
  #Before we simulate shock time, we make sure each series has enough points 
  #following the shock time
  max_of_shock_lengths <- max(level_shock_length, vol_shock_length)
  warnings('Each series has shock ', 5)
  
  # Simulate shock times
  if ( is.null(shock_time_vec) == TRUE)
  {
        shock_time_vec <- c()
        for (i in 1:(n+1))
        {
          #Note: the T* must be at least 'max_of_shock_lengths' after T*
          #Also, shock must come from point 30 onward.
          shock_time_vec[i] <- rdunif(1, 30, Tee[i]-max_of_shock_lengths) 
        }
  }
  
  ############ Simulate Structure of Covariates ############
  
  # Now generate the covariates.  These will be correlated GARCH processes, ideally. 
  # Since multivariate GARCH processes take take technical care to simulate, we first use VAR.
  
  # https://math.stackexchange.com/questions/1529000/how-to-create-a-random-matrix-whose-spectral-radius-1
  
  #Random parameters for the VAR
  param_matrix_entries <- runif(p**2, min = -1/p, max = 1/p)
  simVAR_params <- matrix(param_matrix_entries, nrow = p, byrow = T)
  
  ############ Simulate all n+1 series   ############ 
  
  #Create null lists for depvar and indepvar output
  Y <- vector(mode = "list", length = n+1)
  X <- vector(mode = "list", length = n+1)
  level_shock_vec <- c()
  vol_shock_vec <- c()
  T_star_plus_1_return_vec <- c()
  
  # For each of n+1 series...
  for (i in 1:(n+1)){
    
    #Epsilon vector for the VAR
    innovations_matrix_entries <- rnorm(Tee[i] * p, sd = sigma_x)
    sim_VAR_innovations <- matrix(innovations_matrix_entries, ncol = p, byrow = T)
    
    VAR_process <- VAR.sim(B = simVAR_params, 
                           lag = 1, 
                           include = "none", 
                           n = Tee[i],
                           innov = sim_VAR_innovations) + 1
    
    #Level model
    if (level_model == 'M1'){
                      #Level Shock
                      level_shock_vec[i] <- mu_eps_star + # This is the non-stochastic term
                                            rgnorm(1, 
                                            mu = 0, 
                                            alpha = level_GED_alpha, 
                                            beta = level_GED_beta) # This is the stochastic term
          
          level_shock_mean <- mu_eps_star             
          level_shock_var <- ((level_GED_alpha)**2) * gamma(3/level_GED_beta) / (gamma(1/level_GED_beta)) # https://search.r-project.org/CRAN/refmans/gnorm/html/gnorm.html
                      
                      } 
    else if (level_model == 'M2') { 
          level_shock_vec[i] <- mu_eps_star + 
                                as.numeric(as.matrix(VAR_process[shock_time_vec[i],])) %*% rnorm(p,M2_mu_eps_star,sigma_eps_star) + 
                                rgnorm(1, 
                                       mu = 0, 
                                       alpha = level_GED_alpha, 
                                       beta = level_GED_beta) #What's the variance of this sum?
          
          level_shock_mean <- mu_eps_star 
          level_shock_var <- ((level_GED_alpha)**2) * gamma(3/level_GED_beta) / (gamma(1/level_GED_beta)) + # https://search.r-project.org/CRAN/refmans/gnorm/html/gnorm.html
            (sigma_x**2) * (sigma_eps_star**2)
    }
    else {level_shock_vec[i] <- rnorm(1, 0, sigma_GARCH_innov); level_shock_mean <- 0; level_shock_var <- sigma_GARCH_innov**2}
    
    #Vol model
    if (vol_model == 'M1'){
      
      #Create volatility shock w*
      vol_shock_vec[i] <- rnorm(1, mu_omega_star, vol_shock_sd)
      vol_shock_mean <- vol_shock_multiplier * mu_omega_star 
      vol_shock_var <- (vol_shock_multiplier**2) * vol_shock_sd**2
      
      shock_indicator <- c(
        rep(0, shock_time_vec[i]), 
        rep(vol_shock_vec[i], vol_shock_length), 
        rep(0, Tee[i] - shock_time_vec[i] - vol_shock_length))
      
      #Now add the design matrix to the list X
      X[[i]] <- cbind(VAR_process, shock_indicator)
      
      #Create GARCH model with shock(s)
      GARCH_innov_vec <- c(
        rnorm(shock_time_vec[i], 0, sigma_GARCH_innov), 
        level_shock_vec[i],
        rnorm(Tee[i] - shock_time_vec[i] - 1, 0, sigma_GARCH_innov))
      
      Y[[i]] <- garchxSim(Tee[i], arch = arch_param, garch = garch_param, 
                          xreg =  as.matrix( X[[i]][,(p+1)] ),
                          innovations = GARCH_innov_vec, verbose = TRUE) 
      
    } 
    
    else if (vol_model == 'M2') { 
      vol_shock_vec[i] <- rnorm(1, mu_omega_star, vol_shock_sd)
        as.numeric(as.matrix(VAR_process[shock_time_vec[i],])) %*% rnorm(p,M2_mu_omega_star,sigma_eps_star) 
           #What's the variance of this sum?
      
      vol_shock_mean <- vol_shock_multiplier * mu_omega_star 
      vol_shock_var <- (vol_shock_multiplier**2) * vol_shock_sd**2 + (sigma_x**2) * (sigma_eps_star**2)
      
      shock_indicator <- c(
        rep(0, shock_time_vec[i]), 
        rep(vol_shock_vec[i], vol_shock_length), 
        rep(0, Tee[i] - shock_time_vec[i] - vol_shock_length))
      
      #Now add the design matrix to the list X
      X[[i]] <- cbind(VAR_process, shock_indicator)
      
      #Create GARCH model with shock(s)
      GARCH_innov_vec <- c(
        rnorm(shock_time_vec[i], 0, sigma_GARCH_innov), 
        level_shock_vec[i],
        rnorm(Tee[i] - shock_time_vec[i] - 1, 0, sigma_GARCH_innov))
      
      Y[[i]] <- garchxSim(Tee[i], arch = arch_param, garch = garch_param, 
                          xreg =  as.matrix( X[[i]][,(p+1)] ),
                          innovations = GARCH_innov_vec, verbose = TRUE) 
    }
    
    else { 

      #Create volatility shock w*
      vol_shock_vec[i] <- rnorm(1, 0, sigma_GARCH_innov)
      vol_shock_mean <- 0
      vol_shock_var <- 0
        
        #Now add the design matrix to the list X
        X[[i]] <- cbind(VAR_process)
        
        #Create GARCH model with shock(s)
        GARCH_innov_vec <- rnorm(Tee[i], 0, sigma_GARCH_innov)
        
        Y[[i]] <- garchxSim(Tee[i], arch = arch_param, garch = garch_param, 
                            innovations = GARCH_innov_vec, verbose = TRUE) 
    } #end conditionals that create vol shocks
    
    T_star_plus_1_return_vec[i] <- Y[[i]][,1][shock_time_vec[1]+1,]
    
  } #end loop for n+1 series
  
  ## Compute summary statistics for output
  level_shock_kurtosis <- gamma(5/level_GED_beta)*gamma(1/level_GED_beta)/( (gamma(3/level_GED_beta))**2 ) - 3 #https://en.wikipedia.org/wiki/Generalized_normal_distribution
  vol_shock_kurtosis <- 6 / omega_shape
  
  T_star_sigma <- Y[[1]][,3][shock_time_vec[1],]
  T_star_plus_1_sigma <- Y[[1]][,3][shock_time_vec[1]+1,]
  T_star_plus_2_sigma <- Y[[1]][,3][shock_time_vec[1]+2,]
  T_star_plus_3_sigma <- Y[[1]][,3][shock_time_vec[1]+3,]
  
  ##Output
  cat('Simulation Summary Data','\n',
      '-------------------------------------------------------------\n',
      'Donors:', n, '\n',
      'Series lengths:', Tee, '\n',
      'Shock times:', shock_time_vec, '\n',
      'Level Shock at T*+1:', round(level_shock_vec,2), '\n', 
      'T*+1: Return:', round(T_star_plus_1_return_vec,3), '\n', 
      'Volatility Shock at T*+1', round(vol_shock_vec,2), '\n',
      '\n',
      'Volatility of Time Series under Study', '\n',
      '-------------------------------------------------------------\n',
      'Sigma^2 at T*:', round(T_star_sigma,2), '\n', 
      'Sigma^2 at T*+1:', round(T_star_plus_1_sigma,2), '\n', 
      'Sigma^2 at T*+2:', round(T_star_plus_2_sigma,2), '\n', 
      'Sigma^2 at T*+3:', round(T_star_plus_3_sigma,2), '\n', 
      '\n',
      'Level Shock Moments', '\n',
      '-------------------------------------------------------------\n',
      'Level Shock mean:', round(level_shock_mean,4), '(equivalent to a', round(100*level_shock_mean,2), '% daily move).', ' \n',
      'Level Shock variance:', round(level_shock_var,4), '\n',
      'Level Signal to Noise:', abs(round(level_shock_mean / sqrt(level_shock_var),2)) , '\n',
      'Level Shock excess kurtosis:', round(level_shock_kurtosis, 2) , '\n',
      
      '\n',
      'Vol Shock Moments', '\n',
      '-------------------------------------------------------------\n',
      'Vol Shock mean:', round(vol_shock_mean,2), ' \n',
      'Vol Shock variance:', round(vol_shock_var,4), '\n',
      'Vol Signal to Noise:', abs(round(vol_shock_mean / sqrt(vol_shock_var),3)) , '\n',
      'Vol Shock excess kurtosis:', round(vol_shock_kurtosis, 2)
      )
  
  #Plot the donors
  par(mfrow = c(ceiling(sqrt(2*n + 2)), ceiling(sqrt(2*n + 2))))
  for (i in 1:(n+1))
  {
    plot.ts(X[[i]][-c(1:20),], main = paste('Covariates of Donor ', i, sep = ''))
  }
  
  #Plot the time series
  par(mfrow = c(ceiling(sqrt(n+1)), ceiling(sqrt(n+1))))
  for (i in 1:(n+1))
  {
    
    plot.ts(Y[[i]][,1], ylim = c(min(Y[[i]][,1])*1.2, max(Y[[i]][,1])*1.2), 
            main = paste('y_', i, ", GARCH(",arch_param,",",garch_param,")",
                                     "\n level shock = ", 
                                     round( level_shock_vec[i],2), 
                                     ", vol shock = ", 
                                     round(vol_shock_vec[i],2),
                                     sep = ''), ylab = 'Daily Log-Return')
    abline(v = shock_time_vec[i] + 1, col = 'red')

  }
  
  #Plot the volatility series
  par(mfrow = c(ceiling(sqrt(n+1)), ceiling(sqrt(n+1))))
  for (i in 1:(n+1))
  {
    plot.ts(Y[[i]][-c(1:20),3], xlim=c(21, Tee[i]), main = paste('Volatility Series of y_', i, 
                                     ", GARCH(",arch_param,",",garch_param,")",
                                     "\n level shock = ", 
                                     round( level_shock_vec[i],2), 
                                     ", vol shock = ", 
                                     round(vol_shock_vec[i],2),
                                     sep = ''), ylab = 'Sigma^2')
    abline(v = shock_time_vec[i] + 1, col = 'red')
  }
  
  #Items to return in a list
  return(list(X,Y,Tee,shock_time_vec))
}

output <- synth_vol_sim(n = 8, 
                        p = 2, 
                        arch_param = c(.2),
                        garch_param = c(.7),
                        level_model = c('M1','M2','none')[3],
                        vol_model = c('M1','M2','none')[2],
                        sigma_GARCH_innov = (.007), # this is the sd that goes into rnorm
                        sigma_x = .008, 
                        shock_time_vec = NULL, 
                        level_shock_length = 1,
                        vol_shock_length = 10,
                        a = 90, 
                        b = 150, 
                        mu_eps_star = -.0825,
                        M2_mu_eps_star = .3, 
                        sigma_eps_star = .001,
                        mu_omega_star = .8,
                        M2_mu_omega_star = 2,
                        vol_shock_multiplier = 1,
                        vol_shock_sd = .2,
                        omega_shape = .2, #these are to be deleted
                        omega_rate = 2, #these are to be deleted
                        level_GED_alpha = .05 * sqrt(2), 
                        level_GED_beta = 1.8)

#What is the shock time for y_1?
output[[4]][[1]]

output[[2]][[1]]

## Let's hit the time series under study with a GARCH(1,1)
mod <- garchx(output[[2]][[1]][,1], order = c(1,0))
mod
plot.ts(mod$fitted)
#lines(output[[2]][[1]][,3] * (output[[2]][[1]][,5]**2), col = 'red')
plot.ts(mod$residuals)

#Let's look at the model performance ending immediately after T*+1
length_of_y1 <- output[[3]][1]
shock_time_of_y1 <- output[[4]][1]
covariate_mat_2 <-  as.matrix(c(rep(0,shock_time_of_y1), 1, rep(0, length_of_y1 - shock_time_of_y1 - 1) ))
covariate_mat_3 <- as.matrix(cbind(output[[1]][[1]][,-2], covariate_mat_2))
mod2 <- garchx(output[[2]][[1]][,1][c(1:(shock_time_of_y1+1)),], order = c(1,0), xreg = covariate_mat_3[c(1:(shock_time_of_y1+1)),])
mod2
mod3 <- garchx(output[[2]][[1]][,1][c(1:(shock_time_of_y1+1)),], order = c(1,1), xreg = covariate_mat_2[c(1:(shock_time_of_y1+1)),])
mod3

k <- 2
covariate_mat_4 <- as.matrix(c(rep(0,shock_time_of_y1), rep(1,k), rep(0, length_of_y1 - shock_time_of_y1 - k) ))
mod4 <- garchx(output[[2]][[1]][,1][c(1:(shock_time_of_y1+k)),], order = c(0,1), xreg = covariate_mat_4[c(1:(shock_time_of_y1+k)),])
mod4

#Look at covariates of time series of interest
plot.ts(output[[1]][[1]])

#Look at GARCH process of time series of interest
plot.ts(output[[2]][[1]])

#What is the mean of the squared y for time series of interest?
mean(output[[2]][[1]][,1]**2)

# Objective2: estimation function that takes (n+1)*(p+1) time series as input and 
# 1) calculates weight vector w, 
# 2) calculates fixed effects estimate vector omega*, 
# 3) calculates the adjusted estimate for the
# volatility of time series of interest at T*+1 (i.e. the prediction)
# 4) calculates estimate of volatility on T*+1 for each series 
# using each of three families, and 
# 5) calculates the squared-error loss of the prediction
# Estimation/control options
# --Allow user to enter series of unequal lengths
# --Allow user to enter a vector of integers corresponding to the number of days
# the shock effect lasts for each outcome series
# --Allow user to pick a uniform model for each series (e.g. GARCH(1,1)) OR a BIC-minimizing
# model for each series (or mix and match).
# --Allow user to pick error distribution - see ugarchspec


# this function returns the W^* estimated by synthetic control method (SCM)
scm <- function(X, Tstar, scale = FALSE) { # https://github.com/DEck13/synthetic_prediction/blob/master/prevalence_testing/numerical_studies/COP.R
  # X is a list of covariates for the time series
  # X[[1]] should be the covariate of the time series to predict
  # X[[p]] for p = 2,...,n+1 are covariates for donors
  
  # T^* is a vector of shock-effects time points
  # shock effect point must be > 2
  
  # number of time series for pool
  n <- length(X) - 1
  
  # covariate for time series for prediction
  X1 <- X[[1]][Tstar[1] + 1, , drop = FALSE] # we get only 1 row
  
  # covariates for time series pool
  X0 <- c()
  for (i in 1:n) {
    X0[[i]] <- X[[i + 1]][Tstar[i + 1] + 1, , drop = FALSE] #get 1 row from each donor
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
    norm <- as.numeric(crossprod(matrix(X1 - XW)))
    return(norm)
  } #end objective function
  
  # constraint for W
  Wcons <- function(W) sum(W) - 1
  
  # optimization
  return(solnp(par = rep(1/n, n), fun = weightedX0, eqfun = Wcons, eqB = 0, # will use inequality constraints later
        LB = rep(0, n), UB = rep(1, n), control = list(trace = 0)))
  
  #I added the return statement because an implicit return is bad coding form
} #end SCM function





synth_vol_fit <- function()
{
  
}



