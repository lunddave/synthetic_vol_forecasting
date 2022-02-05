# David Lundquist
# Simulations for Synthetic Prediction GARCH

library(quantmod)
library(garchx)
library(lmtest)
library(extraDistr)
library(gnorm)
library(tsDyn)
library(Rsolnp)
library(RColorBrewer)
library(DescTools)

options(scipen = 7)

#####################################################################

########### A function that produces exactly one simulation case for synthetic volatility

#####################################################################

synth_vol_sim <- function(n, 
                          p, 
                          arch_param, 
                          garch_param,
                          asymmetry_param,
                          level_model, 
                          vol_model,
                          t_dist_v_param, 
                          sigma_x, 
                          earliest_shock_time = 100,
                          shock_time_vec, 
                          level_shock_length,
                          vol_shock_length,
                          a, 
                          b, 
                          mu_eps_star, 
                          M22_mu_eps_star, 
                          sigma_eps_star,
                          mu_omega_star, 
                          M22_mu_omega_star,
                          vol_shock_sd,
                          level_GED_alpha,
                          level_GED_beta,
                          ...){
  
  ## Doc String
  
  # synth_vol_sim: function that simulates (n+1)*(p+1) time series: 
  # a response series and p covariate series for each of the n donors 
  # and for the time series under study, as well.  The series
  # must experience an exogenous shock at exactly one discrete time point
  # in the series.
  # --Input:
  #   --n - number of donors (scalar)
  #   --p - number of covariates (scalar)
  #   --arch parameters (vector of length 0 or more)
  #   --garch parameters (vector of length 0 or more)
  #   --asymmetry_param (vector of length 0 or more)
  #   --level_model - model for level of the shock (string)
  #   --vol_model - model for volatility of the shock (string)
  #   --t_dist_v_param (scalar)
  #   --sigma_x - sigma of innovations in covariates (scalar)
  #   --shock_time_vec - optional input to force shock times (vector of integers)
  #   --level_shock_length (scalar)
  #   --vol_shock_length (scalar)
  #   --a - minumum series length (scalar)
  #   --b - maximum series length (scalar)
  #   --mu_eps_star - intercept of level shock for M1 model
  #   --M22_mu_eps_star - intercept of level shock for M22 model
  #   --sigma_eps_star - variance of level shock for all level shocks
  #   --mu_omega_star - intercept of vol shock for M1 model
  #   --M22_mu_omega_star - intercept of vol shock for M22 model
  #   --vol_shock_sd - variance of level shock for all vol shocks
  #   --level_GED_alpha - alpha parameter for level shock stochastic term
  #   --level_GED_beta - beta parameter for level shock stochastic term

  #Simulate series lengths
  Tee <- rdunif(n+1, a, b)
  
  #Before we simulate shock time, we make sure each series has enough points 
  #following the shock time
  max_of_shock_lengths <- max(level_shock_length, vol_shock_length)

  # Simulate shock times
  if ( is.null(shock_time_vec) == TRUE)
  {
        shock_time_vec <- c()
        for (i in 1:(n+1))
        {
          #Note: the T* must be at least 'max_of_shock_lengths' after T*
          #Also, shock must come from point 30 onward.
          shock_time_vec[i] <- rdunif(1, earliest_shock_time, Tee[i]-max_of_shock_lengths) 
        }
  }
  
  ############ Simulate Structure of Covariates ############
  
  # Now generate the covariates.  These will be correlated GARCH processes, ideally. 
  # Since multivariate GARCH processes take take technical care to simulate, we first use VAR.
  
  # https://math.stackexchange.com/questions/1529000/how-to-create-a-random-matrix-whose-spectral-radius-1
  
  #Random parameters for the VAR
  param_matrix_entries <- runif(p**2, min = -1/p, max = 1/p)
  simVAR_params <- matrix(param_matrix_entries, nrow = p, byrow = T)
    #Note: In our meeting on Jan 12, 2022, Dan Eck and I discussed the justification for fixing 
    #a common Phi matrix for all of the n+1 series.  By fixing a common Phi matrix, 
    #we're basically saying that that each covariate has the same structure across all donors.
    #For example, WTI and VIX behave the same way for Lehmann 2008 as they do for
    #Russia-Saudi-OPEC 2014.  We might want to relax this at some point by partitioning the n
    #donors into k sets and saying that the n donors are temporally clustered that way,
    #Doing it this way would imitate the way that Lin and Eck 2021 COP donors fall into
    #the sets Spring 2008, Fall 2008, Spring 2014.
  
    # An additional note on covariates:
    #   
    # Here I am simulating the covariates to be mean = 1 and to have the same innovation variance.
    # In an application, however, it would be good to follow Lin and Eck (2021) in centering and scaling
    # the covariates so that no covariate dominates in the process of a getting a convex combination w.
  
  ############ Simulate all n+1 series   ############ 
  
  #Create null lists for depvar and indepvar output
  Y <- vector(mode = "list", length = n+1)
  X <- vector(mode = "list", length = n+1)
  level_shock_vec <- c()
  vol_shock_vec <- c()
  T_star_plus_1_return_vec <- c()
  xreg <- c()
  
  #Create M21 level and M21 vol cross donor random effects vectors
  M21_vol_cross_donor_random_effect <- rnorm(p, 0, 1) # tk
  M21_vol_cross_donor_random_effect <- rnorm(p, 0, 1) # tk
  
  # For each of n+1 series...
  for (i in 1:(n+1)){
    
    #Epsilon vector for the VAR
    innovations_matrix_entries <- rnorm(Tee[i] * p, sd = sigma_x)
    sim_VAR_innovations <- matrix(innovations_matrix_entries, ncol = p, byrow = T)
    
    VAR_process <- VAR.sim(B = simVAR_params, 
                           lag = 1, 
                           include = "none", 
                           n = Tee[i],
                           innov = sim_VAR_innovations) + 1 #we add this constant to make the
                                                            #covariates positive
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
    
    else if (level_model == 'M21') { 
      level_shock_vec[i] <- mu_eps_star + 
        as.numeric(as.matrix(VAR_process[shock_time_vec[i],])) %*% M21_level_cross_donor_random_effect + 
        rgnorm(1, 
               mu = 0, 
               alpha = level_GED_alpha, 
               beta = level_GED_beta) #What's the variance of this sum?
      
      level_shock_mean <- mu_eps_star 
      level_shock_var <- ((level_GED_alpha)**2) * gamma(3/level_GED_beta) / (gamma(1/level_GED_beta)) + # https://search.r-project.org/CRAN/refmans/gnorm/html/gnorm.html
        (sigma_x**2) * (sigma_eps_star**2)
    }
    
    else if (level_model == 'M22') { 
          level_shock_vec[i] <- mu_eps_star + 
                                as.numeric(as.matrix(VAR_process[shock_time_vec[i],])) %*% rnorm(p,M22_mu_eps_star,sigma_eps_star) + 
                                rgnorm(1, 
                                       mu = 0, 
                                       alpha = level_GED_alpha, 
                                       beta = level_GED_beta) #What's the variance of this sum?
          
          level_shock_mean <- mu_eps_star 
          level_shock_var <- ((level_GED_alpha)**2) * gamma(3/level_GED_beta) / (gamma(1/level_GED_beta)) + # https://search.r-project.org/CRAN/refmans/gnorm/html/gnorm.html
            (sigma_x**2) * (sigma_eps_star**2)
    }
    
    else {level_shock_vec[i] <- rt(1,t_dist_v_param); level_shock_mean <- 0; level_shock_var <- t_dist_v_param/(t_dist_v_param-2)}
    
    #Vol model
    if (vol_model == 'M1'){
      
      #Create volatility shock w*
      vol_shock_vec[i] <- rnorm(1, mu_omega_star, vol_shock_sd)
      vol_shock_mean <- mu_omega_star 
      vol_shock_var <- vol_shock_sd**2
      
      shock_indicator <- c(
                          rep(0, shock_time_vec[i]), 
                          rep(vol_shock_vec[i], vol_shock_length), 
                          rep(0, Tee[i] - shock_time_vec[i] - vol_shock_length))
        
      #Now add the design matrix to the list X
      X[[i]] <- VAR_process
      
      #Create GARCH model with shock(s)
      GARCH_innov_vec <- c(
                            rt(shock_time_vec[i], t_dist_v_param), 
                            level_shock_vec[i],
                            rt(Tee[i] - shock_time_vec[i] - 1, t_dist_v_param))
    
      Y[[i]] <- garchxSim(Tee[i], arch = arch_param, garch = garch_param, asym = asymmetry_param, 
                          xreg =  as.matrix(shock_indicator),
                          innovations = GARCH_innov_vec, verbose = TRUE) 
    } 
    
    else if (vol_model == 'M21') { #tk
      vol_shock_vec[i] <- rnorm(1, mu_omega_star, vol_shock_sd)
      as.numeric(as.matrix(VAR_process[shock_time_vec[i],])) %*% M21_vol_cross_donor_random_effect 
      #What's the variance of this sum?
      
      vol_shock_mean <- mu_omega_star 
      vol_shock_var <- vol_shock_sd**2 + (sigma_x**2) * (sigma_eps_star**2)
      
      shock_indicator <- c(
        rep(0, shock_time_vec[i]), 
        rep(vol_shock_vec[i], vol_shock_length), 
        rep(0, Tee[i] - shock_time_vec[i] - vol_shock_length))
      
      #Now add the design matrix to the list X
      X[[i]] <- VAR_process
      
      #Create GARCH model with shock(s)
      GARCH_innov_vec <- c(
        rt(shock_time_vec[i], t_dist_v_param), 
        level_shock_vec[i],
        rt(Tee[i] - shock_time_vec[i] - 1, t_dist_v_param))
      Y[[i]] <- garchxSim(Tee[i], arch = arch_param, garch = garch_param, 
                          xreg =  as.matrix(shock_indicator),
                          innovations = GARCH_innov_vec, verbose = TRUE) 
    }
    
    else if (vol_model == 'M22') { 
      vol_shock_vec[i] <- rnorm(1, mu_omega_star, vol_shock_sd)
        as.numeric(as.matrix(VAR_process[shock_time_vec[i],])) %*% rnorm(p,M22_mu_omega_star,sigma_eps_star) 
           #What's the variance of this sum?
      
      vol_shock_mean <- mu_omega_star 
      vol_shock_var <- vol_shock_sd**2 + (sigma_x**2) * (sigma_eps_star**2)
      
      shock_indicator <- c(
                        rep(0, shock_time_vec[i]), 
                        rep(vol_shock_vec[i], vol_shock_length), 
                        rep(0, Tee[i] - shock_time_vec[i] - vol_shock_length))
      
      #Now add the design matrix to the list X
      X[[i]] <- VAR_process
      
      #Create GARCH model with shock(s)
      GARCH_innov_vec <- c(
                          rt(shock_time_vec[i], t_dist_v_param), 
                          level_shock_vec[i],
                          rt(Tee[i] - shock_time_vec[i] - 1, t_dist_v_param))
      Y[[i]] <- garchxSim(Tee[i], arch = arch_param, garch = garch_param, 
                          xreg =  as.matrix(shock_indicator),
                          innovations = GARCH_innov_vec, verbose = TRUE) 
    }
    
    else { 

      #Create volatility shock w*
      vol_shock_vec[i] <- rt(1, t_dist_v_param)
      vol_shock_mean <- 0
      vol_shock_var <- 0
        
        #Now add the design matrix to the list X
        X[[i]] <- VAR_process
        
        #Create GARCH model with shock(s)
        GARCH_innov_vec <- rt(Tee[i], t_dist_v_param)
        Y[[i]] <- garchxSim(Tee[i], arch = arch_param, garch = garch_param, 
                            innovations = GARCH_innov_vec, verbose = TRUE) 
    } #end conditionals that create vol shocks
    
    T_star_plus_1_return_vec[i] <- Y[[i]][,1][shock_time_vec[i]+1,]
    
    #Now we calculate the p-value for the volatility spike of length k
    
    #What follows in the next 6-10 lines is for fitting the series using the garchx function
    # indicator_vec <- as.matrix(c(rep(0,shock_time_vec[i]), rep(1,k)))
    # garch_1_1 <- garchx(Y[[i]][1:(shock_time_vec[i]+k),1], 
    #                     order = c(1,1), 
    #                     xreg = indicator_vec[1:(shock_time_vec[i]+k)],
    #                     control = list(eval.max = 1000, iter.max = 1000, rel.tol = 10^(-12), x.tol = 10^(-12)))
    # xreg_est <- round(coeftest(garch_1_1)[ dim(coeftest(garch_1_1))[1], 1],5)
    # xreg_p_value <- round(coeftest(garch_1_1)[ dim(coeftest(garch_1_1))[1], dim(coeftest(garch_1_1))[2]],5)
    # xreg <- c(xreg, c(xreg_est, xreg_p_value))
    
    #What follows in the next 6-10 lines is for fitting the series using the rugrarch function
    indicator_vec <- as.matrix(c(rep(0,shock_time_vec[i]), rep(1,k)))
    # ru_garch_spec <- ugarchspec(
    #                       list(model = "sGARCH", model = "apARCH", garchOrder = c(1, 1, 1),
    #                        external.regressors = indicator_vec[1:(shock_time_vec[i]+k)]),
    #                       
    #                       mean.model = list(armaOrder=c(0,0),include.mean=F), 
    #                       distribution.model="std") 
    
    ru_garch_spec <- ugarchspec(
                          list(model = "sGARCH", model = "apARCH", garchOrder = c(1, 1, 1),
                          external.regressors = indicator_vec),
                          mean.model = list(armaOrder=c(0,0),include.mean=F),
                          distribution.model="std")
    
    ctrl = list(RHO = 1,DELTA = 1e-8,MAJIT = 15000,MINIT = 15000,TOL = 1e-6)
    garch_1_1 <- ugarchfit(ru_garch_spec, Y[[i]][1:(shock_time_vec[i]+k),1],
                           solver = "solnp", solver.control = ctrl)
    
    xreg_est <- round( garch_1_1@fit$matcoef[ dim(garch_1_1@fit$matcoef)[1], 1] , 5)
    xreg_p_value <- round( garch_1_1@fit$matcoef[ dim(garch_1_1@fit$matcoef)[1],  dim(garch_1_1@fit$matcoef)[2]] , 5)
    xreg <- c(xreg, c(xreg_est, xreg_p_value))
    
    # xreg_est <- 99
    # xreg_p_value <- 99
    # xreg <- c(xreg, c(xreg_est, xreg_p_value))
    
  } #end loop for n+1 series
  
  #Now make xreg into a dataframe
  xreg <- data.frame(matrix(xreg, nrow = n+1, byrow = TRUE))
  
  ## Compute summary statistics for output
  level_shock_kurtosis <- gamma(5/level_GED_beta)*gamma(1/level_GED_beta)/( (gamma(3/level_GED_beta))**2 ) - 3 #https://en.wikipedia.org/wiki/Generalized_normal_distribution
  vol_shock_kurtosis <- -9999
  
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
                                     '\n shock est = ', round(xreg[i,1],3), ', pval = ',round(xreg[i,2],3),
                                     sep = ''), ylab = '100 * Daily Log-Return')
    abline(v = shock_time_vec[i] + 1, col = 'red')
    abline(h = 0, col = 'green')

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
                                     '\n shock est = ', round(xreg[i,1],3), ', pval = ',round(xreg[i,2],3),
                                     sep = ''), ylab = 'Sigma^2')
    abline(v = shock_time_vec[i] + 1, col = 'red')
  }
  
  #Items to return in a list
  return(list(X, Y, Tee,shock_time_vec, xreg))
}

# Here is the length of the vol shock we will use
k <- 15

output <- synth_vol_sim(n = 8, 
                        p = 6, 
                        arch_param = c(.23),
                        garch_param = c(.34),
                        asymmetry_param = c(.02),
                        level_model = c('M1','M21','M22','none')[4],
                        vol_model = c('M1','M21','M22','none')[3],
                        t_dist_v_param = 8, # affects variance and kurtosis of returns
                        sigma_x = 100 * .01, 
                        shock_time_vec = NULL, 
                        level_shock_length = 1,
                        vol_shock_length = k,
                        a = 190, 
                        b = 250, 
                        mu_eps_star = 100 * -.0925,
                        M22_mu_eps_star = 100 * .003, 
                        sigma_eps_star = 100 * .005,
                        mu_omega_star = 100 * .02,
                        M22_mu_omega_star = 100 * .02,
                        vol_shock_sd = 100 * .050,
                        level_GED_alpha = .05 * sqrt(2), 
                        level_GED_beta = 1.8)

# #Let's look at estimates and pvalues
# output[[5]]
# 
# #What is the shock time for y_1?
# output[[4]][1]
# 
# shock_time_of_y1 <- output[[4]][1]
# length_of_y1 <- output[[3]][[1]]
# 
# #Let us look at the output for the time series of interest
# output[[2]][[1]]
# 
# #Let us compare the model without the shock indicator and the model with the indicator
# no_indicator <- garchx(output[[2]][[1]][,1][c(1:(shock_time_of_y1+k)),], order = c(1,1))
# no_indicator
# coeftest(no_indicator)
# AIC(no_indicator)
# plot.ts(no_indicator$fitted)
# predict(no_indicator, n.ahead = 1, newxreg = matrix(1))
# 
# covariate_mat <- as.matrix(c(rep(0,shock_time_of_y1), rep(1,k), rep(0, length_of_y1 - shock_time_of_y1 - k) ))
# yes_indicator <- garchx(output[[2]][[1]][,1][c(1:(shock_time_of_y1+k)),], order = c(1,1), xreg = covariate_mat[c(1:(shock_time_of_y1+k)),])
# yes_indicator
# coeftest(yes_indicator)
# AIC(yes_indicator)
# plot.ts(yes_indicator$fitted)
# predict(yes_indicator, n.ahead = 1, newxreg = matrix(1), verbose = TRUE)
# 
# #let us look at the model prior to the shock
# mod5 <- garchx(output[[2]][[1]][,1][c(1:(shock_time_of_y1)),], order = c(1,1))
# coeftest(mod5)
# AIC(mod5)
# plot.ts(fitted(mod5))
# 
# #Look at covariates of time series of interest
# plot.ts(output[[1]][[1]])
# 
# #Look at GARCH process of time series of interest
# plot.ts(output[[2]][[1]])
# 
# #What is the mean of the squared y for time series of interest?
# mean(output[[2]][[1]][,1]**2)


#####################################################################

########### Now the optimization and fitting

#####################################################################

#First, an distance-based weighting method function written by Jilei Lin (PhD Student at GWU)
# this function returns the W^* estimated by synthetic control method (SCM)
dbw <- function(X, 
                Tstar, 
                scale = FALSE,
                sum_to_1 = TRUE,
                nonneg = TRUE,
                bounded_above = TRUE) { # https://github.com/DEck13/synthetic_prediction/blob/master/prevalence_testing/numerical_studies/COP.R
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
  
  # optimization and return statement
  
  # I have added two features
  # 1) The option to remove the sum-to-1 constraint
  # 2) The option to remove the nonnegativity constraint.
  
  #Thus I need if statements to implement these...
  if (sum_to_1 == TRUE) 
  {
    eq_constraint = function(W) sum(W) - 1 
    }
  else{
    eq_constraint = NULL
  }
  
  if (nonneg == TRUE) 
  {
    lower_bound = rep(0, n) 
  }
  else{
    lower_bound = NULL
  }
  

  
  if (bounded_above == TRUE) 
  {
    upper_bound = rep(1, n) 
  }
  else{
    upper_bound = NULL
  }
  
  object_to_return <- solnp(par = rep(1/n, n), 
                            fun = weightedX0, 
                            eqfun = eq_constraint, 
                            eqB = 0, # will use inequality constraints later
                            LB = lower_bound, UB = upper_bound, control = list(trace = 0))
  return(object_to_return$pars)
  
  #I added the return statement because an implicit return is bad coding form
} #end dbw function
#############################################################################

synth_vol_fit <- function(X,
                          Y,
                          T_star,
                          shock_est_vec,
                          shock_lengths)
{
  ## Doc String
  
  # synth_vol_fit: function that takes (n+1)*(p+1) time series AND a vector of 
  # shock times as input and outputs
  # 1) calculates a single weight vector w, 
  # 2) calculates a single fixed effects estimate vector omega*, 
  # 3) calculates the adjustment estimator vector \hat omega* for time series of interest
  # 4) calculates the volatility of time series of interest at T*+1,T*+2,...,T*+k (i.e. the prediction)
  # 5) calculates estimate of volatility on T*+1 for each series using each of three families, and 
  # 6) calculates the squared-error loss of the prediction
  
  # Estimation/control options
  # --Allow user to enter series of unequal lengths
  # --Allow user to enter a vector of integers corresponding to the number of days
  # the shock effect lasts for each outcome series
  # --Allow user to pick a uniform model for each series (e.g. GARCH(1,1)) OR a BIC-minimizing
  # model for each series (or mix and match).
  # --Allow user to pick error distribution - see ugarchspec
  
  ##Input
  # Y, a list of length n+1, with each entry containing a time series
  # X, a list of length n+1, with each entry containing a dataframe of dimension y_i x p
  # shock_time_vec, a vector of length n+1 containing shock time of each series
  # shock_time_lengths, a vector of length n+1 containing shock time length of each series
  
  #First, we get the vectors w for all 7 sensible methods
  w <- list() #initialize
  matrix_of_specs <- matrix(c(rep(TRUE,4), 
                              rep(FALSE,4), 
                              rep(c(TRUE,TRUE,FALSE,FALSE),2), 
                              rep(c(TRUE,FALSE), 4)), 
                              byrow = FALSE, nrow = 8)
  
  #We drop the second row because it's functionally no different from the first
  matrix_of_specs <- matrix_of_specs[-2,]
  
  for (i in 1:nrow(matrix_of_specs))
    {
  w[[i]] <- dbw(X, 
                 T_star, 
                 scale = TRUE,
                 sum_to_1 = matrix_of_specs[i,1],
                 nonneg = matrix_of_specs[i,2],
                 bounded_above = matrix_of_specs[i,3])
  }
  
  # Now we place these linear combinations into a matrix
  w_mat <- matrix(unlist(w), nrow = nrow(matrix_of_specs), byrow = TRUE)
  
  #Second, we calculate omega_star_hat, which is the dot product of w and the estimated shock effects
  omega_star_hat_vec <- as.numeric(w_mat %*% shock_est_vec[-1])
  
  #Third, we get a prediction to T*_+1 
  data_up_through_T_star <- Y[[1]][,1][1:T_star[1],1]
  sigma2_up_through_T_star <- Y[[1]][,3][1:T_star[1],1]
  garch_1_1 <- garchx(data_up_through_T_star, order = c(1,1))
  pred <- as.numeric(predict(garch_1_1, n.ahead = shock_lengths[1]))
  adjusted_pred_vec <- pred + omega_star_hat_vec
  
  #If any predictions are negative, 0 them out
  adjusted_pred_vec <- pmax(adjusted_pred_vec,0)
  
  print(predict(garch_1_1, n.ahead = shock_lengths[1], verbose = TRUE))
  
  #Fourth, we calculate the ground truth of vol in the k-length period of time series of interest
  ground_truth_T_star_plus_1 <- as.numeric(Y[[1]][,3][(T_star[1] + 1),])
  
  #Last, we calculate MSE for first method
  MSE_adjusted <- (ground_truth_T_star_plus_1 - adjusted_pred_vec)**2
  MSE_unadjusted <- (ground_truth_T_star_plus_1 - pred)**2
  alternative_wins <- MSE_adjusted < MSE_unadjusted
  
  #We now make a vector with the names of each of the 7 sensible linear combinations
  linear_comb_names <- c('Convex Hull',
                        'Drop Bounded Below',
                        'Affine Hull',
                        'Drop Sum-to-1',
                        'Conic Hull',
                        'Bounded Above',
                        'Unrestricted')
  
  #Plot the donor pool weights
  par(mfrow=c(3,3))
  for (i in 1:nrow(w_mat))
    {
    barplot(w_mat[i,], 
            main = paste('Donor Pool Weights:\n', linear_comb_names[i]), names.arg = 2:(length(T_star)),
            ylim = c(min(w_mat[i,]),max(w_mat[i,])))
     }

  #Now let's plot the adjustment
  par(mfrow=c(1,2))
  
  trimmed_prediction_vec_for_plotting <- Winsorize(adjusted_pred_vec, probs = c(0, 0.6))
  
  plot(sigma2_up_through_T_star, 
       main = 'GARCH Prediction versus \nAdjusted Predictions versus Actual',
       ylab = '',
       xlab = "Time",
       xlim = c(0, length(data_up_through_T_star) + 5),
       ylim = c(min(0,trimmed_prediction_vec_for_plotting),  max(pred, trimmed_prediction_vec_for_plotting, data_up_through_T_star, Y[[1]][,3][T_star[1]+1],1) ) )
  title(ylab = expression(sigma^2), line = 2.05, cex.lab = 1.99)            # Add y-axis text
  
  lines(y = c(sigma2_up_through_T_star[T_star[1]],  Y[[1]][,3][T_star[1]+1,1]) , 
        x = c(T_star[1], T_star[1] + 1),  lty=2, lwd=2,
        ylim = c(min(0,trimmed_prediction_vec_for_plotting),  max(pred, trimmed_prediction_vec_for_plotting, data_up_through_T_star, Y[[1]][,3][T_star[1]+1],1) ) )

  colors_for_adjusted_pred <- c('black',
                                brewer.pal(length(adjusted_pred_vec) + 1,'Set1'))
  
  points(y = ground_truth_T_star_plus_1, 
         x = T_star[1] + 1, 
         col = colors_for_adjusted_pred[1], 
         cex = 1.3, pch = 16)
  
  points(y = pred, 
         x = (T_star[1] + 1), 
         col = colors_for_adjusted_pred[2], 
         cex = 1.3, pch = 15)

  for (i in 1:(length(adjusted_pred_vec)))
    {
           points(y = adjusted_pred_vec[i], x = (T_star[1] + 1), 
           col = colors_for_adjusted_pred[i+2], cex = 1.9, pch = 10)
  }
  
  labels_for_legend <- c('Actual','GARCH',linear_comb_names)
  
  legend(x = "topleft",  # Coordinates (x also accepts keywords)
         legend = labels_for_legend,
         1:length(labels_for_legend), # Vector with the name of each group
         colors_for_adjusted_pred,   # Creates boxes in the legend with the specified colors
         title = 'Prediction Method',      # Legend title,
         cex = 1.1
  )
  
  plot.ts(fitted(garch_1_1), 
          main = 'Pre-shock GARCH fitted values (green) \nversus Actual (black)',
          ylab = '', col = 'green',
          cex.lab = 3.99)
  lines(sigma2_up_through_T_star, col = 'black')
  title(ylab = expression(sigma^2), line = 2.05, cex.lab = 1.99)            # Add y-axis text
  
  
  return(list(w = round(w_mat,3), 
              omega_star_hat = round(omega_star_hat_vec, 3),
              adjusted_pred = round(adjusted_pred_vec,3),
              garch_pred = round(pred,3),
              ground_truth = round(ground_truth_T_star_plus_1,3), 
              MSE_adjusted = round(MSE_adjusted,3),
              MSE_unadjusted = round(MSE_unadjusted,3),
              alternative_wins = alternative_wins))
}

# Let's now use the function
X_demo <- output[[1]]
Y_demo <- output[[2]]
T_star_demo <- output[[4]]
shock_effect_vec_demo <- output[[5]][,1]
shock_time_lengths_demo <- c(1,rep(0, length(shock_effect_vec_demo) - 1))

synth_vol_fit(X_demo, 
              Y_demo, 
              T_star_demo, 
              shock_effect_vec_demo,
              shock_time_lengths_demo)
