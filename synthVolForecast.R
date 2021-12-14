# David Lundquist
# Simulations for Synthetic Prediction GARCH

library(quantmod)
library(garchx)
library(lmtest)
library(mlVAR)
library(extraDistr)
library(gnorm)

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

synth_vol_sim <- function(n, p, arch_param, garch_param, model, 
                          sigma_GARCH_innov, sigma_x, 
                          shock_time_vec, 
                          length_of_shock,
                          a, b, 
                          mu_eps_star, sigma_eps_star,
                          level_GED_alpha,
                          level_GED_beta,
                          vol_GED_alpha,
                          vol_GED_beta,
                          ...){
  # We will first generate the covariates.  These will be correlated GARCH processes, ideally. 
  # Since multivariate GARCH processes take take technical care to simulate, we first use VAR.
  
  # https://stats.stackexchange.com/questions/71540/how-to-simulate-two-correlated-ar1-time-series
  # https://cran.r-project.org/web/packages/mlVAR/mlVAR.pdf
  
  # Simulate series lengths
  Tee <- rdunif(n+1, a, b)
  
  # Simulate shock times
  if ( is.null(shock_time_vec) == TRUE)
  {
    shock_time_vec <- c()
    for (i in 1:(n+1))
    {
      shock_time_vec[i] <- rdunif(1, 20, Tee[i]- 1)
    }
  }
  
  ############simulate all n+1 series
  
  #Create null lists for depvar and indepvar output
  Y <- vector(mode = "list", length = n+1)
  X <- vector(mode = "list", length = n+1)
  level_shock_vec <- c()
  vol_shock_vec <- c()
  
  # For each of n+1 series...
  for (i in 1:(n+1)){
    #simulate covariates
    matrix_entries <- runif(p**2, min = -.05, max = .05)
    simVAR_params <- matrix(matrix_entries, nrow = p, byrow = T)
    #https://cran.r-project.org/web/packages/mlVAR/mlVAR.pdf
    VAR_process <- simulateVAR(simVAR_params, means = 0, lags = 1, Nt = Tee[i], residuals = sigma_x)
    vol_shock_vec[i] <- rgnorm(1, mu = 0, alpha = vol_GED_alpha, beta = vol_GED_beta)**2
    shock_indicator <- c(
                          rep(0, shock_time_vec[i]), 
                          vol_shock_vec[i], 
                          rep(0, Tee[i] - shock_time_vec[i] - length_of_shock)
                        )
    X[[i]] <- cbind(VAR_process, shock_indicator)
    
    #Now create shocks
    if (model == 'M1'){
                      level_shock_vec[i] <- mu_eps_star + rgnorm(1, mu = 0, alpha = level_GED_alpha, beta = level_GED_beta)
                      } 
    else {
          level_shock_vec[i] <- mu_eps_star + 
                          as.matrix(X[[1]][shock_time_vec[1],]) %*% rnorm(p,0,sigma_eps_star) + 
                          rgnorm(1, mu = 0, alpha = level_GED_alpha, beta = level_GED_beta)
          }
    
    #Create GARCH model with one shock
    GARCH_innov_vec <- c(
                         rnorm(shock_time_vec[i], 0, sigma_GARCH_innov), 
                         level_shock_vec[i],
                         rnorm(Tee[i] - shock_time_vec[i] - 1, 0, sigma_GARCH_innov))
    
    Y[[i]] <- garchxSim(Tee[i], arch = arch_param, garch = garch_param, xreg = as.matrix(VAR_process),
                        innovations = GARCH_innov_vec, verbose = TRUE) #First column
  }
  
  ## Compute summary statistics for output
  shock_mean <- mu_eps_star 
  shock_var <- ((level_GED_alpha)**2) * gamma(3/level_GED_beta) / (gamma(1/level_GED_beta)) # https://search.r-project.org/CRAN/refmans/gnorm/html/gnorm.html
  shock_kurtosis <- gamma(5/level_GED_beta)*gamma(1/level_GED_beta)/( (gamma(3/level_GED_beta))**2 ) - 3 #https://en.wikipedia.org/wiki/Generalized_normal_distribution
  
  T_star_sigma <- Y[[1]][,3][shock_time_vec[1],]
  T_star_plus_1_sigma <- Y[[1]][,3][shock_time_vec[1]+1,]
  T_star_plus_2_sigma <- Y[[1]][,3][shock_time_vec[1]+2,]
  
  ##Output
  cat('Simulation Summary Data','\n',
      '-------------------------------------------------------------\n',
      'Donors:', n, '\n',
      'Series lengths:', Tee, '\n',
      'Shock times:', shock_time_vec, '\n',
      'Level Shock at T*+1:', round(level_shock_vec,3), '\n', 
      'Volatility Shock at T*+1', round(vol_shock_vec,3), '\n',
      'Sigma^2 at T*:', round(T_star_sigma,3), '\n', 
      'Sigma^2 at T*+1:', round(T_star_plus_1_sigma,3), '\n', 
      'Sigma^2 at T*+2:', round(T_star_plus_2_sigma,3), '\n', 
      '\n',
      'Shock Moments', '\n',
      '-------------------------------------------------------------\n',
      'Shock mean:', round(shock_mean,4), '(equivalent to a', round(100*shock_mean,2), '% daily move).', ' \n',
      'Shock variance:', round(shock_var,4), '\n',
      'Coefficient of Variation:', round(shock_mean/ sqrt(shock_var),3) , '\n',
      'Shock excess kurtosis', round(shock_kurtosis, 3)
      )
  
  #Plot time series under study and donors
  par(mfrow = c(ceiling(sqrt(n+1)), ceiling(sqrt(n+1))))
  for (i in 1:(n+1))
  {
    plot.ts(Y[[i]][,1], main = paste('y_', i, sep = ''), ylab = 'Daily Log-Return')
    abline(v = shock_time_vec[i] + 1, col = 'red')
  }
  
  par(mfrow = c(ceiling(sqrt(n+1)), ceiling(sqrt(n+1))))
  for (i in 1:(n+1))
  {
    plot.ts(Y[[i]][,3], main = paste('Volatility Series of y_', i, sep = ''), ylab = 'Sigma^2')
    abline(v = shock_time_vec[i] + 1, col = 'red')
  }
  
  #Items to return in a list
  return(list(X,Y,Tee))
}

output <- synth_vol_sim(n = 6, 
                        p = 3, 
                        arch_param = .55,
                        garch_param = c(.19),
                        model = c('M1','M2')[1],
                        sigma_GARCH_innov = .006,
                        sigma_x = .005, 
                        shock_time_vec = NULL, 
                        length_of_shock = 1,
                        a = 90, 
                        b = 150, 
                        mu_eps_star = -.055,
                        sigma_eps_star = .0005,
                        level_GED_alpha = .18, 
                        level_GED_beta = 1.4,
                        vol_GED_alpha = .4,
                        vol_GED_beta = 1.3)

plot.ts(output[[2]][[1]])

# Objective2: estimation function that takes (n+1)*(p+1) time series as input and 
# 1) calculates weight vector w, 
# 2) calculates fixed effects estimate vector omega*, 
# 3) calculates the adjusted estimate for the
# volatility of time series of interest at T*+1, 
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



