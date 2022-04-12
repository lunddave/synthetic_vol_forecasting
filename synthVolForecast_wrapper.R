# David Lundquist
# Wrapper for Simulations for Synthetic Prediction GARCH

#We will need the two functions here.
source('Desktop/synthetic_vol_forecasting/synthVolForecast.R')

library(quantmod)
library(garchx)
library(lmtest)
library(extraDistr)
library(gnorm)
library(tsDyn)
library(Rsolnp)
library(RColorBrewer)
library(DescTools)

simulate_and_analyze <- function(n = 19, 
                                 p = 9, 
                                 #model = c(1,1,1),
                                 arch_param = c(.26),
                                 garch_param = c(.6),
                                 #asymmetry_param = c(.15),
                                 
                                 level_model = c('M1','M21','M22','none')[4],
                                 vol_model = c('M1','M21','M22','none')[2],
                                 
                                 sigma_GARCH_innov = 1, # the sd that goes into rnorm
                                 sigma_x = 1, # the sd that goes into the covariates
                                 
                                 min_shock_time = 5,
                                 shock_time_vec = NULL, 
                                 
                                 level_shock_length = 1,
                                 vol_shock_length = 4,
                                 extra_measurement_days = 2,
                                 
                                 a = 3*252, 
                                 b = 10*252, 
                                 
                                 mu_eps_star = -4.25,
                                 level_GED_alpha = sqrt(2), # note: beta = 2, alpha = sqrt(2) is N(0,1)
                                 level_GED_beta = 2, # note: beta = 2, alpha = sqrt(2) is N(0,1))
                                 
                                 M21_M22_level_mu_delta = .2, 
                                 M21_M22_level_sd_delta = .1,
                                 
                                 mu_omega_star = .005,
                                 vol_shock_sd = .005,
                                 
                                 M21_M22_vol_mu_delta = .04,
                                 M21_M22_vol_sd_delta = .01, 
                                 
                                 plot_sim = TRUE,
                                 
                                 plot_fit = TRUE,
                                 
                                 # And now the only inputs for the fitting function
                                 inputted_vol_shock_length = rep(4, n+1),
                                 normchoice = 'l1'
                                 ) 
{
  ## Doc String
  
  # simulate_and_analyze: function wraps two functions: 
  # (1) synth_vol_sim
  # (2) synth_vol_fit
  
  # --Input:
  #   --n - number of donors (scalar)
  #   --model - model order for a parameters to be simulated (optional).  
  #     If specified, then arch, garch parameter values overridden.
  #   --p - number of covariates (scalar)
  #   --arch parameters (vector of length 0 or more)
  #   --garch parameters (vector of length 0 or more)
  #   --asymmetry_param (vector of length 0 or more)
  #   --level_model - model for level of the shock (string)
  #   --vol_model - model for volatility of the shock (string)
  #   --sigma_GARCH_innov (scalar)
  #   --sigma_x - sigma of innovations in covariates (scalar)
  #   --min_shock_time - the minimum number of points after the min series length 'a' that the shock can occur at
  #   --shock_time_vec - optional input to force shock times (vector of integers)
  #   --level_shock_length (scalar)
  #   --vol_shock_length (scalar)
  #   --a - minumum series length (scalar)
  #   --b - maximum series length (scalar)
  
  #   --mu_eps_star - intercept for each of the level shock models
  #   --M21_M22_mu_delta - mean of delta for M21, M22 level models
  #   --M21_M22_level_sd_delta - sd of the vector delta in M21 and M22 level models
  
  #   --mu_omega_star - intercept of vol shock for M1 model
  #   --M21_M22_mu_omega_star - mean of delta for M21, M22 vol models
  #   --vol_shock_sd - variance of the error in all volatility models
  #   --M21_M22_vol_sd_delta - sd of the vector delta in M21 and M22 vol models
  
  #   --level_GED_alpha - alpha parameter for level shock stochastic term
  #   --level_GED_beta - beta parameter for level shock stochastic term
  
  #   --inputted_vol_shock_length - vector of length n+1 referring to the shocks lengths to be used in the estimation process
  
  sim_output <- synth_vol_sim(n = n, 
                          p = p, 
                          #model = c(1,1,1),
                          arch_param = arch_param,
                          garch_param = garch_param,
                          #asymmetry_param = c(.15),
                          
                          level_model = level_model,
                          vol_model = vol_model,
                          
                          sigma_GARCH_innov = sigma_GARCH_innov, # the sd that goes into rnorm
                          sigma_x = sigma_x, # the sd that goes into the covariates
                          min_shock_time = min_shock_time,
                          shock_time_vec = shock_time_vec, 
                          level_shock_length = level_shock_length,
                          vol_shock_length = vol_shock_length,
                          a = a, 
                          b = b, 
                          
                          mu_eps_star = mu_eps_star,
                          level_GED_alpha = level_GED_alpha, # note: beta = 2, alpha = sqrt(2) is N(0,1)
                          level_GED_beta = level_GED_beta, # note: beta = 2, alpha = sqrt(2) is N(0,1)
                          
                          M21_M22_level_mu_delta = M21_M22_level_mu_delta, 
                          M21_M22_level_sd_delta = M21_M22_level_sd_delta,
                          
                          mu_omega_star = mu_omega_star,
                          vol_shock_sd = vol_shock_sd, 
                          M21_M22_vol_mu_delta = M21_M22_vol_mu_delta,
                          M21_M22_vol_sd_delta = M21_M22_vol_sd_delta, 

                          plot = plot_sim)
  
  # Let's now use the fitting function
  X_demo <- sim_output[[1]]
  Y_demo <- sim_output[[2]]
  T_star_demo <- sim_output[[4]]
  shock_effect_vec_demo <- sim_output[[5]][,1]
  garch_order_of_simulation <- sapply(sim_output[[6]], length) 
  
  fitting_output <- synth_vol_fit(X = X_demo,
                                  Y = Y_demo,
                                  T_star = T_star_demo,
                                  shock_est_vec = shock_effect_vec_demo,
                                  shock_lengths = inputted_vol_shock_length,
                                  garch_order_of_simulation[1],
                                  garch_order_of_simulation[2],
                                  garch_order_of_simulation[3],
                                  normchoice = normchoice,
                                  plots = plot_fit
                                  )
  return(fitting_output)
}

simulate_and_analyze()
