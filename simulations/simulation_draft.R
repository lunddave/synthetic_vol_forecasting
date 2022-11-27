
########### Items to investigate ########
# 1) using foreach just once to smash through grid+nsim 
# https://cran.r-project.org/web/packages/foreach/vignettes/foreach.html
# 2) benchmarking the runtime https://www.alexejgossmann.com/benchmarking_r/
# 3) this has a grid: https://www.blasbenito.com/post/02_parallelizing_loops_with_r/
# 4 observe here the nested foreach structure https://www.r-bloggers.com/2013/06/grid-search-for-free-parameters-with-parallel-computing/
# Discuss this with AM.  There is a decision to be made between a nested for each setup AND a replication column
# 6 Use chunking in this explanatory page, as well as nested for each discussion
# https://cran.r-project.org/web/packages/foreach/vignettes/nested.html

# Links only:
#   https://cran.r-project.org/web/packages/foreach/vignettes/foreach.html
#   https://www.alexejgossmann.com/benchmarking_r/
#   https://www.blasbenito.com/post/02_parallelizing_loops_with_r/
#   https://www.r-bloggers.com/2013/06/grid-search-for-free-parameters-with-parallel-computing/
#   https://cran.r-project.org/web/packages/foreach/vignettes/nested.html

# MC
library("parallel")
library("doParallel")
library("foreach")
source("/home/david/Desktop/synthetic_vol_forecasting/synthVolForecast_wrapper.R",
       echo = FALSE,
       verbose = FALSE)

registerDoParallel(cores = detectCores() - 3)
set.seed(13) #tk do we want to vary this?
RNGkind("L'Ecuyer-CMRG")

nsim <- 1

############ We build our parameter grid ############ 
donor_pool_size <- c(3,6)
p <- c(3,6)
alpha <- c(.2,.4)
beta <- c(.2, .4)
vol_model <- c('M1','M21','M22')
level_model <- c('M1','M21','M22','none')[4]
vol_shock_length <- c(1,2,3)
level_shock_length <- c(1,2)[-1]
extra_measurement_days <- c(1,2)
replication_number <- seq(1, nsim, 1)
optimization_norm <- c('l1','l2')
#mu_eps_star <- c(-6, -10)
#level_GED_alpha <- c(sqrt(2), sqrt(5)) # note: beta = 2, alpha = sqrt(2) is N(0,1
#level_GED_beta <- c(.7, 2) # note: beta = 2, alpha = sqrt(2) is N(0,1))
#M21_M22_level_mu_delta <- c(.6, .9)
#M21_M22_level_sd_delta <- c(.4, .6)
mu_omega_star <- c(.05, .15, .33)
vol_shock_sd <- c(.01, .02, .05)
M21_M22_vol_mu_delta <- c( .06, .09)
M21_M22_vol_sd_delta <- c( .04, .06)


list_of_vars <- list(donor_pool_size
                    , p
                    , alpha
                    , beta
                    , vol_model
                    , level_model
                    , vol_shock_length
                    #, level_shock_length
                    , extra_measurement_days
                    , optimization_norm
                    #, mu_eps_star
                    #, level_GED_alpha
                    #, level_GED_beta
                    #, M21_M22_level_mu_delta
                    #, M21_M22_level_sd_delta
                    , mu_omega_star
                    , vol_shock_sd
                    , M21_M22_vol_mu_delta
                    , M21_M22_vol_sd_delta)

names(list_of_vars) <- list('donor_pool_size'
                           , 'p'
                           , 'alpha'
                           , 'beta'
                           , 'vol_model'
                           , 'level_model'
                           , 'vol_shock_length'
                           #, 'level_shock_length'
                           , 'extra_measurement_days'
                           , 'optimization_norm'
                           #, 'mu_eps_star'
                           #, 'level_GED_alpha'
                           #, 'level_GED_beta'
                           #, 'M21_M22_level_mu_delta'
                           #, 'M21_M22_level_sd_delta'
                           , 'mu_omega_star'
                           , 'vol_shock_sd'
                           , 'M21_M22_vol_mu_delta'
                           , 'M21_M22_vol_sd_delta')

# Because of the many many combinations, we have to subset even at the time of creation
gridd <- expand.grid(list_of_vars)

# Now we get rid of arrangements where alpha + beta >= 1
gridd_subset <- gridd[gridd$alpha + gridd$beta < 1,]

# Get rid of models where level shock is long than vol shock
#gridd_subset <- gridd_subset[gridd_subset$level_shock_length <= gridd_subset$vol_shock_length,]

# # Now we get rid of unnecessary rows owing to M2vol and M2level
# keep_condition_1 <- (gridd_subset$vol_model == 'M1') & 
#   (gridd_subset$M21_M22_level_mu_delta == min(gridd_subset$M21_M22_level_mu_delta)) & 
#   (gridd_subset$M21_M22_level_sd_delta == min(gridd_subset$M21_M22_level_sd_delta)) & 
#   (gridd_subset$M21_M22_vol_mu_delta == min(gridd_subset$M21_M22_vol_mu_delta)) & 
#   (gridd_subset$M21_M22_vol_sd_delta == min(gridd_subset$M21_M22_vol_sd_delta))
# 
# keep_condition_2 <- (gridd_subset$level_model == 'M1') & 
#   (gridd_subset$M21_M22_level_mu_delta == min(gridd_subset$M21_M22_level_mu_delta)) & 
#   (gridd_subset$M21_M22_level_sd_delta == min(gridd_subset$M21_M22_level_sd_delta)) & 
#   (gridd_subset$M21_M22_vol_mu_delta == min(gridd_subset$M21_M22_vol_mu_delta)) & 
#   (gridd_subset$M21_M22_vol_sd_delta == min(gridd_subset$M21_M22_vol_sd_delta))
# 
# gridd_subset <- gridd_subset[keep_condition_1 & keep_condition_2,]

nrow(gridd_subset) 
head(gridd_subset, n = 1)

sim_params <- gridd_subset

############ end of parameter grid construction ############ 


########################### Begin parallel architecture ###############################

# simulation time
system.time( 
  
              output_n_sim_1 <- foreach(
                                        n = sim_params$donor_pool_size, 
                                        p = sim_params$p, 
                                        #model = c(1,1,1),
                                        arch_param = sim_params$alpha,
                                        garch_param = sim_params$beta,
                                        #asymmetry_param = c(.15),
                                        
                                        level_model = sim_params$level_model,
                                        vol_model = sim_params$vol_model,
                                        
                                        level_shock_length = sim_params$level_shock_length,
                                        vol_shock_length = sim_params$vol_shock_length ,
                                        
                                        mu_omega_star = sim_params$mu_omega_star,
                                        vol_shock_sd = sim_params$vol_shock_sd, 
                                        M21_M22_vol_mu_delta = sim_params$M21_M22_vol_mu_delta,
                                        M21_M22_vol_sd_delta = sim_params$M21_M22_vol_sd_delta,
                                        
                                        #Now we choose how we want foreach to combine the output of each sim
                                        .combine = 'rbind',
                                        .errorhandling = "remove" #pass is another option
                            
                            ) %dopar% {
            
                            foreach(1:nsim, .combine='rbind') %do% { #begin inner loop 
    
                            to_return <- simulate_and_analyze(n = n,
                                                              p = p,
                                                              #model = c(1,1,1),
                                                              arch_param = arch_param,
                                                              garch_param = garch_param,
                                                              #asymmetry_param = c(.15),
                
                                                              level_model = level_model,
                                                              vol_model = vol_model,
                
                                                              level_shock_length = level_shock_length,
                                                              vol_shock_length = vol_shock_length,

                                                              mu_omega_star = mu_omega_star,
                                                              vol_shock_sd = vol_shock_sd,
                                                              M21_M22_vol_mu_delta = M21_M22_vol_mu_delta,
                                                              M21_M22_vol_sd_delta = M21_M22_vol_sd_delta
                                                              )
            
                              return(to_return)
            
                                                              } # end inner loop
            
              } # end outer loop
) #end system.time

########################### End parallel architecture ###############################

# Save output
save(output_n_sim_1, file = "/home/david/Desktop/synthetic_vol_forecasting/simulation_results/output_n_sim_1")

save(output_n_sim_1, file = paste("simulation_results/simcount_",nsim,"_savetime_",Sys.time(),sep="") )


stopImplicitCluster()

View(output_n_sim_1)
