
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
library(dplyr)
library(parallel)
library(doParallel)
library(foreach)
source("/home/david/Desktop/synthetic_vol_forecasting/synthVolForecast_wrapper.R",
       echo = FALSE,
       verbose = FALSE)

registerDoParallel(6)
set.seed(13) #tk do we want to vary this?
RNGkind("L'Ecuyer-CMRG")

start_time <- Sys.time()
nsim <- 1
permutation_shift <- 0

############ We build our parameter grid ############ 
donor_pool_size <- c(2)
p <- c(2,4)
alpha <- c(.1,.2,.3)
beta <- c(.1,.2,.3)
vol_model <- c('M1','M21','M22')[2]
level_model <- c('M1','M21','M22','none')[4]
vol_shock_length <- c(1,2,3)
level_shock_length <- c(0)
extra_measurement_days <- c(1,2)
replication_number <- seq(1, nsim, 1)
optimization_norm <- c('l1','l2')[2]
mu_eps_star <- c(0, -6, -10)
#level_GED_alpha <- c(sqrt(2), sqrt(5)) # note: beta = 2, alpha = sqrt(2) is N(0,1)
#level_GED_beta <- c(.7, 2) # note: beta = 2, alpha = sqrt(2) is N(0,1))
M21_M22_level_mu_delta <- c(0, .6, .9)
M21_M22_level_sd_delta <- c(0, .4, .6)
mu_omega_star <- c(0,.03)
vol_shock_sd <- c(0,.005)
M21_M22_vol_mu_delta <- c(0,  .03)
M21_M22_vol_sd_delta <- c(0, .005)

list_of_vars <- list(donor_pool_size
                    , p
                    , alpha
                    , beta
                    , vol_model
                    , level_model
                    , vol_shock_length
                    , level_shock_length
                    , extra_measurement_days
                    , optimization_norm
                    , mu_eps_star
                    # #, level_GED_alpha
                    # #, level_GED_beta
                    , M21_M22_level_mu_delta
                    , M21_M22_level_sd_delta
                    , mu_omega_star
                    , vol_shock_sd
                    , M21_M22_vol_mu_delta
                    , M21_M22_vol_sd_delta
                    )

names(list_of_vars) <- list('donor_pool_size'
                           , 'p'
                           , 'alpha'
                           , 'beta'
                           , 'vol_model'
                           , 'level_model'
                           , 'vol_shock_length'
                           , 'level_shock_length'
                           , 'extra_measurement_days'
                           , 'optimization_norm'
                           , 'mu_eps_star'
                           # #, 'level_GED_alpha'
                           # #, 'level_GED_beta'
                           , 'M21_M22_level_mu_delta'
                           , 'M21_M22_level_sd_delta'
                           , 'mu_omega_star'
                           , 'vol_shock_sd'
                           , 'M21_M22_vol_mu_delta'
                           , 'M21_M22_vol_sd_delta'
                           )

# Because of the many many combinations, we have to subset even at the time of creation
gridd <- expand.grid(list_of_vars)

## Now we cull the expansion by removing combinations that do not make sense
## or will lead to failures.

# Now we get rid of arrangements where alpha + beta >= 1
gridd_subset <- gridd[gridd$alpha + gridd$beta < 1,]

# Get rid of models where level shock is long than vol shock
gridd_subset <- gridd_subset[gridd_subset$level_shock_length <= gridd_subset$vol_shock_length,]

# Get rid of M1 level models where

# level_shock_length
# , 'M21_M22_level_mu_delta'
# , 'M21_M22_level_sd_delta'

# is not zero:

gridd_subset <- gridd_subset[ (gridd_subset$level_model == 'M1' &
                                 gridd_subset$level_shock_length != 0 &
                                 gridd_subset$mu_eps_star != 0 &
                                 gridd_subset$M21_M22_level_mu_delta == 0 &
                                 gridd_subset$M21_M22_level_sd_delta == 0) |
                                
                                (gridd_subset$level_model %in% c('M21','M22') &
                                   gridd_subset$level_shock_length != 0 &
                                   gridd_subset$mu_eps_star != 0 &
                                   gridd_subset$M21_M22_level_mu_delta != 0 &
                                   gridd_subset$M21_M22_level_sd_delta != 0) |
                                
                                (gridd_subset$level_model == 'none' &
                                   gridd_subset$level_shock_length == 0 &
                                   gridd_subset$mu_eps_star == 0 &
                                   gridd_subset$M21_M22_level_mu_delta == 0 &
                                   gridd_subset$M21_M22_level_sd_delta == 0),]


# Get rid of M1 vol models where

# , 'M21_M22_vol_mu_delta'
# , 'M21_M22_vol_sd_delta'

# is not zero:

gridd_subset <- gridd_subset[ (gridd_subset$vol_model == 'M1' &
                                 gridd_subset$M21_M22_vol_mu_delta == 0 &
                                 gridd_subset$M21_M22_vol_sd_delta == 0) |
                                
                                (gridd_subset$vol_model != 'M1' &
                                   gridd_subset$M21_M22_vol_mu_delta != 0)
                              ,]

## tk TO DO

# a) fix the 1-14 names for geometric sets
# b) fix parameter combinations that don't make sense
# c) look into why some combinations often lead to GARCH estimation failure
# d) 
# 


# Take stock of what we have

sim_params <- gridd_subset


# This last bit of code will help us see how the parameter combos vary
length_unique <- function(x) { return(length(unique(x)))}

sim_params_check <- sim_params %>% group_by(vol_model) %>% 
  summarise(across(everything(), length_unique),
            .groups = 'drop')  %>%
  as.data.frame()

# Take stock of what we have

sim_params <- gridd_subset

grid_row_count <- nrow(sim_params)

############ end of parameter grid construction ############ 


########################### Begin parallel architecture ###############################

# simulation time
system.time( 
  
                            output <- foreach(
                                        n = sim_params$donor_pool_size
                                        , p = sim_params$p
                                        #,model = c(1,1,1),
                                        ,arch_param = sim_params$alpha
                                        ,garch_param = sim_params$beta
                                        #,asymmetry_param = c(.15)

                                        ,level_model = sim_params$level_model
                                        ,vol_model = sim_params$vol_model

                                        ,level_shock_length = sim_params$level_shock_length
                                        ,vol_shock_length = sim_params$vol_shock_length
                                        ,extra_measurement_days = sim_params$extra_measurement_days
                                        ,optimization_norm = sim_params$optimization_norm
                                        ,mu_eps_star = sim_params$mu_eps_star

                                        ,mu_omega_star = sim_params$mu_omega_star
                                        ,vol_shock_sd = sim_params$vol_shock_sd
                                        ,M21_M22_vol_mu_delta = sim_params$M21_M22_vol_mu_delta
                                        ,M21_M22_vol_sd_delta = sim_params$M21_M22_vol_sd_delta

                                        #Now we choose how we want foreach to combine the output of each sim
                                        ,.combine = 'rbind'
                                        ,.errorhandling = "remove" #pass is another option

                            ) %dopar% {
            
                            foreach(1:nsim
                                    ,
                                    #Now we choose how we want foreach to combine the output of each sim
                                    .combine = 'rbind'
                                    ,
                                    .errorhandling = "remove" #pass is another option
                                    ) %do% { #begin inner loop 
    
                                                        to_return <- simulate_and_analyze(n = n
                                                              ,p = p
                                                              #,model = c(1,1,1)
                                                              ,arch_param = arch_param
                                                              ,garch_param = garch_param

                                                              #asymmetry_param = c(.15),

                                                              ,level_model = level_model
                                                              ,vol_model = vol_model

                                                              ,level_shock_length = level_shock_length
                                                              ,vol_shock_length = vol_shock_length
                                                              ,extra_measurement_days = extra_measurement_days
                                                              ,normchoice = optimization_norm
                                                              ,mu_eps_star = mu_eps_star

                                                              ,mu_omega_star = mu_omega_star
                                                              ,vol_shock_sd = vol_shock_sd
                                                              ,M21_M22_vol_mu_delta = M21_M22_vol_mu_delta
                                                              ,M21_M22_vol_sd_delta = M21_M22_vol_sd_delta
                                                              ,permutation_shift = permutation_shift
                                                              )

                              return(to_return)
            
                                                              } # end inner loop
            
              } # end outer loop
) #end system.time

########################### End parallel architecture ###############################

# Save output

recovery_rate <- round( nrow(output) / (grid_row_count * nsim),3)

end_time <- Sys.time()

running_hours <- round(difftime(end_time, start_time, units="hours"),3)

# Save output
save(output, file = paste("/home/david/Desktop/synthetic_vol_forecasting/simulation_results/simcount_",
                          nsim,
                          "_savetime_", format(Sys.time(), "%a%b%d%X%Y"),"_runtime_",running_hours,
                          "_grid_size",grid_row_count,
                          "_recovery_",recovery_rate,
                          "_permute_",permutation_shift,
                          ".Rdata",sep="") )


stopImplicitCluster()

