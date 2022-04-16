

########### Items to investigate ########
# 1) using foreach just once to smash through grid+nsim 
# https://cran.r-project.org/web/packages/foreach/vignettes/foreach.html
# 2) benchmarking the runtime https://www.alexejgossmann.com/benchmarking_r/
# 3) this has a grid: https://www.blasbenito.com/post/02_parallelizing_loops_with_r/
# 4) https://www.r-bloggers.com/2013/06/grid-search-for-free-parameters-with-parallel-computing/
# 5) observe here the nested foreach structure https://www.r-bloggers.com/2013/06/grid-search-for-free-parameters-with-parallel-computing/
# Discuss this with AM.  There is a decision to be made between a nested for each setup aAND a replication column
# 6) Principle: number of cores should be equal to the number of threads https://www.blasbenito.com/post/02_parallelizing_loops_with_r/
# 7) Use chunking in this explanatory page, as well as nested for each discussion
# https://cran.r-project.org/web/packages/foreach/vignettes/nested.html

# MC
library("parallel")
library("doParallel")
library("foreach")
source("sim_source_data_and_fitting.R")

registerDoParallel(cores = ncores)
set.seed(13)
RNGkind("L'Ecuyer-CMRG")

nsim <- 50

############ We build our parameter grid ############ 
donor_pool_size <- c(10,20)
p <- c(5,10, 15)
alpha <- c(0, .15,.65, .75)
beta <- c(0, .15,.65, .75) 
vol_model <- c('M1','M21','M22')
level_model <- c('M1','M21','M22','none')
vol_shock_length <- c(1,2,3)
level_shock_length <- c(1,2)
extra_measurement_days <- c(1,2,3)
replication_number <- seq(1, nsim, 1)
optimization_norm <- c('l1','l2')
mu_eps_star <- c(-6, -10)
level_GED_alpha <- c(sqrt(2), sqrt(5)) # note: beta = 2, alpha = sqrt(2) is N(0,1
level_GED_beta <- c(.7, 2) # note: beta = 2, alpha = sqrt(2) is N(0,1))
M21_M22_level_mu_delta <- c(.6, .9)
M21_M22_level_sd_delta <- c(.4, .6)
mu_omega_star <- c(.5, 1.5, 3)
vol_shock_sd <- c(.2, .4, .6)
M21_M22_vol_mu_delta <- c(.3, .6, .9)
M21_M22_vol_sd_delta <- c(.2, .4, .6)

list_of_vars <- list(donor_pool_size
                    , p
                    , alpha
                    , beta
                    , vol_model
                    , level_model
                    , vol_shock_length
                    , level_shock_length
                    , extra_measurement_days
                    , replication_number
                    , optimization_norm
                    , mu_eps_star
                    , level_GED_alpha
                    , level_GED_beta
                    , M21_M22_level_mu_delta
                    , M21_M22_level_sd_delta
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
                           , 'level_shock_length'
                           , 'extra_measurement_days'
                           , 'replication_number'
                           , 'mu_eps_star'
                           , 'level_GED_alpha'
                           , 'level_GED_beta'
                           , 'M21_M22_level_mu_delta'
                           , 'M21_M22_level_sd_delta'
                           , 'mu_omega_star'
                           , 'vol_shock_sd'
                           , 'M21_M22_vol_mu_delta'
                           , 'M21_M22_vol_sd_delta')

# Because of the many many combinations, we have to subset even at the time of creation
gridd <- expand.grid(list_of_vars)

# Now we get rid of arrangements where alpha + beta >= 1
gridd_subset <- gridd[gridd$alpha + gridd$beta < 1,]

# Get rid of models where level shock is long than vol shock
gridd_subset <- gridd_subset[gridd_subset$level_shock_length <= gridd_subset$vol_shock_length,]

# Get rid of arrangement with no vol shock and no level shock
gridd_subset <- gridd_subset[gridd_subset$vol_model < 4 || gridd_subset$level_model < 4,]

# Now we get rid of unnecessary rows owing to M2vol and M2level
keep_condition_1 <- (gridd_subset$vol_model == 1) & 
  (gridd_subset$M21_M22_level_mu_delta == min(gridd_subset$M21_M22_level_mu_delta)) & 
  (gridd_subset$M21_M22_level_sd_delta == min(gridd_subset$M21_M22_level_sd_delta)) & 
  (gridd_subset$M21_M22_vol_mu_delta == min(gridd_subset$M21_M22_vol_mu_delta)) & 
  (gridd_subset$M21_M22_vol_sd_delta == min(gridd_subset$M21_M22_vol_sd_delta))

keep_condition_2 <- (gridd_subset$level_model == 1) & 
  (gridd_subset$M21_M22_level_mu_delta == min(gridd_subset$M21_M22_level_mu_delta)) & 
  (gridd_subset$M21_M22_level_sd_delta == min(gridd_subset$M21_M22_level_sd_delta)) & 
  (gridd_subset$M21_M22_vol_mu_delta == min(gridd_subset$M21_M22_vol_mu_delta)) & 
  (gridd_subset$M21_M22_vol_sd_delta == min(gridd_subset$M21_M22_vol_sd_delta))

gridd_subset <- gridd_subset[keep_condition_1 & keep_condition_2,]

nrow(gridd_subset) 
head(gridd_subset)

sim_params <- gridd_subset

############ end of parameter grid construction ############ 


########################### Begin parallel architecture ###############################

# simulation time
system.time(
  output_T_100 <- foreach(
    n = sim_params$donor_pool_size,
    p = sim_params$p,
    alpha = sim_params$alpha,
    beta = sim_params$beta,
    vol_model = sim_params$vol_model,
    level_model = sim_params$level_model,
    vol_shock_length = sim_params$vol_shock_length,
    level_shock_length = sim_params$level_model,
    .combine = 'rbind'
  ) %dopar% {
    
  #function simulate_and_analyze goes here
    
    #returning prediction error as percentage
    return(m.i$prediction.error * 100)
    
  }
) #end system.time

# simulation time
system.time(
  output_T_100 <- lapply(1:nrow(sim_params), FUN = function(j) {
    
    # parameters to populate
    shape.K.T <- sim_params[j, 1]
    n <- sim_params[j, 2]
    
    # %do% evaluates sequentially
    # %dopar% evaluates in parallel
    # .combine results
    out <- foreach(k = 1:nsim, .combine = rbind) %dopar% {
      # result
      
      # We force 
      Sys.sleep(30)
      
      study <- sim.normal.gammaX.decay(mu.gamma.delta = 2, 
                                       mu.alpha = 10, sigma = 0.1, 
                                       sigma.alpha = 0.05, 
                                       sigma.delta.gamma = 0.1, 
                                       p = 13, B = 1000, scale = 2, 
                                       n = n, H = 8, ell = 4,
                                       Kshape = shape.K.T, Tshape = shape.K.T)
      return(study)
    }
    # return results
    out
  })
)
########################### End parallel architecture ###############################

# Save output
save(output_T_100, file = "output_decay_T_100.RData")

stopImplicitCluster()


