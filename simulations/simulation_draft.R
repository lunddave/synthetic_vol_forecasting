

# MC
library("parallel")
library("doParallel")
library("foreach")
source("sim_source_data_and_fitting.R")

registerDoParallel(cores = ncores)
set.seed(13)
RNGkind("L'Ecuyer-CMRG")

nsim <- 100

############ We build our parameter grid ############ 
donor_pool_size <- c(10,20,30)
p <- c(5,10,15)
alpha <- c(.15,.3,.65, .75)
beta <- c(.15,.3,.65, .75)
vol_model <- seq(1,4,1)
level_model <- seq(1,4,1)
vol_shock_length <- c(1,2,3)
level_shock_length <- c(1,2)

list_of_vars <- list(donor_pool_size
                    , p
                    , alpha
                    , beta
                    , vol_model
                    , level_model
                    , vol_shock_length
                    , level_shock_length)

names(list_of_vars) <- list('donor_pool_size'
                           , 'p'
                           , 'alpha'
                           , 'beta'
                           , 'vol_model'
                           , 'level_model'
                           , 'vol_shock_length'
                           , 'level_shock_length')

gridd <- expand.grid(list_of_vars)

# Now we get rid of arrangements where alpha + beta >= 1
gridd_subset <- gridd[gridd$alpha + gridd$beta < 1,]

# Get rid of arrangement with no vol shock and no level shock
gridd_subset <- gridd_subset[gridd_subset$vol_model < 4 || gridd_subset$level_model < 4,]


nrow(gridd_subset) * sum(100, 200, 400)
head(gridd_subset)

sim_params <- gridd_subset

############ end of parameter grid construction ############ 

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


save(output_T_100, file = "output_decay_T_100.RData")





