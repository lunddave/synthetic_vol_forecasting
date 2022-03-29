

########### Items to investigate ########
# 1) using foreach just once to smash through grid+nsim 
# https://cran.r-project.org/web/packages/foreach/vignettes/foreach.html
# 2) benchmarking the runtime https://www.alexejgossmann.com/benchmarking_r/
# 3) this has a grid: https://www.blasbenito.com/post/02_parallelizing_loops_with_r/
# 4) https://www.r-bloggers.com/2013/06/grid-search-for-free-parameters-with-parallel-computing/
# 5) observe here the nested foreach structure https://www.r-bloggers.com/2013/06/grid-search-for-free-parameters-with-parallel-computing/
# Discuss this with AM.  There is a decision to be made between a nested for each setup aAND a replication column
# 6) Principle: number of cores should be equal to the number of threads https://www.blasbenito.com/post/02_parallelizing_loops_with_r/

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
alpha <- c(0, .15,.3,.65, .75)
beta <- c(0, .15,.3,.65, .75) 
vol_model <- seq(1,4,1)
level_model <- seq(1,4,1)
vol_shock_length <- c(1,2,3)
level_shock_length <- c(1,2)
replication_number <- seq(1, nsim, 1)

list_of_vars <- list(donor_pool_size
                    , p
                    , alpha
                    , beta
                    , vol_model
                    , level_model
                    , vol_shock_length
                    , level_shock_length
                    , replication_number)

names(list_of_vars) <- list('donor_pool_size'
                           , 'p'
                           , 'alpha'
                           , 'beta'
                           , 'vol_model'
                           , 'level_model'
                           , 'vol_shock_length'
                           , 'level_shock_length'
                           , 'replication_number')

gridd <- expand.grid(list_of_vars)

# Now we get rid of arrangements where alpha + beta >= 1
gridd_subset <- gridd[gridd$alpha + gridd$beta < 1,]

# Get rid of arrangement with no vol shock and no level shock
gridd_subset <- gridd_subset[gridd_subset$vol_model < 4 || gridd_subset$level_model < 4,]


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


