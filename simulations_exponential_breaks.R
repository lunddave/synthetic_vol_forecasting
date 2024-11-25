
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
library(parallel)
library(doParallel)
library(foreach)
library(doRNG)
library(minpack.lm)

source("/Users/davidlundquist/Desktop/PhD/synthetic_vol_forecasting/exponential_breaks.R",
       echo = FALSE,
       verbose = FALSE)

command_args <- commandArgs(trailingOnly = TRUE)

registerDoParallel(6)
#set.seed(13) #tk do we want to vary this?
#RNGkind("L'Ecuyer-CMRG")

start_time <- Sys.time()
#nsim <- as.numeric(command_args[1])
nsim <- 1

############ We build our parameter grid ############
n <- c(10)
p <- c(15)
covariate_sigma <- c(1)
alpha <- c(5)
eta <- c(-5)
a <- 3*252
b <- 10*252
shock_sd <- 1
mu_delta <- c(.05)

list_of_vars <- list(n
                     ,p
                     ,covariate_sigma
                     ,alpha
                     ,eta
                     ,a
                     ,b
                     ,shock_sd
                     ,mu_delta
)

names(list_of_vars) <- list('n'
                            ,'p'
                            ,'covariate_sigma'
                            ,'alpha'
                            ,'eta'
                            ,'a'
                            ,'b'
                            ,'shock_sd'
                            ,'mu_delta'
)

# Because of the many many combinations, we have to subset even at the time of creation
gridd <- expand.grid(list_of_vars)

sim_params <- gridd

# This last bit of code will help us see how the parameter combos vary
length_unique <- function(x) { return(length(unique(x)))}

# sim_params_check <- sim_params %>% group_by(vol_model) %>% summarise(across(everything(), length_unique),.groups = 'drop')  %>% as.data.frame()
#
# sim_params_check

grid_row_count <- nrow(sim_params)

global_seed <- 1986
rng <- RNGseq(nsim * grid_row_count, global_seed)

############ end of parameter grid construction ############


########################### Begin parallel architecture ###############################

# simulation time
system.time(

  output <- foreach(
    n = sim_params$n
    , p = sim_params$p
    ,covariate_sigma = sim_params$covariate_sigma
    ,alpha = sim_params$alpha
    ,eta = sim_params$eta
    ,a = sim_params$a
    ,b = sim_params$b
    ,shock_sd = sim_params$shock_sd
    ,mu_delta = sim_params$mu_delta

    ,outer_loop_counter = icount()

    #Now we choose how we want foreach to combine the output of each sim
    ,.combine = 'rbind'
    ,.errorhandling = "remove" #pass is another option

  ) %dopar% {

    foreach(1:nsim
            , inner_loop_counter = icount()
            , selected_seed = rng[(outer_loop_counter-1)*nsim + 1:nsim]
            #Now we choose how we want foreach to combine the output of each sim
            ,.combine = 'rbind'
            ,
            .errorhandling = "remove" #pass is another option
    ) %do% { #begin inner loop

      setRNG(selected_seed)

      to_return <- exp_break_maker(n 
                                 ,p 
                                 ,covariate_sigma 
                                 ,alpha 
                                 ,eta 
                                 ,a 
                                 ,b 
                                 ,shock_sd 
                                 ,mu_delta 
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
save(output, file = paste("~/Desktop/PhD/synthetic_vol_forecasting/exponential_sims/simcount_",
                          nsim,
                          "_savetime_", format(Sys.time(), "%a%b%d%X%Y"),"_runtime_",running_hours,"_hr_",
                          "_grid_size",grid_row_count,
                          "_recovery_",recovery_rate,
                          "_seed_", global_seed,
                          ".Rdata",sep="") )

stopImplicitCluster()
