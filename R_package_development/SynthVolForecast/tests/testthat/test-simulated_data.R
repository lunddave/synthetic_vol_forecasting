test_that("Test dbw using simulated data.",
          {
            
            ### BEGIN EXAMPLE WITH SIMULATED DATA
            
            Tee <- 384 #Specify length of time series
            n <- 4 #Specify number of donors
            shock_time_vec <- rep(Tee/2, n+1)
            p <- 6
            Y <- list()
            X <- list()
            
            for (i in 1:(n+1)){
              #we populate the outcome series list
              set.seed(i)
              Y[[i]] <- rnorm(Tee)
              
              #we populate the covariate list
              set.seed(-2 * i)
              X[[i]] <- matrix(rnorm(p * Tee), ncol = p, byrow = FALSE)
            }
            
            expect_output(dbw(X
                              ,dbw_indices
                              ,shock_time_vec
                              ,scale = FALSE
                              ,center = FALSE
                              ,sum_to_1 = 1
                              ,bounded_below_by = 0
                              ,bounded_above_by = 1
                              ,princ_comp_count = min(length(shock_time_vec), ncol(X[[1]]) )
                              ,normchoice = c('l1', 'l2')[2]
                              ,penalty_normchoice = c('l1', 'l2')[1]
                              ,penalty_lambda = 0
                              ,Y = NULL
                              ,Y_lookback_indices = rep(list(c(1,2)),dbw_indices)
                              ,X_lookback_indices = rep(c(1,2),length(dbw_indices))
            )
            ) 
            
          }
)

test_that("Test SynthPrediction using simulated data.",
          {
            
            ### BEGIN EXAMPLE WITH SIMULATED DATA
            
            Tee <- 384 #Specify length of time series
            n <- 4 #Specify number of donors
            shock_time_vec <- rep(Tee/2, n+1)
            p <- 6
            Y <- list()
            X_sp <- list()
            
            for (i in 1:(n+1)){
              #we populate the outcome series list
              set.seed(i)
              Y[[i]] <- rnorm(Tee)
              
              #we populate the covariate list
              set.seed(-2 * i)
              X_sp[[i]] <- matrix(rnorm(p * Tee), ncol = p, byrow = FALSE)
            }
            
            expect_output(SynthPrediction(Y
                                          ,X_sp
                                          ,shock_time_vec
                                          ,rep(1, n+1)
                                          ,dbw_indices = NULL
                                          ,covariate_indices = NULL
                                          ,plots = TRUE
                                          ,display_ground_truth_choice = TRUE),
                          cran = FALSE,
                          error = FALSE,
                          transform = NULL,
                          variant = NULL,
                          cnd_class = FALSE
            )
            
          }
)


test_that("Test SynthVolForecast using simulated data.",
          {

            ### BEGIN EXAMPLE WITH SIMULATED DATA

            Tee <- 384 #Specify length of time series
            n <- 4 #Specify number of donors
            shock_time_vec <- rep(Tee/2, n+1)
            p <- 6
            Y <- list()
            X <- list()

            for (i in 1:(n+1)){
              #we populate the outcome series list
              set.seed(i)
              Y[[i]] <- rnorm(Tee)

              #we populate the covariate list
              set.seed(-2 * i)
              X[[i]] <- matrix(rnorm(p * Tee), ncol = p, byrow = FALSE)
            }

            expect_output(SynthVolForecast(Y
                                           ,X
                                           ,shock_time_vec
                                           ,rep(1, n+1)
                                           ,garch_order = c(1,1)
                                           ,plots = TRUE),
                          cran = FALSE,
                          error = FALSE,
                          transform = NULL,
                          variant = NULL,
                          cnd_class = FALSE
            )

          }
          )

test_that("Test SynthPrediction using simulated data.",
          {

            ### BEGIN EXAMPLE WITH SIMULATED DATA

            Tee <- 384 #Specify length of time series
            n <- 4 #Specify number of donors
            shock_time_vec <- rep(Tee/2, n+1)
            p <- 6
            Y <- list()
            X_sp <- list()

            for (i in 1:(n+1)){
              #we populate the outcome series list
              set.seed(i)
              Y[[i]] <- rnorm(Tee)

              #we populate the covariate list
              set.seed(-2 * i)
              X_sp[[i]] <- matrix(rnorm(p * Tee), ncol = p, byrow = FALSE)
            }

            expect_output(SynthPrediction(Y
                                          ,X_sp
                                          ,shock_time_vec
                                          ,rep(1, n+1)
                                          ,dbw_indices = NULL
                                          ,covariate_indices = NULL
                                          ,plots = TRUE
                                          ,display_ground_truth_choice = TRUE),
                          cran = FALSE,
                          error = FALSE,
                          transform = NULL,
                          variant = NULL,
                          cnd_class = FALSE
            )

          }
)
