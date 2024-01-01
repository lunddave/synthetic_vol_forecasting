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
