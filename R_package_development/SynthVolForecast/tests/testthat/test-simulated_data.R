options(digits = 7)

### BEGIN EXAMPLE WITH SIMULATED DATA

Tee <- 384 #Specify length of time series
n <- 4 #Specify number of donors
shock_time_vec <- rep(Tee/2, n+1)
p <- 6

Y <- list()

#we populate the covariate list
for (i in 1:(n+1)){
  set.seed(i)
  Y[[i]] <- rnorm(Tee)
}

X <- list()

#we populate the covariate list
for (i in 1:(n+1)){
  set.seed(-2 * i)
  X[[i]] <- matrix(rnorm(p * Tee), ncol = p, byrow = FALSE)
}


# #save to rds file
# save(X, file = "X.RData")
# save(Y, file = 'Y.RData')

test_that("Test SynthVolForecast using simulated data.",
          {
            # X <- readRDS(test_path("testdata", "helper_X.rds"))
            # Y <- readRDS(test_path("testdata", "helper_Y.rds"))
            # n <- length(X)
            shock_time_vec <- rep(length(Y[[1]])/2, n+1)

            expect_output(

              SynthVolForecast(Y
                               ,X
                               ,shock_time_vec
                               ,rep(1, n+1)
                               ,garch_order = c(1,1)
                               ,plots = TRUE)
              ,
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
            # X <- readRDS(test_path("testdata", "helper_X.rds"))
            # Y <- readRDS(test_path("testdata", "helper_Y.rds"))
            # n <- length(X)
            shock_time_vec <- rep(length(Y[[1]])/2, n+1)

            expect_output(

              SynthPrediction(Y
                              ,X
                              ,shock_time_vec = shock_time_vec
                              ,rep(1, n+1)
                              ,dbw_indices = NULL
                              ,covariate_indices = NULL
                              ,plots = TRUE
                              ,display_ground_truth_choice = TRUE
              )
              ,
              cran = FALSE,
              error = FALSE,
              transform = NULL,
              variant = NULL,
              cnd_class = FALSE
            )
          }
)
