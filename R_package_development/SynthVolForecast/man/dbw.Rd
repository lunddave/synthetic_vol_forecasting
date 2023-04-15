\name{dwb}
\alias{dwb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dwb(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as

dbw <- function(X,
                Tstar,
                scale = FALSE,
                sum_to_1 = 1,
                bounded_below_by = 0,
                bounded_above_by = 1,
                normchoice = c('l1', 'l2')[2],
                penalty_normchoice = c('l1', 'l2')[1],
                penalty_lambda = 0
                ) { # https://github.com/DEck13/synthetic_prediction/blob/master/prevalence_testing/numerical_studies/COP.R
  # X is a list of covariates for the time series
  # X[[1]] should be the covariate of the time series to predict
  # X[[p]] for p = 2,...,n+1 are covariates for donors

  # T^* is a vector of shock-effects time points
  # shock effect point must be > 2

  # number of time series for pool
  n <- length(X) - 1

  # COVARIATE FOR TIME SERIES UNDER STUDY AT TSTAR
  X1 <- X[[1]][Tstar[1], , drop = FALSE] # we get only 1 row

  # LOOP for grab TSTAR covariate vector for each donor
  X0 <- c()
  for (i in 1:n) {
    X0[[i]] <- X[[i + 1]][Tstar[i + 1], , drop = FALSE] #get 1 row from each donor
  }

        if (scale == TRUE) { #begin if statement
          dat <- rbind(X1, do.call('rbind', X0)) # do.call is for cluster computing?
          dat <- apply(dat, 2, function(x) scale(x, center = TRUE, scale = TRUE))
          X1 <- dat[1, , drop = FALSE]
          X0 <- c()
          for (i in 1:n) {
            X0[[i]] <- dat[i + 1, , drop = FALSE] #we are repopulating X0[[i]] with scaled+centered data
          } #end loop
        } #end if statement

  # objective function
  weightedX0 <- function(W) {
    # W is a vector of weight of the same length of X0
    n <- length(W)
    p <- ncol(X1)
    XW <- matrix(0, nrow = 1, ncol = p)
    for (i in 1:n) {
      XW <- XW + W[i] * X0[[i]]
    } #end of loop

    if (normchoice == 'l1') {
      norm <- as.numeric(norm(matrix(X1 - XW), type = "1"))
    }
    else {
      norm <- as.numeric(crossprod(matrix(X1 - XW)))
    }

    #now add penalty
    if (penalty_normchoice == 'l1' & penalty_lambda > 0) {
      norm <- norm + penalty_lambda * norm(as.matrix(W), type = "1")
    }
    else if (penalty_normchoice == 'l2' & penalty_lambda > 0) {
      norm <- norm + penalty_lambda * as.numeric(crossprod(matrix(W)))
    }
    else {norm <- norm}

    return(norm)
  } #end objective function

  # optimization and return statement

  # I have added features
  # 1) The option to remove the sum-to-1 constraint
  # 2) The option to change the lower bound to -1 or NA
  # 3) option to change the upper bound to NA'
  # 4) option to choose l1 or l2 norm as distance function

  #Thus I need if statements to implement these...

  # conditional for sum to 1
  if (is.na(sum_to_1) == FALSE) {eq_constraint <- function(W) sum(W) - 1 }
  else{eq_constraint = NULL}

  # conditional for bounding below
  if (is.na(bounded_below_by) == FALSE)
          {
            lower_bound = rep(bounded_below_by, n)
          }
          else if (is.na(bounded_below_by) == TRUE)  {
            lower_bound = NULL
          }

  #conditional for bounding above
  if (is.na(bounded_above_by) == FALSE)
          {
            upper_bound = rep(1, n)
          }
          else if (is.na(bounded_above_by) == TRUE)  {
            upper_bound = NULL
          }

  object_to_return <- solnp(par = rep(1/n, n),
                            fun = weightedX0,
                            eqfun = eq_constraint,
                            eqB = 0,
                            LB = lower_bound, UB = upper_bound,
                            control = list(trace = 0
                                           , 1.0e-8
                                           , tol = 1e-9
                                           , outer.iter = 10000
                                           , inner.iter = 10000))
  return(object_to_return$pars)

} #END dbw function

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
