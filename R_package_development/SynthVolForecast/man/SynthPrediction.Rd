\name{SynthPrediction}
\alias{SynthPrediction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function that uses synthetic methods to predict observable levels.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SynthPrediction(Y_series_list, covariates_series_list, shock_time_vec, shock_length_vec, dbw_scale = TRUE, dbw_center = TRUE, dbw_indices = NULL, covariate_indices = NULL, geometric_sets = NULL, days_before_shocktime_vec = NULL, arima_order = NULL, user_ic_choice = c("aicc", "aic", "bic")[1], plots = TRUE, display_ground_truth_choice = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y_series_list}{
%%     ~~Describe \code{Y_series_list} here~~
}
  \item{covariates_series_list}{
%%     ~~Describe \code{covariates_series_list} here~~
}
  \item{shock_time_vec}{
%%     ~~Describe \code{shock_time_vec} here~~
}
  \item{shock_length_vec}{
%%     ~~Describe \code{shock_length_vec} here~~
}
  \item{dbw_scale}{
%%     ~~Describe \code{dbw_scale} here~~
}
  \item{dbw_center}{
%%     ~~Describe \code{dbw_center} here~~
}
  \item{dbw_indices}{
%%     ~~Describe \code{dbw_indices} here~~
}
  \item{covariate_indices}{
%%     ~~Describe \code{covariate_indices} here~~
}
  \item{geometric_sets}{
%%     ~~Describe \code{geometric_sets} here~~
}
  \item{days_before_shocktime_vec}{
%%     ~~Describe \code{days_before_shocktime_vec} here~~
}
  \item{arima_order}{
%%     ~~Describe \code{arima_order} here~~
}
  \item{user_ic_choice}{
%%     ~~Describe \code{user_ic_choice} here~~
}
  \item{plots}{
%%     ~~Describe \code{plots} here~~
}
  \item{display_ground_truth_choice}{
%%     ~~Describe \code{display_ground_truth_choice} here~~
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
David Lundquist
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
