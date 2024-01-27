\name{dbw}
\alias{dbw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function that carries out distance-based weighting.
}
\description{
A function that carries out distance-based weighting.}

\usage{
dbw(X, dbw_indices, shock_time_vec, scale = FALSE, center = FALSE, sum_to_1 = 1, bounded_below_by = 0, bounded_above_by = 1, normchoice = c("l1", "l2")[2], penalty_normchoice = c("l1", "l2")[1], penalty_lambda = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
nxp matrix of covariates
}
  \item{dbw_indices}{
Indices of X to be used in distance-based-weighting
}
  \item{shock_time_vec}{
Vector of shock times, with first entry corresponding to time series under study.
}
  \item{scale}{
TRUE/FALSE, whether or not to scale each column to sd = 1}
  \item{center}{
TRUE/FALSE, whether or not to de-mean each column}
  \item{sum_to_1}{
TRUE/FALSE, whether optimization should arrive at weights that sum to 1}
  \item{bounded_below_by}{
TRUE/FALSE, whether weights should no less than 0
}
  \item{bounded_above_by}{
TRUE/FALSE, whether weights should be no greater than 1
}
  \item{normchoice}{
L1 or L2 norm used in objection function}
  \item{penalty_normchoice}{
L1 or L2 norm used in optimization penalty}
}
  \item{penalty_lambda}{
Penalty used in optmization}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
vector of length n
}
\references{
Lin and Eck 2021
}
\author{
David Lundquist}
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
