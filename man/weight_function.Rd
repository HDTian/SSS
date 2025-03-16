\name{weight_function}
\alias{weight_function}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
weight function
}
\description{
obtain the smoothed and re-weighted weight function values over coords (e.g. tpoints)
}
\usage{
weight_function(coords, Z , X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coords}{a vector of the exposure values over which the weight function values should be computed}
  \item{Z}{a vector corresponding to the instrument}
   \item{X}{a vector corresponding to the exposure}

}
\details{
the inputted vector \code{Z}, \code{X}, \code{Y} should be strictly matched. No missing or NA value is allowed.

the weight function is to be positive (so ensure your instrument and exposure is defined so that \eqn{cov(Z,X)>0}), and it is then smoothed according the log-transformation via kernel smoothing

the weight function is standarized so its integral is equal to 1, so ensure that your \code{coords} range is large enough to cover the weight funciton domain; it is suggested to cover the mimial and maximal values of the exposure
}
\value{
it returns a vector of the weight function values over \code{coords}
}
\examples{
set.seed(100)
Z <- rnorm(1000, 0 , 1)
X <- 0.15*Z + rnorm(1000,0,1)
w_res <- weight_function( seq(  min(X), max(X) , length=1000 )  , Z , X)
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
