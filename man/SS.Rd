\name{SS}
\alias{SS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
parametric fitting using specified change-points
}
\description{
this function adopts two 'S' methods (Stratificaiton, Scalar-on-function regression) and parametric fitting with specified change-points to obtain the estimated effect shape function
}
\usage{
SS(Z,X,Y,pos,Ns_used=100,x_baseline_used = NA,precision =100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{a vector corresponding to the instrument}
   \item{X}{a vector corresponding to the exposure}
    \item{Y}{a vector corresponding to the outcome}
    \item{pos}{a vector indicating the user-specified change-point values}
    \item{Ns_used}{a integer indicating the number of strata}
    \item{x_baseline_used}{a value indicating the baseline exposure value used for defining the effect shape function \eqn{h(x)}. Defult is to take the exposure mean as the baseline level}
    \item{precision}{a integer indicating the degree of precision in the exposure quntiles, which is equal to the number of internal uniform grid. Default is 100.}
}
}
\details{
the inputted vector \code{Z}, \code{X}, \code{Y} should be strictly matched. No missing or NA value is allowed.
}
\value{
It returns a list, containing

\item{bhat}{the estimated effect}
\item{Var_bhat}{the covariance matrix of \code{bhat}}
\item{hx_ggdata}{the data frame used for ggplot of the estimated effect shape function}
\item{hx}{plot of the estimated effet shape function}
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
N<-50000
set.seed(100)
Z <- rbinom( N , 1 , 0.5 )
U <- rnorm( N , 0 , 1 )
Ex <- rnorm( N , 0 , 1 )
alphaZ<- 0.15
X <-  alphaZ*Z + U + Ex    # weak instrument
summary( lm(  X  ~ Z )  )$r.squared   #  < 0.01
Ey <- rnorm( N , 0 , 1 )
Y <- 2*(X+1)*(X>-1) -2*(X-1)*(X>1)   + U + Ey  # one change-point with effect size 1.0
SSS_res<- SSS(Z,X,Y,x_baseline_used = 0)

SS_res <- SS(Z,X,Y,x_baseline_used = 0, pos = SSS_res$posterior_mean ) # use the predicted change-point values as the specified change-point in SS
SS_res$hx
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
