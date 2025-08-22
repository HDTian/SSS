\name{SSS}
\alias{SSS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
SSS method for studying the nonlinear effect in instrumental variable analysis
}
\description{
this function adopts thee 'S' methods (Stratificaiton, Scalar-on-function regression, and SuSiE) to obtain the estimated effect shape function
}
\usage{
SSS(Z,X,Y,
    Ns_used=100,
    x_baseline_used = NA,
    tpoints_used = 'quantile',
    boundary_used = NULL,
    precision =100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{a vector corresponding to the instrument}
   \item{X}{a vector corresponding to the exposure}
    \item{Y}{a vector corresponding to the outcome}
    \item{Ns_used}{a integer indicating the number of strata}
    \item{x_baseline_used}{a value indicating the baseline exposure value used for defining the effect shape function \eqn{h(x)}. Defult is to take the exposure mean as the baseline level}
    \item{tpoints_used}{a character indicates the type of exposure candidate points. Default is quantile, means the candidate points will be the sample quantile points. Alternative is the 'uniform', corresponding to the uniform candidate points between min(exposure) and max(exposure)}
    \item{boundary_used}{a length-2 vector indicates the left and right quantiles (when tpoints_used='quantile'), or the left and right exposure values (when tpoints_used='uniform')}
    \item{precision}{a integer indicating the degree of precision in the exposure quntiles, which is equal to the number of internal uniform grid. Default is 100.}
}
\details{
the inputted vector \code{Z}, \code{X}, \code{Y} should be strictly matched. No missing or NA value is allowed.
}
\value{
It returns a list, containing

\item{Covariate}{the covariate matrix of the sparse regression, which is ready to be fitted by SuSiE}
\item{L_used}{the number of non-zero change-points detected}
\item{PIPres}{the posterior inclusion probability matrix; each row corresponds to one detected change-point}
\item{posterior_mode}{the vector of the posterior predcited value for each detected change-point; using posterior mode}
\item{posterior_mean}{the vector of the posterior predcited value for each detected change-point; using posterior mean}
\item{PIP_ggdata}{the data frame used for ggplot of the posterior inclusion probabilities}
\item{PIP_plot}{posterior inclusion probabilities plot for each detected change-point}
\item{hx_ggdata}{the data frame used for ggplot of the estimated effect shape function}
\item{hx}{plot of the estimated effet shape function}

}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
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
Y <- 1.0*(X-1)*(X>1)  + U + Ey  # one change-point with effect size 1.0

SSS_res<- SSS(Z,X,Y,x_baseline_used = 0)

#SSS_res$PIP_plot
#SSS_res$hx
}

