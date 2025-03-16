\name{my_int_function}
\alias{my_int_function}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
integration calculator
}
\description{
calculate the integration using numerical approximation
}
\usage{
my_int_function(xvec, yvec1, yvec2=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xvec}{a vector corresponding the x-axis locations}
   \item{yvec1}{a vector corresponding the y-axis values over \code{xvec}}
   \item{yvec2}{another optional vector corresponding the y-axis values over \code{xvec}}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
it returns a list contains
\item{sum}{the integral}
\item{areas}{the elemnts over the \code{xvec}}
}

\examples{
# calculate pi
xvec<-seq( 0,1, length=10000   )
yvec <- sqrt(1- (1-xvec)^2)
int_res<- my_int_function(xvec, yvec)$sum

#estimated pi value
int_res*4  #3.14159

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
