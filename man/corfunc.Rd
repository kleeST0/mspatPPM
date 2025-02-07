\name{corfunc}
\alias{corfunc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute Correlation Function for Specified Model}
\description{
Calculates the spatial correlation function for a given set of distances and model parameters, with support for different covariance function models such as the exponential (\code{"exp"}) model.
}
\usage{
corfunc(distance, param, model="exp")
}
\arguments{
  \item{distance}{A numeric vector specifying the spatial distances at which the correlation function is to be evaluated.}

  \item{param}{A numeric vector of parameters for the chosen correlation function model. The interpretation of parameters depends on the specified model.}

  \item{model}{A character string specifying the correlation function model. Options include \code{"exp"} (default) for the exponential model. }
}
\value{
A numeric vector of correlation values corresponding to the input \code{distance} vector.
}
\references{
Lee KH, Coull BA, Majumder S, La Riviere P, Mark Welch J, Starr JR. (2025+) A Bayesian Multivariate Spatial Point Pattern Model: Application to Oral Microbiome FISH Image Data. \emph{under review}.
}
\author{
Kyu Ha Lee
}
%\note{
%....
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{spatial}
\keyword{correlation}
\keyword{covariance}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
