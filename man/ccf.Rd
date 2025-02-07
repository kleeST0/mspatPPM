\name{ccf}
\alias{ccf}
\alias{ccf2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cross-Correlation Function (CCF) Computation for Spatial Point Patterns}
\description{
The \code{ccf} and \code{ccf2} functions compute the CCF for multivariate spatial point pattern models under specified covariance models, such as the exponential correlation function.
}
\usage{
ccf(distance, alpha1, alpha2, param.corfunc, model = "exp")
ccf2(distance, alpha1, alpha2, sigSq1, sigSq2, param.corfunc, model = "exp")
}
\arguments{
  \item{distance}{A numeric vector specifying the spatial distances at which the cross-correlation function is to be evaluated.}

  \item{alpha1}{A numeric vector representing the scaling or loading parameter for the first point type.}

  \item{alpha2}{A numeric vector representing the scaling or loading parameter for the second point type.}

  \item{param.corfunc}{A numeric parameter (or vector of parameters) controlling the spatial correlation function (e.g., range or decay parameter).}

  \item{sigSq1}{A numeric value specifying the variance of the first process (used in \code{ccf2}).}

  \item{sigSq2}{A numeric value specifying the variance of the second process (used in \code{ccf2}).}

  \item{model}{A character string specifying the form of the correlation function. The default is \code{"exp"} for the exponential correlation model.}
}
\details{
The \code{ccf} function computes the cross-correlation function for spatial point processes by scaling the covariance using the \code{alpha1} and \code{alpha2} parameters. It assumes that the covariance structure is specified by a parametric model such as the exponential correlation model.

The \code{ccf2} function extends \code{ccf} by incorporating variances (\code{sigSq1} and \code{sigSq2}) of the respective processes in the computation of the cross-correlation function.

The supported correlation model is currently the \code{"exp"} (exponential) model, defined as:
\deqn{\rho(d) = \exp(-d / \theta),}
where \eqn{\theta} is the range parameter specified by \code{param.corfunc}, and \eqn{d} is the spatial distance.
}
\value{
Both functions return a numeric vector of the computed cross-correlation function values evaluated at the specified distances.
}
\references{
Lee KH, Coull BA, Majumder S, La Riviere P, Mark Welch J, Starr JR. (2025+) A Bayesian Multivariate Spatial Point Pattern Model: Application to Oral Microbiome FISH Image Data. \emph{under review}.
}
\author{
Kyu Ha Lee
}
%\note{
% ....
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
# Example usage of ccf
distance <- seq(0, 10, by = 0.5)
alpha1 <- 0.8
alpha2 <- 0.6
param.corfunc <- 2.0

# Compute CCF using exponential model
result_ccf <- ccf(distance, alpha1, alpha2, param.corfunc)
plot(distance, result_ccf, type = "l", main = "Cross-Correlation Function (CCF)",
     xlab = "Distance", ylab = "CCF")
     
# Example usage of ccf2
sigSq1 <- 1.5
sigSq2 <- 1.2

# Compute CCF with variances
result_ccf2 <- ccf2(distance, alpha1, alpha2, sigSq1, sigSq2, param.corfunc)
plot(distance, result_ccf2, type = "l", main = "CCF with Variances",
     xlab = "Distance", ylab = "CCF")
     
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{spatial}
\keyword{correlation}
\keyword{Bayesian}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
