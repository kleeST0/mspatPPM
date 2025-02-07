\name{post_ccf}
\alias{post_ccf}
\alias{post_ccf2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Posterior Cross-Correlation Functions (CCF) for Multivariate LGCP Models}
\description{
Compute the posterior distributions of the cross-correlation function (CCF) for multivariate Log-Gaussian Cox Process (LGCP) models. Two versions of the function (\code{post_ccf} and \code{post_ccf2}) provide alternative methods for estimating the CCF.
}
\usage{
post_ccf(fit, distance_vec, p.ind=NULL)
post_ccf2(fit, distance_vec, p.ind=NULL)
}
\arguments{
  \item{fit}{An object containing the fitted multivariate LGCP model, typically obtained using the \code{mLGCP} function.}

  \item{distance_vec}{A numeric vector of spatial distances at which the CCF is to be evaluated.}

  \item{p.ind}{An optional vector specifying the indices of the posterior samples to be used. If \code{NULL} (default), the CCF is computed for all available posterior simulations.}
}
\details{
The functions \code{post_ccf} and \code{post_ccf2} estimate the posterior distributions of the cross-correlation function (CCF) for spatial point patterns modeled using multivariate LGCPs. The CCF quantifies the interdependence between different types (or marks) of points as a function of spatial distance.

- \code{post_ccf}: Computes the posterior CCF using a direct estimation method based on posterior samples from the fitted model.
- \code{post_ccf2}: Provides an alternative formulation for posterior CCF estimation.

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
\keyword{Bayesian}
\keyword{cross-correlation}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
