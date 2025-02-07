\name{post_PV}
\alias{post_PV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Posterior Proportion of Variance (PV) for Multivariate LGCP Models}
\description{
Compute the posterior distribution of the proportion of variance (PV) for a multivariate Log-Gaussian Cox Process (LGCP) model at specified spatial distances.
}
\usage{
post_PV(fit, distance_vec, p.ind = NULL)
}
\arguments{
  \item{fit}{An object containing the fitted multivariate LGCP model, typically obtained using the \code{mLGCP} function.}

  \item{distance_vec}{A numeric vector specifying the spatial distances at which the proportion of variance (PV) is to be evaluated.}

  \item{p.ind}{An optional vector specifying the indices of the posterior samples to be used. If \code{NULL} (default), the CCF is computed for all available posterior simulations.}
}
\details{
The function \code{post_PV} estimates the posterior distribution of the **proportion of variance (PV)** for multivariate LGCP models at specified distances. PV quantifies the relative contribution of inter-type variance to the total variance, providing insights into the spatial dependence structure of point patterns.

The function uses the posterior samples of the fitted model parameters to compute PV values and summarizes them to provide posterior means and credible intervals.

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
\keyword{variance}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
