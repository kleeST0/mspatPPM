\name{PV}
\alias{PV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Proportion of Variance for Multivariate LGCP Models}
\description{
Calculates the proportion of variance (PV) for a given distance in a multivariate Log-Gaussian Cox Process (LGCP) model, using specified model parameters and covariance function settings.
}
\usage{
PV(distance, alpha_j, sigSq_j, param.corfunc_j, param.corfunc_0, model = "exp") 
}
\arguments{
  \item{distance}{A numeric vector specifying the spatial distances at which the proportion of variance (PV) is to be evaluated.}

  \item{alpha_j}{A numeric vector representing the variance scaling parameter for the latent Gaussian process associated with type j.}

  \item{sigSq_j}{A numeric value representing the variance parameter for the Gaussian process associated with type j.}

  \item{param.corfunc_j}{A numeric vector of parameters for the correlation function of the latent Gaussian process for type j.}

  \item{param.corfunc_0}{A numeric vector of parameters for the shared correlation function across all types.}

  \item{model}{A character string specifying the covariance function model to use. Options include \code{"exp"} (default) for the exponential correlation function. Future implementations may include \code{"matern"}.}
}
\details{
The **proportion of variance (PV)** quantifies the relative contribution of the inter-type covariance to the total variance at a given distance. It is a key summary statistic in multivariate spatial point pattern analysis using LGCP models.
}
\value{
A numeric vector of PV values corresponding to the input distances.
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


\examples{
# Example usage of the PV function
distance <- seq(0, 5, length.out = 50)
alpha_j <- 1.0
sigSq_j <- 0.5
param.corfunc_j <- c(0.2)  # Example correlation parameter for type j
param.corfunc_0 <- c(0.3)  # Example shared correlation parameter

# Compute PV using the exponential model
pv_values <- PV(distance, alpha_j, sigSq_j, param.corfunc_j, param.corfunc_0, model = "exp")

# Plot the results
plot(distance, pv_values, type = "l", main = "Proportion of Variance",
     xlab = "Distance", ylab = "PV")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{spatial}
\keyword{variance}
\keyword{Bayesian}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
