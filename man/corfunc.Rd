\name{corfunc}
\alias{corfunc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute Correlation Function for Specified Model}
\description{
Calculates the spatial correlation function for a given set of distances and model parameters, with support for different covariance function models. The supported correlation models are currently \code{"exp"} (exponential) and \code{"matern"} (\enc{Matérn}{Matern}).
}
\usage{
corfunc(distance, param, model = "exp", nu = NULL)
}
\arguments{
  \item{distance}{A numeric vector specifying the spatial distances at which the correlation function is to be evaluated.}

  \item{param}{A numeric vector of parameters for the chosen correlation function model. The interpretation of parameters depends on the specified model.}

  \item{model}{A character string specifying the covariance function model to use. Options include \code{"exp"} (default) for the exponential correlation function and \code{"matern"} for the \enc{Matérn}{Matern} correlation function.}
  \item{nu}{Fixed \enc{Matérn}{Matern} smoothness parameter. Required when \code{model = "matern"}.}
}

\details{
For \code{model = "exp"}, the correlation function is
\eqn{\exp(-d/\phi)}. For \code{model = "matern"}, the correlation function is
\deqn{
\frac{2^{1-\nu}}{\Gamma(\nu)}\left(\frac{d}{\phi}\right)^\nu
K_\nu\left(\frac{d}{\phi}\right).
}
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


\examples{
distance <- seq(0, 2, by = 0.05)

rho_exp <- corfunc(distance, param = 0.5, model = "exp")
rho_mat <- corfunc(distance, param = 0.5, model = "matern", nu = 1.5)

plot(distance, rho_exp, type = "l", ylim = c(0, 1),
     xlab = "Distance", ylab = "Correlation")
lines(distance, rho_mat, lty = 2)
legend("topright", legend = c("Exponential", "Matern"),
       lty = c(1, 2), bty = "n")
}


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
