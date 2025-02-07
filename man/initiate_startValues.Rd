\name{initiate_startValues}
\alias{initiate_startValues}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Initialize Starting Values for Multivariate LGCP Model}
\description{
Initializes starting values for parameters of the Bayesian multivariate Log-Gaussian Cox Process (LGCP) model.
}
\usage{
initiate_startValues(ppp, K, model="MLGCP", gamma=NULL, m=NULL, sigSq=NULL,
phi=NULL, gamma0=NULL, alpha=NULL, phi0=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ppp}{A \code{ppp} object representing spatial point pattern data, organized using the \code{SettingCompGrid} function.}

  \item{K}{An integer specifying the number of latent processes to initialize for the multivariate LGCP model.}

  \item{model}{A character string specifying the model type. The default is \code{"MLGCP"} for the multivariate LGCP model.}

  \item{gamma}{Optional starting values for \eqn{\gamma}. If \code{NULL}, default initial values are generated.}

  \item{m}{Optional starting values for \eqn{m}. If \code{NULL}, default initial values are generated.}

  \item{sigSq}{Optional starting values for \eqn{\sigma^2} of the latent Gaussian processes. If \code{NULL}, default initial values are generated.}

  \item{phi}{Optional starting values for  \eqn{\phi}. If \code{NULL}, default initial values are generated.}

  \item{gamma0}{Optional starting values for \eqn{\gamma_0}. If \code{NULL}, default initial values are generated.}

  \item{alpha}{Optional starting values for \eqn{\alpha}. If \code{NULL}, default initial values are generated.}

  \item{phi0}{Optional starting values for \eqn{\phi_0}. If \code{NULL}, default initial values are generated.}
}
%\details{
%...
%}
\value{
\code{initiate_startValues} returns a list containing starting values that can be used for \code{mmzipBvs}.
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
\keyword{initialization}
