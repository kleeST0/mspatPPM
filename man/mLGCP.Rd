\name{mLGCP}
\alias{mLGCP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multivariate Bayesian Log-Gaussian Cox Process (LGCP) Model}
\description{
Fit a multivariate Bayesian LGCP model for spatial point pattern analysis of imaging data, incorporating latent processes to capture multivariate dependence.
}
\usage{
mLGCP(data, K, hyperParams, startValues, mcmcParams, distance_vec,
corr_model = "exp", nu = NULL, nu0 = NULL)
}
\arguments{
\item{data}{A \code{ppp} object representing spatial point pattern data, organized using the \code{SettingCompGrid} function.}

\item{K}{An integer specifying the number of latent processes to model the multivariate dependence structure.}

\item{hyperParams}{A list specifying the hyperparameters for the priors.}

\item{startValues}{A list specifying the initial values for the MCMC algorithm. These can be initialized using the \code{initiate_startValues} function}

\item{mcmcParams}{A list specifying MCMC parameters}

\item{distance_vec}{A vector of spatial distances at which the inter-type CCF and proportion of variance PV will be evaluated.}

\item{corr_model}{Character string specifying the correlation model. Current options are \code{"exp"} for the exponential correlation function and \code{"matern"} for the \enc{Matérn}{Matern} correlation function.}
\item{nu}{A numeric vector of length \code{J} giving the fixed \enc{Matérn}{Matern} smoothness parameters for the taxon-specific latent processes. Ignored when \code{corr_model = "exp"}.
}
\item{nu0}{A numeric vector of length \code{K} giving the fixed \enc{Matérn}{Matern} smoothness parameters for the shared latent processes. Ignored when \code{corr_model = "exp"}. Required when \code{corr_model = "matern"} and \code{K > 0}.
}
}
\details{
This function implements a multivariate Bayesian LGCP model for spatial point pattern analysis. The model captures multivariate dependence structures using latent Gaussian processes and estimates posterior distributions of quantities of interest via MCMC sampling. The function allows for flexible specification of priors, covariance functions, and computational parameters. When \code{corr_model = "matern"}, the smoothness parameters \code{nu} and \code{nu0} are treated as fixed and the range parameters are estimated within the MCMC algorithm.

The output includes estimates of inter-type CCF and PV, and posterior distributions of model parameters.
}
\value{
The output includes estimates of inter-type CCF and PV, and posterior distributions of model parameters.
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

%%\seealso{
\examples{
library(spatstat.geom)

set.seed(1)

x <- runif(300)
y <- runif(300)
marks <- factor(sample(c("A", "B", "C"), 300, replace = TRUE))

multi_ppp <- ppp(x, y, window = owin(c(0, 1), c(0, 1)), marks = marks)
comp_grid <- SettingCompGrid(ppp = multi_ppp, MN = c(32, 32), ext = 2)

J <- 3
K <- 2
startValues <- initiate_startValues(ppp = comp_grid, K = K)

p.GAM <- 0.05
p.PHI <- 0.18
p.M <- 0.18
p.ALP <- 0.18
p.PHI0 <- 0.18
p.GAM0 <- 0.05
p.SIG <- 1 - p.GAM - p.GAM0 - p.PHI - p.M - p.ALP - p.PHI0

mcmcParams <- list(
  run = list(numReps = 20, thin = 1, burninPerc = 0, nSave = 10),
  pr.move = list(p.GAM = p.GAM, p.GAM0 = p.GAM0, p.PHI = p.PHI,
                 p.PHI0 = p.PHI0, p.M = p.M, p.SIG = p.SIG, p.ALP = p.ALP),
  L.gam = 5, L.gam0 = 5,
  eps.gam = rep(0.005, J), eps.gam0 = rep(0.005, K),
  M.gam.var = 1, M.gam0.var = 1,
  propVar.phi = 0.001, propVar.phi0 = 0.001
)

hyperParams <- list(
  mu.m = 0, sigSq.m = 1e5, sigSq.alp = 1,
  a.sig = 10, b.sig = 10,
  l.phi = log(0.0001), u.phi = log(0.1),
  l.phi0 = log(0.0001), u.phi0 = log(0.1),
  a.lam = 0.7, b.lam = 0.7
)

distance_vec <- seq(0, 1, by = 0.05)

\donttest{
fit_exp <- mLGCP(
  data = comp_grid, K = K,
  hyperParams = hyperParams,
  startValues = startValues,
  mcmcParams = mcmcParams,
  distance_vec = distance_vec,
  corr_model = "exp"
)

fit_matern <- mLGCP(
  data = comp_grid, K = K,
  hyperParams = hyperParams,
  startValues = startValues,
  mcmcParams = mcmcParams,
  distance_vec = distance_vec,
  corr_model = "matern",
  nu = rep(1.5, J),
  nu0 = rep(1.5, K)
)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{spatial}
\keyword{Bayesian}
\keyword{LGCP}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
