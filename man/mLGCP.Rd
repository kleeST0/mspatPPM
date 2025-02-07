\name{mLGCP}
\alias{mLGCP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multivariate Bayesian Log-Gaussian Cox Process (LGCP) Model}
\description{
Fit a multivariate Bayesian LGCP model for spatial point pattern analysis of imaging data, incorporating latent processes to capture multivariate dependence.
}
\usage{
mLGCP(data, K, hyperParams, startValues, mcmcParams, distance_vec, model="exp")
}
\arguments{
\item{data}{A \code{ppp} object representing spatial point pattern data, organized using the \code{SettingCompGrid} function.}
  
  \item{K}{An integer specifying the number of latent processes to model the multivariate dependence structure.}

  \item{hyperParams}{A list specifying the hyperparameters for the priors.}

  \item{startValues}{A list specifying the initial values for the MCMC algorithm. These can be initialized using the \code{initiate_startValues} function}

  \item{mcmcParams}{A list specifying MCMC parameters}

  \item{distance_vec}{A vector of spatial distances at which the inter-type CCF and proportion of variance PV will be evaluated.}

  \item{model}{A character string specifying the covariance function model. The default is \code{"exp"} for an exponential covariance function.}
}
\details{
This function implements a multivariate Bayesian LGCP model for spatial point pattern analysis. The model captures multivariate dependence structures using latent Gaussian processes and estimates posterior distributions of quantities of interest via MCMC sampling. The function allows for flexible specification of priors, covariance functions, and computational parameters.

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
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
library(spatstat.geom)

# Create multi-type point pattern data
x <- runif(500)  # x-coordinates
y <- runif(500)  # y-coordinates
marks <- factor(sample(c("A", "B", "C"), 500, replace = TRUE))  # Marks

# Create a ppp object
multi_ppp <- ppp(x, y, window = owin(c(0, 1), c(0, 1)), marks = marks)

# Set computational grid
comp_grid <- SettingCompGrid(ppp = multi_ppp, MN = c(32, 32), ext = 2)

# Specify model parameters
J <- 3
K <- 2  # Number of latent processes
startValues <- initiate_startValues(ppp = comp_grid, K = K)

# Set MCMC parameters
p.GAM <- 0.05
p.PHI <- 0.18
p.M <- 0.18
if(K >0)
{
    p.ALP <- 0.18
    p.PHI0 <- 0.18
    p.GAM0 <- 0.05
}else
{
    p.ALP <- 0
    p.PHI0 <- 0
    p.GAM0 <- 0
}
p.SIG <- 1-p.GAM-p.GAM0-p.PHI-p.M-p.ALP-p.PHI0

mcmcParams <- list(
  run = list(numReps = 50, thin = 1, burninPerc = 0, nSave = 10),
  pr.move=list(p.GAM=p.GAM, p.GAM0=p.GAM0, p.PHI=p.PHI,
    p.PHI0=p.PHI0, p.M=p.M, p.SIG=p.SIG, p.ALP=p.ALP),
  L.gam=20, L.gam0=20, eps.gam=rep(0.005, J), eps.gam0=rep(0.005, K),
  M.gam.var=1, M.gam0.var=1, propVar.phi=0.001, propVar.phi0=0.001
)

# Hyperparameters
hyperParams <- list(
  mu.m = 0, sigSq.m = 1e5, sigSq.alp = 1, a.sig = 10, b.sig = 10,
  l.phi = log(0.0001), u.phi = log(0.01), l.phi0 = log(0.0001), u.phi0 = log(0.01),
    a.lam = 0.7, b.lam = 0.7
)

# Define distance vector
distance_vec <- seq(0, 1, by = 0.05)

# Fit the Multivariate LGCP model
fit <- mLGCP(data = comp_grid, K = K, hyperParams = hyperParams,
             startValues = startValues, mcmcParams = mcmcParams, distance_vec = distance_vec)

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
