\name{BMM_meta}
\alias{BMM_meta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bayesian Meta-Analysis for Multi-Level Image Data}
\description{
The \code{BMM_meta} function implements a Bayesian meta-analysis model to synthesize image-level results (e.g., cross-correlation function (CCF) or proportion of variance (PV)) across multiple subjects and levels, providing global and subject-specific estimates with uncertainty quantification.
}
\usage{
BMM_meta(z, var, cluster, measure="ccf", startValues=NULL, hyperParams=NULL, 
numReps, burnin, thin)
}
\arguments{
  \item{z}{A numeric vector containing the transformed image-level summary statistics
  (e.g., Fisher z-transformed CCF or logit-transformed PV).}
  
  \item{var}{A numeric vector of variances corresponding to the summary statistics in \code{z}. e.g. posterior variance}
  
  \item{cluster}{A factor or integer vector indicating the cluster or group membership for each
summary statistic (e.g., subject-level clustering).}
  
  \item{measure}{A character string specifying the measure to analyze; defaults to \code{"ccf"}.
  Options include \code{"ccf"} for cross-correlation function or \code{"pv"} for proportion of variance.}
  
  \item{startValues}{A list specifying initial values for the MCMC algorithm. Defaults to \code{NULL},
  in which case starting values are automatically initialized.}
  
  \item{hyperParams}{A list containing hyperparameters for the Bayesian model priors. Defaults to \code{NULL},
  using standard non-informative priors.}
  
  \item{numReps}{An integer specifying the total number of MCMC iterations.}
  
  \item{burnin}{The proportion of initial MCMC iterations to discard for burn-in.}
  
  \item{thin}{An integer specifying the thinning interval for saving MCMC samples.}
}

\details{
This function employs a Bayesian multilevel modeling framework to synthesize image-specific
quantities (such as CCF or PV) across multiple images and subjects. The model accounts for
hierarchical dependencies and provides full posterior distributions for subject-specific
and global parameters.

The function uses a Markov Chain Monte Carlo (MCMC) algorithm for posterior sampling,
allowing for robust inference and quantification of uncertainty.
}
  %\value{
  %~Describe the value returned
    %If it is a LIST, use
    %\item{comp1 }{Description of 'comp1'}
    %\item{comp2 }{Description of 'comp2'}
   %...
  %}
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
N = 5
Nimages = 20
zz <- rnorm(N*Nimages, 0, 1)
vv <- runif(N*Nimages, 0.5, 1)
clust <- rep(1:N, each =Nimages)

res <- BMM_meta(zz, vv, clust, measure="ccf", numReps=100, burnin=0.5, thin=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{Bayesian}
\keyword{meta-analysis}
\keyword{MCMC}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
