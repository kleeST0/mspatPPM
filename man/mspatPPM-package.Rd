\name{mspatPPM-package}
\alias{mspatPPM-package}
\alias{mspatPPM}
\docType{package}
\title{
Multivariate Spatial Point Process Models for Image Data Analysis
}
\description{
A Bayesian framework for spatial point pattern analysis of imaging data using a flexible multivariate log-Gaussian Cox process (LGCP) model.
}
\details{
The package includes following functions:
\tabular{ll}{
\code{mLGCP} \tab Implements a Bayesian multivariate log-Gaussian Cox process (LGCP) model \cr
\code{SettingCompGrid} \tab Constructs a computational grid \cr
\code{plot_data} \tab Visualizes spatial point patterns mapped onto the computational grid \cr
\code{initiate_startValues} \tab Initializes starting values for the MCMC algorithm \cr
\code{corfunc} \tab Defines the spatial covariance structure of the latent processes \cr
\code{ccf} \tab Calculates the standardized inter-type cross-correlation function (CCF) \cr
\code{ccf2} \tab Provides an alternative formulation of the CCF \cr
\code{PV} \tab Quantifies the proportion of inter-type variance relative to the total variance \cr
\code{post_ccf} \tab Estimates the posterior distribution of the inter-type CCF. \cr
\code{post_ccf2} \tab Computes the posterior distribution of the alternative CCF formulation, \cr
\code{post_PV} \tab Generates the posterior distribution of the proportion of variance \cr
\code{BMM_meta} \tab Implements a hierarchical Bayesian meta-analysis framework \cr
}
\tabular{ll}{
Package: \tab mspatPPM \cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2025-2-15\cr
License: \tab AGPL-3\cr
LazyLoad: \tab yes\cr
}
}
\author{
Kyu Ha Lee, Brent Coull, Suman Majumder, Jacqueline Starr \cr
Maintainer: Kyu Ha Lee <klee15239@gmail.com>
}
\references{
Lee KH, Coull BA, Majumder S, La Riviere P, Mark Welch J, Starr JR. (2025+) A Bayesian Multivariate Spatial Point Pattern Model: Application to Oral Microbiome FISH Image Data. \emph{under review}.
}
%%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
%\seealso{
%Optional links to other man pages, e.g. ~~
%}
%\examples{
%}
