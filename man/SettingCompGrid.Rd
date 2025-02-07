\name{SettingCompGrid}
\alias{SettingCompGrid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Set Computational Grid for Bayesian Inference in Multivariate LGCP Models}
\description{
Organize and construct a computational grid for spatial point pattern data to facilitate Bayesian inference in multivariate LGCP models.
}
\usage{
SettingCompGrid(ppp, MN, ext, object.only = F, exclude=F)
}
\arguments{
  \item{ppp}{A \code{ppp} object representing spatial point pattern data, which serves as the input for grid construction.}

  \item{MN}{A numeric vector of length 2 specifying the number of grid cells in the x and y directions (e.g., \code{c(M, N)}), which determines the resolution of the computational grid.}

  \item{ext}{A numeric value indicating the amount of extension to be added to the observation window.}

  \item{object.only}{If \code{TRUE}, the analysis will focus on the object specified within its esimtated boundary; otherwise, the entire window will be considered for the intensity function. Defaults to \code{FALSE}.}

  \item{exclude}{An optional argument. If \code{TRUE}, the last mark (type) in the \code{ppp} object will be treated as coordinates to be excluded from the computational grid. Defaults to \code{FALSE}.
}
}
\details{
The \code{SettingCompGrid} function creates a computational grid by overlaying the observation window of the provided \code{ppp} object (spatial point pattern) with a grid of specified resolution.

This grid serves as a foundational component for fitting multivariate LGCP models, allowing the spatial point pattern data to be structured appropriately for computational efficiency.
}
%%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 %% \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%%}
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
# Load the spatstat library
library(spatstat.geom)

# Create multi-type point pattern data
x <- runif(10)  # x-coordinates
y <- runif(10)  # y-coordinates
marks <- factor(sample(c("A", "B"), 10, replace = TRUE))  # Marks: "A" or "B"

# Combine into a ppp object
multi_type_ppp <- ppp(x, y, window = owin(c(0, 1), c(0, 1)), marks = marks)

# Set computational grid with 64x64 grid cells and extension
grid <- SettingCompGrid(ppp = multi_type_ppp, MN = c(64, 64), ext = 2)
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
