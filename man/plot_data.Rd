\name{plot_data}
\alias{plot_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot Spatial Point Pattern Data}
\description{
Visualize spatial point pattern data on the computational grid
}
\usage{
plot_data(ppp, type)
}
\arguments{
    \item{ppp}{A \code{ppp} object representing spatial point pattern data, organized using the \code{SettingCompGrid} function.}

  \item{type}{A character value specifying the type (mark) of points to plot.}
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
# Load the spatstat library
library(spatstat.geom)

# Create multi-type point pattern data
x <- runif(500)  # x-coordinates
y <- runif(500)  # y-coordinates
marks <- factor(sample(c("A", "B"), 500, replace = TRUE))  # Marks: "A" or "B"

# Combine into a ppp object
multi_type_ppp <- ppp(x, y, window = owin(c(0, 1), c(0, 1)), marks = marks)

# Set computational grid with 64x64 grid cells and extension
ppp_grid <- SettingCompGrid(ppp = multi_type_ppp, MN = c(64, 64), ext = 2)

plot_data(ppp_grid, type = "A")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{spatial}
\keyword{plot}
\keyword{visualization}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
