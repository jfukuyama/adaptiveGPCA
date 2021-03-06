% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adaptivegpca-package.R
\name{plot.adaptivegpca}
\alias{plot.adaptivegpca}
\title{Plot an adaptivegpca object}
\usage{
\method{plot}{adaptivegpca}(x, type = c("scree", "samples", "variables"),
  axes = c(1, 2), ...)
}
\arguments{
\item{x}{An object of class \code{adaptivegpca}}

\item{type}{What type of plot to make. \code{scree} will make a
scree plot showing the eigenvalues, \code{samples} will plot
the samples, and \code{variables} will plot the variables.}

\item{axes}{Which axes to plot.}

\item{...}{Not used.}
}
\description{
Plots the output from \code{\link{adaptivegpca}}, either a scree
plot, the samples, or the variables.
}
\examples{
data(AntibioticSmall)
out.agpca = adaptivegpca(AntibioticSmall$X, AntibioticSmall$Q, k = 2)
plot(out.agpca)
plot(out.agpca, type = "samples")
plot(out.agpca, type = "variables")
}
