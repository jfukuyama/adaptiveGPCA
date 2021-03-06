% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adaptivegpca-package.R
\name{gpcaEvecs}
\alias{gpcaEvecs}
\title{gPCA using pre-computed eidendecomposition}
\usage{
gpcaEvecs(X, evecs, evals, D = rep(1, nrow(X)), k)
}
\arguments{
\item{X}{Data matrix.}

\item{evecs}{Eigenvectors of \code{Q}, the inner product/similarity
matrix.}

\item{evals}{Eigenvalues of \code{Q}.}

\item{D}{Sample weights}

\item{k}{The number of components to return.}
}
\description{
Performs gPCA with pre-computed eigenvectors and eigenvalues.
}
\keyword{internal}
