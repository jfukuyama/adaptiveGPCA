#' Antibiotic time course experiment.
#'
#' A phyloseq object describing a time course experiment in which
#' three people two courses of cipro and had their gut microbiomes
#' sampled. See Dethlefsen and Relman, PNAS (2010), at
#' https://www.ncbi.nlm.nih.gov/pubmed/20847294 for more details.
#'
#' @format A phyloseq object.
#' @name AntibioticPhyloseq
NULL

#' A subset of the antibiotic data
#'
#' This is a smaller version of the \code{AntibioticPhyloseq} dataset,
#' for use in the examples so that the running time isn't so long. It
#' has the same samples and a randomly selected set of 200 of the
#' taxa. It is stored as a list with three components: the normalized
#' OTU abundances (\code{X}), the similarity matrix for the taxa
#' (\code{Q}), and the diagonal weight matrix (\code{D}, the identity
#' matrix).
#'
#' @format A list with three components.
#' @name AntibioticSmall
NULL
