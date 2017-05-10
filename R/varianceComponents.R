#' Estimate parameters in hierarchical model
#'
#' Estimates the values of \eqn{r} and \eqn{\sigma} in a model \eqn{X \sim N(0, \sigma^2
#' (r Q + (1 - r)I))}.
#'
#' @param X An \eqn{n \times p} data matrix. 
#' @param Q A \eqn{p \times p} matrix giving the prior variance on the
#' rows of \code{X}.
#' @param Qeig If the eigendecomposition of \code{Q} is already computed, it
#' can be included here.
#' 
#' @return A list with \eqn{r} and \eqn{\sigma}. 
#' @examples
#' data(AntibioticSmall)
#' estimateComponents(AntibioticSmall$X, AntibioticSmall$Q)
#' @export
estimateComponents <- function(X, Q, Qeig = NULL) {
    if(is.null(Qeig)) {
        Qeig = eigen(Q, symmetric = TRUE)
    }
    Xtilde = X %*% Qeig$vectors
    f = function(r) -likelihoodR(Xtilde, r, Qeig$values)
    r = stats::optimize(f, c(0,1))$minimum
    sigma2 = sigma2OfR(Xtilde, r, Qeig$values)
    return(list(r = r, sigma = sqrt(sigma2)))
}


#' Derivative of the likelihood
#'
#' Derivative of the likelihood in the hierarchical model as a
#' function of \eqn{r}.
#'
#' @param Xtilde The transformed data
#' @param r r
#' @param D The eigenvalues of Q
#' @keywords internal
gradLik <- function(Xtilde, r, D) {
    n = nrow(Xtilde)
    p = ncol(Xtilde)
    a = -.5 * n * sum((D - 1) / (r * D + 1 - r))
    b = -.5 * n * p * gradSigma2OfR(Xtilde, r, D) / sigma2OfR(Xtilde, r, D)
    return(a + b)
}

#' Derivative of \eqn{\sigma^2}
#'
#' Derivative of \eqn{\sigma^2(r)} in the hierarchical model
#'
#' @inheritParams gradLik
#' @keywords internal
gradSigma2OfR <- function(Xtilde, r, D) {
    p = length(D)
    n = nrow(Xtilde)
    g = -sum(apply(Xtilde, 1, function(x)
        sum(x^2 * (D - 1) / (n * p * (r * D + 1 - r)^2))))
    return(g)
}

#' Value of \eqn{\sigma^2} that maximizes the likelihood for a given value of
#' \eqn{r}
#'
#' @inheritParams gradLik
#' @keywords internal
sigma2OfR <- function(Xtilde, r, D) {
    p = length(D)
    n = nrow(Xtilde)
    sigma2 = sum(apply(Xtilde, 1, function(x) sum(x^2 / (p * n * (r * D + 1 - r)))))
    return(sigma2)
}

#' The likelihood at a given value of \eqn{r} and the maximizing \eqn{\sigma} for
#' that value of \eqn{r}.
#'
#' @inheritParams gradLik
#' @keywords internal
likelihoodR <- function(Xtilde, r, D) {
    sigma = sqrt(sigma2OfR(Xtilde, r, D))
    return(likelihood(Xtilde, sigma, r, D))
}

#' The likelihood at a given value of \eqn{r} and \eqn{\sigma}
#'
#' @inheritParams gradLik
#' @param sigma Overall scaling factor.
#' @keywords internal
likelihood <- function(Xtilde, sigma, r, D) {
    p = ncol(Xtilde)
    n = nrow(Xtilde)
    s1 = n * sum(log(sigma^2 * (r * D + 1 - r)))
    s2 = sum(apply(Xtilde, 1, function(x) sum(x^2 / (sigma^2 * (r * D + 1 - r)))))
    l = -.5 * (s1 + s2) - (p / 2) * log(2 * pi)
    return(l)
}


#' Variance along eigenvectors of Q
#' 
#' Project the sample points stored in the rows of \code{X} along the
#' eigenvectors of \code{Q} and find the variance along each of the
#' projections.
#' 
#' @param X An \eqn{n \times p} data matrix, each row corresponding to a sample.
#' @param Q A \eqn{p \times p} similarity matrix, either as a matrix
#' or as its eigendecomposition (the output from \code{eigen}).
#'
#' @return A vector containing the variance of the samples along each
#' of the eigenvectors of \code{Q}.
#' @examples
#' data(AntibioticSmall)
#' voe = varianceOnEvecs(AntibioticSmall$X, AntibioticSmall$Q)
#' @export
varianceOnEvecs <- function(X, Q) {
    if(is.list(Q) && !is.null(Q$vectors) && !is.null(Q$values)) {
        Qeig = Q
    } else {
        Qeig = eigen(Q, symmetric = TRUE)
    }
    Xtilde = X %*% Qeig$vectors
    vars = apply(Xtilde, 2, stats::var)
    return(vars)
}

#' Likelihood of data in two-parameter model
#'
#' @param r1 Coefficient of Q
#' @param r2 Coefficient of Q^(-1) in the noise part.
#' @param Xtilde The data projected onto the eigenvectors of Q.
#' @param D The eigenvalues of Q
#'
#' @return The marginal likelihood of the data given r1 and r2.
#' @keywords internal
likelihood_two_params <- function(r1, r2, Xtilde, D) {
    n = nrow(Xtilde)
    p = ncol(Xtilde)
    f = r1 * D + (1 - r1) * r2 * D^(-1) + (1 - r1) * (1 - r2)
    xtnorm = apply(Xtilde, 1, function(x) sum(x^2 / f))
    l = -n * .5 * sum(log(f)) - sum(.5 * p * log(xtnorm))
    return(l)
}


#' Estimate variance components
#'
#' Estimate variance components in a two-parameter model where \eqn{X \sim
#' N(0, \sigma^2 (r_1 Q + (1 - r_1) (r_2 Q^(-1) + (1 - r_2) I)))}
#'
#' @param X An n x p data matrix. 
#' @param Q A p x p psd matrix giving the similarity between the
#' variables.
#' 
#' @return A vector with r1 and r2
#' @keywords internal
estimateComponents2 <- function(X, Q) {
    if(is.list(Q) & !is.null(Q$vectors) & !is.null(Q$values)) {
        Qeig = Q
    } else {
        Qeig = eigen(Q, symmetric = TRUE)
    }
    Qeig$values = ncol(X) * Qeig$values / sum(Qeig$values)
    Xtilde = X %*% Qeig$vectors
    f = function(r) -likelihood_two_params(r[1], r[2], Xtilde, Qeig$values)
    r = stats::optim(c(.5, .5), f, lower = c(0,0), upper = c(1,1))$par
    return(r)
}
