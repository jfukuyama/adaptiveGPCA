

#' Estimate parameters in hierarchical model
#'
#' Estimates the values of r and sigma in a model X sim N(0, sigma^2
#' (r Q + (1 - r)I)).
#'
#' @param X A n x p matrix with n samples in p-dimensional space.
#' @param Q The prior variance on the means.
#' @param maxit The maximum number of bisections to try.
#' @param tol Stop if the absolute value of the derivative is less than this.
#' @param Qeig If the eigendecomposition of Q is already computed, it
#' can be included here.
#'
#' 
#' @return A list with r, sigma, and the values of the likelihood on a
#' grid of values between 0 and 1.
#' @export
estimateComponents <- function(X, Q, maxit = 8, tol = 10^(-10),
                               Qeig = NULL) {
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
#' function of r
#'
#' @param Xtilde The transformed data
#' @param r r
#' @param D The eigenvalues of Q
gradLik <- function(Xtilde, r, D) {
    n = nrow(Xtilde)
    p = ncol(Xtilde)
    a = -.5 * n * sum((D - 1) / (r * D + 1 - r))
    b = -.5 * n * p * gradSigma2OfR(Xtilde, r, D) / sigma2OfR(Xtilde, r, D)
    return(a + b)
}

#' Derivative of sigma^2
#'
#' Derivative of sigma^2(r) in the hierarchical model
#'
#' @inheritParams gradLik
gradSigma2OfR <- function(Xtilde, r, D) {
    p = length(D)
    n = nrow(Xtilde)
    g = -sum(apply(Xtilde, 1, function(x)
        sum(x^2 * (D - 1) / (n * p * (r * D + 1 - r)^2))))
    return(g)
}

#' Value of sigma^2 that maximizes the likelihood for a given value of
#' r
#'
#' @inheritParams gradLik
sigma2OfR <- function(Xtilde, r, D) {
    p = length(D)
    n = nrow(Xtilde)
    sigma2 = sum(apply(Xtilde, 1, function(x) sum(x^2 / (p * n * (r * D + 1 - r)))))
    return(sigma2)
}

#' The likelihood at a given value of r and the maximizing sigma for
#' that value of r
#'
#' @inheritParams gradLik
likelihoodR <- function(Xtilde, r, D) {
    sigma = sqrt(sigma2OfR(Xtilde, r, D))
    return(likelihood(Xtilde, sigma, r, D))
}

#' The likelihood at a given value of r and sigma
#'
#' @inheritParams gradLik
#' @param sigma Overall scaling factor. 
likelihood <- function(Xtilde, sigma, r, D) {
    p = ncol(Xtilde)
    n = nrow(Xtilde)
    s1 = n * sum(log(sigma^2 * (r * D + 1 - r)))
    s2 = sum(apply(Xtilde, 1, function(x) sum(x^2 / (sigma^2 * (r * D + 1 - r)))))
    l = -.5 * (s1 + s2) - (p / 2) * log(2 * pi)
    return(l)
}


#' Variance along eigenvectors
#' 
#' Find the variance of the data along each of the eigenvectors
#' @param X The data, each row a sample. 
#' @param Q The inner product matrix, either as a matrix or as its
#' eigendecomposition (the output from eigen).
#'
#' @return A vector containing the variances along each of the
#' eigenvectors.
#' @export
varianceOnEvecs <- function(X, Q) {
    if(is.list(Q) & !is.null(Q$vectors) & !is.null(Q$values)) {
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
#' Estimate variance components in a two-parameter model where X ~
#' N(0, sigma^2 (r1 Q + (1 - r1) (r2 Q^(-1) + (1 - r2) I)))
#'
#' @param X An n x p matrix with rows corresponding to observations.
#' @param Q A p x p psd matrix giving the structure.
#' 
#' @return A vector with r1 and r2
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
