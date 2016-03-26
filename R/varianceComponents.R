

#' Estimate parameters in hierarchical model
#'
#' Estimates the values of r and sigma in a model X sim N(0, sigma^2 (r Q + (1 - r)I)) using bisection. 
#'
#' @param X A n x p matrix with n samples in p-dimensional space.
#' @param Q The prior variance on the means.
#' @param maxit The maximum number of bisections to try.
#' @param tol Stop if the absolute value of the derivative is less than this. œ
#'
#' 
#' @return A list with r and sigma.
#' @export
estimateComponents <- function(X, Q, maxit = 8, tol = 10^(-10), plot = FALSE, Qeig = NULL) {
    if(is.null(Qeig)) {
        Qeig = eigen(Q, symmetric = TRUE)
    }
    Xtilde = X %*% Qeig$vectors
    # first search along a grid to get an approximation of where the max is
    rvec = seq(0,1,length.out = 100)
    l = sapply(rvec, function(r) likelihoodR(Xtilde, r, Qeig$values))
    # then do bisection in the area where this max is
    start = which.max(l)
    loidx = max(start - 1, 1)
    hiidx = min(start + 1, 100)
    lo = rvec[loidx]; hi = rvec[hiidx]
    for(i in 1:maxit) {
        r = (lo + hi) / 2
        g = gradLik(Xtilde, r, Qeig$values)
        if(g > 0) {
            lo = r
        } else {
            hi = r
        }
        if(abs(g) <= tol) break
    }
    sigma2 = sigma2OfR(Xtilde, r, Qeig$values)
    if(plot){
        plot(l ~ rvec, type = 'l')
        
    }
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
    #s1 = n * p * 2 * log(sigma) + n * sum(log(r * D + 1 - r))
    s1 = n * sum(log(sigma^2 * (r * D + 1 - r)))
    #prec = solve(sigma^2 * (r * Q + (1 - r) * diag(rep(1,p))))
    #s2 = sum(apply(X, 1, function(x) x %*% prec %*% x))
    s2 = sum(apply(Xtilde, 1, function(x) sum(x^2 / (sigma^2 * (r * D + 1 - r)))))
    l = -.5 * (s1 + s2) - (p / 2) * log(2 * pi)
    return(l)
}
