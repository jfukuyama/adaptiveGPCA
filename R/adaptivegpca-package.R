

#' Adaptive gPCA
#'
#' Performs adaptive generalized PCA
#'
#' @param X A data matrix
#' @param Q The inner product matrix on the rows of X, can be given as
#' an eigendecomposition (formatted as the output from eigen()).
#' @param k The number of components to return. 
#'
#' @export
adaptivegPCA <- function(X, Q, k = 2) {
    if(is.matrix(Q)) {
        Qeig = eigen(Q, symmetric = TRUE)
    } else if(is.list(Q) & !is.null(Q$vectors) & !is.null(Q$values)) {
        Qeig = Q
    } else {
        stop("Q is not formatted correctly")
    }
    evecs = Qeig$vectors
    auto = estimateComponents(X, Q, maxit = 30, Qeig = Qeig)
    r = auto$r
    evals = (rep(1/(1-r), ncol(X)) + r^(-1) * Qeig$values^(-1))^(-1)
    J = X %*% evecs
    J = sweep(J, 2, STATS = sqrt(evals), FUN = "*")
    Jsvd = svd(J, nv = k, nu = k)
    V = evecs %*% sweep(Jsvd$v, 1, STATS = evals^(-.5), FUN = "*")
    U = Jsvd$u
    QV = evecs %*% (sweep(t(evecs), 1, STATS = evals, FUN = "*") %*% as.matrix(V))
    return(list(V = V, U = U, QV = QV, vars = Jsvd$d^2 / sum(Jsvd$d^2), r = r))
}



#' Make a sequence of ordinations
#'
#' Gives a list of ordinations going from structured to PCA.
#' 
#' @param X A data matrix.
#' @param Q An inner product on the rows of X. 
#' @param k The number of components to compute for each ordination.
#' @param rvec The values of r for which to make the ordinations.
#' @param findReflections Whether or not flip the axes so as to make
#' neighboring ordinations as close as possible. If k is very large
#' this should be false since all possible axis combinations are
#' searched over.
#' 
#' @return A list with three components: one holding the location
#' points, one holding the species points, and one holding the
#' variance fractions. Each list is itself a list of data frames
#' (location/species points) or of vectors (for the variances).
#' @export
gpcaFullFamily <- function(X, Q, k = 2, rvec = (0:100)/100, findReflections = TRUE) {
    if(is.matrix(Q)) {
        Qeig = eigen(Q, symmetric = TRUE)
    } else if(is.list(Q) & !is.null(Q$vectors) & !is.null(Q$values)) {
        Qeig = Q
    } else {
        stop("Q is not formatted correctly")
    }
    evecs = Qeig$vectors
    locationList = list()
    speciesList = list()
    varsList = list()
    for(r in rvec) {
        if(r == 0) {
            evals = Qeig$values
        } else if (r == 1) {
            evals = rep(1, ncol(X))
        } else {
            evals = (rep(1/(1-r), ncol(X)) + r^(-1) * Qeig$values^(-1))^(-1)
        }
        gpcaout = gpcaEvecs(X, evecs, evals = evals, k = 2)
        locationList[[as.character(r)]] = gpcaout$U
        speciesList[[as.character(r)]] = gpcaout$sp
        varsList[[as.character(r)]] = gpcaout$vars
    }
    if(findReflections) {
        for(i in 2:length(rvec)) {
            s = findReflection(locationList[[i-1]], locationList[[i]])
            for(j in 1:length(s)) {
                locationList[[i]][,j] = locationList[[i]][,j] * s[j]
                speciesList[[i]][,j] = speciesList[[i]][,j] * s[j]
            }
        }
    }
    return(list(locations = locationList, species = speciesList, vars = varsList))
}

#' Find reflection
#'
#' Find a reflection of one data frame so that it most closely matches
#' another.
#'
#' @param df1 The base data frame.
#' @param df2 The data frame that will be reflected across either the
#' x-axis or the y-axis (or both or neither) so that the points in it
#' most closely match df1.
#'
#' @return A vector of length 2: Multiplying the first column of df2
#' by the first element and multiplying the second column of df2 by
#' the second element gives the optimal reflection.
findReflection <- function(df1, df2) {
    l = list()
    for(i in 1:ncol(df1)) {
        l[[i]] = c(-1,1)
    }
    options = expand.grid(l)
    distances = apply(options, 1, function(x) {
        df2mod = df2
        for(i in 1:ncol(df2)) {
            df2mod[,i] = df2mod[,i] * x[i]
        }
        return(sum((df1 - df2mod)^2))
    })
    return(as.numeric(options[which.min(distances),]))
}


#' Performs gPCA with pre-computed eigenvectors and eigenvalues
#' @param X Data matrix.
#' @param evecs Eigenvectors of Q
#' @param evals Eigenvalues of Q
#' @param D Inner product matrix for the columns.
#' @param k The number of components to return. 
gpcaEvecs <- function(X, evecs, evals, D = rep(1, nrow(X)), k) {
    J = X %*% evecs
    J = sweep(J, 2, STATS = sqrt(evals), FUN = "*")
    J = sweep(J, 1, STATS = D, FUN = "*")
    Jsvd = svd(J, nv = k, nu = k)
    evalsginvsqrt = evals^(-.5)
    evalsginvsqrt[evals == 0]= 0
    V = evecs %*% sweep(Jsvd$v, 1, STATS = evalsginvsqrt, FUN = "*")
    U = sweep(Jsvd$u, 1, D^(-.5), FUN = "*")
    sp = evecs %*% (sweep(t(evecs), 1, STATS = evals, FUN = "*") %*% as.matrix(V))
    return(list(V = V, U = U, sp = sp, vars = Jsvd$d^2 / sum(Jsvd$d^2)))
}

#' gPCA
#'
#' Performs standard gPCA
#'
#' @param X A data matrix
#' @param Q An inner product matrix for the rows.
#' @param D An inner product matrix for the columns.
#' @param k The number of components to return.
#'
#' @return A list with principal components, scores, and variable positions.
#' @export
gPCA <- function(X, Q, D = rep(1, nrow(X)), k) {
    Qeig = eigen(Q)
    return(gpcaEvecs(X, Qeig$vectors, Qeig$values, D, k))
}

#' Normalizes a matrix.
#'
#' Normalizes a count matrix and returns the corresponding metric.
#'
#' @param X The matrix to be normalized.
#'
#' @return A list with the normalized matrix and the corresponding metric. 
normalizeMatrix <- function (X) 
{
    wL = rowSums(X) / sum(X)
    X = diag((1 / wL)) %*% X / sum(X)
    Xtilde = (diag(1, nrow(X)) - matrix(1, nrow = nrow(X), ncol = 1) %*% 
        matrix(wL, nrow = 1)) %*% X
    D = diag(wL)
    D = D/max(D)
    return(list(Xtilde = Xtilde, D = D))
}


#' Make the input matrices for PCA
#'
#' Takes a phyloseq object and (optionally) some constraints and
#' creates the matrices necessary to do a generalized PCA.
#'
#' @param physeq A phyloseq object.
#'
#' @importFrom ape vcv
#'
#' @return A list of the matrix to perform gPCA on (X), the norm for
#' the rows (sigma), and the norm for the columns (D).
#' 
#' @export
processPhyseq <- function (physeq) 
{
    if (taxa_are_rows(physeq)) 
        X = t(otu_table(physeq))
    else X = otu_table(physeq)
    Q = vcv(phy_tree(physeq), scale = FALSE)
    Q = Q / sum(diag(Q)) * nrow(Q)
    nm = normalizeMatrix(X)
    Xtilde = nm$Xtilde
    D = nm$D
    return(list(X = Xtilde, Q = Q, D = D))
}


