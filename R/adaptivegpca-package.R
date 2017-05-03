#' Adaptive gPCA
#'
#' Performs adaptive generalized PCA, a dimensionality-reduction
#' method which takes into account similarities between the
#' variables. See Fukuyama (2017) at https://arxiv.org/abs/1702.00501
#' for more details.
#'
#' @param X A data matrix
#' @param Q The inner product matrix on the rows of X, can be given as
#' an eigendecomposition (formatted as the output from eigen()).
#' @param weights Data weights for the rows of X. 
#' @param k The number of components to return.
#'
#' @export
adaptivegPCA <- function(X, Q, weights = rep(1, nrow(X)), k = 2) {
    if(is.matrix(Q)) {
        Qeig = eigen(Q, symmetric = TRUE)
    } else if(is.list(Q) & !is.null(Q$vectors) & !is.null(Q$values)) {
        Qeig = Q
    } else {
        stop("Q is not formatted correctly")
    }
    # normalize so that the trace of Q is the same as the trace of the identity matrix
    Qeig$values = ncol(X) * Qeig$values / sum(Qeig$values)
    evecs = Qeig$vectors
    auto = estimateComponents(X, Q, maxit = 30, Qeig = Qeig)
    r = auto$r
    sigma2 = auto$sigma^2
    evals = (rep(1/(sigma2 * (1-r)), ncol(X)) + (sigma2 * r)^(-1) * Qeig$values^(-1))^(-1)
    out.gpca = gpcaEvecs(X, evecs, evals, weights, k)
    return(list(V = out.gpca$V, U = out.gpca$U, QV = out.gpca$QV,
                lambda = out.gpca$lambda, vars = out.gpca$vars,
                r = r, evals = evals, sig = auto$sigma))
}

#' Make a sequence of ordinations
#'
#' Creates a sequence of gPCA data representations. One end of the
#' sequence doesn't take into account the variable structure at all
#' (PCA), and the other weights the structure of the variables very
#' highly.
#' 
#' @param X A data matrix.
#' @param Q An inner product on the rows of X.
#' @param weights A vector of weights for the rows of X. 
#' @param k The number of components to compute for each ordination.
#' @param rvec The values of r for which to make the ordinations.
#' @param findReflections Whether or not flip the axes so as to make
#' neighboring ordinations as close as possible. If k is very large
#' this should be false since all possible axis combinations are
#' searched over.
#' @param returnLong Return a long data frame with the
#' samples/variables instead of a list of data frames.
#' @param sampledata Extra sample data to be included along with the
#' sample scores.
#' @param variabledata Extra variable data to be included along with
#' the variable loadings.
#' 
#' @return A list with three components: one holding the location
#' points, one holding the species points, and one holding the
#' variance fractions. Each list is itself a list of data frames
#' (location/species points) or of vectors (for the variances).
#' @export
gpcaFullFamily <- function(X, Q, weights = rep(1, nrow(X)), k = 2,
        rvec = (0:100)/100, findReflections = TRUE,
        returnLong = FALSE, sampledata = NULL, variabledata = NULL) {
    if(is.matrix(Q)) {
        Qeig = eigen(Q, symmetric = TRUE)
    } else if(is.list(Q) & !is.null(Q$vectors) & !is.null(Q$values)) {
        Qeig = Q
    } else {
        stop("Q is not formatted correctly")
    }
    Qeig$values = Qeig$values / sum(Qeig$values) * ncol(X)
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
        gpcaout = gpcaEvecs(X, evecs, evals = evals, D = weights, k = 2)
        locationList[[as.character(r)]] = gpcaout$U
        speciesList[[as.character(r)]] = gpcaout$QV
        colnames(locationList[[as.character(r)]]) = paste("Axis", 1:ncol(gpcaout$U), sep = "")
        colnames(speciesList[[as.character(r)]]) = paste("Axis", 1:ncol(gpcaout$QV), sep = "")
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
    if(returnLong) {
        locationsAugmented = lapply(1:length(rvec), function(i) {
            df = data.frame(locationList[[i]], r = rvec[i])
            if(!is.null(sampledata))
                df = cbind(df, sampledata)
            return(df)
        })
        speciesAugmented = lapply(1:length(rvec), function(i) {
            df = data.frame(speciesList[[i]], r = rvec[i])
            if(!is.null(variabledata))
                df = cbind(df, variabledata)
            return(df)
        })
        locationsfull = Reduce(rbind, locationsAugmented)
        speciesfull = Reduce(rbind, speciesAugmented)
        return(list(locations = locationsfull, species = speciesfull, vars = varsList,
                    X = X, Qeig = Qeig))
    }
    return(list(locations = locationList, species = speciesList, vars = varsList,
                X = X, Qeig = Qeig))
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

#' gPCA using pre-computed eidendecomposition
#' 
#' Performs gPCA with pre-computed eigenvectors and eigenvalues.
#' 
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
    QV = evecs %*% (sweep(t(evecs), 1, STATS = evals, FUN = "*") %*% as.matrix(V))
    colnames(U) = colnames(V) = colnames(QV) = paste("Axis", 1:ncol(U), sep = "")
    return(list(V = V, U = U, QV = QV, lambda = Jsvd$d, vars = Jsvd$d^2 / sum(Jsvd$d^2)))
}

#' gPCA
#'
#' Performs standard gPCA with k components on a data matrix X with
#' row inner product Q and column inner product D.
#'
#' @param X A data matrix
#' @param Q An inner product matrix for the rows, either as a matrix
#' or an eigendecomposition. 
#' @param D An inner product matrix for the columns.
#' @param k The number of components to return.
#'
#' @return A list with principal components, scores, and variable positions.
#' @export
gPCA <- function(X, Q, D = rep(1, nrow(X)), k) {
    if(is.matrix(Q)) {
        Qeig = eigen(Q, symmetric = TRUE)
    } else if(is.list(Q) & !is.null(Q$vectors) & !is.null(Q$values)) {
        Qeig = Q
    } else {
        stop("Q is not formatted correctly")
    }
    return(gpcaEvecs(X, Qeig$vectors, Qeig$values, D, k))
}

#' Normalizes a matrix.
#'
#' Normalizes a count matrix X for correspondence analysis and returns
#' the corresponding metric. The normalization for X is as follows:
#' First the row sums of X are computed, giving the weights for each
#' sample. These weights are stored in a matrix D, which defines an
#' inner product on the columns of X. Then the vectors of counts
#' stored in the rows of X are replaced with proportions, and the
#' resulting matrix is centered according to the inner product defined
#' by D. Both the centered data matrix and D are returned to the user. 
#' 
#' @param X The matrix to be normalized.
#'
#' @return A list with the normalized matrix (Xtilde) and the
#' corresponding metric (D).
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
#' Takes a phyloseq object and creates the matrices necessary to do a
#' generalized PCA.
#'
#' @param physeq A phyloseq object.
#' @param ca If TRUE, do the normalization as for correspondence
#' analysis (transform counts to relative abundances, compute sample
#' weights, center the relative abundances according to the sample
#' weights). Otherwise, simply center the data. 
#'
#' @return A list of the matrix to perform gPCA on (X), the norm for
#' the rows (Q), and the sample weights (weights).
#' @importFrom phyloseq otu_table taxa_are_rows phy_tree
#' @importFrom ape vcv
#' @export
processPhyloseq <- function (physeq, ca = FALSE) 
{
    if (taxa_are_rows(physeq)) 
        X = as(t(otu_table(physeq)), "matrix")
    else X = as(otu_table(physeq), "matrix")
    Q = vcv(phy_tree(physeq), scale = FALSE)
    Q = Q / sum(diag(Q)) * nrow(Q)
    if(ca) {
        nm = normalizeMatrix(X)
        Xtilde = nm$Xtilde
        weights = diag(nm$D)
    } else {
        Xtilde = scale(X, scale = FALSE)
        weights = rep(1, nrow(X))
    }
    return(list(X = Xtilde, Q = Q, weights = weights))
}
