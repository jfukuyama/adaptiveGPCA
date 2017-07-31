#' Adaptive gPCA
#'
#' Performs adaptive generalized PCA, a dimensionality-reduction
#' method which takes into account similarities between the
#' variables. See \href{https://arxiv.org/abs/1702.00501}{Fukuyama,
#' J. (2017)} for more details.
#'
#' @param X A \eqn{n \times p} data matrix. 
#' @param Q A \eqn{p \times p} similarity matrix on the variables defining
#' an inner product on the rows of \code{X}, can also be given as an
#' eigendecomposition (formatted as the output from \code{eigen}).
#' @param k The number of components to return.
#' @param weights A vector of length \eqn{n} containing weights for
#' the rows of \code{X}.
#' @return A list containing the row/sample scores (\code{U}), the
#' variable loadings (\code{QV}), the proportion of variance explained
#' by each of the principal components (\code{vars}), the value of
#' \eqn{r} that was used (\code{r}).
#' @examples
#' data(AntibioticSmall)
#' out.agpca = adaptivegpca(AntibioticSmall$X, AntibioticSmall$Q, k = 2)
#' @export
adaptivegpca <- function(X, Q, k = 2, weights = rep(1, nrow(X))) {
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
    auto = estimateComponents(X, Q, Qeig = Qeig)
    r = auto$r
    sigma2 = auto$sigma^2
    evals = (rep(1/(sigma2 * (1-r)), ncol(X)) + (sigma2 * r)^(-1) * Qeig$values^(-1))^(-1)
    out.gpca = gpcaEvecs(X, evecs, evals, weights, k)
    out = list(V = out.gpca$V, U = out.gpca$U, QV = out.gpca$QV,
                lambda = out.gpca$lambda, vars = out.gpca$vars,
               r = r, evals = evals, sig = auto$sigma)
    class(out) = "adaptivegpca"
    return(out)
}


#' Print an adaptivegpca object
#'
#' @param x \code{adaptivegpca} object.
#' @param ... Not used.
#' @method print adaptivegpca
#' @export
print.adaptivegpca <- function(x, ...) {
    cat("An object of class adaptivegpca\n")
    cat("-------------------------------\n")
    cat(paste("Number of axes:", ncol(x$V), "\n"))
    cat(paste("Value of r chosen:", round(x$r, digits = 3), "\n"))
    cat(paste("Fraction of variance explained\nby first", ncol(x$V), "axes:\n"))
    cat(round(x$vars[1:ncol(x$V)], digits = 3), "\n")
}

#' Plot an adaptivegpca object
#'
#' Plots the output from \code{\link{adaptivegpca}}, either a scree
#' plot, the samples, or the variables.
#' 
#' @param x An object of class \code{adaptivegpca}
#' @param type What type of plot to make. \code{scree} will make a
#'     scree plot showing the eigenvalues, \code{samples} will plot
#'     the samples, and \code{variables} will plot the variables.
#' @param axes Which axes to plot. 
#' @param ... Not used.
#' @export
#' @importFrom ggplot2 ggplot aes_string geom_point ylab xlab
#' @examples
#' data(AntibioticSmall)
#' out.agpca = adaptivegpca(AntibioticSmall$X, AntibioticSmall$Q, k = 2)
#' plot(out.agpca)
#' plot(out.agpca, type = "samples")
#' plot(out.agpca, type = "variables")
#' @method plot adaptivegpca
plot.adaptivegpca <- function(x, type = c("scree", "samples", "variables"), axes = c(1,2), ...) {
    type = match.arg(type)
    if(type == "scree") {
        p = ggplot(data.frame(index = seq_along(x$vars), evals = x$vars), aes_string(x = "index", y = "evals")) +
            geom_point() + ylab("Fraction of variance explained") + xlab("Axis")
    } else {
        check_axes(axes, x)
        axis.names = paste("Axis", axes, sep = "")
        axis.labels = paste("Axis ", axes, ": ", round(x$vars[axes] * 100, digits = 1), "%", sep = "")
        if(type == "samples") {
            p = ggplot(data.frame(x$U[,axes])) +
                geom_point(aes_string(x = axis.names[1], y = axis.names[2])) +
                xlab(axis.labels[1]) + ylab(axis.labels[2])
        } else if(type == "variables") {
            p = ggplot(data.frame(x$QV[,axes])) +
                geom_point(aes_string(x = axis.names[1], y = axis.names[2])) +
                xlab(axis.labels[1]) + ylab(axis.labels[2])
        } else {
            stop("Not a valid plotting type.")
        }
    }
    return(p)
}


#' Check that the axes specified are valid
#' @param axes A set of user-specified axes.
#' @param x Object of class \code{adaptivegpca}.
#' @keywords internal
check_axes <- function(axes, x) {
    ## check that the axes are all integers
    if(any(axes != floor(axes)))
        stop("Please specify axes as whole numbers.")
    if(length(axes) != 2)
        stop("Please specify two axes to plot.")
    if(max(axes) > ncol(x$V) | min(axes) < 1)
        stop(paste("Please specify axes within 1 and", ncol(x$V)))
    return(TRUE)
}
    

#' Make a sequence of ordinations
#'
#' Creates a sequence of gPCA data representations. One end of the
#' sequence (\eqn{r = 0}) doesn't do any regularization according to
#' the variable structure (and so is just standard PCA), and the other
#' (\eqn{r = 1}) does a maximal amount of regularization according to
#' the variable structure.
#' 
#' @param X A data matrix of size \eqn{n \times p}. 
#' @param Q A \eqn{p \times p} similarity matrix defining an inner
#' product on the rows of \code{X}.
#' @param weights A vector of weights for the rows of \code{X}. 
#' @param k The number of components to compute for each ordination.
#' @param rvec The values of \eqn{r} for which to make the ordinations.
#' @param findReflections Whether or not flip the axes so as to make
#' neighboring ordinations as close as possible. If \code{k} is very
#' large this should be false since all possible axis combinations are
#' searched over.
#' @param returnLong Return a long data frame with the
#' samples/variables instead of a list of data frames.
#' @param sampledata Extra sample data to be included along with the
#' sample scores.
#' @param variabledata Extra variable data to be included along with
#' the variable loadings.
#' 
#' @return A list containing elements for the sample points
#' (\code{locationList}), the species points (\code{speciesList}), and
#' the variance fractions (\code{varsList}). Each element is itself a
#' list of data frames (location/species points) or of vectors (for
#' the variances).
#' @examples
#' data(AntibioticSmall)
#' out.ff = gpcaFullFamily(AntibioticSmall$X, AntibioticSmall$Q, k = 2)
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
#' @return A vector of length \code{ncol(df1)}: Multiplying the first
#' column of df2 by the first element and multiplying the second
#' column of df2 by the second element and so on gives the optimal
#' reflection.
#' @keywords internal
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
#' @param evecs Eigenvectors of \code{Q}, the inner product/similarity
#' matrix.
#' @param evals Eigenvalues of \code{Q}. 
#' @param D Sample weights
#' @param k The number of components to return.
#' @keywords internal
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
#' Performs standard gPCA with \code{k} components on a data matrix \code{X} with
#' row inner product \code{Q} and weights \code{D}.
#'
#' @param X A data matrix of size \eqn{n \times p}. 
#' @param Q An inner product matrix for the rows, either as a \eqn{p
#' \times p} matrix or an eigendecomposition of such a matrix.
#' @param D Sample weights, a vector of length \eqn{n}. 
#' @param k The number of components to return.
#'
#' @return A list with variable loadings on the principal axes
#' (\code{QV}), sample/row scores (\code{U}), the fraction of the
#' variance explained by each of the axes (\code{vars}).
#' @examples
#' data(AntibioticSmall)
#' out.gpca = gpca(AntibioticSmall$X, AntibioticSmall$Q, k = 2)
#' @export
gpca <- function(X, Q, D = rep(1, nrow(X)), k) {
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
#' @return A list with the normalized matrix (\code{Xtilde}) and the
#' row weights (\code{D}).
#' @keywords internal
normalizeMatrix <- function (X) 
{
    wL = rowSums(X) / sum(X)
    X = diag((1 / wL)) %*% X / sum(X)
    Xtilde = (diag(1, nrow(X)) - matrix(1, nrow = nrow(X), ncol = 1) %*% 
        matrix(wL, nrow = 1)) %*% X
    D = wL/max(wL)
    return(list(Xtilde = Xtilde, D = D))
}


#' Make the input matrices for adaptive gPCA
#'
#' Takes a phyloseq object and creates the matrices necessary to do
#' adaptive gPCA.
#'
#' @param physeq A \code{\link[phyloseq]{phyloseq}} object, from the
#' phyloseq package.
#' @param ca If TRUE, do the normalization as for correspondence
#' analysis (transform counts to relative abundances, compute sample
#' weights, center the relative abundances according to the sample
#' weights). Otherwise, simply center the data. 
#'
#' @return A list of the matrix to perform adaptive gPCA on
#' (\code{X}), the species similarity matrix (\code{Q}), and the
#' sample weights (\code{weights}).
#' @importFrom phyloseq otu_table taxa_are_rows phy_tree
#' @importFrom ape vcv
#' @examples
#' data(AntibioticPhyloseq)
#' pp = processPhyloseq(AntibioticPhyloseq)
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
        weights = nm$D
    } else {
        Xtilde = scale(X, scale = FALSE)
        weights = rep(1, nrow(X))
    }
    return(list(X = Xtilde, Q = Q, weights = weights))
}
