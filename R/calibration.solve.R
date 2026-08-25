`.calibration.solve` <-
function(system, rhs, tol = 256 * .Machine$double.eps)
##########################################################################
# Solve the small symmetric systems arising during calibration.          #
#                                                                        #
# A Cholesky factorization is substantially faster than a full SVD for   #
# the usual full-rank positive-definite case.  Ill-conditioned or rank-  #
# deficient systems fall back to the same truncated-SVD principle used   #
# by MASS::ginv, preserving the robust behaviour needed by calibration.  #
##########################################################################
{
    system <- as.matrix(system)
    rhs.is.vector <- is.null(dim(rhs))
    rhs.names <- names(rhs)
    rhs.dimnames <- dimnames(rhs)
    rhs <- if (rhs.is.vector) matrix(as.numeric(rhs), ncol = 1L) else as.matrix(rhs)

    if (NROW(system) != NCOL(system))
        stop("Calibration system must be a square matrix.")
    if (NROW(rhs) != NROW(system))
        stop("Calibration system and right-hand side have incompatible dimensions.")
    if (!all(is.finite(system)) || !all(is.finite(rhs)))
        stop("Calibration system contains non-finite values.")

    n <- NROW(system)
    if (n == 0L) {
        solution <- if (rhs.is.vector) numeric(0) else matrix(numeric(0), 0L, NCOL(rhs))
        attr(solution, "solver") <- "empty"
        attr(solution, "rank") <- 0L
        return(solution)
    }

    # The Hessian is theoretically symmetric.  Symmetrising removes tiny
    # floating-point asymmetries that can otherwise make chol() reject it.
    system <- (system + t(system)) / 2
    factor <- try(chol(system), silent = TRUE)

    if (!inherits(factor, "try-error") &&
        all(is.finite(factor)) &&
        is.finite(rcond(factor)) &&
        rcond(factor)^2 > tol) {
        solution <- backsolve(factor, forwardsolve(t(factor), rhs))
        residual <- system %*% solution - rhs
        rel.residual <- max(abs(residual) / (1 + abs(rhs)))
        if (all(is.finite(solution)) && rel.residual <= sqrt(tol)) {
            if (rhs.is.vector) {
                solution <- drop(solution)
                names(solution) <- rhs.names
            }
            else {
                dimnames(solution) <- list(rownames(system), rhs.dimnames[[2L]])
            }
            attr(solution, "solver") <- "chol"
            attr(solution, "rank") <- n
            return(solution)
        }
    }

    # Truncated SVD: algebraically equivalent to MASS::ginv(system, tol) %*%
    # rhs, but also records the numerical rank for diagnostics.
    decomposition <- svd(system, nu = n, nv = n)
    positive <- decomposition$d > max(tol * decomposition$d[1L], 0)
    rank <- sum(positive)
    if (rank == 0L) {
        solution <- matrix(0, nrow = n, ncol = NCOL(rhs))
    }
    else {
        solution <- decomposition$v[, positive, drop = FALSE] %*%
            (crossprod(decomposition$u[, positive, drop = FALSE], rhs) /
             decomposition$d[positive])
    }
    if (rhs.is.vector) {
        solution <- drop(solution)
        names(solution) <- rhs.names
    }
    else {
        dimnames(solution) <- list(rownames(system), rhs.dimnames[[2L]])
    }
    attr(solution, "solver") <- "svd"
    attr(solution, "rank") <- rank
    solution
}
