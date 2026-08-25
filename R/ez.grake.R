`ez.grake` <-
function (sample.total, mm, ww, calfun, eta = rep(0, NCOL(mm)), bounds, population, 
    epsilon, maxit, sigma2)
###################################################################
#  Versione modificata della funzione grake del package survey.   #
#  NOTA: Per motivi di efficienza sono state eliminate alcune     #
#        funzionalita' originali NON NECESSARIE per e.calibrate   #
#        (require di MASS e attr(g,"eta")).                       #
###################################################################
{
    if (!inherits(calfun, "calfun")) 
        stop("'calfun' must be of class 'calfun'")
    Fm1 <- calfun$Fm1
    dF <- calfun$dF
    xeta <- drop(mm %*% eta)
    fm1 <- Fm1(xeta, bounds, sigma2)
    g <- 1 + fm1
    deriv <- dF(xeta, bounds, sigma2)
    iter <- 0L
    step.halvings <- 0L
    solvers <- character(0)
    repeat ({
        iter <- iter + 1L
        Tmat <- crossprod(mm * ww * deriv, mm)
        misfit <- population - sample.total - colSums(mm * ww * fm1)
        current.error <- max(abs(misfit)/(1 + abs(population)))
        deta <- .calibration.solve(Tmat, misfit,
                                   tol = 256 * .Machine$double.eps)
        solvers <- c(solvers, attr(deta, "solver"))
        deta <- as.numeric(deta)

        # Damped Newton step.  Full Newton steps remain the default; a step is
        # halved only if it produces non-finite values or increases the largest
        # normalized calibration residual.  This prevents avoidable overflows
        # in raking/logit problems with population totals far from the sample.
        step <- 1
        accepted <- FALSE
        for (backtrack in 0:30) {
            eta.new <- eta + step * deta
            xeta.new <- drop(mm %*% eta.new)
            fm1.new <- Fm1(xeta.new, bounds, sigma2)
            g.new <- 1 + fm1.new
            if (all(is.finite(g.new), is.finite(fm1.new))) {
                misfit.new <- population - sample.total -
                    colSums(mm * ww * fm1.new)
                new.error <- max(abs(misfit.new)/(1 + abs(population)))
                if (is.finite(new.error) &&
                    (new.error <= current.error || step <= 2^-30)) {
                    accepted <- TRUE
                    break
                }
            }
            step <- step / 2
            step.halvings <- step.halvings + 1L
        }

        if (!accepted) {
            eta.new <- eta
            xeta.new <- xeta
            fm1.new <- fm1
            g.new <- g
            misfit.new <- misfit
        }

        eta <- eta.new
        xeta <- xeta.new
        fm1 <- fm1.new
        g <- g.new
        misfit <- misfit.new
        deriv <- dF(xeta, bounds, sigma2)
        if (all(abs(misfit)/(1 + abs(population)) < epsilon)) 
            break
        if (iter >= maxit) {
            achieved <- (abs(misfit)/(1 + abs(population)))
            worst.ind <- which.max(achieved)
            worst.achieved <- achieved[worst.ind]
            warning("Failed to converge: worst achieved epsilon= ", worst.achieved, " in ", 
                iter, " iterations (variable ",names(worst.achieved),"), see the 'ecal.status' attribute of the returned design.")
            # Build a diagnostic structure to be sent to ecal.status
            ko.ind <- which(achieved >= epsilon)
            reldiff <- -(misfit/(1 + population))
            ko.reldiff <- reldiff[ko.ind]
            ko.names <- names(ko.reldiff)
            ko.pop <- population[ko.ind]
            ko.est <- (population - misfit)[ko.ind]
            ko.diff <- ko.est - ko.pop
            diagnosys <- data.frame(`Variable` = ko.names,
                                    `Population.Total` = ko.pop,
                                    `Achieved.Estimate` = ko.est,
                                    `Difference` = ko.diff,
                                    `Relative.Difference` = ko.reldiff,
                                    row.names = NULL)
            # Set missed constraints indices as useful rownames, i.e. the
            # positions of non convergence cells in df.population for each
            # partition (if any):
            rownames(diagnosys) <- ko.ind
            diagnosys <- diagnosys[order(abs(ko.reldiff), decreasing = TRUE), , drop = FALSE]
            # Attach diagnostics to g as an attribute
            attr(g, "failed") <- achieved
            attr(g, "fail.diagnostics") <- diagnosys
            break
        }
    })
    attr(g, "iterations") <- iter
    attr(g, "step.halvings") <- step.halvings
    attr(g, "solvers") <- unique(solvers)
    g
}
