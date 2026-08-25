# Reproducible benchmark for the CAE-INS calibration solver improvements.
# Run from the package root with:
#   Rscript tools/benchmarks/benchmark-calibration.R

set.seed(20260825)
devtools::load_all(".", quiet = TRUE)

legacy.grake <- function(sample.total, mm, ww, calfun,
                         eta = rep(0, NCOL(mm)), bounds, population,
                         epsilon, maxit, sigma2) {
    Fm1 <- calfun$Fm1
    dF <- calfun$dF
    xeta <- drop(mm %*% eta)
    g <- 1 + Fm1(xeta, bounds, sigma2)
    iter <- 1
    repeat {
        Tmat <- crossprod(mm * ww * dF(xeta, bounds, sigma2), mm)
        misfit <- population - sample.total -
            colSums(mm * ww * Fm1(xeta, bounds, sigma2))
        deta <- MASS::ginv(Tmat, tol = 256 * .Machine$double.eps) %*%
            misfit
        eta <- eta + deta
        xeta <- drop(mm %*% eta)
        g <- 1 + Fm1(xeta, bounds, sigma2)
        misfit <- population - sample.total -
            colSums(mm * ww * Fm1(xeta, bounds, sigma2))
        if (all(abs(misfit) / (1 + abs(population)) < epsilon)) break
        iter <- iter + 1
        if (iter > maxit) {
            attr(g, "failed") <- abs(misfit) / (1 + abs(population))
            break
        }
    }
    g
}

elapsed.per.call <- function(expr, times) {
    expr <- substitute(expr)
    system.time(for (i in seq_len(times)) eval(expr, parent.frame()))[["elapsed"]] /
        times
}

median.paired <- function(legacy, modern, times = 7L) {
    legacy <- substitute(legacy)
    modern <- substitute(modern)
    timings <- matrix(NA_real_, nrow = times, ncol = 2L,
                      dimnames = list(NULL, c("legacy", "modern")))
    for (i in seq_len(times)) {
        order <- if (i %% 2L) c("legacy", "modern") else c("modern", "legacy")
        for (implementation in order) {
            gc(FALSE)
            expression <- if (implementation == "legacy") legacy else modern
            timings[i, implementation] <- system.time(
                eval(expression, parent.frame())
            )[["elapsed"]]
        }
    }
    apply(timings, 2L, median)
}

# Solver-only benchmark: this isolates the factorization improvement.
factor <- matrix(rnorm(120 * 120), nrow = 120)
system <- crossprod(factor) + diag(120)
rhs <- rnorm(120)
legacy.solver <- elapsed.per.call(
    MASS::ginv(system, tol = 256 * .Machine$double.eps) %*% rhs,
    times = 30
)
modern.solver <- elapsed.per.call(
    ReGenesees:::.calibration.solve(system, rhs),
    times = 30
)

# End-to-end generalized-raking benchmark with a known feasible solution.
n <- 40000
p <- 60
mm <- cbind(`(Intercept)` = 1,
            matrix(rnorm(n * (p - 1)), nrow = n, ncol = p - 1))
ww <- runif(n, 0.5, 2)
sigma2 <- rep(1, n)
sample.total <- colSums(mm * ww)
eta.true <- c(0.02, rnorm(p - 1, sd = 0.002))
g.true <- exp(drop(mm %*% eta.true))
population <- colSums(mm * ww * g.true)

paired.time <- median.paired(
    legacy.grake(
        sample.total, mm, ww, ReGenesees:::cal.raking,
        bounds = c(-Inf, Inf), population = population,
        epsilon = 1e-9, maxit = 50, sigma2 = sigma2
    ),
    ReGenesees:::ez.grake(
        sample.total, mm, ww, ReGenesees:::cal.raking,
        bounds = c(-Inf, Inf), population = population,
        epsilon = 1e-9, maxit = 50, sigma2 = sigma2
    ),
    times = 7L
)
legacy.time <- paired.time[["legacy"]]
modern.time <- paired.time[["modern"]]

legacy.g <- legacy.grake(
    sample.total, mm, ww, ReGenesees:::cal.raking,
    bounds = c(-Inf, Inf), population = population,
    epsilon = 1e-9, maxit = 50, sigma2 = sigma2
)
modern.g <- ReGenesees:::ez.grake(
    sample.total, mm, ww, ReGenesees:::cal.raking,
    bounds = c(-Inf, Inf), population = population,
    epsilon = 1e-9, maxit = 50, sigma2 = sigma2
)
modern.error <- max(abs(colSums(mm * ww * modern.g) - population) /
                    (1 + abs(population)))

# Stress case: the undamped first Newton step overflows exp(), whereas the
# damped implementation backtracks to a finite, convergent step.
stress.mm <- matrix(1, nrow = 10, ncol = 1)
stress.population <- 1e6
legacy.stress <- try(suppressWarnings(legacy.grake(
    colSums(stress.mm), stress.mm, rep(1, 10), ReGenesees:::cal.raking,
    bounds = c(-Inf, Inf), population = stress.population,
    epsilon = 1e-9, maxit = 50, sigma2 = rep(1, 10)
)), silent = TRUE)
modern.stress <- ReGenesees:::ez.grake(
    colSums(stress.mm), stress.mm, rep(1, 10), ReGenesees:::cal.raking,
    bounds = c(-Inf, Inf), population = stress.population,
    epsilon = 1e-9, maxit = 50, sigma2 = rep(1, 10)
)

results <- data.frame(
    metric = c(
        "solver_legacy_ms", "solver_modern_ms", "solver_speedup",
        "grake_legacy_s", "grake_modern_s", "grake_speedup",
        "grake_max_abs_weight_difference", "grake_max_relative_constraint_error",
        "stress_legacy_converged", "stress_modern_converged",
        "stress_modern_step_halvings"
    ),
    value = c(
        1000 * legacy.solver, 1000 * modern.solver,
        legacy.solver / modern.solver,
        legacy.time, modern.time, legacy.time / modern.time,
        max(abs(legacy.g - modern.g)), modern.error,
        as.numeric(!inherits(legacy.stress, "try-error") &&
                   is.null(attr(legacy.stress, "failed")) &&
                   all(is.finite(legacy.stress))),
        as.numeric(is.null(attr(modern.stress, "failed")) &&
                   all(is.finite(modern.stress))),
        attr(modern.stress, "step.halvings")
    )
)

dir.create("quality_reports", showWarnings = FALSE)
print(results, row.names = FALSE)
write.result <- try(suppressWarnings(write.csv(
    results, "quality_reports/calibration_benchmark.csv", row.names = FALSE
)), silent = TRUE)
if (inherits(write.result, "try-error"))
    message("Benchmark results could not be written; printed results remain authoritative.")
