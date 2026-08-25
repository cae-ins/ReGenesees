test_that("calibration system solver uses a fast and accurate full-rank path", {
    system <- crossprod(matrix(c(1, 2, 0, 1, 3, 1), ncol = 2))
    rhs <- c(2, -1)

    solution <- ReGenesees:::.calibration.solve(system, rhs)

    expect_identical(attr(solution, "solver"), "chol")
    expect_identical(attr(solution, "rank"), 2L)
    expect_equal(as.numeric(solution), as.numeric(solve(system, rhs)),
                 tolerance = 1e-13)
})

test_that("calibration system solver preserves generalized-inverse behaviour", {
    system <- matrix(c(1, 1, 1, 1), nrow = 2)
    rhs <- c(2, 2)

    solution <- ReGenesees:::.calibration.solve(
        system, rhs, tol = sqrt(.Machine$double.eps)
    )
    reference <- c(1, 1)

    expect_identical(attr(solution, "solver"), "svd")
    expect_identical(attr(solution, "rank"), 1L)
    expect_equal(as.numeric(solution), as.numeric(reference), tolerance = 1e-13)
})

test_that("calibration system solver handles multiple right-hand sides", {
    system <- crossprod(matrix(c(1, 2, 0, 1, 3, 1), ncol = 2))
    rhs <- diag(2)

    solution <- ReGenesees:::.calibration.solve(system, rhs)

    expect_identical(attr(solution, "solver"), "chol")
    expect_equal(as.vector(solution), as.vector(solve(system)), tolerance = 1e-13)
})

test_that("damped raking converges after an otherwise overflowing Newton step", {
    mm <- matrix(1, nrow = 10, ncol = 1,
                 dimnames = list(NULL, "(Intercept)"))
    population <- c(`(Intercept)` = 1e6)

    g <- ReGenesees:::ez.grake(
        sample.total = colSums(mm), mm = mm, ww = rep(1, 10),
        calfun = ReGenesees:::cal.raking, bounds = c(-Inf, Inf),
        population = population, epsilon = 1e-9, maxit = 50,
        sigma2 = rep(1, 10)
    )

    expect_null(attr(g, "failed"))
    expect_true(all(is.finite(g)))
    expect_gt(attr(g, "step.halvings"), 0L)
    expect_equal(sum(g), unname(population), tolerance = 1e-3)
})

test_that("public calibration methods retain upstream reference results", {
    data("data.examples", package = "ReGenesees")
    design <- e.svydesign(
        data = example, ids = ~towcod + famcod,
        strata = ~SUPERSTRATUM, weights = ~weight
    )
    model <- ~marstat:sex - 1
    model.matrix <- model.matrix(model, design$variables)
    expected.first <- c(
        486.22401127050148, 483.31816693082862, 483.31816693082862,
        483.31816693082862, 483.31816693082862, 483.31816693082862,
        483.27659214114203, 480.43333563908106, 480.43333563908106,
        483.31816693082862
    )

    old.status <- if (exists("ecal.status", envir = .GlobalEnv,
                             inherits = FALSE)) {
        get("ecal.status", envir = .GlobalEnv)
    } else NULL
    if (exists("ecal.status", envir = .GlobalEnv, inherits = FALSE))
        rm("ecal.status", envir = .GlobalEnv)
    on.exit({
        if (is.null(old.status)) {
            if (exists("ecal.status", envir = .GlobalEnv, inherits = FALSE))
                rm("ecal.status", envir = .GlobalEnv)
        } else {
            assign("ecal.status", old.status, envir = .GlobalEnv)
        }
    }, add = TRUE)

    for (method in c("linear", "raking", "logit")) {
        args <- list(
            design = design, df.population = pop03,
            calmodel = model, calfun = method
        )
        if (identical(method, "logit")) args$bounds <- bounds
        calibrated <- do.call(e.calibrate, args)
        expect_false(exists("ecal.status", envir = .GlobalEnv,
                            inherits = FALSE))
        expect_true(is.list(attr(calibrated, "ecal.status")))
        calibrated.weights <- weights(calibrated)
        achieved <- colSums(model.matrix * calibrated.weights)
        target <- as.numeric(pop03[1, ])
        max.relative.error <- max(
            abs(achieved - target) / (1 + abs(target))
        )

        expect_lt(max.relative.error, 1e-7)
        expect_equal(
            calibrated.weights[1:10], expected.first,
            tolerance = 3e-9
        )
    }
})
