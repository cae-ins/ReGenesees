test_that("rowsum PSU aggregation matches aggregate.data.frame", {
    data <- data.frame(
        strata = factor(c("B", "A", "B", "A", "B")),
        id = c(2, 1, 2, 1, 3),
        y = c(1, 2, 3, 4, 5),
        w = c(2, 3, 4, 5, 6)
    )
    fast <- ReGenesees:::.delta.rowsum(
        data, c("y", "w"), "id", "strata"
    )
    names(fast)[names(fast) == "strata"] <- "strata"
    reference <- aggregate(
        data[, c("y", "w")],
        by = list(strata = data$strata, id = data$id), FUN = sum
    )
    fast <- fast[order(fast$strata, fast$id), ]
    reference <- reference[order(reference$strata, reference$id), ]
    row.names(fast) <- row.names(reference) <- NULL

    expect_equal(fast, reference)
})

test_that("vectorized svyDelta domain path retains upstream results", {
    data("Delta.el", package = "ReGenesees")
    s1$domain <- factor(ifelse(s1$x > median(s1$x), "high", "low"))
    s2$domain <- factor(ifelse(s2$x > median(s2$x), "high", "low"))
    design1 <- e.svydesign(
        ids = ~id, strata = ~strata, weights = ~w, data = s1
    )
    design2 <- e.svydesign(
        ids = ~id, strata = ~strata, weights = ~w, data = s2
    )

    result <- svyDelta(
        expression(y.2 - y.1), design1, design2, by = ~domain,
        rho.STRAT = "noJump", vartype = c("se", "var"), conf.int = TRUE
    )

    expected <- matrix(c(
        -176.105649991724, 472.675347485683, -1102.53230744362,
        750.321007460170, 223421.984120711,
        660.343764779205, 465.828719848298, -252.663749087857,
        1573.35127864627, 216996.396235504
    ), nrow = 2, byrow = TRUE)
    expect_equal(
        unname(as.matrix(as.data.frame(result)[, -1, drop = FALSE])),
        expected, tolerance = 1e-12
    )
    expect_identical(as.character(result$domain), c("high", "low"))
})

test_that("vectorized cluster svyDelta retains covariance diagnostics", {
    data("Delta.clus", package = "ReGenesees")
    design1 <- e.svydesign(
        ids = ~id, strata = ~strata, weights = ~w, data = sclus1
    )
    design2 <- e.svydesign(
        ids = ~id, strata = ~strata, weights = ~w, data = sclus2
    )

    result <- svyDelta(
        expression(y.2 - y.1), design1, design2,
        rho.STRAT = "noJump", vartype = "var"
    )
    diagnostics <- attr(result, "details")

    expect_equal(result$Delta, 1067.58256979977, tolerance = 1e-12)
    expect_equal(result$VAR, 72934.3642899499, tolerance = 1e-12)
    expect_equal(diagnostics$rho, 0.922468516503005, tolerance = 1e-12)
    expect_equal(diagnostics$CoV, 404208.26555342, tolerance = 1e-12)
    expect_identical(diagnostics$nc, 3L)
})
