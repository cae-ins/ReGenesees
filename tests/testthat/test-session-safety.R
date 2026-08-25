test_that("package-local contrast handling leaves global contrasts unchanged", {
    before <- getOption("contrasts")
    data <- data.frame(group = ordered(c("a", "b", "c")))

    actual <- ReGenesees:::.rg.model.matrix(~group, data)
    expected <- stats::model.matrix(
        ~group, data,
        contrasts.arg = list(group = "contr.treatment")
    )

    expect_identical(getOption("contrasts"), before)
    expect_equal(unname(actual), unname(expected))
})

test_that("GVF database mutates without unlocking the namespace", {
    original <- GVF.db$get(verbose = FALSE)
    on.exit(GVF.db$assign(original, verbose = FALSE), add = TRUE)

    expect_true(is.environment(GVF.db))
    GVF.db$insert("CV ~ I(1/Y^3)", Estimator.kind = "Total",
                  Resp.to.CV = "resp", verbose = FALSE)
    updated <- GVF.db$get(verbose = FALSE)

    expect_equal(NROW(updated), NROW(original) + 1L)
    expect_identical(tail(updated$GVF.model, 1L), "CV ~ I(1/Y^3)")
})
