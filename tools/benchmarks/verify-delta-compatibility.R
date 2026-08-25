# Compare svyDelta results and performance with the upstream source tree.
# Usage from either package tree:
#   Rscript verify-delta-compatibility.R write baseline.rds
#   Rscript verify-delta-compatibility.R compare baseline.rds

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L || !arguments[1L] %in% c("write", "compare"))
    stop("Usage: verify-delta-compatibility.R <write|compare> <baseline.rds>")

mode <- arguments[1L]
baseline.path <- arguments[2L]
source.root <- Sys.getenv("REGENESEES_SOURCE", unset = ".")
devtools::load_all(source.root, quiet = TRUE)

extract <- function(result) {
    list(
        values = unname(as.matrix(as.data.frame(result))),
        names = names(result),
        row.names = row.names(result),
        details = unname(as.matrix(attr(result, "details"))),
        variance = unname(as.matrix(vcov(result)))
    )
}

data("Delta.el", package = "ReGenesees")
s1$domain <- factor(ifelse(s1$x > median(s1$x), "high", "low"))
s2$domain <- factor(ifelse(s2$x > median(s2$x), "high", "low"))
des1 <- e.svydesign(ids = ~id, strata = ~strata, weights = ~w, data = s1)
des2 <- e.svydesign(ids = ~id, strata = ~strata, weights = ~w, data = s2)

element.cases <- list(
    element_full = svyDelta(expression(y.2 - y.1), des1, des2,
                            vartype = c("se", "cv", "var"), conf.int = TRUE),
    element_no_jump = svyDelta(expression(y.2 - y.1), des1, des2,
                               rho.STRAT = "noJump", vartype = "var"),
    element_no_strata = svyDelta(expression(y.2 - y.1), des1, des2,
                                 rho.STRAT = "noStrat", vartype = "var"),
    element_independent = svyDelta(expression(y.2 - y.1), des1, des2,
                                   des.INDEP = TRUE, vartype = "var"),
    element_ratio = svyDelta(
        expression((y.2/x.2 - y.1/x.1) / (y.1/x.1)),
        des1, des2, rho.STRAT = "noJump", vartype = "var"
    ),
    element_by_domain = svyDelta(
        expression(y.2 - y.1), des1, des2, by = ~domain,
        rho.STRAT = "noJump", vartype = c("se", "cv", "var"),
        conf.int = TRUE
    )
)

data("Delta.clus", package = "ReGenesees")
sclus1$domain <- factor(ifelse(sclus1$x > median(sclus1$x), "high", "low"))
sclus2$domain <- factor(ifelse(sclus2$x > median(sclus2$x), "high", "low"))
dclus1 <- e.svydesign(ids = ~id, strata = ~strata, weights = ~w,
                      data = sclus1)
dclus2 <- e.svydesign(ids = ~id, strata = ~strata, weights = ~w,
                      data = sclus2)

cluster.cases <- list(
    cluster_full = svyDelta(expression(y.2 - y.1), dclus1, dclus2,
                            vartype = "var"),
    cluster_no_jump = svyDelta(expression(y.2 - y.1), dclus1, dclus2,
                               rho.STRAT = "noJump", vartype = "var"),
    cluster_independent = svyDelta(expression(y.2 - y.1), dclus1, dclus2,
                                   des.INDEP = TRUE, vartype = "var"),
    cluster_by_domain = svyDelta(
        expression(y.2 - y.1), dclus1, dclus2, by = ~domain,
        rho.STRAT = "noJump", vartype = c("se", "var"), conf.int = TRUE
    )
)

results <- lapply(c(element.cases, cluster.cases), extract)

# Seeded domain benchmark: repeat the cluster example with stable, unique PSU
# identifiers.  This stresses the domain-indexing and PSU aggregation paths.
make.large <- function(sample, repetitions = 2000L, domains = 20L) {
    blocks <- lapply(seq_len(repetitions), function(block) {
        value <- sample
        value$id <- paste(block, value$id, sep = "-")
        value$domain <- factor(sprintf("D%02d", 1L + (block - 1L) %% domains),
                               levels = sprintf("D%02d", seq_len(domains)))
        value
    })
    do.call(rbind, blocks)
}
large1 <- make.large(sclus1)
large2 <- make.large(sclus2)
large.des1 <- e.svydesign(ids = ~id, strata = ~strata, weights = ~w,
                          data = large1)
large.des2 <- e.svydesign(ids = ~id, strata = ~strata, weights = ~w,
                          data = large2)
benchmark.times <- numeric(5L)
benchmark.result <- NULL
for (iteration in seq_along(benchmark.times)) {
    timing <- system.time({
        benchmark.result <- svyDelta(
            expression(y.2 - y.1), large.des1, large.des2,
            by = ~domain, rho.STRAT = "noJump", vartype = "var"
        )
    })
    benchmark.times[iteration] <- timing[["elapsed"]]
}
benchmark <- list(
    median.elapsed = median(benchmark.times),
    result = extract(benchmark.result)
)

payload <- list(results = results, benchmark = benchmark)
if (identical(mode, "write")) {
    saveRDS(payload, baseline.path)
    cat("Wrote", length(results), "svyDelta reference cases; median benchmark",
        benchmark$median.elapsed, "seconds\n")
} else {
    baseline <- readRDS(baseline.path)
    comparisons <- do.call(rbind, lapply(names(results), function(case) {
        current <- results[[case]]
        reference <- baseline$results[[case]]
        numeric.differences <- c(
            abs(as.numeric(current$values) - as.numeric(reference$values)),
            abs(as.numeric(current$variance) - as.numeric(reference$variance))
        )
        data.frame(
            case = case,
            max_abs_numeric_difference = max(numeric.differences, na.rm = TRUE),
            structure_identical = identical(
                current[c("names", "row.names")],
                reference[c("names", "row.names")]
            )
        )
    }))
    benchmark.difference <- max(abs(
        as.numeric(benchmark$result$values) -
        as.numeric(baseline$benchmark$result$values)
    ), na.rm = TRUE)
    speedup <- baseline$benchmark$median.elapsed / benchmark$median.elapsed
    print(comparisons, row.names = FALSE)
    cat("Upstream median seconds:", baseline$benchmark$median.elapsed, "\n")
    cat("Candidate median seconds:", benchmark$median.elapsed, "\n")
    cat("Speedup:", speedup, "\n")
    cat("Benchmark max absolute result difference:", benchmark.difference, "\n")
    if (any(comparisons$max_abs_numeric_difference > 1e-8) ||
        any(!comparisons$structure_identical) ||
        benchmark.difference > 1e-8)
        stop("Candidate svyDelta results are not compatible with upstream.")
}
