# Compare public calibration results from the upstream and CAE-INS trees.
#
# Usage from either package tree:
#   Rscript verify-upstream-compatibility.R write baseline.rds
#   Rscript verify-upstream-compatibility.R compare baseline.rds

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L || !arguments[1L] %in% c("write", "compare"))
    stop("Usage: verify-upstream-compatibility.R <write|compare> <baseline.rds>")

mode <- arguments[1L]
baseline.path <- arguments[2L]
source.root <- Sys.getenv("REGENESEES_SOURCE", unset = ".")
devtools::load_all(source.root, quiet = TRUE)
data("data.examples", package = "ReGenesees")

design <- e.svydesign(
    data = example, ids = ~towcod + famcod,
    strata = ~SUPERSTRATUM, weights = ~weight
)

cases <- list(
    global_count_aggregate = list(
        df.population = pop01, calmodel = ~1, calfun = "logit",
        bounds = bounds, aggregate.stage = 2
    ),
    margins_aggregate = list(
        df.population = pop02, calmodel = ~sex + marstat - 1,
        calfun = "logit", bounds = bounds, aggregate.stage = 2
    ),
    joint_global = list(
        df.population = pop03, calmodel = ~marstat:sex - 1,
        calfun = "logit", bounds = bounds
    ),
    joint_partitioned = list(
        df.population = pop03p, calmodel = ~marstat - 1,
        partition = ~sex, calfun = "logit", bounds = bounds
    ),
    continuous_global_aggregate = list(
        df.population = pop04, calmodel = ~(x1 + x2 + x3):regcod - 1,
        calfun = "logit", bounds = bounds, aggregate.stage = 2
    ),
    continuous_partitioned_aggregate = list(
        df.population = pop04p, calmodel = ~x1 + x2 + x3 - 1,
        partition = ~regcod, calfun = "logit", bounds = bounds,
        aggregate.stage = 2
    ),
    high_dimensional_global = list(
        df.population = pop05,
        calmodel = ~(age5c + x1):sex:marstat - 1,
        calfun = "logit", bounds = bounds
    ),
    high_dimensional_partitioned = list(
        df.population = pop05p, calmodel = ~age5c + x1 - 1,
        partition = ~sex:marstat, calfun = "logit", bounds = bounds
    ),
    unbounded_linear = list(
        df.population = pop03, calmodel = ~marstat:sex - 1,
        calfun = "linear"
    ),
    raking = list(
        df.population = pop03, calmodel = ~marstat:sex - 1,
        calfun = "raking"
    )
)

results <- lapply(cases, function(arguments) {
    calibrated <- suppressWarnings(do.call(
        e.calibrate, c(list(design = design), arguments)
    ))
    list(
        weights = unname(weights(calibrated)),
        return.code = unname(attr(calibrated, "ecal.status")$return.code)
    )
})

if (identical(mode, "write")) {
    saveRDS(results, baseline.path)
    cat("Wrote", length(results), "upstream calibration cases to", baseline.path, "\n")
} else {
    baseline <- readRDS(baseline.path)
    if (!identical(names(results), names(baseline)))
        stop("Upstream and candidate case names differ.")
    comparison <- do.call(rbind, lapply(names(results), function(case) {
        current <- results[[case]]
        reference <- baseline[[case]]
        data.frame(
            case = case,
            max_abs_weight_difference = max(abs(
                current$weights - reference$weights
            )),
            max_relative_weight_difference = max(abs(
                current$weights - reference$weights
            ) / (1 + abs(reference$weights))),
            return_code_identical = identical(
                current$return.code, reference$return.code
            )
        )
    }))
    print(comparison, row.names = FALSE)
    if (any(comparison$max_relative_weight_difference > 1e-8) ||
        any(!comparison$return_code_identical))
        stop("Candidate calibration results are not compatible with upstream.")
}
