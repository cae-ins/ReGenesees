`.rg.option` <- function(name, default) getOption(name, default)

`.rg.model.matrix` <- function(object, data, contrasts.arg = NULL, ...)
###########################################################################
# Build model matrices with package-local contrast defaults.  ReGenesees   #
# historically changed options("contrasts") at namespace load time.  This  #
# wrapper preserves its factor encoding without changing the user's global  #
# R session.                                                               #
###########################################################################
{
    if (is.null(contrasts.arg)) {
        factor.columns <- vapply(data, is.factor, logical(1L))
        if (any(factor.columns)) {
            setup <- .rg.option(
                "RG.contrasts",
                c(unordered = "contr.treatment", ordered = "contr.treatment")
            )
            factor.data <- data[factor.columns]
            contrasts.arg <- as.list(vapply(
                factor.data,
                function(x) if (is.ordered(x)) setup[["ordered"]] else setup[["unordered"]],
                character(1L)
            ))
        }
    }
    stats::model.matrix(object, data = data, contrasts.arg = contrasts.arg, ...)
}
