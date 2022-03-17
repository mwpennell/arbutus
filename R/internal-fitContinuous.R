## Internal functions for pulling apart fitContinuous objects.  These
## need to be kept in close alignment with geiger, because they use
## undocumented (and therefore potentially unstable) features of
## geiger.

## Some of these could be merged into geiger, but we're going to need
## similar features for every model type and people won't do this for
## us.  So instead, we'll write some unit tests that check that things
## are behaving as expected.  Then we should at least get fairly early
## warning when things break.  Later on, we can look at rolling some
## interface into the different packages if people will let us.


model_type.gfit <- function(fit, ...) {
  ret <- attr(fit$lik, "model")
  if (is.null(ret) || !is.character(ret) || length(ret) != 1)
    stop("Failed to extract model type from geiger fit")
  ret
}

model_data.gfit <- function(fit, ...) {
  e <- environment(fit$lik)
  list(phy=get("phy", e), data=get("dat", e))
}

model_pars.gfit <- function(fit, ...) {
  model <- model_type(fit)
  pars <- as.list(coef(fit))
  if (!("SE" %in% names(pars)))
    pars$SE <- unique(attr(fit$lik, "cache")$SE)

  # The second condition here works around some odd behaviour with
  # geiger where a vector of SE values is given, but some are also
  # being estimated.  This probably affects basically nobody.
  if (length(pars$SE) > 1 || length(unique(attr(fit$lik, "cache")$SE)) > 1)
    stop("Variable SE not yet implemented")

  pars[c(setdiff(names(pars), "SE"), "SE")]
}

#' @method model_info gfit
#' @export
model_info.gfit <- function(fit, ...) {
  m <- list(data=model_data(fit),
            pars=model_pars(fit),
            type=model_type(fit))
  class(m) <- "fitC"
  m
}


model_type.OUwie <- function(fit, ...){
  ret <- fit$model
  if (is.null(ret) || !is.character(ret) || length(ret) != 1)
    stop("Failed to extract model type from OUwie fit")
  ret
}

model_data.OUwie <- function(fit, ...) {
  list(phy=fit$phy, data=fit$data)
}

model_pars.OUwie <- function(fit, ...) {
  list(sigsq1 = fit$solution[2], sigsq2 = fit$solution[4], SE = fit$solution.se)
}

#' @method model_info gfit
#' @export
model_info.OUwie <- function(fit, ...) {
  m <- list(data=model_data(fit),
            pars=model_pars(fit),
            type=model_type(fit))
  class(m) <- "fitOU"
  m
}