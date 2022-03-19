#internal-OUwie
#Functions for using OUwie objects

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

#' @method model_info OUwie
#' @export
model_info.OUwie <- function(fit, ...) {
  m <- list(data=model_data(fit),
            pars=model_pars(fit),
            type=model_type(fit))
  class(m) <- "fitOU"
  m
}