#internal-OUwie
#Functions for using OUwie objects

model_type.OUwie <- function(fit, ...){
  ret <- fit$model
  if (is.null(ret) || !is.character(ret) || length(ret) != 1)
    stop("Failed to extract model type from OUwie fit")
  ret
}

model_data.OUwie <- function(fit, ...) {
  temp <- fit$data[[2]]
  names(temp) <- rownames(fit$data)
  ret <- list(phy=fit$phy, data=temp)
  ret
}

model_pars.OUwie <- function(fit, ...) {
  list(sigsq = fit$solution["sigma.sq",], SE = fit$solution.se)
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