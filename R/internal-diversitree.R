## See internal-fitContinuous.R for general information about what is
## going on here.

## TODO: check issue #9.

model_type.fit.mle <- function(fit, ...) {
  toupper(sub("^fit\\.mle\\.", "", class(fit)[[1]]))
}

model_data.fit.mle <- function(fit, ...) {
  cache <- get.cache(get.likelihood(fit))
  list(phy=cache$info$phy, data=cache$states)
}

model_pars.fit.mle <- function(fit, ...) {
  loadNamespace("diversitree")
  model <- model_type(fit)
  pars <- as.list(coef(fit, full=TRUE))
  pars <- diversitree.to.arbutus.recast(pars, model)
  pars <- diversitree.to.arbutus.se(pars, get.likelihood(fit))
  pars
}

#' @method model_info fit.mle
#' @export
model_info.fit.mle <- function(fit, ...) {
  m <- list(data=model_data(fit),
            pars=model_pars(fit),
            type=model_type(fit))
  class(m) <- "fitC"
  m
}

## Similar, but for mcmcsamples:
model_type.mcmcsamples <- function(fit, ...) {
  lik <- get.likelihood(fit)
  if (is.constrained(lik)) {
    toupper(class(get.likelihood(lik))[[1]])
  } else {
    toupper(class(lik)[[1]])
  }
}

model_data.mcmcsamples <- function(fit, ...) {
  cache <- get.cache(get.likelihood(fit))
  list(phy=cache$info$phy, data=cache$states)
}

model_pars.mcmcsamples <- function(fit, burnin=NA, thin=NA, sample=NA, ...) {
  model <- model_type(fit)
  lik <- get.likelihood(fit)
  pars <- as.data.frame(coef(fit, burnin=burnin, thin=thin,
                             sample=sample, full=TRUE, lik=lik))
  pars <- diversitree.to.arbutus.recast(pars, model)
  pars <- diversitree.to.arbutus.se(pars, lik)
  pars
}

#' @method model_info mcmcsamples
#' @export
model_info.mcmcsamples <- function(fit, burnin=NA, thin=NA, sample=NA, ...) {
  m <- list(data=model_data(fit),
            pars=model_pars(fit, burnin, thin, sample),
            type=model_type(fit))
  class(m) <- c("multiFitC", "fitC")
  m
}

diversitree.to.arbutus.recast <- function(pars, model) {
  if (!(model %in% c("BM", "OU", "EB")))
    stop(sprintf("Cannot extract from models of model %s", model))

  tr <- c(sigsq="s2", z0="theta")
  idx <- match(names(pars), tr)
  i <- !is.na(idx)
  names(pars)[i] <- names(tr)[idx][i]

  if (model == "OU") {
    # NOTE: This is a workaround for the with.optimum=FALSE OU model.
    # Ideally we'd get the root state, but that requires one extra
    # likelihood evaluation because we don't save it; it's not even
    # always available anyway under the vcv method (but that should
    # change)
    #   attr(get.likelihood(fit)(coef(fit), intermediates=TRUE),
    #        "vals")[[1]]
    # but also assumes a bunch of internal things about how the model
    # works.
    if (length(pars) == 2) {
      pars <- pars[c("alpha", "sigsq")]
      pars$z0 <- NA
    } else {
      pars <- pars[c("alpha", "sigsq", "z0")]
    }
  } else if (model == "EB") {
    pars <- pars[c("a", "sigsq")]
    pars$z0 <- NA
  } else if (model == "BM") {
    ## At this point we don't have access to the root state.  Given
    ## that different root treatments are possible in diversitree,
    ## this is not actually very straightforward, either.
    pars$z0 <- NA
  }

  pars
}

diversitree.to.arbutus.se <- function(pars, lik) {
  ## Pull the standard error out.  Also ugly.
  cache <- get.cache(lik)
  sd <- cache$states.sd
  if (is.null(sd)) {# pgls
    sd <- 0
  }
  if (length(unique(sd)) == 1) {
    pars$SE <- sd[[1]]
  } else {
    stop("Variable length states.sd -- cannot deal with this (yet)")
  }
  pars
}

## Two internal diversitree functions copied here
is.constrained <- function(x) {
  inherits(x, "constrained")
}
get.cache <- function(x) {
  if (inherits(x, "big.brother") || is.constrained(x)) {
    get.cache(attr(x, "func"))
  } else {
    environment(x)$cache
  }
}
get.likelihood <- function(...) {
  loadNamespace("diversitree")
  diversitree::get.likelihood(...)
}
