## For now, only BM models are implemented in diversitree.
model_type.fit.mle.pgls <- function(fit, ...) {
  "BM"
}

model_type.mcmcsamples.pgls <- function(fit, ...) {
  "BM"
}

model_data.fit.mle.pgls <- function(fit, ...) {
  cache <- get.cache(get.likelihood(fit))
  list(phy=cache$info$phy, data=drop(residuals(fit)))
}

## NOTE: this passes the coef.mcmcsamples parameters (burnin, thin,
## sample) along to residuals.mcmcsamples().
model_data.mcmcsamples.pgls <- function(fit, ...) {
  cache <- get.cache(get.likelihood(fit))
  list(phy=cache$info$phy, data=residuals(fit, ...))
}

model_pars.fit.mle.pgls <- function(fit, ...) {
  ## TODO: Will only work with BM models; update once diversitree can
  ## do more than this.
  if (model_type(fit) != "BM")
    stop("Hard coded to deal with BM models only")
  s2 <- coef(fit, full=TRUE)[["s2"]]
  as.list(c(sigsq=s2, SE=0, z0=NA))
}

model_pars.mcmcsamples.pgls <- function(fit, ...) {
  ## TODO: Will only work with BM models; update once diversitree can
  ## do more than this.
  if (model_type(fit) != "BM")
    stop("Hard coded to deal with BM models only")
  ## Pass through here to deal with options to model_pars.mcmcsamples
  ## that take care of chain thinning, etc.
  ans <- NextMethod("model_pars")
  ans[c("sigsq", "SE", "z0")]
}
