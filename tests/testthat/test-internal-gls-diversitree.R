library(diversitree)
library(nlme)
library(testthat)

context("GLS (diversitree) (internals)")

model_type <- arbutus:::model_type
model_data <- arbutus:::model_data
model_pars <- arbutus:::model_pars
model_info <- arbutus:::model_info

model_type.fit.mle <- arbutus:::model_type.fit.mle
model_data.fit.mle <- arbutus:::model_data.fit.mle
model_pars.fit.mle <- arbutus:::model_pars.fit.mle
model_info.fit.mle <- arbutus:::model_info.fit.mle

model_type.fit.mle.pgls <- arbutus:::model_type.fit.mle.pgls
model_data.fit.mle.pgls <- arbutus:::model_data.fit.mle.pgls
model_pars.fit.mle.pgls <- arbutus:::model_pars.fit.mle.pgls

model_type.mcmcsamples <- arbutus:::model_type.mcmcsamples
model_data.mcmcsamples <- arbutus:::model_data.mcmcsamples
model_pars.mcmcsamples <- arbutus:::model_pars.mcmcsamples

model_type.mcmcsamples.pgls <- arbutus:::model_type.mcmcsamples.pgls
model_data.mcmcsamples.pgls <- arbutus:::model_data.mcmcsamples.pgls
model_pars.mcmcsamples.pgls <- arbutus:::model_pars.mcmcsamples.pgls

set.seed(1)
phy <- tree.bd(pars=c(1,0), max.taxa=100)
tx <- sim.character(phy, 1)
ty <- sim.character(phy, 1) + 3
data <- data.frame(x=tx, y=ty, row.names=names(tx))

fit.gls.bm.ml   <- gls(y ~ x, data, corBrownian(1, phy), method="ML")
fit.gls.bm.reml <- gls(y ~ x, data, corBrownian(1, phy), method="REML")

lik.pgls.bm.vcv <- make.pgls(phy, y ~ x, data,
                             control=list(method="vcv"))
lik.pgls.bm.con <- make.pgls(phy, y ~ x, data,
                             control=list(method="contrasts"))

s2.ml <- arbutus:::estimate.sigma2.gls(fit.gls.bm.ml)
p.ml <- c(coef(fit.gls.bm.ml), s2=s2.ml)

fit.pgls.bm.vcv <- find.mle(lik.pgls.bm.vcv, c(0, 0, 1))
fit.pgls.bm.con <- find.mle(lik.pgls.bm.con, c(0, 0, 1))

## Draw some mcmc samples, straight from the likelihood.
set.seed(1)
samples <- mcmc(lik.pgls.bm.con, p.ml, 100, w=1, print.every=0)

test_that("The fits are sane", {
  expect_that(fit.pgls.bm.vcv, equals(fit.pgls.bm.con))
})

test_that("Model types are correct", {
  expect_that(model_type(fit.pgls.bm.vcv), is_identical_to("BM"))
  expect_that(model_type(fit.pgls.bm.con), is_identical_to("BM"))
  expect_that(model_type(samples),         is_identical_to("BM"))
})

test_that("Model parameters are correct", {
  cmp <- as.list(c(sigsq=coef(fit.pgls.bm.vcv)[["s2"]], SE=0, z0=NA))
  expect_that(model_pars(fit.pgls.bm.vcv), equals(cmp))
  expect_that(model_pars(fit.pgls.bm.con), equals(cmp))

  cmp <- data.frame(sigsq=coef(samples)[,"s2"], SE=0, z0=NA)
  expect_that(model_pars(samples), equals(cmp))

  nb <- 30
  cmp <- data.frame(sigsq=coef(samples, burnin=nb)[,"s2"], SE=0, z0=NA)
  expect_that(model_pars(samples, burnin=nb), equals(cmp))

  sample <- 10
  thin <- 2
  set.seed(1)
  cmp <- data.frame(sigsq=coef(samples, burnin=nb,
                      thin=thin, sample=sample)[,"s2"], SE=0, z0=NA)
  set.seed(1)
  expect_that(model_pars(samples, burnin=nb, thin=thin, sample=sample),
              equals(cmp))
})

# Testing this after pars, because this is tricky for mcmcsamples.pgls
test_that("Models return their source data", {
  drop.attributes <- function(x) {
    attributes(x) <- attributes(x)["names"]
    x
  }
  cmp.bm <- list(phy=phy, data=drop.attributes(resid(fit.gls.bm.ml)))
  expect_that(model_data(fit.pgls.bm.vcv), equals(cmp.bm, tolerance=1e-5))
  expect_that(model_data(fit.pgls.bm.con), equals(cmp.bm, tolerance=1e-5))

  ## Update the fits to use the same values as gls() and we'll get
  ## better agreement.
  fit.pgls.bm.vcv$par <- p.ml
  fit.pgls.bm.con$par <- p.ml
  expect_that(model_data(fit.pgls.bm.vcv), equals(cmp.bm, tolerance=1e-13))
  expect_that(model_data(fit.pgls.bm.con), equals(cmp.bm, tolerance=1e-13))

  cmp <- list(phy=phy, data=resid(samples))
  expect_that(model_data(samples), equals(cmp, tolerance=1e-13))

  ## Pass along parameters:
  burnin <- 10
  thin <- 2
  sample <- 10
  set.seed(1)
  cmp <- list(phy=phy,
              data=resid(samples, burnin=burnin, thin=thin, sample=sample))
  set.seed(1)
  expect_that(model_data(samples, burnin=burnin, thin=thin, sample=sample),
              equals(cmp, tolerance=1e-13))  
})

test_that("Overall processed object looks legit", {
  obj.bm.vcv <- model_info(fit.pgls.bm.vcv)
  obj.bm.con <- model_info(fit.pgls.bm.con)

  obj.names <- c("data", "pars", "type")
  expect_that(names(obj.bm.vcv), is_identical_to(obj.names))
  expect_that(names(obj.bm.con), is_identical_to(obj.names))
})
