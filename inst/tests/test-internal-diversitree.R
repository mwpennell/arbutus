## Testing internal functions and testing functions in other packages
## are both considered to be poor form.  So this is something that we
## might re-address before release, but it will be good for now.
source("helper-arbutus.R")

context("diversitree (internals)")

data(geospiza)
dat <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))

phy <- dat$phy
states <- dat$dat[,"wingL"]

lik.bm <- make.bm(phy, states)
lik.ou <- make.ou(phy, states)
lik.eb <- make.eb(phy, states)
lik.ou.opt <- make.ou(phy, states, with.optimum=TRUE,
                      control=list(method="pruning"))

fit.bm <- find.mle(lik.bm, .1)
fit.ou <- find.mle(lik.ou, c(.1, .1))
fit.eb <- find.mle(lik.eb, c(.1, 0))
## TODO: Potential problem when dealing with renamed names; won't
## translate properly.  Potentially, find.mle should have dealt with
## this for us.
fit.ou.opt <- find.mle(lik.ou.opt, c(coef(fit.ou), theta=mean(states)))

set.seed(1)
samples.bm <- mcmc(lik.bm, coef(fit.bm), 100, w=0.1, print.every=0)
set.seed(1)
samples.ou <- mcmc(lik.ou, coef(fit.ou), 30,
                   w=c(0.1, 10),
                   lower=c(0, 0),
                   upper=c(Inf, 100),
                   print.every=0)
set.seed(1)
samples.eb <- mcmc(lik.eb, coef(fit.eb), 30,
                   w=c(.2, 10), lower=c(0, -Inf), upper=c(Inf, 0),
                   print.every=0)
set.seed(1)
samples.ou.opt <- mcmc(lik.ou.opt, coef(fit.ou.opt), 30,
                   w=c(0.1, 10, 10),
                   lower=c(0, 0, -10),
                   upper=c(Inf, 100, 10),
                   print.every=0)

## Testing internal code; need some extra imports:
model_type <- arbutus:::model_type
model_data <- arbutus:::model_data
model_pars <- arbutus:::model_pars
model_info <- arbutus:::model_info
model_type.fit.mle <- arbutus:::model_type.fit.mle
model_data.fit.mle <- arbutus:::model_data.fit.mle
model_pars.fit.mle <- arbutus:::model_pars.fit.mle
model_info.fit.mle <- arbutus:::model_info.fit.mle
model_type.mcmcsamples <- arbutus:::model_type.mcmcsamples
model_data.mcmcsamples <- arbutus:::model_data.mcmcsamples
model_pars.mcmcsamples <- arbutus:::model_pars.mcmcsamples
model_info.mcmcsamples <- arbutus:::model_info.mcmcsamples

test_that("Model types are correct", {
  expect_that(model_type(fit.bm), is_identical_to("BM"))
  expect_that(model_type(fit.ou), is_identical_to("OU"))
  expect_that(model_type(fit.eb), is_identical_to("EB"))
  expect_that(model_type(fit.ou.opt), is_identical_to("OU"))
  expect_that(model_type(samples.bm), is_identical_to("BM"))
  expect_that(model_type(samples.ou), is_identical_to("OU"))
  expect_that(model_type(samples.eb), is_identical_to("EB"))
  expect_that(model_type(samples.ou.opt), is_identical_to("OU"))
})

test_that("Models return their source data", {
  cmp <- list(phy=phy, data=states[phy$tip.label])
  expect_that(model_data(fit.bm), equals(cmp))
  expect_that(model_data(fit.ou), equals(cmp))
  expect_that(model_data(fit.eb), equals(cmp))
  expect_that(model_data(fit.ou.opt), equals(cmp))

  expect_that(model_data(samples.bm), equals(cmp))
  expect_that(model_data(samples.ou), equals(cmp))
  expect_that(model_data(samples.eb), equals(cmp))
  expect_that(model_data(samples.ou.opt), equals(cmp))
})

test_that("Processed coefficient names are as expected", {
  expect_that(names(model_pars(fit.bm)),
              is_identical_to(arbutus:::parnames.bm()))
  expect_that(names(model_pars(fit.ou)),
              is_identical_to(arbutus:::parnames.ou()))
  expect_that(names(model_pars(fit.eb)),
              is_identical_to(arbutus:::parnames.eb()))
  expect_that(names(model_pars(fit.ou.opt)),
              is_identical_to(arbutus:::parnames.ou()))

  expect_that(model_pars(fit.bm)$SE, equals(0))
  expect_that(model_pars(fit.ou)$SE, equals(0))
  expect_that(model_pars(fit.eb)$SE, equals(0))

  ## MCMC:
  expect_that(names(model_pars(samples.bm)),
              is_identical_to(arbutus:::parnames.bm()))
  expect_that(names(model_pars(samples.ou)),
              is_identical_to(arbutus:::parnames.ou()))
  expect_that(names(model_pars(samples.eb)),
              is_identical_to(arbutus:::parnames.eb()))
  expect_that(names(model_pars(samples.ou.opt)),
              is_identical_to(arbutus:::parnames.ou()))

  expect_that(model_pars(samples.bm)$SE,
              equals(rep(0, nrow(samples.bm))))
  expect_that(model_pars(samples.ou)$SE,
              equals(rep(0, nrow(samples.ou))))
  expect_that(model_pars(samples.eb)$SE,
              equals(rep(0, nrow(samples.eb))))
})

test_that("MCMC extraction options are passed from model_pars", {
  n <- nrow(samples.bm)
  expect_that(nrow(model_pars(samples.bm)), equals(n))

  nb <- 10
  expect_that(nrow(model_pars(samples.bm, burnin=nb)),
              equals(nrow(coef(samples.bm, burnin=nb))))
  pb <- 0.14234
  expect_that(nrow(model_pars(samples.bm, burnin=pb)),
              equals(nrow(coef(samples.bm, burnin=pb))))
  thin <- 7
  expect_that(nrow(model_pars(samples.bm, thin=thin)),
              equals(nrow(coef(samples.bm, thin=thin))))

  sample <- 10
  set.seed(1)
  p1 <- model_pars(samples.bm, burnin=nb, thin=thin, sample=sample)
  set.seed(1)
  p2 <- coef(samples.bm, burnin=nb, thin=thin, sample=sample)
  expect_that(nrow(p1), equals(nrow(p2)))
  expect_that(p1$sigsq, is_identical_to(as.numeric(p2[,"s2"])))
})

test_that("Overall processed object looks legit", {
  obj.bm <- model_info(fit.bm)
  obj.ou <- model_info(fit.ou)
  obj.eb <- model_info(fit.eb)

  obj.names <- c("data", "pars", "type")
  expect_that(names(obj.bm), is_identical_to(obj.names))
  expect_that(names(obj.ou), is_identical_to(obj.names))
  expect_that(names(obj.eb), is_identical_to(obj.names))

  obj.bm <- model_info(samples.bm)
  obj.ou <- model_info(samples.ou)
  obj.eb <- model_info(samples.eb)

  expect_that(names(obj.bm), is_identical_to(obj.names))
  expect_that(names(obj.ou), is_identical_to(obj.names))
  expect_that(names(obj.eb), is_identical_to(obj.names))
})

test_that("MCMC extraction options are passed from model_info", {
  n <- nrow(samples.bm)
  
  expect_that(nrow(model_info(samples.bm)$pars), equals(n))

  nb <- 10
  pb <- 0.14234
  thin <- 7
  sample <- 10

  expect_that(nrow(model_info(samples.bm, burnin=nb)$pars),
              equals(nrow(coef(samples.bm, burnin=nb))))
  expect_that(nrow(model_info(samples.bm, burnin=pb)$pars),
              equals(nrow(coef(samples.bm, burnin=pb))))
  expect_that(nrow(model_info(samples.bm, thin=thin)$pars),
              equals(nrow(coef(samples.bm, thin=thin))))
  set.seed(1)
  p1 <- model_info(samples.bm, burnin=nb, thin=thin, sample=sample)$pars
  set.seed(1)
  p2 <- coef(samples.bm, burnin=nb, thin=thin, sample=sample)
  expect_that(nrow(p1), equals(nrow(p2)))
  expect_that(p1$sigsq, is_identical_to(as.numeric(p2[,"s2"])))
})

## Constrained functions
test_that("Constrained functions work", {
  lik.ou.c <- constrain(lik.ou.opt, theta ~ 0)
  set.seed(1)
  fit.ou.c <- find.mle(lik.ou.c, c(.1, .1))
  samples.ou.c <- mcmc(lik.ou.c, coef(fit.ou.c), 30,
                       w=c(0.1, 10),
                       lower=c(0, 0),
                       upper=c(Inf, 100),
                       print.every=0)

  ## Can get the type out:
  expect_that(model_type(fit.ou.c), is_identical_to("OU"))
  expect_that(model_type(samples.ou.c),
              is_identical_to("OU"))

  ## Underlying data:
  cmp <- list(phy=phy, data=states[phy$tip.label])
  expect_that(model_data(fit.ou.c), equals(cmp))
  expect_that(model_data(samples.ou.c), equals(cmp))

  ## Parameters -- these need to include all parameters:
  p <- model_pars(fit.ou.c)
  expect_that(names(p), is_identical_to(arbutus:::parnames.ou()))
  expect_that(p$SE, equals(0))
  expect_that(p$z0, equals(0)) # z0 is theta, translated

  p <- model_pars(samples.ou.c)
  expect_that(names(p), is_identical_to(arbutus:::parnames.ou()))
  expect_that(p$SE, equals(rep(0, nrow(samples.ou))))
  expect_that(p$z0, equals(rep(0, nrow(samples.ou))))
})
