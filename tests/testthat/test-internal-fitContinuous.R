## Testing internal functions and testing functions in other packages
## are both considered to be poor form.  So this is something that we
## might re-address before release, but it will be good for now.
context("fitContinuous (internals)")

data(geospiza)
dat <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))

phy <- dat$phy
states <- dat$dat[,"wingL"]
control <- list(niter=10)

set.seed(1)
fit.bm <- fitContinuousQuiet(phy, states, SE=NA, control=control,ncores=1,
                             model="BM")
fit.ou <- fitContinuousQuiet(phy, states, SE=NA, control=control,ncores=1,
                             model="OU")
fit.eb <- fitContinuousQuiet(phy, states, SE=NA, control=control,ncores=1,
                             model="EB")
fit.lam <- fitContinuousQuiet(phy, states, SE=NA, control=control,ncores=1,
                              model="lambda")
fit.del <- fitContinuousQuiet(phy, states, SE=NA, control=control,ncores=1,
                              model="delta")
fit.kap <- fitContinuousQuiet(phy, states, SE=NA, control=control,ncores=1,
                              model="kappa")
fit.wh <- fitContinuousQuiet(phy, states, SE=NA, control=control,ncores=1,
                             model="white")

## And without SE:
fit.bm.noSE <- fitContinuousQuiet(phy, states, control=control,ncores=1,
                                  model="BM")
fit.ou.noSE <- fitContinuousQuiet(phy, states, control=control,ncores=1,
                                  model="OU")
fit.eb.noSE <- fitContinuousQuiet(phy, states, control=control,ncores=1,
                                  model="EB")
fit.lam.noSE <- fitContinuousQuiet(phy, states, control=control,ncores=1,
                                   model="lambda")
fit.del.noSE <- fitContinuousQuiet(phy, states, control=control,ncores=1,
                                   model="delta")
fit.kap.noSE <- fitContinuousQuiet(phy, states, control=control,ncores=1,
                                   model="kappa")
fit.wh.noSE <- fitContinuousQuiet(phy, states, control=control,ncores=1,
                                  model="white")

test_that("Model types are correct", {
  expect_that(model_type(fit.bm), is_identical_to("BM"))
  expect_that(model_type(fit.ou), is_identical_to("OU"))
  expect_that(model_type(fit.eb), is_identical_to("EB"))
  expect_that(model_type(fit.lam), is_identical_to("lambda"))
  expect_that(model_type(fit.del), is_identical_to("delta"))
  expect_that(model_type(fit.kap), is_identical_to("kappa"))
  expect_that(model_type(fit.wh), is_identical_to("white"))

  expect_that(model_type(fit.bm.noSE), is_identical_to("BM"))
  expect_that(model_type(fit.ou.noSE), is_identical_to("OU"))
  expect_that(model_type(fit.eb.noSE), is_identical_to("EB"))
  expect_that(model_type(fit.lam.noSE), is_identical_to("lambda"))
  expect_that(model_type(fit.del.noSE), is_identical_to("delta"))
  expect_that(model_type(fit.kap.noSE), is_identical_to("kappa"))
  expect_that(model_type(fit.wh.noSE), is_identical_to("white"))
})

test_that("Models return their source data", {
  cmp <- list(phy=phy, data=states)
  expect_that(model_data(fit.bm), equals(cmp))
  expect_that(model_data(fit.ou), equals(cmp))
  expect_that(model_data(fit.eb), equals(cmp))
  expect_that(model_data(fit.lam), equals(cmp))
  expect_that(model_data(fit.del), equals(cmp))
  expect_that(model_data(fit.kap), equals(cmp))
  expect_that(model_data(fit.wh), equals(cmp))

  expect_that(model_data(fit.bm.noSE), equals(cmp))
  expect_that(model_data(fit.ou.noSE), equals(cmp))
  expect_that(model_data(fit.eb.noSE), equals(cmp))
  expect_that(model_data(fit.lam.noSE), equals(cmp))
  expect_that(model_data(fit.del.noSE), equals(cmp))
  expect_that(model_data(fit.kap.noSE), equals(cmp))
  expect_that(model_data(fit.wh.noSE), equals(cmp))


})

test_that("Coefficient names are as expected", {
  pars.bm <- c("sigsq", "SE", "z0")
  pars.ou <- c("alpha", "sigsq", "SE", "z0")
  pars.eb <- c("a", "sigsq", "SE", "z0")
  pars.lam <- c("lambda", "sigsq", "SE", "z0")
  pars.kap <- c("kappa", "sigsq", "SE", "z0")
  pars.del <- c("delta", "sigsq", "SE", "z0")
  pars.wh <- c("sigsq", "SE", "z0")

  expect_that(names(coef(fit.bm)), is_identical_to(pars.bm))
  expect_that(names(coef(fit.ou)), is_identical_to(pars.ou))
  expect_that(names(coef(fit.eb)), is_identical_to(pars.eb))
  expect_that(names(coef(fit.lam)), is_identical_to(pars.lam))
  expect_that(names(coef(fit.del)), is_identical_to(pars.del))
  expect_that(names(coef(fit.kap)), is_identical_to(pars.kap))
  expect_that(names(coef(fit.wh)), is_identical_to(pars.wh))

  expect_that(names(coef(fit.bm.noSE)),
              is_identical_to(setdiff(pars.bm, "SE")))
  expect_that(names(coef(fit.ou.noSE)),
              is_identical_to(setdiff(pars.ou, "SE")))
  expect_that(names(coef(fit.eb.noSE)),
              is_identical_to(setdiff(pars.eb, "SE")))
  expect_that(names(coef(fit.lam.noSE)),
              is_identical_to(setdiff(pars.lam, "SE")))
  expect_that(names(coef(fit.del.noSE)),
              is_identical_to(setdiff(pars.del, "SE")))
  expect_that(names(coef(fit.kap.noSE)),
              is_identical_to(setdiff(pars.kap, "SE")))
  expect_that(names(coef(fit.wh.noSE)),
              is_identical_to(setdiff(pars.wh, "SE")))
})

test_that("Processed coefficient names are as expected", {
  expect_that(names(model_pars(fit.bm)),
              is_identical_to(parnames.bm()))
  expect_that(names(model_pars(fit.ou)),
              is_identical_to(parnames.ou()))
  expect_that(names(model_pars(fit.eb)),
              is_identical_to(parnames.eb()))
  expect_that(names(model_pars(fit.lam)),
              is_identical_to(parnames.lambda()))
  expect_that(names(model_pars(fit.del)),
              is_identical_to(parnames.delta()))
  expect_that(names(model_pars(fit.kap)),
              is_identical_to(parnames.kappa()))
  expect_that(names(model_pars(fit.wh)),
              is_identical_to(parnames.white()))


  expect_that(model_pars(fit.bm)$SE, is_greater_than(0))
  expect_that(model_pars(fit.ou)$SE, is_greater_than(0))
  expect_that(model_pars(fit.eb)$SE, is_greater_than(0))
  expect_that(model_pars(fit.lam)$SE, is_greater_than(0))
  expect_that(model_pars(fit.del)$SE, is_greater_than(0))
  expect_that(model_pars(fit.kap)$SE, is_greater_than(0))
  expect_that(model_pars(fit.wh)$SE, is_greater_than(0))

  expect_that(names(model_pars(fit.bm.noSE)),
              is_identical_to(parnames.bm()))
  expect_that(names(model_pars(fit.ou.noSE)),
              is_identical_to(parnames.ou()))
  expect_that(names(model_pars(fit.eb.noSE)),
              is_identical_to(parnames.eb()))
  expect_that(names(model_pars(fit.lam.noSE)),
              is_identical_to(parnames.lambda()))
  expect_that(names(model_pars(fit.del.noSE)),
              is_identical_to(parnames.delta()))
  expect_that(names(model_pars(fit.kap.noSE)),
              is_identical_to(parnames.kappa()))
  expect_that(names(model_pars(fit.wh.noSE)),
              is_identical_to(parnames.white()))


  expect_that(model_pars(fit.bm.noSE)$SE, equals(0))
  expect_that(model_pars(fit.ou.noSE)$SE, equals(0))
  expect_that(model_pars(fit.eb.noSE)$SE, equals(0))
  expect_that(model_pars(fit.lam.noSE)$SE, equals(0))
  expect_that(model_pars(fit.del.noSE)$SE, equals(0))
  expect_that(model_pars(fit.kap.noSE)$SE, equals(0))
  expect_that(model_pars(fit.wh.noSE)$SE, equals(0))

})

test_that("Overall processed object looks legit", {
  obj.bm <- model_info(fit.bm)
  obj.ou <- model_info(fit.ou)
  obj.eb <- model_info(fit.eb)
  obj.lam <- model_info(fit.lam)
  obj.del <- model_info(fit.del)
  obj.kap <- model_info(fit.kap)
  obj.wh <- model_info(fit.wh)

  obj.bm.noSE <- model_info(fit.bm.noSE)
  obj.ou.noSE <- model_info(fit.ou.noSE)
  obj.eb.noSE <- model_info(fit.eb.noSE)
  obj.lam.noSE <- model_info(fit.lam.noSE)
  obj.kap.noSE <- model_info(fit.kap.noSE)
  obj.del.noSE <- model_info(fit.del.noSE)
  obj.wh.noSE <- model_info(fit.wh.noSE)

  obj.names <- c("data", "pars", "type")
  expect_that(names(obj.bm), is_identical_to(obj.names))
  expect_that(names(obj.ou), is_identical_to(obj.names))
  expect_that(names(obj.eb), is_identical_to(obj.names))
  expect_that(names(obj.lam), is_identical_to(obj.names))
  expect_that(names(obj.del), is_identical_to(obj.names))
  expect_that(names(obj.kap), is_identical_to(obj.names))
  expect_that(names(obj.wh), is_identical_to(obj.names))

  expect_that(names(obj.bm.noSE), is_identical_to(obj.names))
  expect_that(names(obj.ou.noSE), is_identical_to(obj.names))
  expect_that(names(obj.eb.noSE), is_identical_to(obj.names))
  expect_that(names(obj.lam.noSE), is_identical_to(obj.names))
  expect_that(names(obj.del.noSE), is_identical_to(obj.names))
  expect_that(names(obj.kap.noSE), is_identical_to(obj.names))
  expect_that(names(obj.wh.noSE), is_identical_to(obj.names))
})

test_that("Nonsense cases should fail", {
  expect_that(model_info(NULL), throws_error())
  expect_that(model_info(unclass(fit.bm)), throws_error())
})

test_that("Variable length SE throws error", {
  set.seed(1)
  se <- runif(length(states))
  names(se) <- names(states)
  fit.bm <- fitContinuous(phy, states, SE=se, control=control,ncores=1,
                          model="BM")
  expect_that(model_pars(fit.bm), throws_error())

  se[1] <- NA
  fit.bm <- fitContinuousQuiet(phy, states, SE=se, control=control,ncores=1,
                               model="BM")
  expect_that(model_pars(fit.bm), throws_error())
})
