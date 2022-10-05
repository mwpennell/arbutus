context("OUwie internals")

data(geospiza)
dat <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))

phy <- dat$phy
phy$node.label <- c(1,1,1,2,2,2,2,2,2,2,2,2)

control <- list(niter=10)

set.seed(1)

df <- data.frame(Genus_species = phy$tip.label, Reg = c(2,2,2,2,2,2,2,2,2,2,1,1,1))
sim_dat <- OUwie::OUwie.sim(phy, df, alpha = c(1e-10,1e-10), sigma.sq = c(0.45,0.90), theta0 = 1.0, theta = c(0,0))
test_df <- sim_dat[[3]]
names(test_df) <- sim_dat[[1]]

fit.bms <- OUwie::OUwie(phy, sim_dat, model = "BMS")

test_that("Model types are correct", {
  expect_that(model_type(fit.bms), is_identical_to("BMS"))
})

test_that("Models return source data",{
  cmp <- list(phy = phy, data = test_df)
  expect_that(model_data(fit.bms), equals(cmp))
})

test_that("Coefficient names are expected", {
  pars.bms <- c("sigsq", "SE")
  expect_that(names(model_pars(fit.bms)), is_identical_to(parnames.bms()))
  expect_that(names(model_pars(fit.bms)), is_identical_to(pars.bms))
  
})


test_that("Overall processed object looks legit", {
  obj.bms <- model_info(fit.bms)
  obj.names <- c("data", "pars", "type")
})

test_that("Nonsense cases should fail", {
  expect_that(model_info(NULL), throws_error())
  expect_that(model_info(unclass(fit.bms)), throws_error())
})
