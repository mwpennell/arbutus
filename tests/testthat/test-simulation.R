# Test the BM simulation code
context("Simulation")

test_that("Character simulation gives correct distribution", {
  data(geospiza)
  phy <- geospiza$phy

  set.seed(1)

  ## Using geiger as a benchmark
  dat.g <- geiger::sim.char(phy, par=1, nsim=100, model="BM")[,,]

  ## And diversitree:
  dat.d <- replicate(100, diversitree::sim.character(phy, 1, x0=1))

  ## And the new arbutus functions
  dat.a <- sim.char.std.bm(phy, 100, x0=1)

  ## Then, do some fits:
  do.fit <- function(x) {
    lik <- make.bm(phy, x, control=list(method="pruning"))
    unname(coef(find.mle(lik, 1, "optimize",
                         control=list(interval=c(0.1, 10)))))
  }

  ## These are fairly slow, which is annoying.
  est.g <- apply(dat.g, 2, do.fit)
  est.d <- apply(dat.d, 2, do.fit)
  est.a <- apply(dat.a, 2, do.fit)

  expect_that(suppressWarnings(ks.test(est.a, est.g))$p.value,
              is_greater_than(0.1))
  expect_that(suppressWarnings(ks.test(est.a, est.d))$p.value,
              is_greater_than(0.1))
})

test_that("Multifurcations are not handled", {
  data(geospiza)
  phy <- geospiza$phy

  set.seed(1)

  ## Now, test on a tree with multifurcations:
  phy.m <- di2multi(phy, 0.05)
  expect_that(is.binary.tree(phy.m), is_false())

  dat.g <- geiger::sim.char(phy.m, par=1, nsim=100, model="BM")[,,]
  dat.a <- sim.char.std.bm(phy.m, 100, x0=1)

  ## Note that this fits on the binary tree, not the multifurcating
  ## tree.
  do.fit <- function(x) {
    lik <- make.bm(phy, x, control=list(method="pruning"))
    unname(coef(find.mle(lik, 1, "optimize",
                         control=list(interval=c(0.1, 10)))))
  }

  est.g <- apply(dat.g, 2, do.fit)
  est.a <- apply(dat.a, 2, do.fit)

  expect_that(suppressWarnings(ks.test(est.a, est.g))$p.value,
              is_greater_than(0.1))
})

test_that("High level interface works", {
  data(geospiza)
  td <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))
  phy <- td$phy
  dat <- td$data[phy$tip.label,"wingL"]
  unit.tree <- make_unit_tree(phy, dat)

  nsim <- 2
  sims <- simulate_char_unit(unit.tree, nsim=nsim)
  expect_that(length(sims), equals(nsim))
  expect_that(sims, is_a("list"))
  expect_that(all(sapply(sims, inherits, "unit.tree")),
              is_true())

  sims <- simulate_char_unit(list(unit.tree, unit.tree), nsim=10)
  expect_that(length(sims), equals(nsim))
  expect_that(sims, is_a("list"))
  expect_that(all(sapply(sims, inherits, "unit.tree")),
              is_true())
})
