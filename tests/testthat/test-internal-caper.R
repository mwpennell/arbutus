## TODO: NEED A PROPER CHECK OF SIGSQ ESTIMATION!!
context("caper (internals)")

set.seed(1)
phy <- tree.bd(pars=c(1,0), max.taxa = 100)
tx <- sim.character(phy, 1)
ty <- sim.character(phy, 1) + 3
data <- data.frame(x=tx, y=ty, phy$tip.label, row.names=names(tx))
colnames(data)[3] <- "species"

## make comparative data for caper
cdat <- caper::comparative.data(phy, data, species)

## fit models
fit.pgls.bm <- caper::pgls(y~x, cdat)
fit.pgls.la <- caper::pgls(y~x, cdat, lambda="ML")
fit.pgls.de <- caper::pgls(y~x, cdat, delta="ML")
fit.pgls.ka <- caper::pgls(y~x, cdat, kappa="ML")

test_that("Model types are correct", {
    expect_that(model_type(fit.pgls.bm), is_identical_to("BM"))
    expect_that(model_type(fit.pgls.la), is_identical_to("lambda"))
    expect_that(model_type(fit.pgls.de), is_identical_to("delta"))
    expect_that(model_type(fit.pgls.ka), is_identical_to("kappa"))
})

test_that("Multiple models throws error", {
    fit.pgls.multi <- caper::pgls(y~x, cdat, lambda="ML", kappa="ML")
    expect_that(model_type(fit.pgls.multi), throws_error())
})

test_that("Models return their source data", {
    pgls.resid <- function(x){
        y <- as.numeric(resid(x))
        names(y) <- rownames(resid(x))
        y
    }
    cmp.bm <- list(phy=phy, data=pgls.resid(fit.pgls.bm))
    cmp.la <- list(phy=phy, data=pgls.resid(fit.pgls.la))
    cmp.de <- list(phy=phy, data=pgls.resid(fit.pgls.de))
    cmp.ka <- list(phy=phy, data=pgls.resid(fit.pgls.ka))

    expect_that(model_data(fit.pgls.bm), equals(cmp.bm))
    expect_that(model_data(fit.pgls.la), equals(cmp.la))
    expect_that(model_data(fit.pgls.de), equals(cmp.de))
    expect_that(model_data(fit.pgls.ka), equals(cmp.ka))
})

test_that("Coefficient names are as expected", {
    pars.bm <- c("sigsq", "SE", "z0")
    pars.la <- c("lambda", "sigsq", "SE", "z0")
    pars.de <- c("delta", "sigsq", "SE", "z0")
    pars.ka <- c("kappa", "sigsq", "SE", "z0")

    expect_that(names(model_pars(fit.pgls.bm)), is_identical_to(pars.bm))
    expect_that(names(model_pars(fit.pgls.la)), is_identical_to(pars.la))
    expect_that(names(model_pars(fit.pgls.de)), is_identical_to(pars.de))
    expect_that(names(model_pars(fit.pgls.ka)), is_identical_to(pars.ka))
})

test_that("Overall processed object looks legit", {
    obj.bm <- model_info(fit.pgls.bm)
    obj.la <- model_info(fit.pgls.la)
    obj.de <- model_info(fit.pgls.de)
    obj.ka <- model_info(fit.pgls.ka)

    obj.names <- c("data", "pars", "type")
    expect_that(names(obj.bm), is_identical_to(obj.names))
    expect_that(names(obj.la), is_identical_to(obj.names))
    expect_that(names(obj.de), is_identical_to(obj.names))
    expect_that(names(obj.ka), is_identical_to(obj.names))

    expect_that(class(obj.bm), is_identical_to("fitC"))
    expect_that(class(obj.la), is_identical_to("fitC"))
    expect_that(class(obj.de), is_identical_to("fitC"))
    expect_that(class(obj.ka), is_identical_to("fitC"))
})

