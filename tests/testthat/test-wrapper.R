context("wrapper (arbutus)")

data(finch)
phy <- finch$phy
dat1 <- finch$data[,"wingL"]
data <- as.data.frame(finch$data)
multiphy <- list(phy, phy)
class(multiphy) <- "multiPhylo"


## fit bm using fitContinuous
fitc.res <- fitContinuous(phy, dat=dat1, SE=0, control=list(niter=10), ncores=1)


## fit bm using diversitree
bm.lik <- make.bm(phy, dat1, control=list(method="pruning", backend="C"))

## mle
bm.mle.res <- find.mle(bm.lik, x.init = 1)

## mcmc
bm.mcmc.res <- mcmc(bm.lik, x.init = 1, nsteps=10, w=1, print.every=0)


## fit correlation with gls
gls.res <- gls(wingL~tarsusL, data=data, correlation = corBrownian(phy=phy))


## fit correlation with caper
dd <- cbind(data, phy$tip.label)
colnames(dd)[ncol(dd)] <- "species"
cd <- caper::comparative.data(phy, dd, names.col = species)

caper.res <- caper::pgls(wingL~tarsusL, data=cd)

## fit correlation with phylolm
phylolm.res <- phylolm(wingL~tarsusL, data=data, phy=phy, model="BM")


## fit correlation with pgls diversitree
pgls.lik <- make.pgls(phy, wingL~tarsusL, data=data)

## find mle for pgls model
pgls.mle.res <- find.mle(pgls.lik, x.init = c(0.5, 0.5, 0.5))

pgls.mcmc.res <- mcmc(pgls.lik, x.init = c(0.5, 0.5, 0.5), nsteps=10, w=1, print.every = 0)


test_that("arbutus() runs properly", {
    mc.fitc <- arbutus(fitc.res, nsim=10)
    mc.bm.mle <- arbutus(bm.mle.res, nsim=10)
    mc.bm.mcmc <- arbutus(bm.mcmc.res, nsim=10)
    mc.gls <- arbutus(gls.res, nsim=10)
    mc.caper <- arbutus(caper.res, nsim=10)
    mc.phylolm <- arbutus(phylolm.res, nsim=10, phy=phy)
    mc.pgls.mle <- arbutus(pgls.mle.res, nsim=10)
    mc.pgls.mcmc <- arbutus(pgls.mcmc.res, nsim=10)
    mc.phy <- arbutus(phy, data=dat1, nsim=10)
    mc.multiphy <- arbutus(multiphy, data=dat1, nsim=10)

    expect_that(mc.fitc, is_a("arbutus"))
    expect_that(mc.bm.mle, is_a("arbutus"))
    expect_that(mc.bm.mcmc, is_a("arbutus"))
    expect_that(mc.gls, is_a("arbutus"))
    expect_that(mc.caper, is_a("arbutus"))
    expect_that(mc.phylolm, is_a("arbutus"))
    expect_that(mc.pgls.mle, is_a("arbutus"))
    expect_that(mc.pgls.mcmc, is_a("arbutus"))
    expect_that(mc.phy, is_a("arbutus"))
    expect_that(mc.multiphy, is_a("arbutus"))

    expect_that(nrow(mc.fitc$sim), equals(10))
    expect_that(nrow(mc.bm.mle$sim), equals(10))
    expect_that(nrow(mc.bm.mcmc$sim), equals(10))
    expect_that(nrow(mc.gls$sim), equals(10))
    expect_that(nrow(mc.caper$sim), equals(10))
    expect_that(nrow(mc.phylolm$sim), equals(10))
    expect_that(nrow(mc.pgls.mle$sim), equals(10))
    expect_that(nrow(mc.pgls.mcmc$sim), equals(10))
    expect_that(nrow(mc.phy$sim), equals(10))
    expect_that(nrow(mc.multiphy$sim), equals(2))

    expect_that(nrow(mc.fitc$obs), equals(1))
    expect_that(nrow(mc.bm.mle$obs), equals(1))
    expect_that(nrow(mc.bm.mcmc$obs), equals(10))
    expect_that(nrow(mc.gls$obs), equals(1))
    expect_that(nrow(mc.caper$obs), equals(1))
    expect_that(nrow(mc.phylolm$obs), equals(1))
    expect_that(nrow(mc.pgls.mle$obs), equals(1))
    expect_that(nrow(mc.pgls.mcmc$obs), equals(10))
    expect_that(nrow(mc.phy$obs), equals(1))
    expect_that(nrow(mc.multiphy$obs), equals(2))

    stats <- unlist(names(default_pic_stat()))

    st.fitc <- names(pvalue_arbutus(mc.fitc))
    st.bm.mle <- names(pvalue_arbutus(mc.bm.mle))
    st.bm.mcmc <- names(pvalue_arbutus(mc.bm.mcmc))
    st.gls <- names(pvalue_arbutus(mc.gls))
    st.caper <- names(pvalue_arbutus(mc.caper))
    st.phylolm <- names(pvalue_arbutus(mc.phylolm))
    st.pgls.mle <- names(pvalue_arbutus(mc.pgls.mle))
    st.pgls.mcmc <- names(pvalue_arbutus(mc.pgls.mcmc))
    st.phy <- names(pvalue_arbutus(mc.phy))
    st.multiphy <- names(pvalue_arbutus(mc.multiphy))

    expect_that(st.fitc, equals(stats))
    expect_that(st.bm.mle, equals(stats))
    expect_that(st.bm.mcmc, equals(stats))
    expect_that(st.gls, equals(stats))
    expect_that(st.caper, equals(stats))
    expect_that(st.phylolm, equals(stats))
    expect_that(st.pgls.mle, equals(stats))
    expect_that(st.pgls.mcmc, equals(stats))
    expect_that(st.phy, equals(stats))
    expect_that(st.multiphy, equals(stats))

})


test_that("if incomplete data is supplied, method fails",{

    ## phylolm without phylogeny
    expect_that(arbutus(phylolm.res), throws_error())

    ## phylogeny without data
    expect_that(arbutus(phy), throws_error())
    expect_that(arbutus(multiphy), throws_error())
})


test_that("stats are passed from wrapper to internals", {
    ## define some function
    foo <- function(x)
        mean(x$pics[,"contrasts"])


    ## create a named list with function
    new.stat <- list(mean=foo)

    mc.fitc <- arbutus(fitc.res, nsim=10, stats=new.stat)
    mc.bm.mle <- arbutus(bm.mle.res, nsim=10, stats=new.stat)
    mc.bm.mcmc <- arbutus(bm.mcmc.res, nsim=10, stats=new.stat)
    mc.gls <- arbutus(gls.res, nsim=10, stats=new.stat)
    mc.caper <- arbutus(caper.res, nsim=10, stats=new.stat)
    mc.phylolm <- arbutus(phylolm.res, nsim=10, stats=new.stat, phy=phy)
    mc.pgls.mle <- arbutus(pgls.mle.res, nsim=10, stats=new.stat)
    mc.pgls.mcmc <- arbutus(pgls.mcmc.res, nsim=10, stats=new.stat)
    mc.phy <- arbutus(phy, data=dat1, nsim=10, stats=new.stat)
    mc.multiphy <- arbutus(multiphy, data=dat1, nsim=10, stats=new.stat)

    st.fitc <- names(pvalue_arbutus(mc.fitc))
    st.bm.mle <- names(pvalue_arbutus(mc.bm.mle))
    st.bm.mcmc <- names(pvalue_arbutus(mc.bm.mcmc))
    st.gls <- names(pvalue_arbutus(mc.gls))
    st.caper <- names(pvalue_arbutus(mc.caper))
    st.phylolm <- names(pvalue_arbutus(mc.phylolm))
    st.pgls.mle <- names(pvalue_arbutus(mc.pgls.mle))
    st.pgls.mcmc <- names(pvalue_arbutus(mc.pgls.mcmc))
    st.phy <- names(pvalue_arbutus(mc.phy))
    st.multiphy <- names(pvalue_arbutus(mc.multiphy))

    stats <- "mean"

    expect_that(st.fitc, equals(stats))
    expect_that(st.bm.mle, equals(stats))
    expect_that(st.bm.mcmc, equals(stats))
    expect_that(st.gls, equals(stats))
    expect_that(st.caper, equals(stats))
    expect_that(st.phylolm, equals(stats))
    expect_that(st.pgls.mle, equals(stats))
    expect_that(st.pgls.mcmc, equals(stats))
    expect_that(st.phy, equals(stats))
    expect_that(st.multiphy, equals(stats))
})
