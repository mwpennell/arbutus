#' @title Assess the adequacy of continuous trait models
#'
#' @description Use test statistics to assess model adequacy for
#' phylogenetic models of continuous character evolution. This function is
#' a simply a wrapper which does all the steps involved in evaluating model
#' adequacy using the approach outlined in Pennell, FitzJohn, Cornwell and Harmon (in review).
#'
#' @param x a fitted model object or a \code{phylo}/\code{multiPhylo} (see Details).
#'
#' @param nsim the number of datasets to simulate. This is passed as an argument to the
#' function \code{\link{simulate_char_unit}}.
#'
#' @param stats a named list of test statistics to calculate on observed and simulated
#' datasets. See \code{\link{calculate_pic_stat}} for details. If nothing is supplied,
#' the function uses the default test statistics (see \code{\link{default_pic_stat}}
#'
#' @param ... additional arguments to be passed to \code{\link{make_unit_tree}}
#'
#' @details \code{arbutus} was developed to evaluate the adequacy (or absolute goodness-of-fit)
#' of phylogenetic models of continuous character evolution. The basic principle
#' underlying the approach is that if the generating model is a Brownian motion
#' process, the contrasts (sensu Felsenstein, 1985) will be independent and
#' identically distributed (I.I.D.). We can evaluate this condition by calculating
#' a set of test statistics on the contrasts of our observed data, then simulating datasets under
#' a Brownian motion process along the phylogeny and calculating the test statistics
#' on the contrasts of each simulated data set. We can then compare our observed test statistics to the
#' simulated distribution of summary statistics.
#'
#' The I.I.D. property of the contrasts will hold only for the case of Brownian motion.
#' However, the approach can be extended to any arbitrarily complex model of continuous character
#' evolution (as long as the model is based on a multivariate normal distribution) by first
#' rescaling the phylogeny based on the fitted parameter estimates, creating what we refer to
#' as a 'unit.tree'. Once the rescaling is done, the contrasts will be I.I.D. if the model
#' that was fit is adequate.
#'
#' P-values represent the two-tailed probability that the observed summary statistic
#' (or distribution of summary statistics) came from the same distribution as the simulated
#' summary statistics. Low p-values provide evidence that the model is inadequate.
#'
#' The function \code{arbutus} wraps several other functions:
#' \enumerate{
#'  \item{\code{\link{make_unit_tree}}: }{rescales phylogeny based on fitted parameter values}
#'  \item{\code{\link{calculate_pic_stat}}: }{calculates test statistics on observed data}
#'  \item{\code{\link{simulate_char_unit}}: }{simulates datasets under Brownian motion with rate 1}
#'  \item{\code{\link{calculate_pic_stat}}: }{calculates test statistics on simulated data}
#'  \item{\code{\link{compare_pic_stat}}: }{compares observed to simulated test statistics}
#' }
#'
#' \code{arbutus} interacts with objects produced by fitting evolutionary models using a variety of packages.
#' Currently supported objects are as follows:
#'
#' \itemize{
#'  \item a \code{gfit} object returned from fitting a model of continuous character evolution using
#'   \code{fitContinuous} in the \code{geiger} package.
#'
#'  \item a \code{fit.mle} object returned from fitting a model of continuous character evolution
#'   using \code{find.mle} in the \code{diversitree} package.
#'
#'  \item a \code{mcmcsamples} object returned from fitting a model of continuous character evolution
#'   using MCMC methods in the \code{diversitree} package.
#' \code{\link{make_unit_tree}} will apply the same
#'   trait dataset to a set of unit trees based on sampled parameters. By default this will create a
#'   unit tree for every sample in the mcmc chain. To modify this, additional arguments can be use.
#'   \code{burnin} specifies how many samples to remove from the beginning of the chain.
#'   \code{thin} specifies the thinning interval (e.g. if \code{thin=5}, the function will create a unit
#'   tree from every fifth parameter set sampled.
#'   \code{sample} specifies how many samples to draw from the MCMC run.
#'
#'  \item a \code{gls} object returned from fitting a phylogenetic generalized least squares model
#'   of character correlation using \code{gls} in the \code{nlme} package.
#'
#'  \item a \code{pgls} object returned from fitting a phylogenetic generalized least squares model
#'   of character correlation using \code{pgls} in the \code{caper} package.
#'
#'  \item a \code{phylolm} object returned from fitting a phylogenetic generalized linear model of
#'   character correlation using \code{phylolm} in the \code{phylolm} package. As the phylogeny is not
#'   returned with the \code{phylolm} object, a \code{phy} argument must also be specified.
#'
#'  \item a \code{phylo} object. If a \code{phylo} object is supplied, the tree is assumed to have been
#'   rescaled previously. A second argument \code{data} must also be provided included the trait
#'   data as a named vector with names equal to the tip.labels of the phylogeny.
#'
#'  \item a \code{multiPhylo object}. If a \code{multiPhylo} object is supplied, the tree is assumed to have been
#'   rescaled previously. A second argument \code{data} must also be provided included the trait
#'   data as a named vector with names equal to the tip.labels of the phylogenies. Note that this
#'   function will append the same data set to every tree in the \code{multiPhylo} object.
#' }
#'
#'
#' @export arbutus
#'
#' @seealso \code{\link{make_unit_tree}}, \code{\link{calculate_pic_stat}}, \code{\link{simulate_char_unit}}, \code{\link{compare_pic_stat}}
#'
#'
#' @examples
#' ## finch data
#' data(finch)
#' phy <- finch$phy
#' data <- finch$data[,"wingL"]
#'
#' \dontrun{
#' require(geiger)
#' ## fit Brownian motion model
#' ## using geiger's fitContinuous function
#' fit.bm <- fitContinuous(phy=phy, dat=data, model="BM",
#'                                  control=list(niter=10), ncores=1)
#'
#' ## check adequacy of BM model
#' ## get p-values for default summary statistics
#' modelad.bm <- arbutus(fit.bm, nsim=10)
#' modelad.bm
#'
#'
#' ## fit Ornstein-Uhlenbeck model
#' ## again, using geiger's fitContinuous function
#' fit.ou <- suppressWarnings(fitContinuous(phy=phy, dat=data,
#'                model="OU", control=list(niter=10), ncores=1))
#'
#' ## check adequacy of OU model
#' modelad.ou <- arbutus(fit.ou, nsim=10)
#'
#' require(diversitree)
#' ## fit Brownian motion model using ML
#' ## using diversitree's find.mle function
#'
#' bmlik <- make.bm(phy, data)
#' fit.bm.dt <- find.mle(bmlik, 0.1)
#'
#' ## this creates a 'fit.mle' object which can be used
#' ## in 'arbutus'
#' modelad.bm.dt <- arbutus(fit.bm.dt)
#'
#' ## fit Brownian motion model using MCMC
#' mcmc.bm.dt <- mcmc(bmlik, x.init=1, nsteps=1000, w=1)
#'
#' ## construct a unit tree object from mcmcsamples
#' ## removing 100 samples as burnin
#' ## and sampling 10 parameter sets
#' modelad.bm.dt.mcmc <- arbutus(mcmc.bm.dt,
#'                        burnin=100, samples=10)
#'
#'
#' require(nlme)
#' ## Use pgls to look for a correlation between two traits
#'
#' t1 <- data
#' t2 <- td$data[,"tarsusL"]
#' dd <- cbind.data.frame(t1, t2)
#'
#' ## fit gls model with corPagel correlation structure
#' fit.gls <- gls(t1~t2, data=dd, correlation=corPagel(phy=phy, value=1))
#'
#' ## this creates a 'gls' object which can be used
#' ## in 'arbutus'
#' modelad.gls <- arbutus(fit.gls)
#'
#' }
#'
arbutus <- function(x, nsim=1000, stats=NULL, ...){

    ## create unit tree
    unit.tree <- make_unit_tree(x, ...)

    ## compute test statistics on observed data
    obs <- calculate_pic_stat(unit.tree, stats=stats)

    ## simulate n datasets
    sim.dat <- simulate_char_unit(unit.tree, nsim=nsim)

    ## compute test statistics on simulated data
    sim <- calculate_pic_stat(sim.dat, stats=stats)

    ## compare test statistics between observed and simulated
    res <- compare_pic_stat(obs, sim)

    res
}
