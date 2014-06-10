#' @title Compare observed and simulated test statistics
#'
#' @description Calculates two--tailed p-values from comparing
#' observed to simulated test statistics in order to evaluate
#' model adequacy for phylogenetic models of continuous character evolution.
#'
#' @param obs a data.frame containing the observed test statistics
#' with the column names denoting the names of the statistics used.
#'
#' @param sim a data.frame containing the simulated test statistics
#' with the column names denoting the names of the statistics used.
#'
#' @details The column names must be the same for the observed and simulated
#' test statistics (though not necessarily in the same order). Both
#' \code{obs} and \code{sim} are designed to be obtained
#' from the function \code{\link{calculate_pic_stat}}.
#'
#' If there is only 1 row in \code{obs} (such as the case if ML estimate
#' of model is used), then the function will compare the observed test statistic(s) to
#' the distribution of simulated simulated test statistics.
#'
#' If there is more than 1 row in \code{obs} (such as the case if we have
#' a posterior distribution of parameter estimates or the test statistics are calculated
#' across many trees), then the function requires that the number of rows in
#' \code{obs} and \code{sim} be equal. This is because the function
#' calculates the p-values by performing pairwise comparisons. the function will automatically
#' check that this condition is met and produce an error if it is not.
#'
#' In both cases, the p-values that are returned are two--tailed p-values as
#' the function considers the probability of the observed statistics being larger
#' or smaller than the simulated statistics.
#'
#' The distribution of simulated test statistics along with the observed test
#' statistic using a generic \code{\link{plot}} function (see \code{\link{plot.arbutus}})
#' for more details. For the purposes of plotting, the observed and simulated test
#' statistics used as arguments are returned with the p-values.
#'
#' The p-values can be extracted from the output of \code{\link{compare_pic_stat}} using the
#' function \code{\link{pvalue_arbutus}}
#'
#' @return an object of class \code{arbutus} which is a list containing the following objects:
#' \describe{
#'  \item{p.values}{Two--tailed p-values for all test statistics}
#'  \item{obs}{The input data frame of observed test statistics}
#'  \item{sim}{The input data frame of simulated test statistics}
#' }
#'
#' @seealso \code{\link{calculate_pic_stat}}, \code{\link{plot.arbutus}}, \code{\link{pvalue_arbutus}},
#' \code{\link{arbutus}}
#'
#' @export compare_pic_stat
#'
#' @examples
#' data(finch)
#' phy <- finch$phy
#' dat <- finch$data[,"wingL"]
#' unit.tree <- make_unit_tree(phy, data=dat)
#'
#' ## calculate default test stats on observed data
#' obs <- calculate_pic_stat(unit.tree, stats=NULL)
#'
#' ## simulate data on unit.tree
#' simdat <- simulate_char_unit(unit.tree, nsim=10)
#'
#' ## calculate default test stats on simulated data
#' sim <- calculate_pic_stat(simdat, stats=NULL)
#'
#' ## compare simulated to observed test statistics
#' res <- compare_pic_stat(obs, sim)
#'
#' ## extract p-values
#' pvalue_arbutus(res)
#'
#' ## plot results
#' plot(res)
#' 
compare_pic_stat <- function(obs, sim){
    ## check to make sure names are the same
    if (!all(names(obs) %in% names(sim)))
        stop("Column names must match between observed and simulated")
    ## make sure names are in correct order
    if (ncol(obs) > 1)
        obs <- obs[,names(sim)]
    ## determine if single/multiple obs test stat supplied
    tmp <- nrow(obs) == 1
    p.values <- vector(length=ncol(obs))
    if (tmp){
        for (i in seq_len(ncol(obs))){
            o <- obs[,i]
            s <- sim[,i]
            pr <- length(s[s > o])/(length(s) + 1)
            pl <- length(s[s <= o])/(length(s) + 1)
            p <- min(pr, pl)
            p.values[i] <- p*2
        }
    } else {
        ## check to make sure number of rows is the same
        if(nrow(obs) != nrow(sim))
            stop("If inputting multiple observed test stats, the data.frame must be the same size as that of the simulated test stats")
        
        for (i in seq_len(ncol(obs))){
            ## pairwise comparison
            pr <- vector(length=nrow(obs))
            pl <- vector(length=nrow(obs))
            for (j in seq_len(nrow(obs))){
                pr[j] <- sim[j,i] > obs[j,i]
                pl[j] <- sim[j,i] <= obs[j,i]
            }
            p.values[i] <- min(length(which(pr)), length(which(pl)))/nrow(obs)*2
        }
    }
    names(p.values) <- colnames(obs)
    res <- list(p.values=p.values, obs=obs, sim=sim)
    class(res) <- c("arbutus", "list")
    res
}
               
    









#' @title Print arbutus object
#'
#' @description Generic function for printing 'arbutus' objects
#' from \code{\link{compare_pic_stat}}
#'
#' @param x an \code{arbutus} object
#' @param ... additional arguments to be passed to print
#'
#' @return prints p.values to screen
#'
#'
#' @method print arbutus
#' @export
#'
#' @examples
#' data(finch)
#' phy <- finch$phy
#' dat <- finch$data[,"wingL"]
#'
#' ## use arbutus() to assess model adequacy
#' ## pretend as if tree has already been rescaled
#' p <- arbutus(phy, data=dat, nsim=10)
#'
#' ## print output
#' p
#' 
print.arbutus <- function(x, ...){
   cat("Two-tailed p-values for test statistics: \n")
   cat("\n")
   print(unclass(x$p.values))
}






#' @title Extract p--values for test statistics
#'
#' @description Utility function for extracting p-values from the output
#' of \code{\link{compare_pic_stat}}
#'
#' @param x an \code{arbutus} object from the function \code{\link{compare_pic_stat}}
#'
#' @return a named vector of two-tailed p-values
#'
#' @export pvalue_arbutus
#'
#' @seealso \code{\link{compare_pic_stat}}
#'
#'
#' @examples
#' data(finch)
#' phy <- finch$phy
#' dat <- finch$data[,"wingL"]
#' unit.tree <- make_unit_tree(phy, data=dat)
#'
#' ## calculate default test stats on observed data
#' obs <- calculate_pic_stat(unit.tree, stats=NULL)
#'
#' ## simulate data on unit.tree
#' sim.dat <- simulate_char_unit(unit.tree, nsim=10)
#'
#' ## calculate default test stats on simulated data
#' sim <- calculate_pic_stat(sim.dat, stats=NULL)
#'
#' ## compare simulated to observed test statistics
#' res <- compare_pic_stat(obs, sim)
#'
#' ## get p-values
#' pvalue_arbutus(res)
#'
#' ## note these are returned by default with print.arbutus
#' res
#'
pvalue_arbutus <- function(x){
    if (!inherits(x, "arbutus"))
        stop("'x' must be an arbutus object")
    
    x$p.values
}





#' @title Multivariate measure of model adequacy
#'
#' @description Computes Mahalanobis distance between the observed test statistics and
#' the simulated test statistics as a multivariate measure of model fit
#'
#' @param x an \code{arbutus} object inherited from \code{\link{compare_pic_stat}}
#'
#' @details This function computes the Mahalanobis distance between the observed and simulated
#' test statistics. The Mahalanobis distance (see \code{\link[stats]{mahalanobis}} is a unit-less and
#' scale-invariant of the distance
#' between a single data point (our observed test statistics) and a common point (here the mean of the
#' simulated test statistics), taking into account the covariance between the test statistics from
#' the simulated data. It assumes that the distribution of test statistics is multivariate normal.
#' For the default test statistics (see \code{\link{default_pic_stat}}, this condition should be met -- the
#' exception being the \code{\link{pic_stat_dcdf}} statistic, which as it is bounded at 0, will not be. As a result,
#' if \code{\link{pic_stat_dcdf}} is included in the set of test statistic, this function will take the
#' natural log of the values before computing the Mahalanobis distance. All other test statistics will
#' be assumed to be normally distributed and used as is.
#'
#' While the Mahalanobis distance may be a useful test measure in some circumstances, we recommend checking
#' the test statistics individually for a number of reasons. First, our procedure for calculating the p-values
#' for the test statistics is general and does not depend on assumptions regarding the distribution of values.
#' Second, the interpretation from the individual p-values is much more clear from the perspective of either
#' posterior predictive or parametric bootstrapping theory. Third, and most importantly, the fact that some test
#' statistics capture the variation in the data better than others provides useful information as to how and why the
#' model is inadequate.
#'
#' If only one set of observed
#' test statistics are available (e.g. from fitting a model using maximum likelihood to a single tree), a single distance
#' will be returned. If multiple sets of test statistics are available (e.g. from fitting a model using a Bayesian MCMC),
#' the function will return a distribution of distances.
#'
#' @return the Mahalanobis distance between the observed and simulated test statistic. 
#'
#' @export mahalanobis_arbutus
#'
#' @seealso \code{\link{compare_pic_stat}}, \code{\link[stats]{mahalanobis}}
#'
#' @examples
#' data(finch)
#' phy <- finch$phy
#' dat <- finch$data[,"wingL"]
#' unit.tree <- make_unit_tree(phy, data=dat)
#'
#' ## calculate default test stats on observed data
#' obs <- calculate_pic_stat(unit.tree, stats=NULL)
#'
#' ## simulate data on unit.tree
#' sim.dat <- simulate_char_unit(unit.tree, nsim=10)
#'
#' ## calculate default test stats on simulated data
#' sim <- calculate_pic_stat(sim.dat, stats=NULL)
#'
#' ## compare simulated to observed test statistics
#' res <- compare_pic_stat(obs, sim)
#'
#' ## calculate Mahalanobis distance
#' mahalanobis_arbutus(res)
#'
mahalanobis_arbutus <- function(x){
    if (!inherits(x, "arbutus"))
        stop("x must be a 'arbutus' object")
    ## log ks d-statistic
    if ("d.cdf" %in% names(pvalue_arbutus(x))){
        x$obs[,"d.cdf"] <- log(x$obs[,"d.cdf"])
        x$sim[,"d.cdf"] <- log(x$sim[,"d.cdf"])
    }
    obs <- as.matrix(x$obs)
    sim <- as.matrix(x$sim)
    cv <- cov(sim)

    m <- mahalanobis(x=obs, center=colMeans(sim), cov=cv)
    m
}
