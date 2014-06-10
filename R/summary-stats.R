#' @title Calculate test statistics on a unit tree
#'
#' @description Calculates a set of test statistics on the contrasts included with a 'unit.tree' object.
#'
#' @param unit.tree a \code{unit.tree} object or list of unit.trees
#' @param stats a named list of test statistics to calculate on the unit.tree.
#'   If no \code{stats} argument supplied, default test statistics are used.
#'
#' @details This function can be applied to either a single 'unit.tree' of object or a list of 'unit.tree' objects.
#' If \code{stats=NULL} default test statistics are used. The default test statistics are the following:
#'  \enumerate{
#'   \item{m.sig: }{The mean of the squared contrasts. This is equivalent to the REML estimate of sigsq.}
#'
#'   \item{c.var: }{The coefficient of variation of the absolute value of the contrasts.}
#'
#'   \item{s.var: }{The slope of a linear model fit between the contrasts and their expected variances.}
#'
#'   \item{s.asr: }{The slope of a linear model fit between the contrasts and their inferred ancestral state.}
#'
#'   \item{s.hgt: }{The slope of a linear model fit between the contrasts and the node height at which they were calculated.}
#'
#'   \item{d.cdf: }{The D-statistic from a KS test comparing the distribution of the contrasts to a normal distribution  with mean 0 and variance equal to the square root of the squared mean of the contrasts.}
#'  }
#' User defined test statistics can be supplied as a named list of functions (see examples). The functions supplied
#' must take a unit.tree as argument and perform an operation on at least one of the elements of the object
#' (see \code{\link{make_unit_tree}} for details).
#'
#' @return A data.frame with the calculated test statistics across all unit.trees provided.
#'
#'
#' @export calculate_pic_stat
#'
#' @seealso \code{\link{default_pic_stat}}, \code{\link{pic_stat_msig}}, \code{\link{pic_stat_cvar}}, \code{\link{pic_stat_svar}}, \code{\link{pic_stat_shgt}}, \code{\link{pic_stat_sasr}}, \code{\link{pic_stat_dcdf}}
#'
#'
#' @examples
#' data(finch)
#' phy <- finch$phy
#' dat <- finch$data[,"wingL"]
#' unit.tree <- make_unit_tree(phy, data=dat)
#'
#' ## use default statistics
#' test.stat <- calculate_pic_stat(unit.tree, stats=NULL)
#' test.stat
#'
#' ## user defined statistics
#' mean.con <- function(x) mean(x$pics[,"contrasts"])
#' max.con <- function(x) max(x$pics[,"contrasts"])
#'
#' test.stat.user <- calculate_pic_stat(unit.tree,
#'                     stats=list(mean = mean.con, max = max.con))
#' test.stat.user
#'  
calculate_pic_stat <- function(unit.tree, stats=NULL){
  stats <- check_pic_stat(stats)
  if (inherits(unit.tree, "unit.tree")){
    # One unit.tree
    res <- do.call(cbind.data.frame,
                   lapply(stats, function(x) x(unit.tree)))
  } else if (is.list(unit.tree) && length(unit.tree) > 0) {
    # Multiple unit.trees
    res <- do.call(rbind, lapply(unit.tree, calculate_pic_stat, stats))
  } else {
    # Anything else
    stop("'unit.tree' must be a single unit.tree, or a list of them")
  }
  res
}









#' @title REML estimate of BM parameter
#'
#' @description Calculates the mean of the squared contrasts,
#' which is equal to the REML estimate of the BM rate parameter.
#' Used as a test statistic to evaluate model adequacy
#'
#' @param unit.tree a \code{unit.tree} object
#'
#' @details This test statistic is used to assess the ability of the model
#' to adequate estimate the rate of evolution. On datasets simulated on a unit.tree
#' (see \code{\link{simulate_char_unit}}), the expected value is equal to 1.
#'
#' The REML estimate of sigma2 is included as a default
#' test statistic in the function \code{\link{calculate_pic_stat}}.
#'
#' The test statistic computed from a single \code{unit.tree} can be visualized
#' with the function \code{\link{pic_stat_msig_plot}}.
#'
#' @return \code{m.sig} estimate of sigma2
#'
#' @export pic_stat_msig
#'
#' @seealso \code{\link{calculate_pic_stat}}, \code{\link{default_pic_stat}}, \code{\link{pic_stat_msig_plot}}
#'
#'
#' @examples
#' data(finch)
#' phy <- finch$phy
#' dat <- finch$data[,"wingL"]
#' unit.tree <- make_unit_tree(phy, data=dat)
#'
#' ## Estimate sigsq from unit.tree
#'
#' pic_stat_msig(unit.tree)
#'
#' ## Visualization
#'
#' pic_stat_msig_plot(unit.tree)
#' 
pic_stat_msig <- function(unit.tree){
    ## make sure the unit.tree is of class 'unit.tree'
    assert.is.unit.tree(unit.tree)

    mean(unit.tree$pics[,"contrasts"]^2)
}








#' @title Evaluate normality of contrasts
#'
#' @description Performs Kolmogorov-Smirnov test comparing distribution
#' of contrasts with that of a normal distribution with mean 0 and variance
#' equal to the square root of the mean of squared contrasts
#'
#' @param unit.tree a \code{unit.tree} object
#'
#' @details This test statistics is used to evaluate whether the assumption of multivariate normailty
#' is appropriate. If the model which generated the data is the fitted model, we expect the
#' square root ofthe mean of squared contrasts to be equal to 1. The empirical estimate is used
#' rather than assume a variance of 1 to reduce the overlap between this test statistic and the
#' REML estimate of sigma2 (see \code{\link{pic_stat_msig}}). The Kolmogorov-Smirnov (KS) test is a
#' non-parameteric test which computes the maximum distance \eqn{D} between two cumulative distribution functions.
#' Running the test multiple times on the same data will produce slightly different values due to the fact
#' that the null distribution is produced by randomly drawing from a normal distribution.
#'
#' The KS-D statistic is included as a default test statistic
#' in the function \code{\link{calculate_pic_stat}}.
#'
#' The test statistic computed from a single \code{unit.tree} can be visualized
#' with the function \code{\link{pic_stat_dcdf_plot}}.
#'
#' @return \code{d.cdf} the D-statistic from a KS-test
#'
#' @export pic_stat_dcdf
#'
#' @seealso \code{\link{calculate_pic_stat}}, \code{\link{default_pic_stat}}, \code{\link[stats]{ks.test}}
#'
#'
#' @examples
#' data(finch)
#' phy <- finch$phy
#' dat <- finch$data[,"wingL"]
#' unit.tree <- make_unit_tree(phy, data=dat)
#'
#' ## KS-D statistic
#' pic_stat_dcdf(unit.tree)
#'
#' ## Visualization
#' pic_stat_dcdf_plot(unit.tree)
#'
pic_stat_dcdf <- function(unit.tree){
    ## make sure the unit.tree is of class 'unit.tree'
    assert.is.unit.tree(unit.tree)
        
    ## simulate a null distribution
    sd <- sqrt(mean(unit.tree$pics[,"contrasts"]^2))
    nd <- rnorm(10000, mean=0, sd=sd)

    ## KS test
    ksbm <- suppressWarnings(ks.test(unit.tree$pics[,"contrasts"], nd)$statistic)

    unname(ksbm)
}











#' @title Relationship between contrasts and their variances
#'
#' @description Fits a linear model between the absolute value of constrasts
#' and their expected variances
#'
#' @param unit.tree a \code{unit.tree} object
#'
#' @details This test statistic estimates the slope from fitting
#' a linear model between the contrasts and their expected variance. It is used
#' to evaluate whether the model is adequately capturing variation relative to branch
#' lengths. If the model is adequate, we expect there not to be a relationship
#' between the contrasts and their variances (i.e. slope ~ 0).
#'
#' The slope of the contrasts and their variances is included as a default test statistic
#' in the function \code{\link{calculate_pic_stat}}
#'
#' The test statistic computed from a single \code{unit.tree} can be visualized
#' with the function \code{\link{pic_stat_svar_plot}}.
#'
#' @return \code{s.var} the estimated slope paramter
#'
#' @export pic_stat_svar
#'
#' @seealso \code{\link{calculate_pic_stat}}, \code{\link{default_pic_stat}},
#' \code{\link{pic_stat_svar_plot}}, \code{\link[stats]{lm}}
#'
#'
#' @examples
#' data(finch)
#' phy <- finch$phy
#' dat <- finch$data[,"wingL"]
#' unit.tree <- make_unit_tree(phy, data=dat)
#'
#' ## estimate slope
#' pic_stat_svar(unit.tree)
#'
#' ## visualization
#' pic_stat_svar_plot(unit.tree)
#' 
pic_stat_svar <- function(unit.tree){
    ## make sure the unit.tree is of class 'unit.tree'
    assert.is.unit.tree(unit.tree)
    
    ## sd of the pics
    sd.pic <- sqrt(unit.tree$pics[,"variance"])

    ## absolute value of the pics
    abs.pic <- abs(unit.tree$pics[,"contrasts"])

    ## fit a linear model
    cv <- lm(abs.pic ~ sd.pic)

    unname(coef(cv)["sd.pic"])
}









#' @title Relationship between contrasts and node height
#'
#' @description Fits a linear model between the absolute value of the contrasts
#' and the node height at which they were calculated
#'
#' @param unit.tree a \code{unit.tree} object
#'
#' @details This test statistic is the estimate of the slope from fitting
#' a linear model between the contrasts and the node height (i.e. tree depth)
#' at which they were calculated. It is used
#' to evaluate whether the model is adequately capturing variation relative to tree depth.
#' If the model is adequate, we expect there not to be a relationship
#' between the contrasts and their node height (i.e. slope ~ 0). Note that the node heights
#' being used here are the absolute node ages from the original phylogeny and not the adjusted
#' node ages obtained from computing contrasts (see \code{\link{pic}}).
#'
#' The slope of the contrasts and their variances is included as a default test statistic
#' in the function \code{\link{calculate_pic_stat}}.
#'
#' The test statistic computed from a single \code{unit.tree} can be visualized
#' with the function \code{\link{pic_stat_svar_plot}}.
#' 
#' @return \code{s.hgt} the estimated slope parameter
#'
#' @seealso \code{\link{calculate_pic_stat}}, \code{\link{default_pic_stat}}, \code{\link{pic_stat_svar_plot}}, \code{\link[stats]{lm}}
#'
#' @export pic_stat_shgt
#'
#'
#' @examples
#' data(finch)
#' phy <- finch$phy
#' dat <- finch$data[,"wingL"]
#' unit.tree <- make_unit_tree(phy, data=dat)
#'
#' ## estimate slope
#' pic_stat_shgt(unit.tree)
#'
#' ## visualizaiton
#' pic_stat_shgt_plot(unit.tree)
#' 
pic_stat_shgt <- function(unit.tree){
    ## make sure the unit.tree is of class 'unit.tree'
    assert.is.unit.tree(unit.tree)

    ## node heights 
    ht <- edge.height(unit.tree$phy)
    N <- Ntip(unit.tree$phy)
    nh <- ht$end[c((N + 1):nrow(ht))]
    names(nh) <- rownames(ht[c((N + 1):nrow(ht)),])

    ## absolute value of the pics
    abs.pic <- abs(unit.tree$pics[,"contrasts"])

    ## fit a linear model

    ch <- lm(abs.pic ~ nh)

    unname(coef(ch)["nh"])
}








#' @title Relationship between contrasts and inferred ancestral state
#'
#' @description Fits a linear model between the absolute value of the contrasts
#' and the ancestral state inferred at the node at which the contrast was calculated
#'
#' @param unit.tree a \code{unit.tree} object
#'
#' @details This test statistic is the estimate of the slope from fitting
#' a linear model between the contrasts and the inferred ancestral state at the node
#' at which the contrasts were calculated. It is used
#' to evaluate whether the model is adequately capturing variation relative to ancestral state.
#' The ancestral states are calculated using the contrast based approach and are not technically
#' ancestral state estimates (see Felsenstein 1985 -- American Naturalist); this is conceptually and statistically
#' different from the maximum likelihood approach to reconstructing ancestral states.
#' The contrast based approach is
#' used here as there is a clear expectation of the relationship between the the contrasts and
#' trait value. If the model is adequate, we expect there not to be a relationship
#' between the contrasts and their node height (i.e. slope ~ 0).
#'
#' The slope of the contrasts and their variances is included as a default test statistic
#' in the function \code{\link{calculate_pic_stat}}.
#'
#' The test statistic computed from a single \code{unit.tree} can be visualized
#' with the function \code{\link{pic_stat_sasr_plot}}.
#' 
#' @return \code{s.asr} the estimated slope parameter
#'
#' @seealso \code{\link{calculate_pic_stat}}, \code{\link{default_pic_stat}}, \code{\link{pic_stat_sasr_plot}},
#' \code{\link[stats]{lm}}, \code{\link[ape]{ace}}
#'
#'
#' @export pic_stat_sasr
#'
#' @examples
#' data(finch)
#' phy <- finch$phy
#' dat <- finch$data[,"wingL"]
#' unit.tree <- make_unit_tree(phy, data=dat)
#'
#' ## estimate slope
#' pic_stat_sasr(unit.tree)
#'
#' ## visualization
#' pic_stat_sasr_plot(unit.tree)
#' 
pic_stat_sasr <- function(unit.tree){
    ## make sure the unit.tree is of class 'unit.tree'
    assert.is.unit.tree(unit.tree)
    
    ## reconstruct ancestors using pic method
    asr <- ace(unit.tree$data, unit.tree$phy, method="pic")$ace

    ## absolute value of pics
    abs.pic <- abs(unit.tree$pics[,"contrasts"])

    ## fit linear model
    ca <- lm(abs.pic ~ asr)

    unname(coef(ca)["asr"])
}






#' @title Coefficient of variation of contrasts
#'
#' @description Calculates the coefficient of variation on a set of contrasts
#'
#' @param unit.tree a \code{unit.tree} object
#'
#' @details This test statistic is used to evaluate whether the model is sufficiently
#' capturing the variation in rates across the tree. The coefficient of variation is the sd/mean.
#'
#' The coefficient of variation of the contrasts is included as a default test statistic
#' in the function \code{\link{calculate_pic_stat}}.
#'
#' The test statistic computed from a single \code{unit.tree} can be visualized
#' with the function \code{\link{pic_stat_sasr_plot}}.
#'
#' @return \code{c.var} the coefficient of variation of the contrasts
#'
#' @export pic_stat_cvar
#'
#' @seealso \code{\link{calculate_pic_stat}}, \code{\link{default_pic_stat}}, \code{\link{pic_stat_cvar_plot}}
#'
#'
#' @examples
#' data(finch)
#' phy <- finch$phy
#' dat <- finch$data[,"wingL"]
#' unit.tree <- make_unit_tree(phy, data=dat)
#'
#' ## estimate variance of contrasts
#' pic_stat_cvar(unit.tree)
#'
#' ## visualization
#' pic_stat_cvar_plot(unit.tree)
#'
pic_stat_cvar <- function(unit.tree){
    ## make sure the unit.tree is of class 'unit.tree'
    assert.is.unit.tree(unit.tree)

    con <- abs(unit.tree$pics[,"contrasts"])
    sd(con)/mean(con)
}







#' @title Internal function for getting default test statistics
#'
#' @description Creates a list of the default test statistics to be used
#' to asses model adequacy
#'
#' @details The following test statistics are produced by this function:
#'  \enumerate{
#'   \item{m.sig: }{The mean of the squared contrasts. This is equivalent to the REML estimate of sigsq.}
#'
#'   \item{c.var: }{The coefficient of variation of the absolute value of the contrasts.}
#'
#'   \item{s.var: }{The slope of a linear model fit between the contrasts and their expected variances.}
#'
#'   \item{s.anc: }{The slope of a linear model fit between the contrasts and their inferred ancestral state.}
#'
#'   \item{s.hgt: }{The slope of a linear model fit between the contrasts and the node height at which they were calculated.}
#'
#'   \item{d.cdf: }{The D-statistic from a KS test comparing the distribution of the contrasts to a normal distribution  with mean 0 and variance equal to the square root of the squared mean of the contrasts.}
#'  }
#'
#' @return a named list of functions
#'
#' @export default_pic_stat
#'
#' @keywords internal
#'
#' @seealso \code{\link{calculate_pic_stat}}, \code{\link{pic_stat_msig}}, \code{\link{pic_stat_cvar}},
#' \code{\link{pic_stat_svar}}, \code{\link{pic_stat_sasr}}, \code{\link{pic_stat_shgt}},
#' \code{\link{pic_stat_dcdf}}
#'
#'
#' @examples
#' ## get default test statistics
#' stats <- default_pic_stat()
#' stats
default_pic_stat <- function()
    list("m.sig"=pic_stat_msig, "c.var"=pic_stat_cvar, "s.var"=pic_stat_svar, "s.asr"=pic_stat_sasr,
         "s.hgt"=pic_stat_shgt, "d.cdf"=pic_stat_dcdf)










#' @title Internal function to check test statistics
#'
#' @description Makes sure test statistics sent to \code{\link{calculate_pic_stat}} are in the correct format
#'
#' @param stats a named list of test statistics
#'
#' @return named list of test statistics after being checked.
#'
#' @seealso \code{\link{calculate_pic_stat}}
#'
#' @export check_pic_stat
#'
#' @keywords internal
#'
#'
#' @examples
#' ## produce list of test stats (using default test statistics)
#' check_pic_stat(default_pic_stat())
#'
#' ## use custom list of test statistics
#' foo <- function(x) mean(x$pics[,"contrasts"])
#'
#' my.stats <- list(mean.of.contrasts=foo)
#' check_pic_stat(my.stats)
#
check_pic_stat <- function(stats) {
  if (is.null(stats))
    stats <- default_pic_stat()
  else if (!is.list(stats) || is.null(names(stats)) ||
           any(names(stats) == "" || length(stats) == 0))
    ## (note that the second check does not always enforce
    ## length(stats) > 0)
    stop("'stats' must be a named list with at least one element")
  else if (!all(sapply(stats, is.function)))
    stop("All elements of stats must be functions")
  stats
}
