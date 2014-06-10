## plotting functions for arbutus test statistics


#' @title Plot contrasts versus their variances
#'
#' @description Plot the contrasts derived from a rescaled 'unit.tree' against
#' their expected variances.
#'
#' @param x a \code{unit.tree} object (see \code{\link{make_unit_tree}}).
#'
#' @param xlab x-axis label. Defaults to "Variance".
#'
#' @param ylab y-axis label. Defaults to "Contrast".
#'
#' @param las plot parameter (see \code{\link{par}}). Defaults to 1.
#'
#' @param col plot colours. The first colour is the colour of the points, the second, the colour of the regression line. If no argument, uses default colours.
#'
#' @param pch plot parameter (see \code{\link{par}}). Defaults to 19.
#'
#' @param cex plot parameter (see \code{\link{par}}). Defaults to 1.
#'
#' @param lwd plot parameter (see \code{\link{par}}). Defaults to 2.
#'
#' @param ... additional arguments to be passed to \code{\link{plot}}.
#'
#' @details This function generates a plot of the observed contrasts versus their
#' expected variances (see \code{\link[ape]{pic}}).
#' The slope of this plot is equal to
#' the test statistic \code{\link{pic_stat_svar}}
#' and is used as a default test statistic in
#' \code{\link{calculate_pic_stat}}.
#'
#' If the model is a good fit to the data (i.e. an adequate model),
#' the slope ~ 0. The purpose of the plot is to get a quick visual glance
#' at the data and to assess whether or not the model inadequacy may be driven by outliers.
#'
#' The function can currently only take a single unit.tree
#' (does not integrate across multiple trees).
#'
#' @seealso \code{\link{make_unit_tree}}, \code{\link{calculate_pic_stat}},
#' \code{\link{pic_stat_svar}}
#'
#' @export pic_stat_svar_plot
#'
#' @examples
#' ## finch data
#' data(finch)
#' phy <- finch$phy
#' data <- finch$data[,"wingL"]
#'
#' u <- make_unit_tree(phy, data=data)
#'
#' pic_stat_svar_plot(u)
#' 
pic_stat_svar_plot <- function(x, xlab="Variance", ylab="Contrast",
                               las=1, col=NA, pch=19,
                               cex=1, lwd=2, ...){
    if (is.na(col))
        col <- default_colours_arbutus()

    if (length(col) == 1)
        col <- rep(col, 2)

    co <- abs(x$pics[,"contrasts"])
    va <- x$pics[,"variance"]
    res <- lm(co~va)
    plot(x=va, y=co, las=las, xlab=xlab, ylab=ylab,
         col=col[1], pch=pch, cex=cex, ...)
    abline(coef=coef(res), col=col[2], lwd=lwd)
}





#' @title Plot contrasts versus their inferred ancestral states
#'
#' @description Plot the contrasts derived from a rescaled 'unit.tree' against
#' their inferred ancestral state.
#'
#' @param x a \code{unit.tree} object (see \code{\link{make_unit_tree}}).
#'
#' @param xlab x-axis label. Defaults to "Ancestral state".
#'
#' @param ylab y-axis label. Defaults to "Contrast".
#'
#' @param las plot parameter (see \code{\link{par}}). Defaults to 1.
#'
#' @param col plot colours. The first colour is the colour of the points, the second, the colour of the regression line. If no argument, uses default colours.
#'
#' @param pch plot parameter (see \code{\link{par}}). Defaults to 19.
#'
#' @param cex plot parameter (see \code{\link{par}}). Defaults to 1.
#'
#' @param lwd plot parameter (see \code{\link{par}}). Defaults to 2.
#'
#' @param ... additional arguments to be passed to \code{\link{plot}}.
#'
#' @details This function generates a plot of the observed contrasts versus their
#' inferred ancestral state (see \code{\link[ape]{pic}}, \code{\link[ape]{ace}}).
#' Note that the ancestral states are estimated using the 'pic' method, which
#' differs somewhat from the ML estimate of the ancestral states.
#' The slope of this plot is equal to
#' the test statistic \code{\link{pic_stat_sasr}}
#' and is used as a default test statistic in
#' \code{\link{calculate_pic_stat}}. 
#'
#' If the model is a good fit to the data (i.e. an adequate model),
#' the slope ~ 0. The purpose of the plot is to get a quick visual glance
#' at the data and to assess whether or not the model inadequacy may be driven by outliers.
#'
#' The function can currently only take a single unit.tree
#' (does not integrate across multiple trees).
#'
#' @seealso \code{\link{make_unit_tree}}, \code{\link{calculate_pic_stat}},
#' \code{\link{pic_stat_sasr}}
#'
#' @export pic_stat_sasr_plot
#'
#' @examples
#' ## finch data
#' data(finch)
#' phy <- finch$phy
#' data <- finch$data[,"wingL"]
#'
#' u <- make_unit_tree(phy, data=data)
#'
#' pic_stat_sasr_plot(u)
#' 
pic_stat_sasr_plot <- function(x, xlab="Ancestral state",
                               ylab="Contrast",las=1, col=NA,pch=19,
                               cex=1, lwd=2, ...){

    assert.is.unit.tree(x)
    
    if (is.na(col))
        col <- default_colours_arbutus()

    if (length(col) == 1)
        col <- rep(col, 2)
    
    co <- abs(x$pics[,"contrasts"])
    an <- ace(x$data, x$phy, method="pic")$ace
    res <- lm(co~an)
    plot(x=an, y=co, las=las, xlab=xlab, ylab=ylab,
         col=col[1], pch=pch, cex=cex, ...)
    abline(coef=coef(res), col=col[2], lwd=lwd)
}












#' @title Plot contrasts versus their node height
#'
#' @description Plot the contrasts derived from a rescaled 'unit.tree' against
#' the node height (i.e. depth) at which they are computed
#'
#' @param x a \code{unit.tree} object (see \code{\link{make_unit_tree}}).
#'
#' @param xlab x-axis label. Defaults to "Node height".
#'
#' @param ylab y-axis label. Defaults to "Contrast".
#'
#' @param las plot parameter (see \code{\link{par}}). Defaults to 1.
#'
#' @param col plot colours. The first colour is the colour of the points, the second, the colour of the regression line. If no argument, uses default colours.
#'
#' @param pch plot parameter (see \code{\link{par}}). Defaults to 19.
#'
#' @param cex plot parameter (see \code{\link{par}}). Defaults to 1.
#'
#' @param lwd plot parameter (see \code{\link{par}}). Defaults to 2.
#'
#' @param ... additional arguments to be passed to \code{\link{plot}}.
#'
#' @details This function generates a plot of the observed contrasts versus their
#' the node height at which they are computed (see \code{\link[ape]{pic}}).
#' The slope of this plot is equal to
#' the test statistic \code{\link{pic_stat_shgt}}
#' and is used as a default test statistic in
#' \code{\link{calculate_pic_stat}}. 
#'
#' If the model is a good fit to the data (i.e. an adequate model),
#' the slope ~ 0. The purpose of the plot is to get a quick visual glance
#' at the data and to assess whether or not the model inadequacy may be driven by outliers.
#'
#' The function can currently only take a single unit tree (does not integrate across multiple trees).
#'
#' @seealso \code{\link{make_unit_tree}}, \code{\link{calculate_pic_stat}},
#' \code{\link{pic_stat_shgt}}
#'
#' @export pic_stat_shgt_plot
#'
#' @examples
#' ## finch data
#' data(finch)
#' phy <- finch$phy
#' data <- finch$data[,"wingL"]
#'
#' u <- make_unit_tree(phy, data=data)
#'
#' pic_stat_shgt_plot(u)
#' 
pic_stat_shgt_plot <-  function(x, xlab="Node height",
                               ylab="Contrast",las=1, col=NA, pch=19,
                               cex=1, lwd=2, ...){

    assert.is.unit.tree(x)
    
    if (is.na(col))
        col <- default_colours_arbutus()

    if (length(col) == 1)
        col <- rep(col, 2)

    ## extract pics
    co <- abs(x$pics[,"contrasts"])

    ht <- edge.height(x$phy)
    N <- Ntip(x$phy)
    nh <- ht$end[c((N + 1):nrow(ht))]
    names(nh) <- rownames(ht[c((N + 1):nrow(ht)),])

    res <- lm(co~nh)
    plot(x=nh, y=co, las=las, xlab=xlab, ylab=ylab,
         col=col[1], pch=pch, cex=cex, ...)
    abline(coef=coef(res), col=col[2], lwd=lwd)
}






#' @title Plot empirical cumulative distribution of contrasts versus a normal
#' distribution
#'
#' @description Plot the empirical cumulative distribution of the contrats
#' derived from a rescaled 'unit.tree' against the ecdf of a normal distribution
#' with mean 0 and sd of the observed contrasts.
#'
#' @param x a \code{unit.tree} object (see \code{\link{make_unit_tree}}).
#'
#' @param xlab x-axis label. Defaults to "Sample quantiles".
#'
#' @param ylab y-axis label. Defaults to "expression(hat(F)[n](x))".
#'
#' @param col plot colours. The first colour is the colour of the cdf of the
#' contrats and the second is that of the normal distribution. If no argument,
#' uses default colours.
#'
#' @param las plot parameter (see \code{\link{par}}). Defaults to 1.
#'
#' @param lwd plot parameter (see \code{\link{par}}). Defaults to 2.
#'
#' @param legend logical, whether legend should be plotted.
#'
#' @param legend.obs legend label for distribution of contrasts. Defaults to
#' "Contrats"
#'
#' @param legend.norm legend label for normal distribution. Defaults to
#' "Normal"
#'
#' @param cex plot parameter for legend (see \code{link{par}}). Defaults to 1.
#'
#' @param ... additional arguments to be passed to \code{\link{plot}}.
#'
#' @details This function plots the empirical cumulative distribution function
#' of the contrasts computed on the unit tree against that of a normal
#' distribution. The empirical sd of the contrasts is used rather than the
#' expected (i.e., 1) to distinguish non-normailty vs. inflated/deflated
#' variance.
#'
#' The D-statistics from the KS test (used in \code{\link{pic_stat_dcdf}}) is
#' the maximum distance between the two curves. This is a default test
#' test statistic in \code{\link{calculate_pic_stat}}.
#'
#' If the model is adequate, the two curves should line up.
#'
#' The function can currently only take a single unit.tree (does not
#' integrate across unit.trees)
#'
#' @seealso \code{\link{make_unit_tree}}, \code{\link{calculate_pic_stat}},
#' \code{\link{pic_stat_dcdf}}
#'
#' @export pic_stat_dcdf_plot
#'
#' @examples
#' ## finch data
#' data(finch)
#' phy <- finch$phy
#' data <- finch$data[,1]
#'
#' u <- make_unit_tree(phy, data=data)
#'
#' pic_stat_dcdf_plot(u)
#' 
pic_stat_dcdf_plot <- function(x, xlab="Sample quantiles",
                               ylab=expression(hat(F)[n](x)), col=NA, las=1,
                               lwd=2, legend=TRUE, legend.obs="Contrasts",
                               legend.norm="Normal", cex=1, ...){

    assert.is.unit.tree(x)
    
    if (is.na(col))
        col <- default_colours_arbutus()

    if (length(col) == 1)
        col <- rep(col,2)

    ## extract pics
    co <- x$pics[,"contrasts"]
    ## ecdf of contrats
    co.ecdf <- sort(co)
    n <- length(co.ecdf)
    ## ecdf of a normal distribution
    nd <- rnorm(10000, mean=0, sd=sqrt(mean(x$pics[,"contrasts"]^2)))
    nd.ecdf <- ecdf(nd)
    plot(co.ecdf, (1:n)/n, type="s", lwd=lwd,
         las=las, xlab=xlab, ylab="", col=col[1], ...)
    lines(nd.ecdf, col=col[2], do.points=FALSE, lwd=lwd)
    mtext(text=expression(hat(F)[n](x)), side=2, line=2.5)
    if (legend)
        legend("bottomright", c(legend.obs, legend.norm), lwd=lwd,
           col=col, bty="n", cex=cex)
}




#' @title Plot histogram of squared contrasts
#'
#' @description This is a diagnostic plot for looking at the distribution
#' of squared contrasts derived from a rescaled 'unit.tree' object.
#'
#' @param x a \code{unit.tree} object (see \code{\link{make_unit_tree}}).
#'
#' @param xlab x-axis label. Defaults to "Squared contrasts".
#'
#' @param ylab y-axis label. Defaults to "Counts".
#'
#' @param col plot colours. The first colour is used for the histogram, the
#' second designates the mean. If no argument supplied, uses default colours.
#'
#' @param breaks number of breaks for histogram. Defaults to the number of
#' contrats/10.
#'
#' @param las plot parameter (see \code{\link{par}}). Defaults to 1.
#'
#' @param lwd plot parameter (see \code{\link{par}}). Defaults to 2.
#'
#' @param ... aditional arguments to be passed to \code{\link{plot}}.
#'
#' @details This function is just intended to be a quick diagnostic for looking
#' at the distribution of squared contrats. The mean value (the vertical line
#' in the plot) is equal to the value of \code{\link{pic_stat_msig}} and is
#' the REML estimate of the Brownian motion rate parameter. It is a default
#' test statistic in \code{\link{calculate_pic_stat}}.
#'
#' If the model is adequate the mean should be close to 1.
#'
#' This function can currently only take a single unit.tree (does not
#' integrate across unit.trees)
#'
#' @seealso \code{\link{make_unit_tree}}, \code{\link{calculate_pic_stat}},
#' \code{\link{pic_stat_msig}}
#'
#' @export pic_stat_msig_plot
#'
#' @examples
#' ## finch data
#' data(finch)
#' phy <- finch$phy
#' data <- finch$data[,1]
#'
#' u <- make_unit_tree(phy, data=data)
#'
#' pic_stat_msig_plot(u)
#' 
pic_stat_msig_plot <- function(x, xlab="Squared contrasts", ylab="Counts",
                               col=NA, breaks=NA, las=1, lwd=2, ...){

    assert.is.unit.tree(x)
    
    if (is.na(col))
        col <- default_colours_arbutus()

    if (length(col) == 1)
        col <- rep(col,2)

    sq <- x$pics[,"contrasts"]^2

    if (is.na(breaks))
        breaks <- round(length(sq)/10)

    xlim <- range(sq)
    hist.sq <- hist(sq, breaks, plot=FALSE)
    ylim <- range(hist.sq$counts)
    plot(NA, las=las, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
    hist_outline(hist.sq, col=col[1], lwd=lwd, density=FALSE)
    abline(v=mean(sq), col=col[2], lwd=lwd)
}



#' @title Plot histogram of contrasts with variance
#'
#' @description This is a diagnostic plot for looking at the distribution
#' of contrasts derived from a rescaled 'unit.tree' object.
#'
#' @param x a \code{unit.tree} object (see \code{\link{make_unit_tree}}).
#'
#' @param xlab x-axis label. Defaults to "Contrasts".
#'
#' @param ylab y-axis label. Defaults to "Counts".
#'
#' @param col plot colours. The first colour is used for the histogram, the
#' second designates the mean. If no argument supplied, uses default colours.
#'
#' @param breaks number of breaks for histogram. Defaults to the number of
#' contrats/10.
#'
#' @param las plot parameter (see \code{\link{par}}). Defaults to 1.
#'
#' @param lwd plot parameter (see \code{\link{par}}). Defaults to 2.
#'
#' @param ... aditional arguments to be passed to \code{\link{plot}}.
#'
#' @details This function is just intended to be a quick diagnostic for looking
#' at the distribution of contrats. The standard deviation (the arrows
#' in the plot)/mean is equal to the value of \code{\link{pic_stat_cvar}}. It
#' is a default test statistic in \code{\link{calculate_pic_stat}}.
#'
#' If the model is adequate the sd should be close to 1.
#'
#' This function can currently only take a single unit.tree (does not
#' integrate across unit.trees)
#'
#' @seealso \code{\link{make_unit_tree}}, \code{\link{calculate_pic_stat}},
#' \code{\link{pic_stat_cvar}}
#'
#' @export pic_stat_cvar_plot
#'
#' @examples
#' ## finch data
#' data(finch)
#' phy <- finch$phy
#' data <- finch$data[,1]
#'
#' u <- make_unit_tree(phy, data=data)
#'
#' pic_stat_cvar_plot(u)
#' 
pic_stat_cvar_plot <- function(x, xlab="Contrasts", ylab="Counts",
                               col=NA, breaks=NA, las=1, lwd=2, ...){

    assert.is.unit.tree(x)

    if (is.na(col))
        col <- default_colours_arbutus()

    if (length(col) == 1)
        col <- rep(col, 2)

    co <- abs(x$pics[,"contrasts"])

    if (is.na(breaks))
        breaks <- round(length(co)/10)

    xlim <- range(co)
    hist.co <- hist(co, breaks, plot=FALSE)
    ylim <- range(hist.co$counts)

    ## arrows
    sdcon <- c((mean(co) - sd(co)), (mean(co) + sd(co)))
    
    plot(NA, las=las, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
    hist_outline(hist.co, col=col[1], lwd=lwd, density=FALSE)
    arrows(sdcon[1], 0, sdcon[2], 0, code=3, angle=90, length=0.02,
           lwd=lwd, col=col[2])
    
}

    

#' @title Plot the simulated test statistics versus observed test statistic
#' (point estimate of parameters)
#'
#' @description Plots a histogram of the simulated values along with
#' the observed values for a single test statistic when a point estimate is used
#'
#' @param x an \code{arbutus} object inherited from
#' \code{\link{compare_pic_stat}}
#'
#' @param stat the name of the test statistic to plot
#'
#' @param xlab x-axis label. Defaults to \code{stat}.
#'
#' @param ylab y-axis label. Defaults to "Counts".
#'
#' @param col colours used in plot. The first colour corresponds to the
#' histogram of simulated test statistics. The second colour to the line
#' representing the observed value of the test statistic. If no argument is
#' supplied, the function uses the default colours.
#'
#' @param breaks the number of breaks in the histogram. Defaults to
#' 20.
#'
#' @param las plot parameter (see \code{\link{par}}).
#'
#' @param lwd plot parameter (see \code{\link{par}}).
#'
#' @param legend logical, specifying whether a legend should be included.
#'
#' @param legend.pos postiion of the legend (see \code{\link{legend}}).
#'
#' @param cex plot parameter for legend (see \code{\link{par}}).
#'
#' @param ... additional arguments to be passed to plot
#'
#' @details This function is designed to give researchers a quick visual
#' inspection of model adequacy by plotting the simulated distribution of
#' a test statistic along with the empirical estimate.
#'
#' If the empirical estimate lies in the tails of the distribution of simulated
#' test statistics, the model can be rejected as inadequate.
#'
#' The p-values corresponding to this plot can be obtained with the function
#' \code{\link{pvalue_arbutus}}.
#'
#' \code{\link{plot.arbutus}} wraps this function and allows results from all
#' test statistics to be viewed in a single plot.
#'
#' @seealso \code{\link{arbutus}}, \code{\link{compare_pic_stat}},
#' \code{\link{plot.arbutus}}
#'
#' @export arbutus_plot_point
#'
#' @examples
#' ## finch data
#' data(finch)
#' phy <- finch$phy
#' data <- finch$data[,1]
#'
#' ## pretend that the tree has already been rescaled
#' modelad <- arbutus(phy, data=data, nsim=20)
#' arbutus_plot_point(modelad, "m.sig")
#'
#'
#' \dontrun{
#' ## use diversitree to fit a BM model with ML
#' require(diversitree)
#' lik <- make.bm(phy, data)
#' f <- find.mle(lik, x.init=1)
#' modelad.bm <- arbutus(f)
#' arbutus_plot_point(modelad.bm, "m.sig")
#' }
#' 
arbutus_plot_point <- function(x, stat, xlab=NA, ylab="Counts", col=NA,
                               breaks=NA, las=1, lwd=3, legend=TRUE,
                               legend.pos="topright",
                               cex=1, ...){

    if (!inherits(x, "arbutus"))
        stop("x needs to be an object of class 'arbutus'")

    sim <- x$sim[,stat]
    obs <- x$obs[,stat]

    ## set defaults
    if (is.na(col))
        col <- default_colours_arbutus()

    if (length(col) != 2)
        stop("If 'col' is input as argument, length must equal 2")

    if (is.na(breaks))
        breaks <- 20

    if (is.na(xlab))
        xlab <- stat

        ## include vlines outside of histogram
    xlim <- c(min(obs, min(sim)), max(obs,max(sim)))
   
    hist.sim <- hist(sim, breaks=breaks, plot=FALSE)

    ## set ylim
    ylim <- range(hist.sim$counts)
    
    plot(NA, las=1, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
    hist_outline(hist.sim, col=col[1], lwd=lwd, density=FALSE)
    abline(v=obs, col=col[2], lwd=lwd)

    if (legend){
        legend(legend.pos, c("Empirical", "Simulated"), lwd=lwd,
               col=rev(col), bty="n", cex=cex)
    }
    
}



#' @title Plot the simulated test statistics versus observed test statistic
#' (distribution of parameter estimates)
#'
#' @description Plots a histogram of the simulated values along a histogram of
#' the observed values for a single test statistic when a distribution of
#' estimates is used
#'
#' @param x an \code{arbutus} object inherited from
#' \code{\link{compare_pic_stat}}
#'
#' @param stat the name of the test statistic to plot
#'
#' @param xlab x-axis label. Defaults to \code{stat}.
#'
#' @param ylab y-axis label. Defaults to "Counts".
#'
#' @param col colours used in plot. The first colour corresponds to the
#' histogram of simulated test statistics. The second colour to the line
#' representing the observed value of the test statistic. If no argument is
#' supplied, the function uses the default colours.
#'
#' @param breaks the number of breaks in the histogram. Defaults to
#' 20.
#'
#' @param las plot parameter (see \code{\link{par}}).
#'
#' @param lwd plot parameter (see \code{\link{par}}).
#'
#' @param legend logical, specifying whether a legend should be included.
#'
#' @param legend.pos postiion of the legend (see \code{\link{legend}}).
#'
#' @param cex plot parameter for legend (see \code{\link{par}}).
#'
#' @param ... additional arguments to be passed to plot
#'
#' @details This function is designed to give researchers a quick visual
#' inspection of model adequacy by plotting the simulated distribution of
#' a test statistic along with the distribution of empirical estimates.
#'
#' If the empirical distribution of the test statistic lies in the tails of the
#' distribution of simulated test statistics, the model can be rejected
#' as inadequate.
#'
#' The p-values corresponding to this plot can be obtained with the function
#' \code{\link{pvalue_arbutus}}.
#'
#' \code{\link{plot.arbutus}} wraps this function and allows results from all
#' test statistics to be viewed in a single plot.
#'
#' @seealso \code{\link{arbutus}}, \code{\link{compare_pic_stat}},
#' \code{\link{plot.arbutus}}
#'
#' @export arbutus_plot_dist
#'
#' @examples
#' ## finch data
#' data(finch)
#' phy <- finch$phy
#' data <- finch$data[,1]
#'
#' ## pretend that the tree has already been rescaled
#' ## and that we have a list of rescaled trees
#' mphy <- lapply(seq_len(50), function(x) {phy})
#' class(mphy) <- "multiPhylo"
#' modelad <- arbutus(mphy, data=data)
#' arbutus_plot_dist(modelad, "m.sig")
#'
#'
#' \dontrun{
#' ## use diversitree to fit a BM model with MCMC
#' require(diversitree)
#' lik <- make.bm(phy, data)
#' f <- mcmc(lik, x.init=1, nsteps=1000, w=1, print.every=0
#' modelad.bm <- arbutus(f, sample=100)
#' arbutus_plot_dist(modelad.bm, "m.sig")
#' }
#' 
arbutus_plot_dist <- function(x, stat, xlab=NA, ylab="Counts", col=NA,
                              breaks=NA, las=1, lwd=3, legend=TRUE,
                              legend.pos="topright",
                              cex=1, ...){

    if (!inherits(x, "arbutus"))
        stop("x needs to be an object of class 'arbutus'")

    sim <- x$sim[,stat]
    obs <- x$obs[,stat]

    ## set defaults
    if (is.na(col))
        col <- default_colours_arbutus()

    if (length(col) != 2)
        stop("If 'col' is input as argument, length must equal 2")

    if (is.na(breaks))
        breaks <- 20 

    if (is.na(xlab))
        xlab <- stat
    
    hist.sim <- hist(sim, breaks=breaks, plot=FALSE)
    hist.obs <- hist(obs, breaks=breaks, plot=FALSE)

    ## set plot limits
    ylim <- range(hist.sim$counts, hist.obs$counts)
    xlim <- range(hist.sim$breaks, hist.obs$breaks)
    
    plot(NA, las=1, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
    hist_outline(hist.sim, col=col[1], lwd=lwd, density=FALSE)
    hist_outline(hist.obs, col=col[2], lwd=lwd, density=FALSE)

    if (legend){
        legend(legend.pos, c("Empirical", "Simulated"), lwd=lwd,
               col=rev(col), bty="n", cex=cex)
    }
}



#' @title Plot the simulated test statistics versus observed test statistics
#'
#' @description Plots a histogram of the simulated values along with
#' the observed values for a single test statistic. Either a point estimate or
#' distribution of observed values may be used. Is the generic \code{\link{plot}}
#' function for objects of class \code{arbutus}.
#'
#' @param x an \code{arbutus} object inherited from
#' \code{\link{compare_pic_stat}}
#'
#' @param xlab x-axis label. Defaults to \code{stat}.
#'
#' @param ylab y-axis label. Defaults to "Counts".
#'
#' @param col colours used in plot. The first colour corresponds to the
#' histogram of simulated test statistics. The second colour to the line, or
#' histogram,
#' representing the observed value of the test statistic. If no argument is
#' supplied, the function uses the default colours.
#'
#' @param breaks the number of breaks in the histogram. Defaults to
#' 20.
#'
#' @param las plot parameter (see \code{\link{par}}).
#'
#' @param lwd plot parameter (see \code{\link{par}}).
#'
#' @param legend logical, specifying whether a legend should be included (see
#' details).
#'
#' @param legend.pos postiion of the legend (see \code{\link{legend}}).
#'
#' @param cex plot parameter for legend (see \code{\link{par}}).
#'
#' @param n.row number of rows to plot when plotting multiple test statistics.
#' Defaults to 2.
#'
#' @param n.col number of columns to plot when plotting multiple test statistics.
#' Defaults to \code{ncol(obs)/2}.
#'
#' @param mar margins for each plot. Defaults to 5,4,2,2.
#'
#' @param ... additional arguments to be passed to plot
#'
#' @details This function is designed to give researchers a quick visual
#' inspection of model adequacy by plotting the simulated distribution of
#' all test statistics versus the emprical values.
#'
#' This function is a simply a wrapper for the functions
#' \code{\link{arbutus_plot_point}} and \code{\link{arbutus_plot_dist}}. It
#' creates plots for all test statistics used in \code{\link{arbutus}} and
#' \code{\link{compare_pic_stat}}. The function can be supplied either a point
#' estimate for the observed test statistics or a distribution of estimates.
#'
#' If the empirical estimate lies in the tails of the distribution of simulated
#' test statistics, the model can be rejected as inadequate.
#'
#' The p-values corresponding to this plot can be obtained with the function
#' \code{\link{pvalue_arbutus}}.
#'
#' We are still working on a general function for plotting legends. For the time
#' being, if \code{legend=TRUE}, a legend will appear with each test statistic.
#' This is obviously not ideal, but we are working on it.
#'
#' If users wish to create their own plots, all of the information is stored
#' in the \code{arbutus} object. \code{x$obs} contains a data.frame of observed
#' test statistics and \code{x$sim}, a data.frame of simulated test statistics.
#' 
#' @seealso \code{\link{arbutus}}, \code{\link{compare_pic_stat}},
#' \code{\link{arbutus_plot_point}}, \code{\link{arbutus_plot_dist}}
#'
#' @method plot arbutus
#' @export
#'
#' @examples
#' ## finch data
#' data(finch)
#' phy <- finch$phy
#' data <- finch$data[,1]
#'
#' ## pretend that the tree has already been rescaled
#' modelad <- arbutus(phy, data=data, nsim=20)
#' plot(modelad)
#'
#'
#' \dontrun{
#' ## use diversitree to fit a BM model with ML
#' require(diversitree)
#' lik <- make.bm(phy, data)
#' f <- find.mle(lik, x.init=1)
#' modelad.bm.ml <- arbutus(f)
#' plot(modelad.bm.ml, legend=FALSE)
#'
#' ## use diversitree to fit a BM model with MCMC
#' m <- mcmc(lik, x.init=1, nsteps=1000, w=1, print.every=0)
#' modelad.bm.mcmc <- arbutus(m, sample=100)
#' plot(modelad.bm.mcmc)
#' }
#' 
plot.arbutus <- function(x, xlab=NA, ylab="Counts", col=NA,
                         breaks=NA, las=1, lwd=3, legend=TRUE,
                         legend.pos="topright",
                         cex=1, n.row=2, n.col=NA, mar=c(5,4,2,2), ...){

    obs <- x$obs
    sim <- x$sim
    stats <- colnames(obs)
    type <- nrow(obs) == 1

    if (is.na(n.col) && ncol(obs) > 1)
        n.col <- round(ncol(obs)/2)

    if (is.na(n.col) && ncol(obs) == 1)
        n.col <- 1
 
    par(mfrow=c(n.row,n.col), mar=mar)

    if (type)
        tmp <- sapply(stats, function(i)
                      arbutus_plot_point(x, i, xlab=xlab, ylab=ylab,
                                         col=col, breaks=breaks, las=las,
                                         legend=legend, legend.pos=legend.pos,
                                         cex=cex,...))

    if (!type)  
        tmp <- sapply(stats, function(i)
                      arbutus_plot_dist(x, i, xlab=xlab, ylab=ylab,
                                         col=col, breaks=breaks, las=las,
                                         legend=legend, legend.pos=legend.pos,
                                         cex=cex,...))

}



## Set default colors for arbutus plotting
default_colours_arbutus <- function()
    c("#252A30", "orangered1")

## functions for making histogram oulines
hist_outline <- function(h, ..., density=TRUE) {
  xy <- hist_xy(h, density)
  lines(xy, ...)
}
hist_fill <- function(h, ..., density=TRUE) {
  xy <- hist_xy(h, density)
  polygon(xy, ...)
}

hist_xy <- function(h, density=TRUE) {
  dx <- diff(h$mids[1:2])
  xx <- rep(with(h, c(mids - dx/2, mids[length(mids)] + 
                      dx/2)), each = 2)
  yy <- c(0, rep(if (density) h$density else h$counts, each = 2), 0)
  list(x=xx, y=yy)
}

    
        

        


  
