## plotting functions for arbutus summary statistics


#' @title Plot contrasts versus their variances
#'
#' @description Plot the contrasts derived from a rescaled 'unit.tree' against
#' their expected variances.
#'
#' @param x a model fitted object or \code{phylo} object that can be passed to \code{\link{as.unit.tree}}.
#'
#' @param col colours used for plotting
#'
#' @param ... additional arguments to be passed to \code{\link{as.unit.tree}}
#'
#' @details This function generates a plot of the observed contrasts versus their
#' expected variances (see \code{\link[ape]{pic}}). The slope of this plot is equal to
#' the summary statistic \code{\link{cor.contrast.var}}
#' and is used as a default summary statistic in \code{\link{summ.stats}}.
#'
#' If the model is a good fit to the data (i.e. an adequate model),
#' the slope ~ 0. The purpose of the plot is to get a quick visual glance
#' at the data and to assess whether or not the model inadequacy may be driven by outliers.
#'
#' The function can currently only take a single dataset (does not integrate across multiple trees).
#'
#' \code{contrast.var.plot} uses the \code{\link[ggplot2]{ggplot2}} library.
#'
#' @seealso \code{\link{as.unit.tree}}, \code{\link{summ.stats}}, \code{\link{cor.contrast.var}}
#'
#' @export contrast.var.plot
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
#'                                  control=list(niter=10))
#'
#' ## plot contrasts versus their variances
#' contrast.var.plot(fit.bm)
#' }
#' 
contrast.var.plot <- function(x, col=NULL, ...){

    ## build unit tree from model object
    ut <- as.unit.tree(x, ...)

    ## create data frame for pics
    dat <- as.data.frame(ut$pics)
    dat$contrasts <- abs(dat$contrasts)

    ## if no colors supply, set defaults
    if (is.null(col))
        col <- def.colours.arbutus()

    if (length(col) > 2)
        warning("Only first two colours used")

    .e <- environment()
    p <- ggplot(dat, aes(variance, contrasts), environment=.e)
    p <- p + stat_smooth(method="lm", fill=col[1], colour=col[1], size=2, alpha=0.5)
    p <- p + geom_point(colour=col[2], fill=col[2])
    p <- p + theme_bw()
    p <- p + theme(plot.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border=element_blank(),
                   axis.line = element_line(color = 'black'))
    p <- p + xlab("E[Var] of contrast")
    p <- p + ylab("Contrast")
    p <- p + ggtitle("Variation relative to branch lengths")

    print(p)
}















#' @title Plot contrasts versus their inferred ancestral states
#'
#' @description Plot the contrasts derived from a rescaled 'unit.tree' against
#' their inferred ancestral state.
#'
#' @param x a model fitted object or \code{phylo} object that can be passed to \code{\link{as.unit.tree}}.
#'
#' @param col colours used for plotting
#'
#' @param ... additional arguments to be passed to \code{\link{as.unit.tree}}
#'
#' @details This function generates a plot of the observed contrasts versus their
#' inferred ancestral state (see \code{\link[ape]{pic}}, \code{\link[ape]{ace}}).
#' Note that the ancestral states are estimated using the 'pic' method, which
#' differs somewhat from the ML estimate of the ancestral states.
#' The slope of this plot is equal to
#' the summary statistic \code{\link{cor.contrast.asr}}
#' and is used as a default summary statistic in \code{\link{summ.stats}}. 
#'
#' If the model is a good fit to the data (i.e. an adequate model),
#' the slope ~ 0. The purpose of the plot is to get a quick visual glance
#' at the data and to assess whether or not the model inadequacy may be driven by outliers.
#'
#' The function can currently only take a single dataset (does not integrate across multiple trees).
#'
#' \code{contrast.asr.plot} uses the \code{\link[ggplot2]{ggplot2}} library.
#'
#' @seealso \code{\link{as.unit.tree}}, \code{\link{summ.stats}}, \code{\link{cor.contrast.asr}}
#'
#' @export contrast.asr.plot
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
#'                                  control=list(niter=10))
#'
#' ## plot contrasts versus their variances
#' contrast.asr.plot(fit.bm)
#' }
#' 
contrast.asr.plot <- function(x, col=NULL, ...){

    ## build unit tree from model object
    ut <- as.unit.tree(x, ...)

    ## extract pics
    pics <- abs(ut$pics[,"contrasts"])

    ## get tree and data, calculate ancestral states
    asr <- ace(ut$data, ut$phy, method="pic")$ace

    ## create a dataframe
    dat <- cbind.data.frame(asr, pics)

    ## if no colors supply, set defaults
    if (is.null(col))
        col <- def.colours.arbutus()

    if (length(col) > 2)
        warning("Only first two colours used")


    .e <- environment()
    p <- ggplot(dat, aes(asr, pics), environment=.e)
    p <- p + stat_smooth(method="lm", fill=col[1], colour=col[1], size=2, alpha=0.5)
    p <- p + geom_point(fill=col[2], colour=col[2])
    p <- p + theme_bw()
    p <- p + theme(plot.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border=element_blank(),
                   axis.line = element_line(color = 'black'))    
    p <- p + xlab("Ancestral state")
    p <- p + ylab("Contrast")
    p <- p + ggtitle("Variation relative to trait value")

    print(p)
}












#' @title Plot contrasts versus their node height
#'
#' @description Plot the contrasts derived from a rescaled 'unit.tree' against
#' the node height (i.e. depth) at which they are computed
#'
#' @param x a model fitted object or \code{phylo} object that can be passed to \code{\link{as.unit.tree}}.
#'
#' @param col colours used for plotting
#'
#' @param ... additional arguments to be passed to \code{\link{as.unit.tree}}
#'
#' @details This function generates a plot of the observed contrasts versus their
#' the node height at which they are computed (see \code{\link[ape]{pic}}).
#' The slope of this plot is equal to
#' the summary statistic \code{\link{cor.contrast.nh}}
#' and is used as a default summary statistic in \code{\link{summ.stats}}. 
#'
#' If the model is a good fit to the data (i.e. an adequate model),
#' the slope ~ 0. The purpose of the plot is to get a quick visual glance
#' at the data and to assess whether or not the model inadequacy may be driven by outliers.
#'
#' The function can currently only take a single dataset (does not integrate across multiple trees).
#'
#' \code{contrast.nh.plot} uses the \code{\link[ggplot2]{ggplot2}} library.
#'
#' @seealso \code{\link{as.unit.tree}}, \code{\link{summ.stats}}, \code{\link{cor.contrast.nh}}
#'
#' @export contrast.nh.plot
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
#'                                  control=list(niter=10))
#'
#' ## plot contrasts versus their variances
#' contrast.nh.plot(fit.bm)
#' }
#' 
contrast.nh.plot <- function(x, col=NULL, ...){

    ## build unit tree from model object
    ut <- as.unit.tree(x, ...)

    ## extract pics
    pics <- abs(ut$pics[,"contrasts"])

    ht <- edge.height(ut$phy)
    N <- Ntip(ut$phy)
    nh <- ht$end[c((N + 1):nrow(ht))]
    names(nh) <- rownames(ht[c((N + 1):nrow(ht)),])

    ## create a dataframe
    dat <- cbind.data.frame(nh, pics)

    ## if no colors supply, set defaults
    if (is.null(col))
        col <- def.colours.arbutus()

    if (length(col) > 2)
        warning("Only first two colours used")

    .e <- environment()
    p <- ggplot(dat, aes(nh, pics), environment=.e)
    p <- p + stat_smooth(method="lm", fill=col[1], colour=col[1], size=2, alpha=0.5)
    p <- p + geom_point(fill=col[2], color=col[2])
    p <- p + theme_bw()
    p <- p + xlab("Node height")
    p <- p + ylab("Contrast")
    p <- p + ggtitle("Variation relative to time")

    print(p)
}




#' @title Plot summary statistics for model adequacy
#'
#' @description Plots the distribution of simulated summary statistics
#' along with the observed summary statistics (either a point estimate or
#' a distribution of values) for visual inspection of model adequacy
#'
#' @param x a 'phy.ss' object inherited from \code{\link{phy.model.check}}
#' (or, alternatively \code{\link{compare.summ.stats}})
#'
#' @param col a vector of two colours to be used in the plot; the first is for the observed summary
#' statistics and the second for the simulated summary statistics.
#'
#' @param ... additional arguments to be passed to plot (not currently functional)
#'
#' @details This function is designed to give researchers a quick visual inspection of
#' model adequacy by plotting the simulated distribution of summary statistics along
#' with the empirical estimate(s). If only a point estimated of the observed summary statistics
#' is input (such as from estimating the parameters on a single tree using maximum likelihood),
#' the function will plot the simulated summary statistics with a vertical line representing the
#' observed summary statistic. If a distribution of observed summary statistics is used, the function
#' will plot both distributions (simulated and observed) on the same plot. This function is not
#' designed to be flexible. If researchers want to create their own plots in a different style, all
#' the necessary information is available in the 'phy.ss' object output by both \code{\link{phy.model.check}}
#' and \code{\link{compare.summ.stats}} and can be indexed by \code{x$summ.stats.obs} and
#' \code{x$summ.stats.sim}.
#'
#' This plotting function using the \code{\link[ggplot2]{ggplot2}} library.
#'
#' @method plot phy.ss
#' @export
#'
#' @seealso \code{\link{phy.model.check}}, \code{\link{compare.summ.stats}}
#'
#' @examples
#' data(finch)
#' phy <- finch$phy
#' dat <- finch$data[,"wingL"]
#' 
#' p <- phy.model.check(phy, data=dat, nsim=20)
#'
#' plot(p)
#' 
plot.phy.ss <- function(x, col=NULL, ...){

    .e <- environment()

    ## if no colors supply, set defaults
    if (is.null(col))
        col <- def.colours.arbutus()

    if (length(col) > 2)
        warning("Only first two colours used")

    ## bind stats together
    type <- rep("observed", (ncol(x$summ.stats.obs) * nrow(x$summ.stats.obs)))
    obs <- arbutus.melt(x$summ.stats.obs)
    obs <- cbind.data.frame(obs,type)
    type <- rep("simulated", (ncol(x$summ.stats.sim) * nrow(x$summ.stats.sim)))
    sim <- arbutus.melt(x$summ.stats.sim)
    sim <- cbind.data.frame(sim, type)
    stat <- rbind(obs, sim)

    ## identify whether the observed are from dist or not
    if (nrow(x$summ.stats.obs) == 1){ ## histogram
        
        p <- ggplot(stat, aes(x=value), environment = .e)
        p <- p + geom_histogram(data=subset(stat, type=="simulated"),
                                fill=col[1])
        p <- p + geom_vline(data=subset(stat, type=="observed"),
                            aes(xintercept=value), color=col[2], size=1)
        p <- p + facet_wrap(~variable, scales="free")
        p <- p + theme_bw()
        p <- p + theme(strip.background=element_rect(fill="white"),
                       plot.background=element_blank(),
                       panel.grid.major=element_blank(),
                       panel.grid.minor=element_blank(),
                       axis.ticks.y=element_blank(),
                       axis.text.y=element_blank())
        p <- p + xlab("Value of test statistic")
        p <- p + ylab("Density")
        
    } else { ## distribution
        p <- ggplot(stat, aes(x=value, fill=type), environment = .e)
        p <- p + geom_histogram(alpha=0.6, position="identity")
        p <- p + scale_fill_manual(values=col)
        p <- p + facet_wrap(~variable, scales="free")
        p <- p + theme_bw()
        p <- p + theme(strip.background=element_rect(fill="white"),
                       plot.background=element_blank(),
                       panel.grid.major=element_blank(),
                       panel.grid.minor=element_blank(),
                       axis.ticks.y=element_blank(),
                       axis.text.y=element_blank())
        p <- p + xlab("Value of test statistic")
        p <- p + ylab("Density")
    }

    suppressMessages(print(p))
}


## Set default colors for arbutus plotting
def.colours.arbutus <- function()
    c("#f5634a","#ed303c")

## Little function for melting phy.ss objects for plotting
arbutus.melt <- function(df){
    ## get the values
    val <- lapply(seq_len(ncol(df)), function(x) return(df[,x]))
    value <- do.call(c, val)

    ## get the variable
    v <- lapply(seq_len(ncol(df)), function(x) {rep(colnames(df)[x], nrow(df))})
    variable <- do.call(c, v)

    ## bind them together
    out <- cbind.data.frame(variable, value)

    ## factor the variable
    out$variable <- factor(out$variable, levels=colnames(df))

    ## return data frame
    out
}

    




                       
    
    
        

        


  
