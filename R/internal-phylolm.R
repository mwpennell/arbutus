## functions for pulling apart phylolm objects

model_type.phylolm <- function(fit, ...){
    model <- fit$model

    ## treat random and fixed root OU model as the same for now
    if (model == "OUrandomRoot" | model == "OUfixedRoot")
        model <- "OU"

    model
}


## It does not appear that the original data and phy is included
## Must include them separately for now...

model_data.phylolm <- function(fit, phy, ...){
    if (missing(phy))
      stop("phylolm fitted objects do not include the phylogeny. This needs to be entered separately with the argument 'phy'.")
    res <- fit$residuals[,1]
    list(phy=phy, data=res)
}


model_pars.phylolm <- function(fit, ...){
    sigsq <- fit$sigma2

    ## match model to parameters
    model <- model_type(fit)

    if (model == "BM"){
        pars <- list(sigsq=sigsq, z0=NA, SE=0)
    } else {
        pars <- list(sigsq=sigsq, par=fit$optpar, z0=NA, SE=0)

        mm <- c("OU", "EB", "lambda", "kappa", "delta")
        par <- c("alpha", "a", "lambda", "kappa", "delta")
        tmp <- par[which(mm == model)]
        names(pars)[names(pars) == "par"] <- tmp
    }
    pars
}

#' @method model_info phylolm
#' @export
model_info.phylolm <- function(fit, phy, ...){
    m <- list(data=model_data(fit, phy),
              pars=model_pars(fit),
              type=model_type(fit))
    class(m) <- "fitC"
    m
}
