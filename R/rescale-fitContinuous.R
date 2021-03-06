## fxns for rescaling tree according to fitContinuous type models
## adopted from rescale.phylo in GEIGER (written by Jon Eastman)

model_phylo_bm <- function(phy, pars){
    if (pars$sigsq < 0)
        stop("Parameters need to be non-negative")
    phy$edge.length <- phy$edge.length * pars$sigsq
    model_phylo_se(phy, pars)
}

model_phylo_ou <- function(phy, pars){
    if (pars$sigsq < 0 || pars$alpha < 0)
        stop("Parameters need to be non-negative")

    ## rescale branch lengths according to alpha...
    alpha <- pars$alpha
    ## Note that this is more similar to the rescaling in diversitree
    ## than to geiger.
    if (alpha > 0) {
      tt <- branching.times.by.depth(phy)
      t0  <- tt$Tmax - tt$t2
      len <- tt$t2 - tt$t1
      bl <- -exp(-2 * alpha * t0) * expm1(-2 * alpha * len) / (2*alpha)
      phy$edge.length <- bl[phy$edge[,2]]
    }

    ## ...and then by sigsq
    phy$edge.length <- phy$edge.length * pars$sigsq
    model_phylo_se(phy, pars)
}

model_phylo_eb <- function(phy, pars){
    if (pars$sigsq < 0)
        stop("Parameters need to be non-negative")

    ## rescale tree according to "a" parameter...
    a <- pars$a
    if (a != 0){
      tt <- branching.times.by.depth(phy)
      bl <- (exp(a * tt$t2) - exp(a * tt$t1)) / a
      phy$edge.length <- bl[phy$edge[,2]]
    }

    ## ...and then by sigsq
    phy$edge.length <- phy$edge.length * pars$sigsq

    ## return phy
    model_phylo_se(phy, pars)
}


## lambda transformation
## NOTE [RGF]: This function seems too complicated; is there a simpler
## way of writing it?
model_phylo_lambda <- function(phy, pars){
    if (pars$sigsq < 0)
        stop("Parameters need to be non-negative")

    heights <- edge.height(phy)

    ## index edge matrix by tips
    N <- Ntip(phy)
    tips <- match(seq_len(N), phy$edge[,2])

    ## get all path lengths
    pp <- heights$start[N + 1] - heights$end[seq_len(N)]

    ## rescale all branches according to lambda value
    bl <- phy$edge.length * pars$lambda

    ## readjust the tip branches
    bl[tips] <- bl[tips] + (pp - (pp * pars$lambda))
    phy$edge.length <- bl

    ## check to make sure no negative branch lengths produced
    if (any(phy$edge.length < 0))
        warning("negative branch lengths in tranformed tree. Lambda may be too large")

    ## rescale branch lengths according to sigsq
    phy$edge.length <- phy$edge.length * pars$sigsq

    model_phylo_se(phy, pars)
}

model_phylo_kappa <- function(phy, pars){
  if (pars$sigsq < 0 || pars$kappa < 0)
    stop("Parameters need to be non-negative")
  phy$edge.length <- (phy$edge.length^pars$kappa) * pars$sigsq
  model_phylo_se(phy, pars)
}

## NOTE [RGF]: I skipped the tidy-up here, because there is a good
## chance that this doesn't work (see issue #58).
model_phylo_delta <- function(phy, pars){
    if (pars$sigsq < 0 || pars$delta < 0)
        stop("Parameters need to be non-negative")

    ## get all heights
    heights <- edge.height(phy)

    ## get tree depth
    ## same as max(branching.times(phy)) but works for non-ultrametric trees
    N <- Ntip(phy)
    Tmax <- heights$start[N+1]

    ## get braching times in terms of tree depth
    ## note: Geiger's delta transformation does not index the branches (like OU)
    ## need to check if this is correct
    t <- Tmax - heights$start
    l <- heights$start - heights$end

    ## rescale the branches according to delta
    delta <- pars$delta

    bl <- (t + l)^delta - t^delta
    phy$edge.length <- bl[phy$edge[,2]]

    ## Rescaling using geiger's method
    scl <- Tmax^delta
    phy$edge.length <- phy$edge.length/scl * Tmax

    ## rescale branch lengths according to sigsq
    sigsq <- pars$sigsq

    phy$edge.length <- phy$edge.length * sigsq

    ## return phy
    model_phylo_se(phy, pars)
}


## See issue #59
model_phylo_trend <- function(phy, pars)
    stop("trend model is not currently implemented")

## NOTE: White noise transformation turns tree into star phylogeny.
model_phylo_white <- function(phy, pars){
  if (pars$sigsq < 0)
    stop("Parameters need to be non-negative")
  is.tip <- phy$edge[,2] <= Ntip(phy)
  phy$edge.length[!is.tip] <- 0
  phy$edge.length[ is.tip] <- pars$sigsq
  model_phylo_se(phy, pars)
}


#' @method make_model_phylo fitC
#' @export
make_model_phylo.fitC <- function(x, ...) {
    ## get model
    model <- x$type

    ## get tree
    phy <- x$data$phy

    ## get parameters
    pars <- x$pars

    ## Translation function; all have argument list (phy, pars)
    tr <- switch(model,
                 BM=model_phylo_bm,
                 OU=model_phylo_ou,
                 EB=model_phylo_eb,
                 lambda=model_phylo_lambda,
                 kappa=model_phylo_kappa,
                 delta=model_phylo_delta,
                 trend=model_phylo_trend,
                 white=model_phylo_white)

    if (is.data.frame(pars)) {
      rphy <- lapply(seq_len(nrow(pars)), function(i)
                     tr(phy, pars[i,]))
      class(rphy) <- "multiPhylo"
    } else {
      rphy <- tr(phy, pars)
    }

    ## return rescaled phylogeny
    rphy
}

model_phylo_rescale <- function(type) {
  switch(type,
         BM=model_phylo_bm,
         OU=model_phylo_ou,
         EB=model_phylo_eb,
         lambda=model_phylo_lambda,
         kappa=model_phylo_kappa,
         delta=model_phylo_delta,
         trend=model_phylo_trend,
         white=model_phylo_white,
         stop("Unknown model type ", type))
}

#' @title Rescale phylogeny based on fitted model parameters
#'
#' @description Rescales phylogeny to form a 'unit.tree' using parameters
#' from fitted model
#'
#' @param x an object inherited from a \code{model_info} function
#' @param ... additional arguments to be passed to \code{make_model_phylo}
#'
#' @details This is a generic function which rescales the phylogeny based on the model
#' specific information. While the class and information may differ, the object must include
#' the 'phylo' object to be rescaled and model specific information such as model type and
#' parameter values.
#'
#' To include additional types of models, researchers will need to build a new \code{\link{model_info}}
#' function for the given model type. The output should be assigned a novel class.
#' A \code{make_model_phylo.class} object will need to be created to perform the rescaling.
#'
#' @return a \code{phylo} object that can be used to from a \code{unit.tree} object
#'
#' @seealso \code{\link{make_unit_tree}}
#'
#' @export make_model_phylo
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
#' ## get model info using internal arbutus function
#' info.bm <- model_info(fit.bm)
#'
#' ## rescale phylogeny based on model info
#' make_model_phylo(info.bm)
#' }
#'
make_model_phylo <- function(x, ...)
    UseMethod("make_model_phylo")


model_phylo_se <- function(phy, pars) {
  if (pars$SE < 0)
    stop("SE must be non-negative")
  tips <- phy$edge[,2] <= Ntip(phy)
  phy$edge.length[tips] <- phy$edge.length[tips] + pars$SE^2
  phy
}

## Branching times in terms of tree depth
branching.times.by.depth <- function(phy) {
  N <- Ntip(phy)
  heights <- edge.height(phy)

  ## Total tree depth; same as as max(branching.times(phy)) but
  ## works for non-ultrametric trees
  Tmax <- heights$start[N+1]

  idx <- match(1:nrow(heights), phy$edge[,2])
  t1 <- Tmax - heights$end[phy$edge[idx,1]]
  t2 <- heights$start - heights$end + t1
  list(t1=t1, t2=t2, Tmax=Tmax)
}

## Total tree depth; same as as max(branching.times(phy)) but
## works for non-ultrametric trees
tree.depth <- function(phy) {
  N <- Ntip(phy)
  heights <- edge.height(phy)

  ## Total tree depth; same as as max(branching.times(phy)) but
  ## works for non-ultrametric trees
  heights$start[N+1]
}
