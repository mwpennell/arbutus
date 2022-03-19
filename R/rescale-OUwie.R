#rescale OUwie

#' @method make_model_phylo fitOU
#' @export
make_model_phylo.fitOU <- function(x, ...){
  ## get model
  model <- x$type
  
  ## get tree
  phy <- x$data$phy
  
  ## get parameters
  pars <- x$pars
  
  ## Translation function; all have argument list (phy, pars)
  tr <- switch(model,
               BMS=model_phylo_bms)
  
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

model_phylo_bms <- function(phy, pars){
  if (any(pars$sigsq < 0))
    stop("Parameters need to be non-negative")
  n = length(pars$sigsq)
  
  for(reg in 1:n){
    for(node in 1:phy$Nnode){
      if(phy$node.label[node] == reg){
        index <- which(phy$edge[,1] == node + length(phy$tip.label))
        phy$edge.length[index] <- phy$edge.length[index]# * pars$sigsq[reg]^2
      }
    }
  }
  
  #phy <- ifelse(!is.null(pars$SE), model_phylo_ouse(phy, pars), phy)
  phy    
}


model_phylo_ouse <- function(phy, pars) {
  if (pars$SE < 0)
    stop("SE must be non-negative")
  tips <- phy$edge[,2] <= Ntip(phy)
  phy$edge.length[tips] <- phy$edge.length[tips] + pars$SE^2
  phy
}