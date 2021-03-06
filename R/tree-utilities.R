## these are utility fxns for extracting information from trees

## get node heights (this is equivalent to the internal geiger fxn heights.phylo() )

edge.height <- function(phy){
	phy <- reorder(phy, "postorder")
	n <- length(phy$tip.label)
	n.node <- phy$Nnode
	xx <- numeric(n + n.node)
	for (i in nrow(phy$edge):1) xx[phy$edge[i, 2]] <- xx[phy$edge[i, 1]] + phy$edge.length[i]
	root <- ifelse(is.null(phy$root.edge), 0, phy$root.edge)
	labs <- c(phy$tip.label, phy$node.label)
	depth <- max(xx)
	tt <- depth - xx
	idx <- 1:length(tt)
	dd <- phy$edge.length[idx]
	mm <- match(1:length(tt), c(phy$edge[, 2], Ntip(phy) + 1))
	dd <- c(phy$edge.length, root)[mm]
	ss <- tt + dd
	res <- cbind(ss, tt)
	rownames(res) <- idx
	colnames(res) <- c("start", "end")
	data.frame(res)
}


check.tree.data <- function(phy, data, sort=FALSE, warnings=TRUE) {
  if (missing(data))
		stop("If a 'phylo' or 'multiPhylo' object is supplied, 'data' must be included as well")
	if (!is.null(dim(data)))
    stop("Multidimensional data")
  if (!identical(names(data), phy$tip.label))
    stop("Trait data and species do not align")
}

## geiger and diversitree's drop.tip function
prune.phylo <- function(phy, tip, trim.internal = TRUE, subtree = FALSE, root.edge = 0, rooted = is.rooted(phy)){


  if(missing(tip)) return(phy)
  if (is.character(tip)) tip <- which(phy$tip.label %in% tip)
  if(!length(tip)) return(phy)

  phy=as.phylo(phy)
  Ntip <- length(phy$tip.label)
  tip=tip[tip%in%c(1:Ntip)]
  if(!length(tip)) return(phy)


  phy <- reorder(phy)
  NEWROOT <- ROOT <- Ntip + 1
  Nnode <- phy$Nnode
  Nedge <- nrow(phy$edge)

  wbl <- !is.null(phy$edge.length)
  edge1 <- phy$edge[, 1]
  edge2 <- phy$edge[, 2]
  keep <- !(edge2 %in% tip)

  ints <- edge2 > Ntip
  repeat {
    sel <- !(edge2 %in% edge1[keep]) & ints & keep
    if (!sum(sel))
      break
    keep[sel] <- FALSE
  }

  phy2 <- phy
  phy2$edge <- phy2$edge[keep, ]
  if (wbl)
    phy2$edge.length <- phy2$edge.length[keep]
  TERMS <- !(phy2$edge[, 2] %in% phy2$edge[, 1])
  oldNo.ofNewTips <- phy2$edge[TERMS, 2]
  n <- length(oldNo.ofNewTips)
  idx.old <- phy2$edge[TERMS, 2]
  phy2$edge[TERMS, 2] <- rank(phy2$edge[TERMS, 2])
  phy2$tip.label <- phy2$tip.label[-tip]
  if (!is.null(phy2$node.label))
    phy2$node.label <-
      phy2$node.label[sort(unique(phy2$edge[, 1])) - Ntip]
  phy2$Nnode <- nrow(phy2$edge) - n + 1L
  i <- phy2$edge > n
  phy2$edge[i] <- match(phy2$edge[i], sort(unique(phy2$edge[i]))) + n
  storage.mode(phy2$edge) <- "integer"
  collapse.singles(phy2)
}




check.names.phylo <- function(phy, data, data.names = NULL) {

	if (is.null(data.names)) {
		if (is.vector(data)) {
			data.names <- names(data);
		} else {
			data.names <- rownames(data);
		}
	}
	t <- phy$tip.label;
	r1 <- t[is.na(match(t, data.names))];
	r2 <- data.names[is.na(match(data.names, t))];

	r <- list(sort(r1), sort(r2));

	names(r) <- cbind("tree_not_data", "data_not_tree")
	if (length(r1) == 0 && length(r2) == 0) {
		return("OK");
	} else {
		return(r);
	}
}
