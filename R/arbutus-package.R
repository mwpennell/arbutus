#' Arbutus: Evaluate model adequacy for phylogenetic models of continuous character evolution
#'
#' @description Arbutus was developed to evaluate the adequacy (or absolute goodness-of-fit)
#' of phylogenetic models of continuous character evolution. The basic principle
#' underlying the approach is that if the generating model is a Brownian motion
#' process, the contrasts (sensu Felsenstein, 1985) will be independent and
#' identically distributed (I.I.D.). We can evaluate this condition by calculating
#' a set of test statistics the contrasts on our observed data, then simulating datasets under
#' a Brownian motion process along the phylogeny and calculating the test statistics
#' on the contrasts of each simulated data set. We can then compare our observed test statistics to the
#' simulated distribution of test statistics.
#'
#' The I.I.D. property of the contrasts will hold only for the case of Brownian motion.
#' However, the approach can be extended to any arbitrarily complex model of continuous character
#' evolution (as long as the model is based on a multivariate normal distribution) by first
#' rescaling the phylogeny based on the fitted parameter estimates, creating what we refer to
#' as a 'unit.tree'. Once the rescaling is done, the contrasts will be I.I.D. if the model
#' that was fit is a adequate one. 
#'
#' The basic functionality of arbutus can be employed using the function \code{\link{arbutus}}.
#' This function wraps several other functions which can be used separately:
#' \enumerate{
#'  \item{\code{\link{make_unit_tree}}: }{rescales phylogeny based on fitted parameter values}
#'  \item{\code{\link{calculate_pic_stat}}: }{calculates test statistics on observed data}
#'  \item{\code{\link{simulate_char_unit}}: }{simulates datasets under Brownian motion}
#'  \item{\code{\link{calculate_pic_stat}}: }{calculates test statistics on simulated data}
#'  \item{\code{\link{compare_pic_stat}}: }{compares observed to simulated statistics}
#' }
#'
#' Arbutus is designed to interact with a number of other packages for fitting trait models. The
#' architecture is flexible such that additional packages can be incorporated relatively easily.
#'
#' @references Matthew W. Pennell, Richard G. FitzJohn, William K. Cornwell and Luke J. Harmon
#' (in review.).  Model adequacy and the macroevolution of angiosperm functional traits.
#'
#' @references Joseph Felsenstein (1985). Phylogenies and the comparative method. The American
#' Naturalist 125:1-15.
#' 
#' @docType package
#' @import ape
#' @name arbutus
#' @aliases arbutus-arbutus arbutus
#'
NULL







#' Galapagos finch dataset
#'
#' An example dataset including a phylogeny of 13 Galapagos finches with
#' trait data matched to the phylogeny.
#'
#' \describe{
#'  \item{phy}{phylogeny of Galapagos finches}
#'  \item{data}{measured traits for each species in the phylogeny}
#' }
#'
#' @docType data
#' @keywords datasets
#' @format a list with two elements
#' @name finch
NULL
