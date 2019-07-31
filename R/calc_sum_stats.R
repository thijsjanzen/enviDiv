#' calculate summary statistics for a focal tree. The following summary
#'   statistics are calculated:
#'   nLTT, gamma, mean branching times and number of lineages
#' @param focal_tree input phylogenetic tree for which to calculate
#'                   summary statistics.
#' @param emp_tree reference empirical tree, only used to calculate the nLTT.
#' @return vector of 11 summary statistics
#' @export
calc_sum_stats <- function(focal_tree, emp_tree = NULL) {
  # trees that are extinct turn up as NULL:
  if (is.null(focal_tree)) {
    return(rep(Inf, 11))
  }

  if(is.null(emp_tree)) {
    emp_tree <- focal_tree
  }

  # trees with only 2 tips turn up with an
  # empty edgelist (they are only a crown)
  if (is.null(focal_tree$edge.length)) {
    return(rep(Inf, 11))
  }

  focal_tree <- ape::multi2di(focal_tree)

  if (min(ape::branching.times(focal_tree), na.rm = T) < 0) {
    cat("ERROR, negative branch lengths!\n")
    return(rep(Inf, 11))
  }
  if (focal_tree$Nnode + 1 < 3) {
    return(rep(Inf, 11))
  }

  output <- c()
  output[1] <- nLTT::nltt_diff_exact(focal_tree, emp_tree)
  output[2] <- ape::gammaStat(focal_tree) # gamma
  output[3] <- mean(focal_tree$edge.length, na.rm = T) # mean branch length
  output[4] <- focal_tree$Nnode + 1 # number of lineages

  beta_stat <- function(tree) {
    if (!ape::is.rooted(tree)) return(NA)
    if (!ape::is.binary.tree(tree)) {
      dichotomousphylogeny <- ape::multi2di(tree, random = TRUE)
      b <- apTreeshape::maxlik.betasplit(dichotomousphylogeny)$max_lik
      return (b)
    }

    b <- apTreeshape::maxlik.betasplit(tree)$max_lik
    return (b)
  }

  colless_stat <- function(tree)  {
    if (!ape::is.rooted(tree)) return(NA)
    tree <- ape::multi2di(tree)
    x <- apTreeshape::colless(apTreeshape::as.treeshape(tree), norm = NULL)
    return(x)
  }

  sackin_stat <- function(tree)  {
    if (!ape::is.rooted(tree)) return(NA)
    tree <- ape::multi2di(tree)
    x <- apTreeshape::sackin(apTreeshape::as.treeshape(tree), norm = NULL)
    return(x)
  }

  output[5] <- beta_stat(focal_tree) # beta
  output[6] <- colless_stat(focal_tree) # colless
  output[7] <- sackin_stat(focal_tree) # sackin
  output[8] <- phyloTop::avgLadder(focal_tree) # ladder
  output[9] <- phyloTop::cherries(focal_tree) # cherries
  output[10] <- phyloTop::ILnumber(focal_tree) # IL number
  output[11] <- phyloTop::pitchforks(focal_tree) # pitch forks

  return(output)
}
