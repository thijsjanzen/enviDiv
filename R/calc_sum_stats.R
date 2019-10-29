#' @keywords internal
beta_stat <- function(tree) {
  b <- apTreeshape::maxlik.betasplit(tree)$max_lik
  return(b)
}

#' @keywords internal
colless_stat <- function(tree)  {
  if (tree$Nnode + 1 < 3) return(Inf) # 2 nodes doesn't work here

  x <- apTreeshape::colless(apTreeshape::as.treeshape(tree), norm = NULL)
  return(x)
}

#' @keywords internal
sackin_stat <- function(tree)  {
  if (tree$Nnode + 1 < 3) return(Inf) # 2 nodes doesn't work here

  x <- apTreeshape::sackin(apTreeshape::as.treeshape(tree), norm = NULL)
  return(x)
}




#' calculate summary statistics for a focal tree. The following summary
#'   statistics are calculated:
#'   nLTT, gamma, mean branching times and number of lineages
#' @param focal_tree input phylogenetic tree for which to calculate
#'                   summary statistics.
#' @param emp_tree reference empirical tree, only used to calculate the nLTT.
#' @return vector of 8 summary statistics
#' @export
calc_sum_stats <- function(focal_tree, emp_tree = NULL) {
  # trees that are extinct turn up as NULL:
  if (is.null(focal_tree) || class(focal_tree) != "phylo") {
    warning("no tree found, returning infinite statistics")
    return(rep(Inf, 15))
  }

  if (is.null(emp_tree)) {
    emp_tree <- focal_tree
  }

  if (!ape::is.rooted(focal_tree)) {
    warning("tree was not rooted, returning infinite statistics")
    return(rep(Inf, 15))
  }

  focal_tree <- ape::multi2di(focal_tree)

  if (min(ape::branching.times(focal_tree), na.rm = TRUE) < 0) {
    deviation <- min(ape::branching.times(focal_tree))
    warning("the root was < 0, error made was: ", deviation,
            " probably a rounding error\n")
    return(rep(Inf, 15))
  }

  output <- c()
  output[1] <- nLTT::nltt_diff_exact(focal_tree, emp_tree)
  output[2] <- ape::gammaStat(focal_tree) # gamma
  output[3] <- mean(focal_tree$edge.length, na.rm = TRUE) # mean branch length
  output[4] <- focal_tree$Nnode + 1 # number of lineages
  output[5] <- beta_stat(focal_tree) # beta
  output[6] <- colless_stat(focal_tree) # colless
  output[7] <- sackin_stat(focal_tree) # sackin
  output[8] <- phyloTop::avgLadder(focal_tree) # ladder
  output[9] <- phyloTop::cherries(focal_tree) # cherries
  output[10] <- phyloTop::ILnumber(focal_tree) # number of internal nodes with a single child
  output[11] <- phyloTop::pitchforks(focal_tree) #pitchforks
  output[12] <- phyloTop::stairs(focal_tree) # stair casedness

  spectral <- RPANDA::spectR(focal_tree)
  output[13] <- spectral$principal_eigenvalue
  output[14] <- spectral$asymmetry
  output[15] <- spectral$peakedness

  return(output)
}
