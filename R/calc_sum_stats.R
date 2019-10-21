#' @keywords internal
beta_stat <- function(tree) {
  if (!ape::is.rooted(tree)) return(Inf)
  if (!ape::is.binary.tree(tree)) {
    tree <- ape::multi2di(tree, random = TRUE)
  }

  b <- apTreeshape::maxlik.betasplit(tree)$max_lik
  return(b)
}

#' @keywords internal
colless_stat <- function(tree)  {
  if (!ape::is.rooted(tree)) return(Inf)
  if (tree$Nnode + 1 < 3) return(Inf) # 2 nodes doesn't work here
  if (!ape::is.binary.tree(tree)) {
    tree <- ape::multi2di(tree, random = TRUE)
  }
  x <- apTreeshape::colless(apTreeshape::as.treeshape(tree), norm = NULL)
  return(x)
}

#' @keywords internal
sackin_stat <- function(tree)  {
  if (!ape::is.rooted(tree)) return(NA)
  if (tree$Nnode + 1 < 3) return(Inf) # 2 nodes doesn't work here
  if (!ape::is.binary.tree(tree)) {
    tree <- ape::multi2di(tree, random = TRUE)
  }
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
    return(rep(Inf, 8))
  }

  if (is.null(emp_tree)) {
    emp_tree <- focal_tree
  }

  focal_tree <- ape::multi2di(focal_tree)

  if (min(ape::branching.times(focal_tree), na.rm = TRUE) < 0) {
    deviation <- min(ape::branching.times(focal_tree))
    warning("the root was < 0, error made was: ", deviation,
            " probably a rounding error\n")
    if (abs(deviation) > 0.01) {
      warning("The error was huge, no attempt tried to correct,
              aborting this tree\n")
      return(rep(Inf, 8))
    }
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

  return(output)
}
