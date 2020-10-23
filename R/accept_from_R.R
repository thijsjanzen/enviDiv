#' calculate summary statistics for a focal tree. The following summary
#'   statistics are calculated:
#'   nLTT, gamma, mean branching times and number of lineages
#' @param emp_stats statistics to compare against
#' @param newick  vector with newick strings
#' @param threshold maximum threshold value
#' @param sd recorded standard deviation in the population
#' @return vector with true or false for each tree
#' @export
accept_from_R <- function(emp_stats,
                          newick,
                          threshold,
                          sd,
                          emp_brts) {

  results <- rep(0, length(newick))
  for(i in 1:length(newick)) {
    phy <- ape::read.tree(text = newick[i])
    results[i] <- accept_this_tree(phy, emp_stats, threshold, sd)
  }
  return(results)
}

#' @keywords internal
accept_this_tree <- function(phy, emp_stats, threshold, sd, emp_brts) {
  phy <- ape::multi2di(phy)

  for (i in 1:length(emp_stats)) {
    phy_stat <- calc_stat(phy, i, emp_brts)

    for(j in 1:length(phy_stat)) {
       fit <- abs(phy_stat[j]  - emp_stats[i + j - 1]) / sd[i + j - 1]

      if (fit > threshold) return(FALSE)
    }
  }
  return(TRUE)
}

#' @keywords internal
calc_stat <- function(focal_tree, index, brts_emp_tree) {

  if (index == 1) {  # gamma
    return(ape::gammaStat(focal_tree))
  }

  if (index == 2) { # cherries
    return(phyloTop::cherries(focal_tree))
  }

  if (index == 3) {   # pitchforks
    return(phyloTop::pitchforks(focal_tree))
  }

  if (index == 4) {  # ladder
    return(phyloTop::avgLadder(focal_tree))
  }

  if(index == 5) { # colless
    return(apTreeshape::colless(apTreeshape::as.treeshape(focal_tree)))
  }

  if (index == 6) { # we assume the emp tree is in the global environment
    # nltt
    brts_focal_tree <- ape::branching.times(focal_tree)
    lineages_emp_tree <- 2:length(brts_emp_tree)
    lineages_focal_tree <- 2:length(brts_focal_tree)

    return(nLTT::nltt_diff_exact_brts(brts_focal_tree, lineages_focal_tree,
                                      brts_emp_tree, lineages_emp_tree))
  }

  if (index == 7) {
    # beta
    return(apTreeshape::maxlik.betasplit(focal_tree)$max_lik)
  }

  if (index == 8) {
    spectral <- RPANDA::spectR(focal_tree)
    output <- c()
    output[1] <- spectral$principal_eigenvalue
    output[2] <- spectral$asymmetry
    output[3] <- spectral$peakedness
    return(output)
  }
}

