#' calculate summary statistics for a focal tree. The following summary
#'   statistics are calculated:
#'   nLTT, gamma, mean branching times and number of lineages
#' @param focal_tree input phylogenetic tree for which to calculate
#'                   summary statistics.
#' @return vector of 70 summary statistics
#' @export
calc_sum_stats <- function(focal_tree) {
  # trees that are extinct turn up as NULL:
  if (is.null(focal_tree) || class(focal_tree) != "phylo") {
    warning("no tree found, returning infinite statistics")
    return(rep(Inf, 70))
  }

  if (!ape::is.rooted(focal_tree)) {
    warning("tree was not rooted, returning infinite statistics")
    return(rep(Inf, 70))
  }

  focal_tree <- ape::multi2di(focal_tree)

  if (min(ape::branching.times(focal_tree), na.rm = TRUE) < 0) {
    deviation <- min(ape::branching.times(focal_tree))
    warning("the root was < 0, error made was: ", deviation,
            " probably a rounding error\n")
    return(rep(Inf, 70))
  }

  output <- unlist(treestats::calc_all_stats(focal_tree))

  return(output)
}
