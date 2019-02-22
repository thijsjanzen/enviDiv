#' Function to create a phy object tree
#' @param params parameter vector
#' @param water_changes vector with parameter changes
#' @param crown_age crown age of tree to be simulated
#' @param seed pseudo random number generator seed
#' @return phy object
#' @export
create_envidiv_tree <- function(params,
                        water_changes,
                        crown_age,
                        seed = NULL) {

  if(length(params) != 4) {
    stop("params should be a vector of four entries")
  }

  if(is.null(seed)) seed = Sys.time()

  local_newick_string <- create_tree_cpp(params,
                                         water_changes,
                                         seed,
                                         crown_age)

  phy_tree <- ape::read.tree(text = local_newick_string)
  return(phy_tree)
}

#' Function to create a newick string
#' @param params parameter vector
#' @param water_changes vector with parameter changes
#' @param crown_age crown age of tree to be simulated
#' @param seed pseudo random number generator seed
#' @return newick string
#' @export
raw_tree <- function(params,
                        water_changes,
                        crown_age,
                        seed = NULL)  {
  local_newick_string <- create_tree_cpp(params,
                                         water_changes,
                                         seed,
                                         crown_age)
  return(local_newick_string)
}
