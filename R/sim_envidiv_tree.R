#' simulate a tree using the environmental diversification model
#' @param params parameters used to simulate:
#' \itemize{
#'   \item{extinction}{per lineage extinction rate}
#'   \item{sympatric speciation rate at high water}{per lineage rate of
#'   speciation when the water level is high}
#'   \item{sympatric speciation rate at low water}{per lineage rate of
#'   speciation when the water level is low}
#'   \item{allopatric speciation rate}{per allopatric pair rate of speciation}
#'   \item{perturbation}{standard deviation of post-hoc
#'                       branching time perturbation}
#'   \item{water model}{Water model: 1) no water level changes, 2) literature
#'                     water level change, 3) extrapolated water level changes}
#' }
#' @param crown_age age of the crown of the tree
#' @param abc (boolean) is the tree simulated in an ABC fitting scheme,
#'                      or not? additional verbal output is provided if not.
#' @param seed random nmber seed
#' @return phy object
#' @rawNamespace useDynLib(enviDiv)
#' @rawNamespace import(Rcpp)
#' @export
sim_envidiv_tree <- function(params,
                             crown_age,
                             abc = FALSE,
                             seed = NULL) {

  if(crown_age < 0) {
    warning("crown age should be larger than zero\n")
    return(NULL)
  }

  if(is.null(seed)) seed <- as.numeric(Sys.time())

  water_changes <- generate_water(params[6], crown_age)
  local_newick_string <- create_tree_cpp(params,
                                         water_changes,
                                         seed,
                                         crown_age)

  if (local_newick_string == "extinction") {
    if (!abc) warning("Tree went extinct, returning NULL\n")
    return(NULL)
  }

  if (local_newick_string == "overflow") {
    if (!abc) warning("Tree too big, returning NULL")
    return(NULL)
  }

  phy_tree <- phytools::read.newick(text = local_newick_string)

  if (is.null(phy_tree))  {
    warning("phy tree is NULL")
    return(NULL)
  }
  if (is.null(phy_tree$edge.length)) {
    warning("phy_tree$edge.length == NULL")
    return(NULL)
  }

  if (!"phylo" %in% class(phy_tree)) {
    warning("phy is not of class phylo")
    return(NULL)
  }

  if (length(phy_tree$tip.label) == 2) {
    if (!abc) warning("tree has only two tips\n")
    if (abc) return(NULL)
  }

  if (length(ape::is.binary(phy_tree)) > 1) {
    new_phy_tree <- ape::collapse.singles(phy_tree)
    if (length(ape::is.binary(phy_tree)) > 1) {
      cat(params, "\n")
      cat(local_newick_string, "\n")
      cat("ERROR, could not generate binary tree\n")
      stop()
    }
    phy_tree <- new_phy_tree
  }

  if (length(ape::branching.times(phy_tree)) !=
     (-1 + length(phy_tree$tip.label))) {
    new_phy_tree <- ape::collapse.singles(phy_tree)
    if (length(ape::branching.times(new_phy_tree)) !=
       (-1 + length(new_phy_tree$tip.label))) {
      cat(params, "\n")
      cat(local_newick_string, "\n")
      cat("ERROR, could not generate tree without singles\n")
      stop()
    }
    phy_tree <- new_phy_tree
  }

  return(phy_tree)
}
