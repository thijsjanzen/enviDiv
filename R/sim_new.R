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
#'   \item{water rate}{rate of water level change if model = 4 (see below)}
#' }
#' @param model Water model: 1) no water level changes, 2) literature water level change,
#' 3) extrapolated water level changes, 4) using random rates, until literature values.
#' @param crown_age age of the crown of the tree
#' @param max_lin maximum number of extant lineages in the tree.
#' @return phy object
#' @export
sim_envidiv_tree_new <- function(params,
                                 model,
                                 crown_age,
                                 max_lin = 500,
                                 seed = -1) {

  if (crown_age < 0) {
    warning("crown age should be larger than zero\n")
    return(NULL)
  }

  sim_result <- sim_envidiv2_cpp(params,
                                 model,
                                 crown_age,
                                 max_lin,
                                 seed)

  error_code <- sim_result$code

  phy_tree <- NULL
  sim_result$Ltable[, 1] <- crown_age - sim_result$Ltable[, 1]
  not_min1 <- which(sim_result$Ltable[, 4] != -1)
  sim_result$Ltable[not_min1, 4] <- crown_age - sim_result$Ltable[not_min1, 4]

  if (error_code == "done") {
    phy_tree <- treestats::l_to_phylo(sim_result$Ltable,
                                      TRUE)
  }

  return(list("phy" = phy_tree,
              "water" = sim_result$water,
              "error_code" = error_code,
              "ltable" = sim_result$Ltable))
}

#' simulate a tree using the environmental diversification model
#' @param model Water model: 1) no water level changes, 2) literature water level change,
#' 3) extrapolated water level changes, 4) using random rates, until literature values.
#' @param crown_age age of the crown of the tree
#' @param min_lin minimum number of extant lineages in the tree
#' @param max_lin maximum number of extant lineages in the tree.
#' @return phy object
#' @export
sim_envidiv_tree_new_cond <- function(model,
                                      crown_age,
                                      min_lin = 4,
                                      max_lin = 500,
                                      num_tries = 10000) {


  sim_result <- sim_new_cond_cpp(model,
                                 crown_age,
                                 min_lin,
                                 max_lin,
                                 num_tries)

  error_code <- sim_result$code

  phy_tree <- NULL
  sim_result$Ltable[, 1] <- crown_age - sim_result$Ltable[, 1]
  not_min1 <- which(sim_result$Ltable[, 4] != -1)
  sim_result$Ltable[not_min1, 4] <- crown_age - sim_result$Ltable[not_min1, 4]

  if (error_code == "done") {
    phy_tree <- treestats::l_to_phylo(sim_result$Ltable,
                                      TRUE)
  }

  return(list("phy" = phy_tree,
              "water" = sim_result$water,
              "error_code" = error_code,
              "ltable" = sim_result$Ltable,
              "params" = sim_result$parameters))
}
