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
#' @param max_lin maximum number of extant lineages in the tree.
#' @param abc (boolean) is the tree simulated in an ABC fitting scheme,
#'                      or not? additional verbal output is provided if not.
#' @param seed random nmber seed
#' @return phy object
#' @export
sim_envidiv_tree <- function(params,
                             crown_age,
                             max_lin = 500,
                             abc = FALSE,
                             seed = NULL) {

  if (crown_age < 0) {
    warning("crown age should be larger than zero\n")
    return(NULL)
  }

  water_changes <- generate_water(params[7], crown_age)

  sim_result <- create_tree_cpp(params,
                                water_changes,
                                crown_age,
                                max_lin)

  error_code <- sim_result$code

  if (error_code == "extinction") {
    if (!abc) warning("Tree went extinct, returning NULL\n")
    return(NULL)
  }

  if (error_code == "overflow") {
    if (!abc) warning("Tree too big, returning NULL")
    return(NULL)
  }

  phy_tree <- enviDiv::sim_table_to_phy(sim_result$Ltable,
                                        crown_age)

  return(phy_tree)
}

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
#' @param max_lin maximum number of extant lineages in the tree.
#' @param abc (boolean) is the tree simulated in an ABC fitting scheme,
#'                      or not? additional verbal output is provided if not.
#' @param seed random nmber seed
#' @return phy object
#' @export
sim_envidiv_tree2 <- function(params,
                             crown_age,
                             max_lin = 500,
                             abc = FALSE,
                             seed = NULL) {

  if (crown_age < 0) {
    warning("crown age should be larger than zero\n")
    return(NULL)
  }

  water_changes <- generate_water(params[7], crown_age)

  sim_result <- create_tree_cpp(params,
                                water_changes,
                                crown_age,
                                max_lin)

  error_code <- sim_result$code

  if (error_code == "extinction") {
    if (!abc) warning("Tree went extinct, returning NULL\n")
    return(NULL)
  }

  if (error_code == "overflow") {
    if (!abc) warning("Tree too big, returning NULL")
    return(NULL)
  }

  phy_tree <- enviDiv::sim_table_to_phy(sim_result$Ltable,
                                        crown_age)

  return(list("phy" = phy_tree,
              "water" = water_changes,
              "error_code" = error_code,
              "ltable" = sim_result$Ltable))
}

#' correct ltable
#' @param input_matrix ltable
#' @param crown_age age
#' @keywords export
correct_ltable_internal <- function(input_matrix, crown_age) {
  local_l_table <- input_matrix
  local_l_table <- local_l_table[order(local_l_table[, 1]), 1:4]
  local_l_table[, 1] <- crown_age - local_l_table[, 1]
  not_min <- which(local_l_table[, 4] != -1)
  local_l_table[not_min, 4] <- crown_age - local_l_table[not_min, 4]
  #local_l_table <-  local_l_table[order(abs(local_l_table[, 3])), 1:4]

  # time to rename things:
  local_ltable2 <- local_l_table
  for (i in 1:nrow(local_l_table)) {
    old_id <- local_l_table[i, 3]
    old_parent <- local_l_table[i, 2]
    if (i == 1) {
      local_ltable2[i, 2] <- 0
      local_ltable2[i, 3] <- -1

      new_daughters <- which(local_l_table[, 2] == old_id)
      new_daughters <- new_daughters[new_daughters > i]
      local_ltable2[new_daughters, 2] <- -1
    } else if (i == 2) {
      local_ltable2[i, 2] <- -1
      local_ltable2[i, 3] <- 2

      new_daughters <- which(local_l_table[, 2] == old_id)
      new_daughters <- new_daughters[new_daughters > i]
      local_ltable2[new_daughters, 2] <- 2
    } else {
      parent <- local_ltable2[i, 2]
      new_id <- i * sign(parent)

      new_daughters <- which(local_l_table[, 2] == old_id)
      new_daughters <- new_daughters[new_daughters > i]
      local_ltable2[new_daughters, 2] <- new_id

      # update self:
      local_ltable2[i, 3] <- new_id
    }
  }
  return(local_ltable2)
}

#' function to transform an Ltable to a phylogeny
#' @description transform an Ltable to a phy object
#' @param input_matrix input matrix
#' @param crown_age crown age
#' @export
sim_table_to_phy <- function(input_matrix,
                             crown_age) {

  local_l_table <- correct_ltable_internal(input_matrix, crown_age)

  #phy_tree <- DDD::L2phylo(local_l_table)
  phy_tree <- treestats::l_to_phylo(local_l_table, drop_extinct = TRUE)


  if (length(ape::is.binary(phy_tree)) > 1) {
    new_phy_tree <- ape::collapse.singles(phy_tree)
    if (length(ape::is.binary(new_phy_tree)) > 1) {
      stop("could not generate binary tree\n")
    }
    phy_tree <- new_phy_tree
  }

  if (length(treestats::branching_times(phy_tree)) !=
      (-1 + length(phy_tree$tip.label))) {
    new_phy_tree <- ape::collapse.singles(phy_tree)
    if (length(treestats::branching_times(new_phy_tree)) !=
        (-1 + length(new_phy_tree$tip.label))) {
      stop("could not generate tree without singles\n")
    }
    phy_tree <- new_phy_tree
  }

  return(phy_tree)
}


#' function to transform an Ltable to a newick string. Calls
#' sim_table_to_phy and then transforms to newick
#' @description transform an Ltable to a newick string
#' @param input_matrix input matrix
#' @param crown_age crown age
#' @export
sim_table_to_newick <- function(input_matrix,
                                crown_age) {
  phy_tree <- sim_table_to_phy(input_matrix, crown_age)
  return(ape::write.tree(phy_tree))
}
