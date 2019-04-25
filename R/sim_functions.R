#' simulate a tree using the environmental diversification model
#' @param params parameters used to simulate
#' @param crown_age age of the crown of the tree
#' @return phy object
#' @rawNamespace useDynLib(enviDiv)
#' @export
sim_envidiv_tree <- function(params,
                             crown_age) {

  seed = as.numeric(Sys.time())

  water_changes <- generate_water(params[6], crown_age)
  local_newick_string <- create_tree_cpp(params,
                                         water_changes,
                                         seed,
                                         crown_age)

  if(local_newick_string == "extinction") {
    return(NULL)
  }

  if(local_newick_string == "overflow") {
    return(NULL)
  }

  phy_tree <- phytools::read.newick(text = local_newick_string)

  valid_tree <- TRUE
  if(is.null(phy_tree)) valid_tree <- FALSE
  if(is.null(phy_tree$edge.length)) valid_tree <- FALSE
  if (!"phylo" %in% class(phy_tree)) valid_tree <- FALSE
  if(length(phy_tree$tip.label) == 2) valid_tree <- FALSE

  if(!valid_tree) return(NULL)

  output_tree <- phy_tree
  graphics::plot(output_tree)
  if(length(geiger::is.extinct(output_tree)) > 0)  {
    output_tree <- geiger::drop.extinct(phy_tree, tol = 0.01)
  }

  #output_tree <- ape::di2multi(output_tree)
  return(output_tree)
}

#' generate a sequence of waterlevel changes, depending on the model
#' @param water_model can be 1) no water level changes, 2) literature water level changes 3) literature water level changes + exponentially distributed water level changes
#' @param maximum_time crown age
#' @return vector of water level changes in [0, maximum_time], where 0 indicates the start of the tree (e.g. the root)
generate_water <- function(water_model,
                           maximum_time) {
  if(water_model == 1) return(c(0, maximum_time * 2))
  if(water_model == 2) {
    output <- c()
    output[1] = (maximum_time - 1.1);			#0
    output[2] = (maximum_time - 0.55);		#1
    output[3] = (maximum_time - 0.393);		#0
    output[4] = (maximum_time - 0.363);		#1
    output[5] = (maximum_time - 0.295);		#0
    output[6] = (maximum_time - 0.262);		#1
    output[7] = (maximum_time - 0.193);		#0
    output[8] = (maximum_time - 0.169);		#1
    output[9] = (maximum_time - 0.04);		#0
    output[10] = (maximum_time - 0.035);	#1
    output[11] = (maximum_time);				  #1
    return(output)
  }
  if(water_model == 3) {
    output2 <- c()
    waterLevel = 1;
    time = 0;
    tempTime = time;
    while(tempTime < (maximum_time - 1.1)) {
      tempTime = tempTime + stats::rexp(1, 10)
      if(tempTime > (maximum_time - 1.1)) break;
      time = tempTime;
      output2 <- c(output2, time)
      waterLevel = 1 - waterLevel;
    }
    if(waterLevel == 0) {
      output2 <- output2[-length(output2)] #water level has to be high
    }

    output <- c()
    output[1] = (maximum_time - 1.1);			#0
    output[2] = (maximum_time - 0.55);		#1
    output[3] = (maximum_time - 0.393);		#0
    output[4] = (maximum_time - 0.363);		#1
    output[5] = (maximum_time - 0.295);		#0
    output[6] = (maximum_time - 0.262);		#1
    output[7] = (maximum_time - 0.193);		#0
    output[8] = (maximum_time - 0.169);		#1
    output[9] = (maximum_time - 0.04);		#0
    output[10] = (maximum_time - 0.035);	#1
    output[11] = (maximum_time);				  #1

    combined_output <- c(output2, output)
    return(combined_output)
  }
}
