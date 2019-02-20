## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end


create_tree <- function(params,
                        water_changes,
                        crown_age,
                        seed = NULL) {

  if(length(params) != 4) {
    stop("params should be a vector of four entries")
  }

  if(is.null(seed)) seed = sys.time()

  crown_age <- 5
  params <- c(0, 1, 0, 0)
  water_changes = c(0, 2*crown_age)


  local_newick_string = create_tree(params,
                                    water_changes,
                                    seed,
                                    crown_age)


}

