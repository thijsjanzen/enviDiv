#' function to perform ABC-SMC
#' @param number_of_trees number of particles used per iteration of
#'                            the SMC algorithm
#' @param min_tips minimum number of tips
#' @param max_tips maximum number of tips
#' @param model used water model
#' @param crown_age crown age
#' @param write_to_file boolean, if TRUE, results are written to file.
#' @param num_threads number of threads
#' @return a tibble containing the results
#' @export
generate_trees_tbb <- function(number_of_trees = 1000,
                                min_tips = 50,
                                max_tips = 150,
                                model = NULL,
                                crown_age = NULL,
                                file_name= NULL,
                                num_threads = -1) {

  if (is.null(crown_age)) {
    stop("Please either provide a reference tree, or provide the crown age")
  }

  sim_result <- enviDiv::create_ref_table_tbb_par(model = model,
                                                    num_repl = number_of_trees,
                                                    crown_age = crown_age,
                                                    min_lin = min_tips,
                                                    max_lin = max_tips,
                                                    num_threads = num_threads)

  convert_to_phylo <- function(newick_string) {
    phylo_tree <- ape::read.tree(text = newick_string)
    return(phylo_tree)
  }

  phylo_trees <- lapply(sim_result$trees, convert_to_phylo)

  # now we calculate stats
  stats <- future.apply::future_lapply(phylo_trees, calc_sum_stats)

  stat_matrix <- matrix(unlist(stats, use.names = FALSE),
                        ncol = 15,
                        byrow = TRUE)

  results <- cbind(sim_result$parameters, stat_matrix)

  colnames(results) <-
    c("extinct", "sym_high", "sym_low", "allo", "jiggle", "model",
      "nltt", "gamma", "mbr", "num_lin",
      "beta", "colless", "sackin", "ladder", "cherries", "ILnumber",
      "pitchforks", "stairs",
      "spectr_eigen", "spectr_asymmetry", "spectr_peakedness")

  results <- tibble::as_tibble(results)

  readr::write_tsv(results, path = file_name,
                   append = TRUE)

}
