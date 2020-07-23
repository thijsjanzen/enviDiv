#' function to perform ABC-SMC
#' @param number_of_trees number of trees to generate
#' @param min_tips minimum number of tips (inclusive)
#' @param max_tips maximum number of tips (inclusive)
#' @param model used water model
#' @param crown_age crown age
#' @param file_name file name to write output to
#' @param num_threads number of threads
#' @return a matrix with parameter values and the associated summary statistics
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

  cat("simulating trees is done\n")

  # now we calculate stats
  cat("calculating summary statistics for all trees...\n")

  indices <- seq_along(phylo_trees)

  progressr::with_progress({
    p <- progressr::progressor(along = phylo_trees)
    stats <- future.apply::future_lapply(indices, function(x, ...) {
      p(sprintf("x=%g", x))
      calc_sum_stats(phylo_trees[[x]])
    })
  })

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

  readr::write_tsv(results, path = file_name)
  cat(paste("reference table written to:", file_name))
  return(results)
}
