#' function to perform ABC-SMC
#' @param number_of_trees number of trees to generate
#' @param min_tips minimum number of tips (inclusive)
#' @param max_tips maximum number of tips (inclusive)
#' @param model used water model
#' @param crown_age crown age
#' @param file_name_trees file name to write trees
#' @param file_name_stats file name to write stats
#' @param num_threads number of threads
#' @param block_size maximum number of trees for which summary statistics are
#' calculated in parallel. Should fit in memory so should not be too big.
#' Default is 10000
#' @return a matrix with parameter values and the associated summary statistics
#' @export
generate_trees_tbb <- function(number_of_trees = 1000,
                               min_tips = 50,
                               max_tips = 150,
                               model = NULL,
                               crown_age = NULL,
                               file_name_trees = "trees.txt",
                               file_name_stats = "stats.txt",
                               num_threads = -1,
                               block_size = 10000) {

  if (is.null(crown_age)) {
    stop("Please either provide a reference tree, or provide the crown age")
  }


  convert_to_phylo <- function(newick_string) {
    phylo_tree <- ape::read.tree(text = newick_string)
    return(phylo_tree)
  }

  num_blocks <- 1 + floor(number_of_trees / block_size)

  num_done <- 0

  start_time <- Sys.time()

  for (i in seq_len(num_blocks)) {
    total_sampled <- i * block_size
    if (total_sampled > number_of_trees) {
      overshoot <- total_sampled - number_of_trees
      block_size <- block_size - overshoot
    }

    sim_result <- enviDiv::create_ref_table_tbb_par(model = model,
                                                    num_repl = block_size,
                                                    crown_age = crown_age,
                                                    min_lin = min_tips,
                                                    max_lin = max_tips,
                                                    num_threads = num_threads)

    phylo_trees <- lapply(sim_result$trees, convert_to_phylo)

    cat("simulating trees is done\n")
    trees_for_writing <- phylo_trees
    class(trees_for_writing) <- "multiPhylo"

    ape::write.tree(trees_for_writing, file_name_trees, append = T)
    rm(trees_for_writing)
    # now we calculate stats
    cat("calculating summary statistics \n")

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

    results <- cbind(sim_result$parameters[indices], stat_matrix)

    colnames(results) <-
      c("extinct", "sym_high", "sym_low", "allo", "jiggle", "model",
        "nltt", "gamma", "mbr", "num_lin",
        "beta", "colless", "sackin", "ladder", "cherries", "ILnumber",
        "pitchforks", "stairs",
        "spectr_eigen", "spectr_asymmetry", "spectr_peakedness")

    results <- tibble::as_tibble(results)

    if(i == 1) {
      readr::write_tsv(results, path = file_name_stats)
    } else {
      readr::write_tsv(results, path = file_name_stats, append = T)
    }
    num_done <- num_done + block_size
    current_time <- Sys.time()
    num_left <- number_of_trees - num_done
    time_per_tree <- (current_time - start_time)[[1]] / num_done
    time_remaining <- time_per_tree * num_left

    cat("now done: ", num_done, "trees\n")
    cat("time remaining: ", time_remaining, "seconds \n")
  }

  cat(paste("reference table written to:", file_name_stats))
  return(results)
}
