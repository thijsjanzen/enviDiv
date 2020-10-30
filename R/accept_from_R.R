#' calculate summary statistics for a focal tree. The following summary
#'   statistics are calculated:
#'   nLTT, gamma, mean branching times and number of lineages
#' @param emp_stats statistics to compare against
#' @param newick  vector with newick strings
#' @param threshold maximum threshold value
#' @param sd recorded standard deviation in the population
#' @param emp_tree empirical tree, necessary for nLTT calculation
#' @return vector with true or false for each tree
#' @export
accept_from_r <- function(emp_stats,
                          newick,
                          threshold,
                          sd,
                          emp_tree,
                          foreach_blocksize = 100,
                          num_threads = -1) {

  `%dopar%` <- foreach::`%dopar%`

  results <- rep(0, length(newick))

  num_cl <- num_threads
  if (num_threads == -1) num_cl <- parallel::detectCores()

  cl <- parallel::makeForkCluster(num_cl)
  doParallel::registerDoParallel(cl)

  # now we split everything up across threads:
  index_matrix <- split_into_blocks(m = length(newick),
                                    block_size = foreach_blocksize)

  index_matrix <- tibble::as_tibble(index_matrix)

  do_analysis <- function(newick_strings, indices_matrix, i) {
    output <- list()
    cnt <- 1
    start <- indices_matrix$lower[[i]]
    end   <- indices_matrix$upper[[i]]
    for (j in start:end) {
      phy <- ape::read.tree(text = newick_strings[j])
      output[[cnt]] <- accept_this_tree(phy, emp_stats, threshold, sd, emp_tree)
      cnt <- cnt + 1
    }
    return(output)
  }

  indices <- seq_along(index_matrix$upper)
  results <- foreach::foreach(i = indices)  %dopar% {
    do_analysis(newick, index_matrix, i)
  }
  parallel::stopCluster(cl)

  return(unlist(results))
}

#' @keywords internal
get_stats_in_order <- function(focal_tree, emp_tree) {
  output_stats <- c()
  for (i in 1:8) {
    if (i < 8) {
      output_stats[i] <- enviDiv::calc_stat(focal_tree, i, emp_tree)
    } else {
      add <- enviDiv::calc_stat(focal_tree, i, emp_tree)
      output_stats <- c(output_stats, add)
    }
  }
  return(output_stats)
}


#' @keywords internal
accept_this_tree <- function(phy, emp_stats, threshold, sd, emp_tree) {
  phy <- ape::multi2di(phy)

  for (i in 1:8) {
    phy_stat <- tryCatch({
                          calc_stat(phy, i, emp_tree)
                        }, error = function(cond) {
                          return(1e20)
                        })

    for (j in seq_len(phy_stat)) {
      fit <- abs(phy_stat[j]  - emp_stats[i + j - 1]) / sd[i + j - 1]

      if (fit > threshold) return(FALSE)
    }
  }
  return(TRUE)
}

#' calculate statistics
#' @param focal_tree focal tree
#' @param index index
#' @param brts_emp_tree branching times emperical tree
#' @export
calc_stat <- function(focal_tree, index, emp_tree) {

  if (index == 1) {  # gamma
    return(ape::gammaStat(focal_tree))
  }

  if (index == 2) { # cherries
    return(phyloTop::cherries(focal_tree))
  }

  if (index == 3) {   # pitchforks
    return(phyloTop::pitchforks(focal_tree))
  }

  if (index == 4) {  # ladder
    return(phyloTop::avgLadder(focal_tree))
  }

  if (index == 5) { # colless
    return(apTreeshape::colless(apTreeshape::as.treeshape(focal_tree)))
  }

  if (index == 6) { # we assume the emp tree is in the global environment
    return(nLTT::nltt_diff(focal_tree, emp_tree))
  }

  if (index == 7) {
    # beta
    return(apTreeshape::maxlik.betasplit(focal_tree)$max_lik)
  }

  if (index == 8) {
    spectral <- RPANDA::spectR(focal_tree)
    output <- c()
    output[1] <- spectral$principal_eigenvalue
    output[2] <- spectral$asymmetry
    output[3] <- spectral$peakedness
    return(output)
  }
}
