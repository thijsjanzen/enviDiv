#' function wrapper around treestats
#' @export
#' @return list of statistics functions
create_statistics_list <- function() {

  stats <- list()
  stats$rquartet <- treestats::rquartet # moved up to trigger early fail
  stats$colless_quad <- treestats::colless_quad # moved up to trigger early fail
  stats$wiener <- treestats::wiener # moved up to trigger early fail
  stats$gamma <- treestats::gamma_statistic
  stats$sackin <- treestats::sackin
  stats$colless <- treestats::colless
  stats$colless_corr <- treestats::colless_corr
  stats$beta <- treestats::beta_statistic
  stats$blum <- treestats::blum
  stats$pigot_rho <- treestats::pigot_rho

  stats$treeness <- treestats::treeness
  stats$nltt_base <- treestats::nLTT_base
  stats$phylogenetic_div <- treestats::phylogenetic_diversity
  stats$avg_ladder <- treestats::avg_ladder
  stats$max_ladder <- treestats::max_ladder
  stats$cherries <- treestats::cherries
  stats$double_cherries <- treestats::double_cherries
  stats$four_prong <- treestats::four_prong
  stats$il_number <- treestats::ILnumber
  stats$pitchforks <- treestats::pitchforks

  stats$stairs <- treestats::stairs
  stats$imbalance_steps <- treestats::imbalance_steps
  stats$j_one <- treestats::j_one
  stats$b1 <- treestats::b1
  stats$b2 <- treestats::b2
  stats$area_per_pair <- treestats::area_per_pair
  stats$average_leaf_depth <- treestats::average_leaf_depth
  stats$i_stat <- treestats::mean_i
  stats$ew_colless <- treestats::ew_colless
  stats$max_del_width <- treestats::max_del_width

  stats$max_depth <- treestats::max_depth
  stats$avg_vert_depth <- treestats::avg_vert_depth
  stats$max_width <- treestats::max_width
  stats$mw_over_md <- treestats::mw_over_md
  stats$tot_path <- treestats::tot_path_length
  stats$tot_internal_path <- treestats::tot_internal_path
  stats$rogers <- treestats::rogers
  stats$stairs2 <- treestats::stairs2
  stats$tot_coph <- treestats::tot_coph
  stats$var_depth <- treestats::var_leaf_depth

  stats$symmetry_nodes <- treestats::sym_nodes
  stats$mpd <- treestats::mean_pair_dist
  stats$psv <- treestats::psv
  stats$vpd <- treestats::var_pair_dist
  stats$mntd <- treestats::mntd
  stats$j_stat <- treestats::entropy_j
  stats$crown_age <- treestats::crown_age
  stats$tree_height <- treestats::tree_height
  stats$max_betweenness <- treestats::max_betweenness
  stats$diameter <- treestats::diameter

  local_closeness <- function(tree, w, n) {
    return(treestats::max_closeness(tree, weight = w, normalization = ifelse(n ==
                                                                               TRUE, "tips", "none")))
  }
  stats$max_closeness <- function(x) {
    local_closeness(x, FALSE, FALSE)
  }
  stats$max_closenessW <- function(x) {
    local_closeness(x, TRUE, FALSE)
  }

  stats$eigen_centrality <- function(x) {
    return(max(treestats::eigen_centrality(x, weight = FALSE)$eigenvector))
  }
  stats$eigen_centralityW <- function(x) {
    return(max(treestats::eigen_centrality(x, weight = TRUE)$eigenvector))
  }
  stats$mean_branch_length <- treestats::mean_branch_length
  stats$var_branch_length <- treestats::var_branch_length
  stats$mean_branch_length_int <- treestats::mean_branch_length_int
  stats$mean_branch_length_ext <- treestats::mean_branch_length_ext
  stats$var_branch_length_int <- treestats::var_branch_length_int
  stats$var_branch_length_ext <- treestats::var_branch_length_ext


  stats$root_imbalance <- treestats::root_imbalance
  stats$number_of_lineages <- treestats::number_of_lineages

  get_minmax_lapl <- function(phylo) {
    out <- list()
    temp_stats <- treestats::minmax_laplace(phylo, TRUE)

    if (length(temp_stats) >= 2) {
      out$min_laplace <- temp_stats$min
      out$max_laplace <- temp_stats$max
    }
    else {
      out$min_laplace <- NA
      out$max_laplace <- NA
    }
    return(out)
  }
  stats$minmax_lapl <- get_minmax_lapl

  get_minmax_adj <- function(phylo) {
    out <- list()
    temp_stats <- treestats::minmax_adj(phylo, TRUE)

    if (length(temp_stats) >= 2) {
      out$min_adj <- temp_stats$min
      out$max_adj <- temp_stats$max
    }
    else {
      out$min_adj <- NA
      out$max_adj <- NA
    }
    return(out)
  }
  stats$minmax_adj <- get_minmax_adj

  get_laplacian_dist <- function(phylo) {
    temp_stats <- tryCatch(expr = {
      treestats::laplacian_spectrum(phylo)
    }, error = function(e) {
      return(NA)
    })

    out <- list()
    if (length(temp_stats) == 5) {
      out$laplace_spectrum_a <- temp_stats$asymmetry
      out$laplace_spectrum_p <- temp_stats$peakedness
      out$laplace_spectrum_e <- log(temp_stats$principal_eigenvalue)
      out$laplace_spectrum_g <- temp_stats$eigengap[[1]]
    }
    else {
      out$laplace_spectrum_a <- NA
      out$laplace_spectrum_p <- NA
      out$laplace_spectrum_e <- NA
      out$laplace_spectrum_g <- NA
    }
    return(out)
  }
  stats$laplacian_d <- get_laplacian_dist

  return(stats)
}

#' return a vector of statistics names used
#' @export
#' @return vector of names
names_statistics_list <- function() {
  basic_names <- create_statistics_list()
  basic_names <- names(basic_names)
  # last three are compound
  av <- which(basic_names == "minmax_lapl")
  basic_names <- basic_names[-av]
  basic_names <- c(basic_names, c("min_lapl", "max_lapl"))

  av <- which(basic_names == "minmax_adj")
  basic_names <- basic_names[-av]
  basic_names <- c(basic_names, c("min_adjl", "max_adjl"))

  av <- which(basic_names == "laplacian_d")
  basic_names <- basic_names[-av]
  basic_names <- c(basic_names, c("laplace_spectrum_a",
                                  "laplace_spectrum_p",
                                  "laplace_spectrum_e",
                                  "laplace_spectrum_g"))
  return(basic_names)
}
