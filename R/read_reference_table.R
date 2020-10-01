#' collect reference table based on previously generated reference tables for
#' each model
#' @param path path of directory where each subdirectory [1, 2, 3] is located
#' @param return_full_table should the function only return summary stats, or
#' also parameters?
#' @param num_rows number of samples to be taken
#' @param random_rows should these rows be sampled randomly?
#' @return reference table
#' @export
read_reference_table <- function(path,
                                 return_full_table = FALSE,
                                 num_rows = NA,
                                 random_rows = TRUE) {

  reference_table <- c()
  for (i in 1:3) {
    file_name <- paste0(path, i, "/reference_table.txt")
    vy <- readr::read_tsv(file_name, col_names = F)
    colnames(vy)  <-
      c("extinct", "sym_high", "sym_low", "allo", "jiggle", "model",
        "nltt", "gamma", "mbr", "num_lin",
        "beta", "colless", "sackin", "ladder", "cherries", "ILnumber",
        "pitchforks", "stairs",
        "spectr_eigen", "spectr_asymmetry", "spectr_peakedness")

    if (!is.na(num_rows)) {
      vy <- sample_subset(vy, num_rows, random_rows)
    }

    reference_table <- rbind(reference_table, vy)
    rm(vy)
  }

  focal_models <- as.factor(reference_table$model)

  if (!return_full_table) {

    ref_table <- dplyr::select(reference_table,
                               c("nltt", "gamma", "mbr", "num_lin",
                                 "beta", "colless", "sackin", "ladder",
                                 "cherries", "ILnumber", "pitchforks", "stairs",
                                 "spectr_eigen", "spectr_asymmetry",
                                 "spectr_peakedness"))

    return(cbind(focal_models, ref_table))
  } else {
    return(reference_table)
  }
}

#' @keywords internal
sample_subset <- function(vy,
                          num_rows,
                          random_rows) {
  output <- c()
  remaining_indices <- seq_along(vy$extinct)
  num_selected <- 0

  while (length(remaining_indices) > 0 && num_selected < num_rows) {
    num_remaining <- num_rows - num_selected
    indices <- sample(remaining_indices, num_remaining)
    remaining_indices <- remaining_indices[-match(indices, remaining_indices)]

    vz <- vy[indices, ]
    to_remove <- c()
    for (i in seq_along(vz$extinct)) {
      a <- as.numeric(vz[i, ])
      if (sum(is.infinite(a)) || sum(is.na(a)) || sum(is.nan(a))) {
        to_remove <- c(to_remove, i)
      }
    }
    if (length(to_remove) > 0) {
      vz <- vz[-to_remove, ]
    }
    output <- rbind(output, vz)
    num_selected <- num_selected + length(vz$extinct)
  }
  return(output)
}
