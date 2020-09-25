#' collect reference table based on previously generated reference tables for
#' each model
#' @param path path of directory where each subdirectory [1, 2, 3] is located
#' @param return_full_table should the function only return summary stats, or
#' also parameters?
#' @return reference table
#' @export
read_reference_table <- function(path,
                                 return_full_table = FALSE) {

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


    reference_table <- rbind(reference_table, vy)
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
