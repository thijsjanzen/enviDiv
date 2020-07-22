#' function to perform ABC-SMC
#' @param simulations_per_model number of particles used per iteration of
#'                            the SMC algorithm
#' @param min_tips minimum number of tips
#' @param max_tips maximum number of tips
#' @param model used water model
#' @param emp_tree phy object holding phylogeny of the tree to be fitted on
#' @param crown_age crown age
#' @param write_to_file boolean, if TRUE, results are written to file.
#' @param file_name file name
#' @param exp_prior use exponential or uniform prior
#' @return a tibble containing the results
#' @export
create_reference_table <- function(simulations_per_model = 1000,
                                min_tips = 50,
                                max_tips = 150,
                                model = NULL,
                                emp_tree = NULL,
                                crown_age = NULL,
                                write_to_file = FALSE,
                                file_name,
                                exp_prior = FALSE) {

  if (!is.null(model)) {
    generate_from_prior(simulations_per_model,
                        min_tips,
                        max_tips,
                        model,
                        emp_tree,
                        crown_age,
                        write_to_file,
                        file_name,
                        exp_prior)
  } else {
    for (focal_model in 1:3) {
      generate_from_prior(simulations_per_model,
                          min_tips,
                          max_tips,
                          focal_model,
                          emp_tree,
                          crown_age,
                          write_to_file,
                          file_name,
                          exp_prior)
    }
  }
}
