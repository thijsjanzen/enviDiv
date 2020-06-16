#' infer parameter values
#' @param reference_table reference table
#' @param focal_model focal model
#' @param to_analyze matrix with summary statistics of tree to infer paramaters for
#' @return tibble with parameter estimates
#' @export
infer_params_abcrf <- function(reference_table,
                               focal_model,
                               to_analyze) {

  model_data <- subset(reference_table,
                       reference_table$model == focal_model)


  sum_stat_names <- c("nltt", "gamma", "mbr", "num_lin",
                      "beta", "colless", "sackin", "ladder", "cherries",
                      "ILnumber", "pitchforks", "stairs",
                      "spectr_eigen", "spectr_asymmetry", "spectr_peakedness")

  data3 <- dplyr::select(model_data, c("sym_high", sum_stat_names))
  sym_high <- data3$sym_high
  pred_params <- abcrf::regAbcrf(sym_high~., data3, paral = TRUE)
  vv <- stats::predict(pred_params, to_analyze, data3, paral = TRUE)
  sym_estimates <- as.vector(vv$expectation)

  data3 <- dplyr::select(model_data, c("extinct", sum_stat_names))
  extinct <- data3$extinct
  pred_params <- abcrf::regAbcrf(extinct~., data3, paral = TRUE)
  vv <- stats::predict(pred_params, to_analyze, data3, paral = TRUE)
  extinct_estimates <- as.vector(vv$expectation)

  data3 <- dplyr::select(model_data, c("sym_low", sum_stat_names))
  sym_low <- data3$sym_low
  pred_params <- abcrf::regAbcrf(sym_low~., data3, paral = TRUE)
  vv <- stats::predict(pred_params, to_analyze, data3, paral = TRUE)
  symlow_estimates <- as.vector(vv$expectation)


  data3 <- dplyr::select(model_data, c("allo", sum_stat_names))
  allo <- data3$allo
  pred_params <- abcrf::regAbcrf(allo~., data3, paral = TRUE)
  vv <- stats::predict(pred_params, to_analyze, data3, paral = TRUE)
  allo_estimates <- as.vector(vv$expectation)


  data3 <- dplyr::select(model_data, c("jiggle", sum_stat_names))
  jiggle <- data3$jiggle
  pred_params <- abcrf::regAbcrf(jiggle~., data3, paral = TRUE)
  vv <- stats::predict(pred_params, to_analyze, data3, paral = TRUE)
  jiggle_estimates <- as.vector(vv$expectation)

  output <- cbind(extinct_estimates, sym_estimates, symlow_estimates,
                  allo_estimates, jiggle_estimates)

  colnames(output) <- c("Extinction",
                        "Sympatric_high_water",
                        "Sympatric_low_water",
                        "Allopatric",
                        "Posterior_perturbation")

  output <- tibble::as_tibble(output)
  return(output)
}
