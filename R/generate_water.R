#' Generate a sequence of waterlevel changes, depending on the model
#' @param water_model can be 1) no water level changes,
#'                           2) literature water level changes
#'                          3) literature water level changes +
#'                              exponentially distributed water level changes
#' @param maximum_time crown age
#' @return vector of water level changes in [0, maximum_time], where 0
#'         indicates the start of the tree (e.g. the root)
#' @export
generate_water <- function(water_model,
                           maximum_time) {

  if (maximum_time < 0) {
    warning("maximum time has to be > 0")
    return()
  }

  output <- get_waterlevel_cpp(water_model, maximum_time)
  return(output)
}
