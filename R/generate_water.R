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
  if(water_model == 1) return(c(0, maximum_time * 2))
  if(water_model == 2) {
    output <- c()
    output[1] = (maximum_time - 1.1);			#0
    output[2] = (maximum_time - 0.55);		#1
    output[3] = (maximum_time - 0.393);		#0
    output[4] = (maximum_time - 0.363);		#1
    output[5] = (maximum_time - 0.295);		#0
    output[6] = (maximum_time - 0.262);		#1
    output[7] = (maximum_time - 0.193);		#0
    output[8] = (maximum_time - 0.169);		#1
    output[9] = (maximum_time - 0.04);		#0
    output[10] = (maximum_time - 0.035);	#1
    output[11] = (maximum_time);				  #1
    return(output)
  }
  if(water_model == 3) {
    output2 <- c()
    waterLevel = 1;
    time = 0;
    tempTime = time;
    while(tempTime < (maximum_time - 1.1)) {
      tempTime = tempTime + stats::rexp(1, 10)
      if(tempTime > (maximum_time - 1.1)) break;
      time = tempTime;
      output2 <- c(output2, time)
      waterLevel = 1 - waterLevel;
    }
    if(waterLevel == 0) {
      output2 <- output2[-length(output2)] #water level has to be high
    }

    output <- c()
    output[1] = (maximum_time - 1.1);			#0
    output[2] = (maximum_time - 0.55);		#1
    output[3] = (maximum_time - 0.393);		#0
    output[4] = (maximum_time - 0.363);		#1
    output[5] = (maximum_time - 0.295);		#0
    output[6] = (maximum_time - 0.262);		#1
    output[7] = (maximum_time - 0.193);		#0
    output[8] = (maximum_time - 0.169);		#1
    output[9] = (maximum_time - 0.04);		#0
    output[10] = (maximum_time - 0.035);	#1
    output[11] = (maximum_time);				  #1

    combined_output <- c(output2, output)
    return(combined_output)
  }
}
