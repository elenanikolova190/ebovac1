#' Sigmoid function
#'
#'
#' @description Sigmoid function to account for the effect of control measures such as the arrival
#' of protective equipment, community awareness of the ebola virus disease, etc.
#'
#' @param time time in days
#' @param alpha the shape of change (default value )
#' @param shift the midpoint of the time of change (default value )
#' @param reduct_tr the reduction in transmission rate (default 1.5)
#'
#' @return the reduction in transmission rate after the introduction of control intervention
#'
#' @examples sigmoid(time = 20, shape = 0.8, shift = 1, reduct_tr = 0.09),
#'
#'

sigmoid <- function (time, shape = 1, shift = 0, reduct_tr = 0.05)
{
  if (length(time) == 0){
    return(c())
  }
  stopifnot(is.numeric(time), is.numeric(shape), is.numeric(shift), is.numeric(reduct_tr))

  shape <- shape[1]
  shift <- shift[1]

  print(exp(-shape * (time - shift)))

  return(reduct_tr/(1 + exp(-shape * (time - shift))))
}


#' Linear function
#'
#'
#' @description The Linear time varying function to account for the effect of control measures such as the arrival
#' of protective equipment, community awareness of the Ebola virus disease, etc.
#'
#' @param time time in days
#' @param slope the shape of change (default value )
#' @param reduct_tr the reduction in transmission rate (default 1.5)
#'
#' @return the reduction in transmission rate after the introduction of control intervention
#'
#' @examples linear_t_var(time = 20, slope = 0.01),
#'

linear_t_var <- function (time, slope = 0.01)
{
  if (length(time) == 0){
    return(c())
  }

  stopifnot(is.numeric(time), is.numeric(slope))

  slope <- slope[1]

  return(1-(slope * time))

}

#' calc_beta()
#'
#' Calculate beta
#'
#' @description Calculates the transmission rate, beta, depending on the selected by the user option:
#'                 * linear  (when "sel" of the input parameters is set to 1)
#'                * sigmoid  (when "sel" of the input parameters is set to 2)
#'                * constant
#'
#' @param time time in days
#' @param theta a list of parameters
#'
#' @return the reduction in transmission rate after the introduction of control intervention according to the
#'          user selection.
#'
#' @examples calc_beta(time, theta)
#'
#' What is important here is to set the "sel" value of the input parameters tp 1, 2 or else in order to select the
#' typr of change in transmision rate.
#' In case sigmoid is selected the values that need to be set are the shape and midpoint, while for a linear change
#' in transmission tate only the shape (alpha) is needed
#'
calc_beta <- function(time, theta){
  if (length(time) == 0){
    return(c())
  }

  sel <- theta["sel"]
  alpha <- theta["alpha"]
  shift <- theta["shift"]
  reduct_tr <- theta["reduct_tr"]

  result <- 1

  if ( sel == 1 ){
    result <- linear_t_var(time, alpha)
  } else if (sel == 2) {
    result <- sigmoid(time, alpha, shift, reduct_tr)
  }

  return(result)
}

