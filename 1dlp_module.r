#' Title One dose leaky protection model functions
#'
#'

one_dose_leaky <- function(time, state, theta) {

  # transform theta to model parameters
  N            <- theta[["N"]]            # Affected population size
  r_lat        <- 1/theta[["d_lat"]]      # Latency period
  r_inf        <- 1/theta[["d_inf"]]      # Infectious period
  VE           <- theta[["VE"]]           # Vaccine efficacy
  vac_start    <- theta[["vac_start"]]    # Starting day of trial (counted from beginning of epidemic, e.g. day 30)
  vac_duration <- theta[["vac_duration"]] # How many days the trial will go on
  doses_num    <- theta[["Nv"]]           # Number of vaccines administered in trial
  beta         <- theta[["R0"]]*r_inf     # Power of infection
  kappa        <- 1/theta[["immunity_onset"]] # Rate of onset of immunity

  # State variables

  #General Population
  S <- state[["S"]]
  E <- state[["E"]]
  I <- state[["I"]]
  R <- state[["R"]]

  # Vaccinated participants who had been exposed to virus prior to beginning of trial
  Evac <- state[["Evac"]]
  Ivac <- state[["Ivac"]]
  Rvac <- state[["Rvac"]]

  # Vaccinated arm of the trial before onset of immunity
  Sv <- state[["Sv"]]
  Ev <- state[["Ev"]]
  Iv <- state[["Iv"]]
  Rv <- state[["Rv"]]

  # Vaccinated arm of the trial after onset of immunity
  Vp <- state[["Vp"]]
  Ep <- state[["Ep"]]
  Ip <- state[["Ip"]]
  Rp <- state[["Rp"]]

  # Control arm of trial before onset of immunity
  C1 <- state[["C1"]]
  Ec1 <- state[["Ec1"]]
  Ic1 <- state[["Ic1"]]
  Rc1 <- state[["Rc1"]]

  # Control arm of the trial after onset of immunity
  C2 <- state[["C2"]]
  Ec2 <- state[["Ec2"]]
  Ic2 <- state[["Ic2"]]
  Rc2 <- state[["Rc2"]]

  ## Specify total population
  N = S + E + I + R + Sv + Ev + Iv + Rv + Evac + Ivac + Rvac + Vp + Ep + Ip + Rp + C1 + Ec1 + Ic1 + Rc1 + C2 + Ec2 + Ic2 + Rc2


  # Number of vaccines administered daily
  # if(is.null(vac_rate ) & vac_duration != 0) vac_rate <- vac_num / vac_duration

  # Calculate end day of trial
  vac_end <- vac_duration + vac_start

  # Calculate number of vaccinated & number of non-vaccinated people from the total population
  N_not_vacc <- S+E+I+R
  num_vac    <- (N - N_not_vacc)/2

  if(time < vac_start)      # Before vaccination trial start
  {
    lambda <- beta*I/N

    # SEIR compartments of general population
    dS <-  - lambda * S
    dE <-  lambda * S - r_lat * E
    dI <-  r_lat * E - r_inf*I
    dR <-  r_inf * I

    #Compartments related to the vaccination trial
    dEvac <- 0
    dIvac <- 0
    dRvac <- 0
    dSv <- 0
    dEv <- 0
    dIv <- 0
    dRv <- 0
    dVp <- 0
    dEp <- 0
    dIp <- 0
    dRp <- 0
    dC1 <- 0
    dEc1 <- 0
    dIc1 <- 0
    dRc1 <- 0
    dC2 <- 0
    dEc2 <- 0
    dIc2 <- 0
    dRc2 <- 0
  }
  else
  {
    if((time >= vac_start) & (time < vac_end )) {     # Period of Trial
      r_vac <- (doses_num /(S+E))/ vac_duration
    }
    else{
      r_vac <- 0.0                                   # After vaccination has ended
    }


    lambda <- beta*(I + Ivac + Iv + Ip + Ic1 + Ic2)/N

    dS <-  - lambda * S - (r_vac * S) * 2
    dE <-  lambda * S - r_lat * E - (r_vac * E) * 2
    dI <-  r_lat * E - r_inf * I
    dR <-  r_inf * I

    # Vaccinated arm of the trial before onset of immunity
    dSv <-  r_vac * S - lambda * Sv - kappa * Sv
    dEv <- lambda * Sv - r_lat * Ev
    dIv <- r_lat * Ev - r_inf * Iv
    dRv <- r_inf * Iv

    # Participants in trial who have been exposed prior vaccination
    dEvac <- r_vac * E - r_lat * Evac
    dIvac <- r_lat * Evac - r_inf * Ivac
    dRvac <-  r_inf * Ivac

    # Vaccinated arm of trial after onset of immunity
    if(VE > 0){
      dVp <- -(1 - VE) * lambda * Vp + kappa * Sv
    }
    else dVp <- 0

    dEp <- (1 - VE) * lambda * Vp - r_lat * Ep
    dIp <- r_lat * Ep - r_inf * Ip
    dRp <- r_inf *Ip

    # Control arm of the trial before onset of immunity
    dC1 <- r_vac * S  - lambda * C1 - kappa * C1
    dEc1 <- lambda * C1 - r_lat * Ec1
    dIc1 <- r_lat * Ec1 - r_inf * Ic1
    dRc1 <- r_inf * Ic1

    # Control arm of the trial after onset of immunity
    dC2 <- - lambda * C2 + kappa * C1
    dEc2 <- lambda * C2 - r_lat * Ec2
    dIc2 <- r_lat * Ec2 - r_inf * Ic2
    dRc2 <- r_inf * Ic2

  }

  return(list(c(dS,dE,dI,dR, dEvac, dIvac, dRvac,dSv,dEv,dIv,dRv,dVp,dEp,dIp,dRp,dC1,dEc1,dIc1,dRc1,dC2,dEc2,dIc2,dRc2)))

}


#'
#' @return
#' @export
#'
#' @examples one_dose_leaky_model(theta, init, times)
#'

one_dose_leaky_model <- function(theta, init.state, times) {

  # put incidence at 0 in init.state
  traj <- as.data.frame(ode(init.state, times, one_dose_leaky, theta, method = "ode45"))

  return(traj)

}

#' Function: One Dose Leaky Protection Model Module input
#'
#'
#'
onedlpUI <- function(id) {

}


onedlpServer <- function(id, theta) {


  moduleServer(
    id,
    function(input, output, session) {

      init.state <- c(
        S = theta[["N"]]  - theta[["init_inf"]] ,
        E = 0,
        Evac = 0,
        Ivac = 0,
        Rvac = 0,
        I = theta[["init_inf"]] ,
        R = 0,
        Sv = 0,
        Ev = 0,
        Iv = 0,
        Rv = 0,
        Vp = 0,
        Ep = 0,
        Ip = 0,
        Rp = 0,
        C1 = 0,
        Ec1 = 0,
        Ic1 = 0,
        Rc1 = 0,
        C2 = 0,
        Ec2 = 0,
        Ic2 = 0,
        Rc2 = 0
      )

      times <- seq(0,18*10,1)

      vac_start    <- theta[["vac_start"]]
      d_inf        <- theta[["d_inf"]]
      N            <- theta[["N"]]


      traj <- one_dose_leaky_model(theta, init.state, times)
      num_rows <- nrow(traj)

      aggr_data <- traj[vac_start:num_rows,] %>%
        subset( select=c( Sv, Ev, Vp, Iv, Rv, Ip, Rp, Ic1, Ic2, Rc1, Rc2, Evac, Ivac, Rvac)) %>%
        mutate(
          Vac = Sv + Evac,
          Vac_exp = Evac,
          Vac_prot = Vp,
          Vac_inf = Ip, # are excluded here as the infected move to R
#          Vac_inf_washout = Iv,   # Vaccinated infected during washout period
          Cont_inf = Ic1 + Ic2,
          Evac_inf = Ivac,
          Rec_exposed  = Ev               # recently exposed vaccinated people
        ) %>%
        select( Vac, Vac_exp, Vac_prot, Vac_inf, Cont_inf, Evac_inf, Rec_exposed) %>% # Vac_inf_washout) %>%
        # To get the incidence divide the prevalence results by duration of infection
        mutate(
          Vac = as.integer(Vac),#as.integer(Vac/num_vac_daily ), #- lag(Vac, default = first(Vac))),
          Vac_exp = as.integer(Vac_exp),
          Vac_prot = as.integer(Vac_prot/theta[["immunity_onset"]]),
#          Vac_inf_washout = as.integer(Vac_inf_washout / d_inf),
          Vac_inf = as.integer(Vac_inf / theta[["d_inf"]]),
          Cont_inf = as.integer(Cont_inf / theta[["d_inf"]]),
          Evac_inf = as.integer(Evac_inf / theta[["d_inf"]])
        )

      aggr_data

    }
  )
}

