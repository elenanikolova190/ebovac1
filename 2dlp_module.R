#' Title Two doses leaky protection model functions
#'
#'

two_dose_leaky <- function(time, state, theta) {

  # transform theta to model parameters
  N            <- theta[["N"]]            # Affected population size
  r_lat        <- 1/theta[["d_lat"]]      # Latency period
  r_inf        <- 1/theta[["d_inf"]]      # Infectious period
#  r_prot       <- 1/theta[["d_prot"]]     # Rate of protection
  VE           <- theta[["VE"]]           # Vaccine efficacy
  vac_start    <- theta[["vac_start"]]    # Starting day of trial (counted from beginning of epidemic, e.g. day 30)
  vac_duration <- theta[["vac_duration"]] # How many days the trial will go on
  doses_num    <- theta[["Nv"]]           # Number of vaccines administered in trial
  beta         <- theta[["R0"]]*r_inf     # Power of infection
  kappa        <- 1/theta[["immunity_onset"]] # Rate of onset of immunity

  # Additional parameters for second dose
  kappa2       <- 1/theta[["imm_onset_2"]]          # rate of onset after booster
  vac_gap      <- theta[["vac_gap"]]                # gap between first and second dose


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

  # Vaccinated arm of the trial after first dose of vaccine
  Sv <- state[["Sv"]]
  Ev <- state[["Ev"]]
  Iv <- state[["Iv"]]
  Rv <- state[["Rv"]]

  # Vaccinated arm of the trial after onset of partial immunity
  Vp <- state[["Vp"]]
  Ep <- state[["Ep"]]
  Ip <- state[["Ip"]]
  Rp <- state[["Rp"]]

  # Vaccinated arm of the trial after booster
  B <- state[["B"]]
  Eb <- state[["Eb"]]
  Ib <- state[["Ib"]]
  Rb <- state[["Rb"]]

  # Vaccinated arm of the trial after onset of partial immunity
  Bp <- state[["Bp"]]
  Ebp <- state[["Ebp"]]
  Ibp <- state[["Ibp"]]
  Rbp <- state[["Rbp"]]

  # Control arm of the trial after first dose
  C1 <- state[["C1"]]
  Ec1 <- state[["Ec1"]]
  Ic1 <- state[["Ic1"]]
  Rc1 <- state[["Rc1"]]

  # Control arm of the trial after onset of immunity
  C2 <- state[["C2"]]
  Ec2 <- state[["Ec2"]]
  Ic2 <- state[["Ic2"]]
  Rc2 <- state[["Rc2"]]

  # Control arm of the trial after booster dose
  Cb <- state[["Cb"]]
  Ecb <- state[["Ecb"]]
  Icb <- state[["Icb"]]
  Rcb <- state[["Rcb"]]

  # Control arm of the trial after booster dose after period k2
  Cb2 <- state[["Cb2"]]
  Eb2 <- state[["Eb2"]]
  Ib2 <- state[["Ib2"]]
  Rb2 <- state[["Rb2"]]

  # Calculate end day of trial for 1st dose and booster
  vac_end_1st_dose  <- vac_duration + vac_start
  vac_start_booster <- vac_gap + vac_start
  vac_end_booster   <- vac_start_booster + vac_duration

  # Rate of infection
  lambda <- beta*(I + Ivac + Ip + Iv + Ibp + +Ib + Icb + Ib2 + Ic1 + Ic2)/N

  if(time < vac_start)      # Before vaccination trial start
  {

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
    dB <- 0
    dEb <- 0
    dIb <- 0
    dRb <- 0
    dBp <- 0
    dEbp <- 0
    dIbp <- 0
    dRbp <- 0
    dC1 <- 0
    dEc1 <- 0
    dIc1 <- 0
    dRc1 <- 0
    dC2 <- 0
    dEc2 <- 0
    dIc2 <- 0
    dRc2 <- 0
    dCb <- 0
    dEcb <- 0
    dIcb <- 0
    dRcb <- 0
    dCb2 <- 0
    dEb2 <- 0
    dIb2 <- 0
    dRb2 <- 0
  }
  else if((time >= vac_start) & (time < vac_start_booster )) # From the moment trial begins until booster
  {
    # To ensure that vaccination ends end right time
    if((time >= vac_start) & (time < vac_end_1st_dose )){
      r_vac <- (doses_num /(S+E)) / vac_duration
      r_vac_b <- 0.0
    } else {
      r_vac <- 0
      r_vac_b <- 0.0
    }

      # General Population
      dS <-  - lambda * S - (r_vac * S) * 2
      dE <-  lambda * S - r_lat * E - r_vac * E
      dI <-  r_lat * E - r_inf*I
      dR <-  r_inf * I

      # Participants in trial who have been exposed prior to vaccination
      dEvac <- r_vac * E - r_lat * Evac
      dIvac <- r_lat * Evac - r_inf * Ivac
      dRvac <-  r_inf * Ivac

      # Vaccinated arm of the trial after first dose
      dSv <- r_vac * S - lambda * Sv - kappa * Sv
      dEv <- lambda * Sv - r_lat * Ev
      dIv <- r_lat * Ev - r_inf * Iv
      dRv <- r_inf * Iv

      # Vaccinated arm of the trial after onset of partial immunity
      dVp <- kappa*Sv - r_vac_b*Vp - (1 - VE/2) * lambda * Vp
      dEp <- (1 - VE/2) * lambda * Vp - r_lat * Ep
      dIp <- r_lat * Ep - r_inf * Ip
      dRp <- r_inf *Ip

      # Control arm of the trial before onset of immunity
      dC1  <- r_vac * S - lambda * C1 - kappa * C1
      dEc1 <- lambda * C1 - r_lat * Ec1
      dIc1 <- r_lat * Ec1 - r_inf * Ic1
      dRc1 <- r_inf * Ic1

      # Control arm of the trial after onset of partial immunity
      dC2  <-  kappa * C1 - r_vac_b * C2 - lambda * C2
      dEc2 <- lambda * C2 - r_lat * Ec2
      dIc2 <- r_lat * Ec2 - r_inf * Ic2
      dRc2 <- r_inf * Ic2

      # All compartments related to booster are set to 0
      dB <- 0
      dEb <- 0
      dIb <- 0
      dRb <- 0
      dBp <- 0
      dEbp <- 0
      dIbp <- 0
      dRbp <- 0
      dCb <- 0
      dEcb <- 0
      dIcb <- 0
      dRcb <- 0
      dCb2 <- 0
      dEb2 <- 0
      dIb2 <- 0
      dRb2 <- 0

    } else {         # Start of Booster administration

      # Rate of booster administration
      if((time > vac_start_booster) & (time < vac_end_booster)){
        r_vac_b <- (doses_num/Vp) / vac_duration
        r_vac   <- 0.0
      } else {
        r_vac_b <- 0.0
        r_vac   <- 0.0
      }

      # General Population
      dS <-  - lambda * S - (r_vac * S) * 2
      dE <-  lambda * S - r_lat * E - r_vac * E
      dI <-  r_lat * E - r_inf*I
      dR <-  r_inf * I

      # Participants in trial who have been exposed prior to vaccination
      dEvac <- r_vac * E - r_lat * Evac
      dIvac <- r_lat * Evac - r_inf * Ivac
      dRvac <-  r_inf * Ivac

      # Vaccinated arm of the trial after first dose
      dSv <- r_vac * S - lambda * Sv - kappa * Sv
      dEv <- lambda * Sv - r_lat * Ev
      dIv <- r_lat * Ev - r_inf * Iv
      dRv <- r_inf * Iv

      # Vaccinated arm of the trial after onset of partial immunity
      dVp <- kappa*Sv - r_vac_b*Vp - (1 - VE/2) * lambda * Vp
      dEp <- (1 - VE/2) * lambda * Vp - r_lat * Ep
      dIp <- r_lat * Ep - r_inf * Ip
      dRp <- r_inf *Ip

      # Vaccinated arm of the trial after booster
      dB  <- r_vac_b * Vp - (1 - VE/2) * lambda * B - kappa2*B
      dEb <- (1 - VE/2) * lambda * B - r_lat * Eb
      dIb <- r_lat * Eb - r_inf * Ib
      dRb <- r_inf *Ib

      # Vaccinated arm of the trial after onset of partial immunity
      dBp  <- kappa2*B - (1-VE)*lambda*Bp
      dEbp <- (1 - VE) * lambda * Bp - r_lat * Ebp
      dIbp <- r_lat * Ebp - r_inf * Ibp
      dRbp <- r_inf *Ibp

      # Control arm of the trial before onset of immunity
      dC1  <- r_vac * S - lambda * C1 - kappa * C1
      dEc1 <- lambda * C1 - r_lat * Ec1
      dIc1 <- r_lat * Ec1 - r_inf * Ic1
      dRc1 <- r_inf * Ic1

      # Control arm of the trial after onset of partial immunity
      dC2  <-  kappa * C1 - r_vac_b * C2 - lambda * C2
      dEc2 <- lambda * C2 - r_lat * Ec2
      dIc2 <- r_lat * Ec2 - r_inf * Ic2
      dRc2 <- r_inf * Ic2

      # Control arm of the trial after booster dose
      dCb  <-  r_vac_b * C2 - lambda * Cb
      dEcb <- lambda * Cb - r_lat * Ecb
      dIcb <- r_lat * Ecb - r_inf * Icb
      dRcb <- r_inf * Icb

      # Control arm of the trial after onset of immunity after booster
      dCb2 <- kappa2 * Cb - lambda * Cb2
      dEb2 <- lambda * Cb2 - r_lat * Eb2
      dIb2 <- r_lat * Eb2 - r_inf * Ib2
      dRb2 <- r_inf * Ib2

    }

  return(list(c(dS, dE, dI, dR, dEvac, dIvac, dRvac,
                dSv, dEv, dIv, dRv, dVp, dEp, dIp, dRp,
                dB, dEb, dIb, dRb, dB, dEb, dIb, dRb,
                dC1, dEc1, dIc1, dRc1, dC2, dEc2, dIc2, dRc2,
                dCb, dEcb, dIcb, dRcb, dCb2, dEb2, dIb2, dRb2)))


}


#  Wrapper function for desolve ode function
#
#' @return
#' @export
#'
#' @examples two_dose_leaky_model(theta, init, times)
#'

two_dose_leaky_model <- function(theta, init.state, times) {

  # put incidence at 0 in init.state
  traj <- as.data.frame(ode(init.state, times, two_dose_leaky, theta, method = "ode45"))

  return(traj)

}

#' Function: Two Dose Leaky Protection Model Module input
#'
#'
#'
twodlpUI <- function(id) {

}


twodlpServer <- function(id, theta) {


  moduleServer(
    id,
    function(input, output, session) {

      init.state <- c(
        S = theta[["N"]]  - theta[["init_inf"]] ,
        E = 0,
        I = theta[["init_inf"]] ,
        R = 0,
        Evac = 0,
        Ivac = 0,
        Rvac = 0,
        Sv  = 0,
        Ev = 0,
        Iv = 0,
        Rv = 0,
        Vp = 0,
        Ep = 0,
        Ip = 0,
        Rp = 0,
        B = 0,
        Eb = 0,
        Ib = 0,
        Rb = 0,
        Bp = 0,
        Ebp = 0,
        Ibp = 0,
        Rbp = 0,
        C1 = 0,
        Ec1 = 0,
        Ic1 = 0,
        Rc1 = 0,
        C2 = 0,
        Ec2 = 0,
        Ic2 = 0,
        Rc2 = 0,
        Cb = 0,
        Ecb = 0,
        Icb = 0,
        Rcb = 0,
        Cb2 = 0,
        Eb2 = 0,
        Ib2 = 0,
        Rb2 = 0
      )

      times <- seq(0,18*10,1)

      vac_start    <- theta[["vac_start"]]
      d_inf        <- theta[["d_inf"]]
      N            <- theta[["N"]]


      traj <- two_dose_leaky_model(theta, init.state, times)
      num_rows <- nrow(traj)

      print(traj)

      aggr_data <- traj[vac_start:num_rows,] %>%
        subset( select=c( Sv, Ev, Vp, Ip, Rp,Iv, B, Ib, Rb, Ic1, Ic2, Rc1, Rc2, Icb, Rcb, Ib2, Rb2, Evac, Ivac, Rvac)) %>%
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
        select( Vac, Vac_exp, Vac_prot, Vac_inf, Cont_inf, Evac_inf, Rec_exposed ) %>%  #, Vac_inf_washout) %>%
        # To get the incidence divide the prevalence results by duration of infection
        mutate(
          Vac = as.integer(Vac),#as.integer(Vac/num_vac_daily ), #- lag(Vac, default = first(Vac))),
          Vac_exp = as.integer(Vac_exp),
          Vac_prot = as.integer(Vac_prot/theta[["immunity_onset"]]),
  #ß        Vac_inf_washout = as.integer(Vac_inf_washout / d_inf),
          Vac_inf = as.integer(Vac_inf / d_inf),
          Cont_inf = as.integer(Cont_inf / d_inf),
          Evac_inf = as.integer(Evac_inf / d_inf)
        )

      aggr_data

    }
  )
}

