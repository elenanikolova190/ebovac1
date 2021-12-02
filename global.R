

# WELCOME TO EBOVAC

#EBOVAC is licensed under the GNU General Public License (GPL) v2.0 (https://github.com/ceefluz/radar/blob/master/LICENSE)

# INSTALL DEPENDENCIES ----------------------------------------------------

#source('dependencies.R')
# load all packages
#lapply(required_packages, require, character.only = TRUE)

# DATA TRANSFORMATION AND NEW VARIABLES -----------------------------------


# HELP & INTRO DATA ---------------------------------------------------------------

steps <- read_csv2("help.csv")


# FLUID DESIGN FUNCTION ---------------------------------------------------

fluid_design <- function(id, w, x, y, z) {
  fluidRow(
    div(
      id = id,
      column(
        width = 6,
        uiOutput(w),
        uiOutput(y)
      ),
      column(
        width = 6,
        uiOutput(x),
        uiOutput(z)
      )
    )
  )
}


# SUMMERISE FUNCTIONS ------------------------------------------------------------------
vaccinated_15days <- function(Vn, vac_duration){

  vac_daily   <- Vn/vac_duration
  vac_current <- 0
  vac_current[1] <- as.integer(vac_daily*15)

  for(i in 2:10) {
    if(vac_current[i-1] < (Vn - vac_current[1]))
    {
      vac_current[i] <- as.integer(vac_daily * 15 + vac_current[i-1])
    } else{
      vac_current[i] <- as.integer(Vn)
    }
  }
  vac_current
}
