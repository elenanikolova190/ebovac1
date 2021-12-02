#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyjs)
library(shinydashboard)
library(shinycssloaders)
library(shinyWidgets)
library(shinyBS)
library(rintrojs)
library(tidyverse)
library(deSolve)
library(dplyr)
library(DT)
library(rsconnect)


source('1dlp_module.R', local = TRUE)
source('2dlp_module.R')
source('1daonp_module.R')
source('global.R')


ui <- dashboardPage(
  skin = "black",
  title = "EBOVAC",

  # HEADER -----------------------------------------------------------------

  dashboardHeader(
    title = "EBOVAC",
    titleWidth = 300,
    dropdownMenu(
      type = "notifications",
      headerText = strong("HELP"),
      icon = icon("question"),
      badgeStatus = NULL,
      notificationItem(
        text = (steps$text[1]),
        icon = icon("spinner")
      ),
      notificationItem(
        text = steps$text[2],
        icon = icon("address-card")
      ),
      notificationItem(
        text = steps$text[3],
        icon = icon("calendar")
      ),
      notificationItem(
        text = steps$text[4],
        icon = icon("user-md")
      ),
      notificationItem(
        text = steps$text[5],
        icon = icon("ambulance")
      ),
      notificationItem(
        text = steps$text[6],
        icon = icon("flask")
      ),
      notificationItem(
        text = strong(steps$text[7]),
        icon = icon("exclamation")
      )
    ),
    tags$li(
      a(
        strong("ABOUT Ebovac"),
        height = 40,
        title = "",
        target = "_blank"
      ),
      class = "dropdown"
    )
  ),

  # SIDEBAR -----------------------------------------------------------------

  dashboardSidebar(
    width = 300,
      sidebarMenu(
              div(id = "sidebar_button",
                  br(),
                        bsButton(inputId = "confirm",
                                     label = "START EBOVAC",
                                     icon = icon("play-circle"),
                                     style = "danger")
                        ),
               div(class = "inlay", style = "height:15px;width:100%;background-color: #ecf0f5;"),
               menuItem("POPULATION", tabName = "population",
                        startExpanded = FALSE,
                        numericInput("N", "Affected population", min = 10, max = 50000000, value = 100000),
                        numericInput("init_inf", "Initial number of case", value = 10)
               ),
               menuItem("VACCINE", tabName = "vaccine",
                        startExpanded = FALSE,
                        sliderInput("VE", "Vaccine Efficacy", min = 0.0, max = 1.0, value = 0.75),
                        numericInput("imm_onset", "Onset of immunity of protection (in days)", value = 14),
                        numericInput("imm_onset_2", "Onset of immunity afret 2nd dose", value = 10)
               ),
               menuItem("TRIAL", tabName = "trial",
                        startExpanded = FALSE,
                        numericInput("start_date", "Starting day of trial", value = 15),
                        numericInput("num_doses", "Number of vaccine doses", value = 10000),
                        numericInput("duration", "Trial duration in days", value = 30),
                        numericInput("vac_gap", "Time between 1st and 2nd doses", value = 30)
               ),
               menuItem("EPIDEMIC CHARACTERISTICS", tabName = "epidem",
                        startExpanded = FALSE,
                        numericInput("R0", "Reproductive number (R0)", value = 1.5),
                        numericInput("d_lat", "Average latent period (in days)", value = 10),
                        numericInput("d_inf", "Average infectious peiod (in days)", value = 12)
               ),
               menuItem("INTERVENTION", tabName = "interven",
                        startExpanded = FALSE,
                        numericInput("alpha", "Shape of change (alpha)", value = 0.01),
                        numericInput("shift", "Midpoint of the time change (shift)", value = 3),
                        numericInput("reduct", "Reduction in transmission rate", value = 0.09),
                        numericInput("intervention_time", "Time of intervention", value = 20)
               )
        )
  ),  # end of sidebar

  # BODY --------------------------------------------------------------------

  dashboardBody(
    tags$head(
      tags$link(
        rel = "stylesheet",
        type = "text/css",
        href = "ebovac_style.css")
    ),
    useShinyjs(),
    introjsUI(),

    fluidRow(
      h4("One Dose Leaky Protection Model Results"), br(),
      plotOutput("plot1dl"), br(),
      dataTableOutput("table1dl")
#      h4("One Dose All or Nothing Protection Model Results"), br(),
#      plotOutput("plot1daon"), br(),
#      dataTableOutput("table1daon")
 #     tabsetPanel(
#        tabPanel("One dose leaky",
 #                fluidRow(
  #                 column(12,
   #                       plotOutput("plot1dl")
    #              ),
     #             column(12,
      #                   ç
#                  )
#                )
 #       ),
  #      tabPanel("One dose all or nothing",
   #              fluidRow(
    #               column(12,
     #                     plotOutput("plot1daon")
      #             ),
       #            column(12,
        #                  tableOutput("table1daon")
#                   )
 #               )
  #      ),
   #     tabPanel("Two dose leaky",
#                 fluidRow(
 #                  column(12,
  #                        plotOutput("plot1dl")
   #                ),
    #               column(12,
     #                     tableOutput("table1dl")
      #             )
       #          )
  #      ),
   #     tabPanel("Two dose all or nothing",
  #               fluidRow(
   #                column(12,
    #                      plotOutput("plot1dl")
     #              ),
      #             column(12,
       #                   tableOutput("table1dl")
        #           )
         #        )
#        )
 #     )
    )

  ) # end of dashboardbody

) #end of ui function




server <- function(input, output) {

## UPDATE INPUT VARIABLES-------------------------------------------------------
  observeEvent(
    input$confirm,
    {

#      print("entered the observe event sequence")

      # The intput variables parsed into a data frame
      inputVars <- reactiveValues(
        N           = input$N,
        init_inf    = input$init_inf,
        VE          = input$VE,
        Nv          = input$num_doses,
        d_lat       = input$d_lat,
        d_inf       = input$d_inf,
        R0          = input$R0,
        vac_duration     = input$duration,    # vaccination campaign duration in days
        vac_start        = input$start_date,     # vaccination start date relative to the start of epidemic
        immunity_onset   = input$imm_onset,
        immunity_onset_2 = input$imm_onset_2,  # For 2 dose models: onset of immunity for second dose
        vac_gap          = input$vac_gap       # For 2 dose models: gap between doses in days
      )

      theta <- reactiveValuesToList(inputVars)

      aggr_data_1dlp <- onedlpServer("1dlp", theta )
      aggr_data_1daonp <- onedaonpServer("1daonp", theta )
      aggr_data_2dlp <- twodlpServer("2dlp", theta)

      # One dose leaky protection table
      weekly_out_1dlp <- aggr_data_1dlp %>%
        mutate(weeks = 1:nrow(aggr_data_1dlp) %/% 7) %>%
        group_by(weeks) %>%
        summarise( Vaccinated_num = sum(Vac),
                   Vprotected = sum(Vac_prot),
                   Vaccinated = sum(Vac_inf),
                   Control = sum(Cont_inf),
                   Evac_inf = sum(Evac_inf)
        )


      df_plot_1dlp <- weekly_out_1dlp %>%
        subset( select=c(weeks, Control, Vaccinated) ) %>%
        gather(state, Infected, -weeks)

      # One dose all or nothing protection table
      weekly_out_1daon <- aggr_data_1daonp %>%
        mutate(weeks = 1:nrow(aggr_data_1daonp) %/% 7) %>%
        group_by(weeks) %>%
        summarise( Vaccinated_num = sum(Vac),
                   Vprotected = sum(Vac_prot),
                   Vaccinated = sum(Vac_inf),
                   Control = sum(Cont_inf),
                   Evac_inf = sum(Evac_inf)
        )


      df_plot_1daon <- weekly_out_1daon %>%
        subset( select=c(weeks, Control, Vaccinated) ) %>%
        gather(state, Infected, -weeks)



      #SUMMARY TABLE FOR ONE DOSE LEAAKY PROTECTION MODEL-----------------------------------------

      vac_duration <- theta[["vac_duration"]]
      vac_total    <- theta[["Nv"]]

      vaccinated_15days_ <- as.matrix(vaccinated_15days(Vn = vac_total, vac_duration = vac_duration))

      summary_15days <- aggr_data_1dlp %>%
        mutate(periods = 1:nrow(aggr_data_1dlp) %/% 15) %>%
        group_by(periods) %>%
        summarise(
          Vprotected_15 = sum(Vac_prot),
          Vaccinated_Inf_15 = sum(Vac_inf),
          Control_Inf_15 = sum(Cont_inf)
        ) %>%
        subset( select=c( Vaccinated_Inf_15, Control_Inf_15)) %>% t()

      my_name_vector = c("Day15", "Day30", "Day45", "Day60", "Day75", "Day90", "Day105", "Day120", "Day135", "Day150");

      outTable = as.data.frame(t(vaccinated_15days_));
      colnames(outTable) <- my_name_vector;
      sum15Table <- summary_15days[1:2, 1:10]

      colnames(sum15Table) <- my_name_vector;
      summaryTable <- rbind(outTable, sum15Table)
      rownames(summaryTable) <- c("Vaccinated","Vaccinated Infected","Control Arm Infected")


      output$plot1dl <- renderPlot({
        ggplot(data = df_plot_1dlp) + ggtitle("Weekly incidence in trial participants") +
          theme(plot.title = element_text(hjust = 0.5, size=14, face="bold.italic")) +
          geom_point(mapping = aes(x = weeks, y = Infected, color = state))
        })
      output$table1dl <- renderDataTable(summaryTable, caption = "Number of infected trial participants")

      output$plot1daon <- renderPlot({
        ggplot(data = df_plot_1daon ) +
          geom_point(mapping = aes(x = weeks, y = Infected, color = state))
      })
      output$table1daon <- renderDataTable(summaryTable)

      output$plot2dl <- renderPlot({
        ggplot(data = df_plot) +
          geom_point(mapping = aes(x = weeks, y = Infected, color = state))
      })
      output$table2dl <- renderDataTable(summaryTable)

      output$plot2daon <- renderPlot({
        ggplot(data = df_plot) +
          geom_point(mapping = aes(x = weeks, y = Infected, color = state))
      })

      output$table2daon <- renderDataTable(summaryTable)

    }
  )

}

shinyApp(ui, server = server)



