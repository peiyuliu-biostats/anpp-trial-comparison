# /R/ui_outputs.R
# ui module3: all output elements (plots and results).

uimod_content_outputs <- function() {
  fluidRow(
    # -- Row 4: result Outputs --
    column(width = 6,
           box(title = "MSP Antibody-based Trial Results", status = "info", solidHeader = TRUE, width = NULL,
               tabsetPanel(
                 id = "msp_tabs",
                 tabPanel("Biting Rate", plotOutput("plot_lambda_B")),
                 tabPanel("Individual Trajectory", plotOutput("plot_indiv_traj")),
                 tabPanel("Population Trajectory", plotOutput("plot_pop_traj")),
                 tabPanel("Sample Size", plotOutput("plot_ss_trial1"))
               )
           )
    ),
    column(width = 6,
           box(title = "Infection-based Trial Results", status = "success", solidHeader = TRUE, width = NULL,
               tabsetPanel(
                 id = "infection_tabs",
                 tabPanel("Infection Hazard", plotOutput("plot_hazard_I")),
                 tabPanel("Cumulative Incidence", plotOutput("plot_ci")),
                 tabPanel("Sample Size", plotOutput("plot_ss_trial2"))
               )
           )
    )
  )
}