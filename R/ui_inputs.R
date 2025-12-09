# R/ui_inputs.R
# ui Module2: all input controls, buttons, and the records table.

# helper function: parameters input describtions with a tooltip  
param_input_factory <- function(inputId, label_symbol, tooltip_text, value, min, max, step) {
  div(style = "display: flex; align-items: center; margin-bottom: 5px;",
      # Left side: Label and info icon
      div(style = "width: 50%;",
          tags$label(HTML(label_symbol), style = "margin-right: 5px;"),
          shiny::icon("info-circle", title = tooltip_text)
      ),
      # Right side: Numeric input box
      div(style = "width: 50%;",
          numericInput(inputId, label = NULL, value = value, min = min, max = max, step = step)
      )
  )
}

uimod_content_inputs <- function() {
  tagList(
    # -- Row 1: parameter Inputs --
    fluidRow(
      box(title = "MSP Antibody-based Trial Parameters", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 4,
          param_input_factory("lambda_max", "λ_max:", "The maximum number of mosquito bites per person per day within a year", 20, 5, 65, 1),
          param_input_factory("lambda_min", "λ_min:", "The minimum number of mosquito bites per person per day within a year", 2, 1, 35, 1),
          param_input_factory("mu_c", "μ_c:", "The instantaneous increase of individual mosquito salivary protein antibodies occurs at the moment of a mosquito bite", 1, 0.1, 30, 0.1),
          param_input_factory("prop_c", "prop_c:", "The ratio of the variance to the mean of the instantaneous increment C", 1.5, 1.5, 3, 0.1),
          param_input_factory("gamma0_hl", "Half-life:", "Half-life (in days) of mosquito saliva protein antibodies in individuals", 45, 15, 145, 0.1),
          param_input_factory("B0", "B0:", "Initial value of individual mosquito salivary protein antibody", 5, 0.1, 150, 0.1),
          param_input_factory("sigma_epsilon0", "σ_ε:", "Standard deviation of measurement error of B0", 5, 0.1, 20, 0.1)
      ),
      box(title = "Infection-based Trial Parameters", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 4,
          param_input_factory("PV", "PV:", "Prevalence of infectious mosquitoes with dengue virus", 0.3, 0.1, 0.7, 0.01),
          param_input_factory("P_hv", "P_hv(%):", "Probability of mosquitoes infectiousness for dengue fever", 30, 1, 95, 0.1),
          param_input_factory("Pc", "P_c:", "The overall probability that a single random bite will result in human Dengue infection.", 0.0005, 0.0001, 0.001, 0.0001),
          param_input_factory("Ps", "Ps(%):", "Proportion of susceptible people to dengue fever", 30, 1, 95, 0.1)
      ),
      box(title = "Common Parameters", status = "warning", solidHeader = TRUE, collapsible = TRUE, width = 4,
          param_input_factory("T", "T:", "Trial duration (days)", 365, 120, 730, 1),
          param_input_factory("theta", "θ:", "Effectiveness of biting intervention.", 0.5, 0.01, 0.99, 0.01),
          param_input_factory("alpha", "α:", "Type I Error", 0.05, 0.01, 0.1, 0.01),
          param_input_factory("power", "1 - β:", "Power", 0.8, 0.8, 0.99, 0.01)
      )
    ),
    
    # -- Row 2: add Record Button --
    fluidRow(
      column(width = 12, align = "center",
             actionButton("add_record", "Add Record", icon = icon("plus"), style = "color: #fff; background-color: #6f42c1; border-color: #6f42c1; margin-bottom: 15px;")
      )
    ),
    
    # -- Row 3: records Table --
    fluidRow(
      box(title = "Parameter Records", status = "primary", solidHeader = TRUE, width = 12,
          actionButton("remove_record", "Remove Selected Record", icon = icon("trash"), style = "margin-bottom: 10px;"),
          DT::dataTableOutput("records_table")
      )
    )
  )
}