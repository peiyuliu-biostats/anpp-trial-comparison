uimod_content_doc <- function(which = c("model_description", "author")) {
  which <- match.arg(which)
  
  if (which == "model_description") {
    
    tabItem(
      tabName = "model_description",
      withMathJax(),
      fluidRow(
        box(
          title = "Model Description for the Trial Comparison App",
          solidHeader = TRUE, status = "primary", width = 12,
          
          p(strong("Author:"), " Peiyu Liu, Department of Biostatistics, University of Florida"),
          hr(),
          
          h3("Overview"),
          p("This application provides an interactive platform to simulate and compare two distinct clinical trial designs for evaluating vector-control interventions. The core objective is to contrast the statistical power and required sample sizes of a novel serological endpoint (based on anti-MSP antibodies) against a traditional clinical endpoint (based on dengue infection) [1]."),
          p("The app is built upon two underlying mechanistic models [1]:"),
          tags$ul(
            tags$li(strong("The MSP Antibody-based Model (Trial 1):"), " A stochastic model that simulates how an individual's antibody levels change over time in response to seasonal mosquito bites [1]."),
            tags$li(strong("The Infection-based Model (Trial 2):"), " A simplified epidemiological model that calculates the risk of clinical infection based on the same seasonal biting pattern [1].")
          ),
          p("By adjusting the parameters in the sidebar, you can explore how different immunological, epidemiological, and trial design assumptions influence the efficiency of each approach [1]."),
          
          hr(),
          h3("1. The MSP Antibody-based Trial Model"),
          p("This model simulates the dynamics of antibodies against mosquito salivary protein (MSP). It assumes that each mosquito bite acts as a small 'immunological boost' and the total antibody level is the cumulative result of all past bites, balanced by natural decay [1]."),
          
          h4("Biting Rate (\\(\\lambda(t)\\))"),
          p("The foundation of the model is a seasonal mosquito biting rate. We assume the number of bites a person receives per day follows a sinusoidal pattern over a year [1]."),
          
          h4("Antibody Impulse Response and Decay"),
          p("The model is based on three key immunological processes [1]:"),
          tags$ul(
            tags$li(strong("Impulse Magnitude (C)"), "[1]"),
            tags$li(strong("Antibody Decay (\\(\\gamma\\))"), "[1]"),
            tags$li(strong("Baseline and Noise"), "[1]")
          ),
          
          h4("How Sample Size is Calculated (Trial 1)"),
          p("The required sample size is calculated using a standard two-sample t-test framework on the log-transformed antibody levels to detect the difference in mean antibody levels with the specified statistical power (1-β) and significance level (α) [1]."),
          
          hr(),
          h3("2. The Infection-based Trial Model"),
          p("This model calculates the sample size for a traditional trial where the endpoint is the first documented clinical infection with dengue [1]."),
          
          h4("From Bites to Infection Hazard"),
          p("The daily risk of infection (the 'hazard rate') is directly proportional to the seasonal biting rate, \\(\\lambda(t)\\), but is modified by several key probabilities that you control [1]."),
          
          h4("Cumulative Incidence"),
          p("The app integrates this daily hazard over the trial duration (T) to calculate the cumulative incidence [1]."),
          
          h4("How Sample Size is Calculated (Trial 2)"),
          p("The sample size calculation for this trial is based on comparing the cumulative infection risk between the control and intervention groups using a standard method for time-to-event data based on the Average Hazard Ratio (AHR) [1]."),
          
          hr(),
          h3("3. Common Parameters"),
          p("These parameters define the overall trial design and are applied to both models to ensure a fair comparison [1].")
        )
      )
    )
    
  } else if (which == "author") {
    
    tabItem(
      tabName = "author",
      fluidRow(
        box(
          title = "About This Application",
          solidHeader = TRUE,
          status = "primary",
          width = 12,
          h3("Author Information"),
          p("This application was developed to compare clinical trial designs based on different endpoints."),
          tags$ul(
            tags$li(HTML("<b>Author:</b> Peiyu Liu")),
            tags$li(HTML("<b>Affiliation:</b> Department of Biostatistics, University of Florida")),
            tags$li(HTML("<b>Contact:</b> <a href='mailto:pyliu0620@outlook.com'>pyliu0620@outlook.com</a>"))
          ),
          hr()
        )
      )
    )
  }
}
