library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(pwr)
library(DT)  

# ==============================================================================
# 1. Helper functions
# ==============================================================================
# Calculate log-normal distribution parameters from mean and variance
find_lognorm_params <- function(target_mean, target_var) {
  if (target_mean <= 0 || target_var <= 0) { return(list(mu = NA, sigma_sq = NA)) }
  cv_sq <- target_var / (target_mean^2)
  if (cv_sq < -1) return(list(mu = NA, sigma_sq = NA))
  sigma_sq <- log(1 + cv_sq)
  mu <- log(target_mean) - sigma_sq / 2
  return(list(mu = mu, sigma_sq = sigma_sq))
}

# Calculate sample size for a two-sample t-test
calculate_sample_size_ttest <- function(mu_c, var_c, mu_i, var_i, power = 0.8, alpha = 0.05) {
  if (any(is.na(c(mu_c, var_c, mu_i, var_i))) || var_c <= 0 || var_i <= 0) return(NA)
  log_mu_c <- log(mu_c) - 0.5 * var_c / mu_c^2
  log_var_c <- log(1 + var_c / mu_c^2)
  log_mu_i <- log(mu_i) - 0.5 * var_i / mu_i^2
  log_var_i <- log(1 + var_i / mu_i^2)
  if (any(is.na(c(log_mu_c, log_var_c, log_mu_i, log_var_i))) || log_var_c <= 0 || log_var_i <= 0) return(NA)
  sd_pooled <- sqrt((log_var_c + log_var_i) / 2)
  d <- abs(log_mu_c - log_mu_i) / sd_pooled
  if (is.na(d) || d < 1e-6) return(NA)
  pwr_result <- try(pwr.t.test(n = NULL, d = d, sig.level = alpha, power = power, type = "two.sample", alternative = "two.sided"), silent = TRUE)
  if (inherits(pwr_result, "try-error")) return(NA)
  return(ceiling(pwr_result$n))
}

# Generate event times for a non-homogeneous Poisson process
generate_event_times <- function(lambda0, lambda_amp, t_end, omega) {
  lambda_max <- lambda0 + lambda_amp
  n_candidate <- rpois(1, lambda_max * t_end)
  if (n_candidate == 0) return(numeric(0))
  candidate_times <- sort(runif(n_candidate, 0, t_end))
  accept_prob <- (lambda0 + lambda_amp * sin(omega * candidate_times)) / lambda_max
  return(candidate_times[runif(n_candidate) <= accept_prob])
}

# Calculate the antibody trajectory for a single individual
calculate_antibody_trajectory <- function(event_times, C_values, B0, epsilon, gamma, times) {
  sapply(times, function(t) {
    relevant_events <- event_times[event_times <= t]
    relevant_Cs <- C_values[1:length(relevant_events)]
    base_A <- if (length(relevant_events) > 0) {
      sum(relevant_Cs * exp(-gamma * (t - relevant_events)))
    } else {0}
    return(base_A + B0 + epsilon)
  })
}

# ==============================================================================
# 2. Shiny UI (User Interface)
# ==============================================================================

# UI helper function: creates a compact parameter input with a tooltip
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

#  Sidebar Settings

ui <- dashboardPage(
  dashboardHeader(title = "Trial Comparison"),
  # -- Sidebar with new structure --
  dashboardSidebar(
    sidebarMenu(id = "tabs",
                tags$li(class = "header", "TRIAL COMPARISON DESIGN"),
                menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
                menuItem("Github", icon = icon("github"), href = "https://github.com/peiyuliu-biostats/anpp-trial-comparison"),
                menuItem("Model Description", icon = icon("book-open"), href = "model_description.html"),
                menuItem("shinyapps.io", icon = icon("cloud-upload-alt"), href = "https://peiyuliu.shinyapps.io/anpp-trial-comparison/"),
                menuItem("Author", tabName = "author", icon = icon("user-circle"))
    )
  ),
  dashboardBody(
    # Add custom CSS to modify the dashboard colors
    tags$head(tags$style(HTML("
      /* Final CSS Rules */

      /* 1. Top header bar to grey */
      .skin-blue .main-header .navbar,
      .skin-blue .main-header .logo {
        background-color: #6c757d !important;
      }

      /* 2. Parameter boxes (and Record box) to purple */
      .skin-blue .box.box-solid.box-primary {
        border-color: #6f42c1 !important;
      }
      .skin-blue .box.box-solid.box-primary > .box-header {
        color: #ffffff !important;
        background-color: #6f42c1 !important;
      }

      /* 3a. MODIFIED: MSP Results box with light purple HEADER ONLY */
      .skin-blue .box.box-solid.box-info {
        border: 1px solid #d3bce3 !important;
        background-color: #ffffff !important; /* White box body */
      }
      .skin-blue .box.box-solid.box-info > .box-header {
        color: #333333 !important;
        background-color: #f2e8f9 !important; /* Light purple header */
      }

      /* 3b. MODIFIED: Infection Results box with light orange HEADER ONLY */
       .skin-blue .box.box-solid.box-success {
        border: 1px solid #f7d3ba !important;
        background-color: #ffffff !important; /* White box body */
      }
      .skin-blue .box.box-solid.box-success > .box-header {
        color: #333333 !important;
        background-color: #fce5d4 !important; /* Light orange header */
      }
      
      /* 4. Sidebar background to light grey with black text */
      .skin-blue .main-sidebar {
        background-color: #f8f9fa !important;
      }
      .skin-blue .sidebar-menu > li > a {
        color: #333333 !important;
      }
      .skin-blue .sidebar-menu > li.header {
         color: #6f42c1 !important;
         background: #f8f9fa !important;
         font-weight: bold;
      }
      .skin-blue .sidebar-menu > li.active > a {
        background: #ffffff !important;
        border-left-color: #6f42c1 !important;
      }
      .skin-blue .sidebar-menu > li:hover > a {
        background: #ffffff !important;
      }
    "))),
    
    # Define the content for each tabItem
    tabItems(
      # Main dashboard content
      tabItem(tabName = "dashboard",
              # -- Row 1: Parameter Inputs --
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
              
              # -- Row 2: Add Record Button --
              fluidRow(
                column(width = 12, align = "center",
                       actionButton("add_record", "Add Record", icon = icon("plus"), style = "color: #fff; background-color: #6f42c1; border-color: #6f42c1; margin-bottom: 15px;")
                )
              ),
              
              # -- Row 3: Records Table --
              fluidRow(
                box(title = "Parameter Records", status = "primary", solidHeader = TRUE, width = 12,
                    actionButton("remove_record", "Remove Selected Record", icon = icon("trash"), style = "margin-bottom: 10px;"),
                    DT::dataTableOutput("records_table")
                )
              ),
              
              # -- Row 4: Result Outputs --
              fluidRow(
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
      ),
      
      # author information tab content
      tabItem(tabName = "author",
              fluidRow(
                box(
                  title = "About This Application",
                  solidHeader = TRUE,
                  status = "primary",
                  width = 12,
                  
                  # You can edit the content below easily
                  h3("Author Information"),
                  p("This application was developed to compare clinical trial designs based on different endpoints."),
                  tags$ul(
                    tags$li(HTML("<b>Author:</b> Peiyu Liu")),
                    tags$li(HTML("<b>Affiliation:</b> Department of Biostatistics, University of Florida")),
                    tags$li(HTML("<b>Contact:</b> <a href='mailto:pyliu0620@outlook.com'>pyliu0620@outlook.com</a>"))
                  ),
                  
                  hr(), # A horizontal line for separation
                  
                  h3("Publication"),
                  p("  "), #TB Filled
                  p("  ") #TBD
                  
                )
              )
      )
    )
  )
)

# ==============================================================================
# 3. Shiny Server (Server Logic)
# ==============================================================================
server <- function(input, output, session) {
  
  # -- a) reactive value to store the records --
  records <- reactiveVal(data.frame())
  
  # -- b) for the "Add Record" button click --
  observeEvent(input$add_record, {
    isolate({
      new_record <- data.frame(
        lambda_max = input$lambda_max, lambda_min = input$lambda_min, mu_c = input$mu_c, prop_c = input$prop_c,
        gamma0_hl = input$gamma0_hl, B0 = input$B0, sigma_epsilon0 = input$sigma_epsilon0, PV = input$PV,
        P_hv = input$P_hv, Pc = input$Pc, Ps = input$Ps, T_val = input$T, theta = input$theta,
        alpha = input$alpha, power = input$power
      )
      records(bind_rows(records(), new_record))
    })
  })
  
  # -- c) Observer for the "Remove Record" button click --
  observeEvent(input$remove_record, {
    selected_rows <- input$records_table_rows_selected
    if (!is.null(selected_rows)) {
      current_records <- records()
      records(current_records[-selected_rows, , drop = FALSE])
    }
  })
  
  # -- d) Render the records table --
  output$records_table <- DT::renderDataTable({
    df <- records()
    if (nrow(df) > 0) {
      colnames(df) <- c("λ_max", "λ_min", "μ_c", "prop_c", "Half-life", "B0", "σ_ε", "PV", "P_hv(%)", "P_c", "Ps(%)", "T(days)", "θ", "α", "1-β")
    }
    DT::datatable(df, selection = 'single', rownames = FALSE, options = list(scrollX = TRUE, pageLength = 5))
  })
  
  # -- e) compute derived parameters from inputs --
  params <- reactive({
    req(input$lambda_max, input$lambda_min)
    lmax <- max(input$lambda_max, input$lambda_min)
    lmin <- min(input$lambda_max, input$lambda_min)
    lambda0 <- (lmax + lmin) / 2
    lambda_amp <- (lmax - lmin) / 2
    omega <- 2 * pi / 365
    gamma <- log(2) / input$gamma0_hl
    var_c <- input$mu_c * input$prop_c
    lognorm_params <- find_lognorm_params(input$mu_c, var_c)
    list(
      lambda0 = lambda0, lambda_amp = lambda_amp, omega = omega, gamma = gamma,
      mu_c = input$mu_c, var_c = var_c, E_C2 = var_c + input$mu_c^2,
      mu_c_ln = lognorm_params$mu, sigma_c_sq_ln = lognorm_params$sigma_sq,
      B0 = input$B0, sigma_epsilon0 = input$sigma_epsilon0, PV = input$PV,
      P_hv = input$P_hv / 100, Pc = input$Pc, Ps = input$Ps / 100, T_val = input$T,
      theta = input$theta, alpha = input$alpha, power = input$power
    )
  })
  
  # -- f) eventReactive to calculate plot data only when "Add Record" is clicked --
  plot_data <- eventReactive(input$add_record, {
    p <- params()
    t_seq <- 1:p$T_val
    lambda_B_t_control <- p$lambda0 + p$lambda_amp * sin(p$omega * t_seq)
    lambda_B_t_interv <- lambda_B_t_control * p$theta
    h_inf <- function(t, lambda0, lambda_amp, gamma, omega) { (lambda0 / gamma) + (lambda_amp / (gamma^2 + omega^2)) * (gamma * sin(omega * t) - omega * cos(omega * t)) }
    g_inf <- function(t, lambda0, lambda_amp, gamma, omega) { (lambda0 / (2 * gamma)) + (lambda_amp / (4 * gamma^2 + omega^2)) * (2 * gamma * sin(omega * t) - omega * cos(omega * t)) }
    mean_control <- p$mu_c * h_inf(t_seq, p$lambda0, p$lambda_amp, p$gamma, p$omega) + p$B0
    mean_interv  <- p$mu_c * h_inf(t_seq, p$lambda0 * p$theta, p$lambda_amp * p$theta, p$gamma, p$omega) + p$B0
    var_control <- p$E_C2 * g_inf(t_seq, p$lambda0, p$lambda_amp, p$gamma, p$omega) + p$var_c * h_inf(t_seq, p$lambda0, p$lambda_amp, p$gamma, p$omega)^2 + p$sigma_epsilon0^2
    var_interv  <- p$E_C2 * g_inf(t_seq, p$lambda0 * p$theta, p$lambda_amp * p$theta, p$gamma, p$omega) + p$var_c * h_inf(t_seq, p$lambda0 * p$theta, p$lambda_amp * p$theta, p$gamma, p$omega)^2 + p$sigma_epsilon0^2
    ss_t1_per_group <- mapply(calculate_sample_size_ttest, mean_control, var_control, mean_interv, var_interv, MoreArgs = list(power = p$power, alpha = p$alpha))
    H_t <- function(t, lambda0, lambda_amp, omega) { lambda0 * t - (lambda_amp / omega) * (cos(omega * t) - 1) }
    cum_hazard_control <- p$PV * p$P_hv * H_t(t_seq, p$lambda0, p$lambda_amp, p$omega) * p$Pc
    cum_hazard_interv <- cum_hazard_control * p$theta
    P_c_incid <- p$Ps * (1 - exp(-cum_hazard_control)) 
    P_i_incid <- p$Ps * (1 - exp(-cum_hazard_interv))
    d <- (2 * (qnorm(1 - p$alpha / 2) + qnorm(p$power))^2) / (log(p$theta))^2
    list(
      df_lambda_B = bind_rows(data.frame(time = t_seq, rate = lambda_B_t_control, group = "Control"), data.frame(time = t_seq, rate = lambda_B_t_interv, group = "Intervention")),
      df_pop_traj = bind_rows(data.frame(time = t_seq, mean = mean_control, group = "Control"), data.frame(time = t_seq, mean = mean_interv, group = "Intervention")),
      df_ss_t1 = data.frame(time = t_seq, sample_size = ss_t1_per_group * 2),
      df_hazard_I = bind_rows(data.frame(time = t_seq, hazard = p$PV * p$P_hv * lambda_B_t_control, group = "Control"), data.frame(time = t_seq, hazard = p$PV * p$P_hv * lambda_B_t_control * p$theta, group = "Intervention")),
      df_ci = bind_rows(data.frame(time = t_seq, ci = P_c_incid, group = "Control"), data.frame(time = t_seq, ci = P_i_incid, group = "Intervention")),
      df_ss_t2 = data.frame(time = t_seq, sample_size = d / (P_c_incid + P_i_incid))
    )
  }, ignoreNULL = FALSE)
  
  # -- g) eventReactive to simulate individual trajectory --
  indiv_traj_data <- eventReactive(input$add_record, {
    p <- params()
    t_end_sim <- 1
    sim_times <- seq(0, t_end_sim, length.out = 500)
    set.seed(123)
    bite_times <- generate_event_times(p$lambda0, p$lambda_amp, t_end_sim, p$omega)
    num_bites <- length(bite_times)
    C_values <- if (num_bites > 0) rlnorm(num_bites, meanlog = p$mu_c_ln, sdlog = sqrt(p$sigma_c_sq_ln)) else numeric(0)
    epsilon <- rnorm(1, 0, p$sigma_epsilon0)
    antibody_values <- calculate_antibody_trajectory(bite_times, C_values, p$B0, epsilon, p$gamma, sim_times)
    traj_df <- data.frame(time = sim_times, antibody = antibody_values)
    bite_df <- if (num_bites > 0) {
      data.frame(time = bite_times, antibody = calculate_antibody_trajectory(bite_times, C_values, p$B0, epsilon, p$gamma, bite_times))
    } else { data.frame(time = numeric(0), antibody = numeric(0)) }
    list(traj_df = traj_df, bite_df = bite_df)
  }, ignoreNULL = FALSE)
  
  # --- h) Render all plots ---
  output$plot_lambda_B <- renderPlot({
    ggplot(plot_data()$df_lambda_B, aes(x = time, y = rate, color = group)) +
      geom_line(size = 1) +
      scale_color_manual(values = c("Control" = "blue", "Intervention" = "orange")) +
      labs(title = "Daily Mosquito Biting Rate", x = "Time (days)", y = "Bites per person per day", color = "Group") +
      theme_minimal() + theme(legend.position = "bottom")
  })
  output$plot_indiv_traj <- renderPlot({ data <- indiv_traj_data(); ggplot(data$traj_df, aes(x = time, y = antibody)) + geom_line(color = "black", size = 1) + geom_point(data = data$bite_df, aes(x = time, y = antibody), color = "red", size = 4, shape = 18) + geom_hline(yintercept = params()$B0, linetype = "dashed", color = "gray50") + annotate("text", x = 0.5, y = params()$B0, label = "B0", vjust = -0.5, color = "gray50") + labs(title = "Antibody Trajectory for One Simulated Individual Over One Day", x = "Time (days)", y = "Antibody Level") + theme_minimal() })
  output$plot_pop_traj <- renderPlot({ ggplot(plot_data()$df_pop_traj, aes(x = time, y = mean, color = group)) + geom_line(size = 1) + scale_color_manual(values = c("Control" = "blue", "Intervention" = "orange")) + labs(title = "Theoretical Population Mean Antibody Trajectory", x = "Time (days)", y = "Mean Antibody Level", color = "Group") + theme_minimal() + theme(legend.position = "bottom") })
  output$plot_ss_trial1 <- renderPlot({ ggplot(plot_data()$df_ss_t1, aes(x = time, y = sample_size)) + geom_line(color = "purple", size = 1) + scale_y_log10() + labs(title = "Total Sample Size (Antibody Endpoint)", subtitle = "Calculated at each time point", x = "Time (days)", y = "Total Sample Size (Log Scale)") + theme_minimal() })
  output$plot_hazard_I <- renderPlot({ ggplot(plot_data()$df_hazard_I, aes(x = time, y = hazard, color = group)) + geom_line(size = 1) + scale_color_manual(values = c("Control" = "blue", "Intervention" = "orange")) + labs(title = "Daily Infection Hazard Rate", x = "Time (days)", y = "Hazard Rate", color = "Group") + theme_minimal() + theme(legend.position = "bottom") })
  output$plot_ci <- renderPlot({ ggplot(plot_data()$df_ci, aes(x = time, y = ci, color = group)) + geom_line(size = 1) + scale_y_continuous(labels = scales::percent) + scale_color_manual(values = c("Control" = "blue", "Intervention" = "orange")) + labs(title = "Cumulative Incidence of Infection", x = "Time (days)", y = "Cumulative Incidence", color = "Group") + theme_minimal() + theme(legend.position = "bottom") })
  output$plot_ss_trial2 <- renderPlot({ ggplot(plot_data()$df_ss_t2, aes(x = time, y = sample_size)) + geom_line(color = "purple", size = 1) + scale_y_log10() + labs(title = "Total Sample Size (Infection Endpoint)", subtitle = "Calculated assuming trial ends at a given time point", x = "Time (days)", y = "Total Sample Size (Log Scale)") + theme_minimal() })
}

shinyApp(ui = ui, server = server)