# R/server_outputs.R
# server module2: the plotting helper functions and the output rendering logic.

# =============================================== 
# 1. plotting and table generation functions
# =============================================== 

#' an interactive DT datatable for parameter records.
#' @param df data frame of records.
#' @return DT datatable 
create_records_table <- function(df) {
  if (nrow(df) > 0) {
    colnames(df) <- c("λ_max", "λ_min", "μ_c", "prop_c", "Half-life", "B0", "σ_ε", "PV", "P_hv(%)", "P_c", "Ps(%)", "T(days)", "θ", "α", "1-β")
  }
  DT::datatable(
    df,
    selection = 'single',
    rownames = FALSE,
    options = list(
      scrollX = TRUE,
      pageLength = 5,
      # sort by the first column (latest record) descending
      order = list(list(0, 'desc'))
    )
  )
}

#' plot for daily mosquito biting rate.
#' @param df Data frame with columns: time, rate, group.
plot_biting_rate <- function(df) {
  ggplot(df, aes(x = time, y = rate, color = group)) +
    geom_line(size = 1) +
    scale_color_manual(values = c("Control" = "blue", "Intervention" = "orange")) +
    labs(title = "Daily Mosquito Biting Rate", x = "Time (days)", y = "Bites per person per day", color = "Group") +
    theme_minimal() +
    theme(legend.position = "bottom")
}

#' a plot for a single individual's antibody trajectory.
#' @param data A list of traj_df (trajectory line) and bite_df (bite points).
#' @param B0 The baseline antibody level for the horizontal line.
plot_individual_trajectory <- function(data, B0) {
  ggplot(data$traj_df, aes(x = time, y = antibody)) +
    geom_line(color = "black", size = 1) +
    geom_point(data = data$bite_df, aes(x = time, y = antibody), color = "red", size = 4, shape = 18) +
    geom_hline(yintercept = B0, linetype = "dashed", color = "gray50") +
    annotate("text", x = 0.5, y = B0, label = "B0", vjust = -0.5, color = "gray50") +
    labs(title = "Antibody Trajectory for One Simulated Individual", x = "Time (days)", y = "Antibody Level") +
    theme_minimal()
}

#' plot for the theoretical population mean antibody trajectory.
#' @param df data frame with columns: time, mean, group.
plot_population_trajectory <- function(df) {
  ggplot(df, aes(x = time, y = mean, color = group)) +
    geom_line(size = 1) +
    scale_color_manual(values = c("Control" = "blue", "Intervention" = "orange")) +
    labs(title = "Theoretical Population Mean Antibody Trajectory", x = "Time (days)", y = "Mean Antibody Level", color = "Group") +
    theme_minimal() +
    theme(legend.position = "bottom")
}

#' plot for the daily infection hazard rate.
#' @param df data frame with columns: time, hazard, group.
plot_infection_hazard <- function(df) {
  ggplot(df, aes(x = time, y = hazard, color = group)) +
    geom_line(size = 1) +
    scale_color_manual(values = c("Control" = "blue", "Intervention" = "orange")) +
    labs(title = "Daily Infection Hazard Rate", x = "Time (days)", y = "Hazard Rate", color = "Group") +
    theme_minimal() +
    theme(legend.position = "bottom")
}

#' plot for the cumulative incidence of infection.
#' @param df data frame with columns: time, ci, group.
plot_cumulative_incidence <- function(df) {
  ggplot(df, aes(x = time, y = ci, color = group)) +
    geom_line(size = 1) +
    scale_y_continuous(labels = scales::percent) +
    scale_color_manual(values = c("Control" = "blue", "Intervention" = "orange")) +
    labs(title = "Cumulative Incidence of Infection", x = "Time (days)", y = "Cumulative Incidence", color = "Group") +
    theme_minimal() +
    theme(legend.position = "bottom")
}

#' plot for total sample size over time (reusable).
#' @param df Data frame with columns: time, sample_size.
plot_sample_size <- function(df, main_title, subtitle_text) {
  ggplot(df, aes(x = time, y = sample_size)) +
    geom_line(color = "purple", size = 1) +
    scale_y_log10() +
    labs(
      title = main_title,
      subtitle = subtitle_text,
      x = "Time (days)",
      y = "Total Sample Size (Log Scale)"
    ) +
    theme_minimal()
}


# ============================= 
# 2. output rendering 
# =============================

servermod_output_render <- function(output, reactives) {
  
  # --- render records Table ---
  output$records_table <- DT::renderDataTable({
    create_records_table(reactives$records())
  })
  
  # --- render MSP Antibody-based Trial Plots ---
  output$plot_lambda_B <- renderPlot({
    plot_biting_rate(reactives$plot_data()$df_lambda_B)
  }) %>% bindCache(reactives$plot_data())
  
  output$plot_indiv_traj <- renderPlot({
    plot_individual_trajectory(
      data = reactives$indiv_traj_data(),
      B0 = reactives$params()$B0
    )
  }) %>% bindCache(reactives$indiv_traj_data(), reactives$params()$B0)
  
  output$plot_pop_traj <- renderPlot({
    plot_population_trajectory(reactives$plot_data()$df_pop_traj)
  }) %>% bindCache(reactives$plot_data())
  
  output$plot_ss_trial1 <- renderPlot({
    plot_sample_size(
      df = reactives$plot_data()$df_ss_t1,
      main_title = "Total Sample Size (Antibody Endpoint)",
      subtitle_text = "Calculated at each time point"
    )
  }) %>% bindCache(reactives$plot_data())
  
  # --- render Infection-based Trial Plots ---
  output$plot_hazard_I <- renderPlot({
    plot_infection_hazard(reactives$plot_data()$df_hazard_I)
  }) %>% bindCache(reactives$plot_data())
  
  output$plot_ci <- renderPlot({
    plot_cumulative_incidence(reactives$plot_data()$df_ci)
  }) %>% bindCache(reactives$plot_data())
  
  output$plot_ss_trial2 <- renderPlot({
    plot_sample_size(
      df = reactives$plot_data()$df_ss_t2,
      main_title = "Total Sample Size (Infection Endpoint)",
      subtitle_text = "Calculated assuming trial ends at a given time point"
    )
  }) %>% bindCache(reactives$plot_data())
}