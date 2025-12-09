# R/server_inputs.R
# server module1: handles inputs, reactive states, and core computations for ploting data.

# ==============================================================================
# 1. input data calculation functions for output data ploting
# ==============================================================================

#------computing functions from original sript-------

find_lognorm_params <- function(target_mean, target_var) {
  if (target_mean <= 0 || target_var <= 0) { return(list(mu = NA, sigma_sq = NA)) }
  cv_sq <- target_var / (target_mean^2)
  if (cv_sq < -1) return(list(mu = NA, sigma_sq = NA))
  sigma_sq <- log(1 + cv_sq)
  mu <- log(target_mean) - sigma_sq / 2
  return(list(mu = mu, sigma_sq = sigma_sq))
}
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
generate_event_times <- function(lambda0, lambda_amp, t_end, omega) {
  lambda_max <- lambda0 + lambda_amp
  n_candidate <- rpois(1, lambda_max * t_end)
  if (n_candidate == 0) return(numeric(0))
  candidate_times <- sort(runif(n_candidate, 0, t_end))
  accept_prob <- (lambda0 + lambda_amp * sin(omega * candidate_times)) / lambda_max
  return(candidate_times[runif(n_candidate) <= accept_prob])
}
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

find_lognorm_params <- function(target_mean, target_var) {
  if (target_mean <= 0 || target_var <= 0) { return(list(mu = NA, sigma_sq = NA)) }
  cv_sq <- target_var / (target_mean^2)
  if (cv_sq < -1) return(list(mu = NA, sigma_sq = NA))
  sigma_sq <- log(1 + cv_sq)
  mu <- log(target_mean) - sigma_sq / 2
  return(list(mu = mu, sigma_sq = sigma_sq))
}

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
generate_event_times <- function(lambda0, lambda_amp, t_end, omega) {
  lambda_max <- lambda0 + lambda_amp
  n_candidate <- rpois(1, lambda_max * t_end)
  if (n_candidate == 0) return(numeric(0))
  candidate_times <- sort(runif(n_candidate, 0, t_end))
  accept_prob <- (lambda0 + lambda_amp * sin(omega * candidate_times)) / lambda_max
  return(candidate_times[runif(n_candidate) <= accept_prob])
}
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


#' calculate outputdata for all main plots based on parameters
#' @param p a list of parameters from the debounced reactive 
#' @return a list of data frames for plotting 
calculate_main_plot_data <- function(p) {
  t_seq <- 1:p$T_val
  
  # biting Rate
  lambda_B_t_control <- p$lambda0 + p$lambda_amp * sin(p$omega * t_seq)
  lambda_B_t_interv <- lambda_B_t_control * p$theta
  
  # functions for antibody calculations
  h_inf <- function(t, lambda0, lambda_amp, gamma, omega) { (lambda0 / gamma) + (lambda_amp / (gamma^2 + omega^2)) * (gamma* sin(omega * t) - omega * cos(omega * t)) }
  g_inf <- function(t, lambda0, lambda_amp, gamma, omega) { (lambda0 / (2*gamma)) + (lambda_amp / (4*gamma^2 + omega^2)) * (2*gamma * sin(omega * t) - omega * cos(omega * t)) }
  
  # Antibody Mean and Variance
  mean_control <- p$mu_c * h_inf(t_seq, p$lambda0, p$lambda_amp, p$gamma, p$omega) + p$B0
  mean_interv  <- p$mu_c * h_inf(t_seq, p$lambda0 * p$theta, p$lambda_amp * p$theta, p$gamma, p$omega) + p$B0
  var_control <- p$E_C2 * g_inf(t_seq, p$lambda0, p$lambda_amp, p$gamma, p$omega) + p$var_c * h_inf(t_seq, p$lambda0, p$lambda_amp, p$gamma, p$omega)^2 + p$sigma_epsilon0^2
  var_interv  <- p$E_C2 * g_inf(t_seq, p$lambda0 * p$theta, p$lambda_amp * p$theta, p$gamma, p$omega) + p$var_c * h_inf(t_seq, p$lambda0 * p$theta, p$lambda_amp * p$theta, p$gamma, p$omega)^2 + p$sigma_epsilon0^2
  
  # sample Size (Trial 1)
  ss_t1_per_group <- mapply(calculate_sample_size_ttest, mean_control, var_control, mean_interv, var_interv, MoreArgs = list(power = p$power, alpha = p$alpha))
  
  # Infection hazard and incidence
  H_t <- function(t, lambda0, lambda_amp, omega) { lambda0 * t - (lambda_amp / omega) * (cos(omega * t) - 1) }
  cum_hazard_control <- p$PV * p$P_hv * H_t(t_seq, p$lambda0, p$lambda_amp, p$omega) * p$Pc
  cum_hazard_interv <- cum_hazard_control * p$theta
  P_c_incid <- p$Ps * (1 - exp(-cum_hazard_control))
  P_i_incid <- p$Ps * (1 - exp(-cum_hazard_interv))
  
  # sample Size (Trial 2)
  d_events <- (2 * (qnorm(1 - p$alpha / 2) + qnorm(p$power))^2) / (log(p$theta))^2
  
  list(
    df_lambda_B = bind_rows(data.frame(time = t_seq, rate = lambda_B_t_control, group = "Control"), data.frame(time = t_seq, rate = lambda_B_t_interv, group = "Intervention")),
    df_pop_traj = bind_rows(data.frame(time = t_seq, mean = mean_control, group = "Control"), data.frame(time = t_seq, mean = mean_interv, group = "Intervention")),
    df_ss_t1 = data.frame(time = t_seq, sample_size = ss_t1_per_group * 2),
    df_hazard_I = bind_rows(data.frame(time = t_seq, hazard = p$PV * p$P_hv * lambda_B_t_control, group = "Control"), data.frame(time = t_seq, hazard = p$PV * p$P_hv * lambda_B_t_control * p$theta, group = "Intervention")),
    df_ci = bind_rows(data.frame(time = t_seq, ci = P_c_incid, group = "Control"), data.frame(time = t_seq, ci = P_i_incid, group = "Intervention")),
    df_ss_t2 = data.frame(time = t_seq, sample_size = d_events / (P_c_incid + P_i_incid))
  )
}

#' antibody trajectory for a single simulated individual
#' @param p ist of parameters from the debounced reactive.
#' @return  list containing two data frames: one for the trajectory line, one for bite points.
calculate_individual_trajectory <- function(p) {
  t_end_sim <- 1 # Simulate over one day
  sim_times <- seq(0, t_end_sim, length.out = 500)
  set.seed(123) # for reproducibility
  
  bite_times <- generate_event_times(p$lambda0, p$lambda_amp, t_end_sim, p$omega)
  num_bites <- length(bite_times)
  
  C_values <- if (num_bites > 0) rlnorm(num_bites, meanlog = p$mu_c_ln, sdlog = sqrt(p$sigma_c_sq_ln)) else numeric(0)
  epsilon <- rnorm(1, 0, p$sigma_epsilon0)
  
  antibody_values <- calculate_antibody_trajectory(bite_times, C_values, p$B0, epsilon, p$gamma, sim_times)
  
  traj_df <- data.frame(time = sim_times, antibody = antibody_values)
  bite_df <- if (num_bites > 0) {
    data.frame(
      time = bite_times, 
      antibody = calculate_antibody_trajectory(bite_times, C_values, p$B0, epsilon, p$gamma, bite_times)
    )
  } else { 
    data.frame(time = numeric(0), antibody = numeric(0)) 
  }
  
  list(traj_df = traj_df, bite_df = bite_df)
}

# ======================================== 
# 2. server input/computing module
# ========================================

servermod_inputs_computing <- function(input) {
  
  # --- 1. reactive State Management ---
  
  # stores the history of parameter records
  records <- reactiveVal(data.frame())
  
  # add a new record to the table when the button is clicked
  observeEvent(input$add_record, {
    # iensures this block only reads input values at the moment of the click
    isolate({
      new_record <- data.frame(
        lambda_max = input$lambda_max, lambda_min = input$lambda_min, mu_c = input$mu_c, prop_c = input$prop_c,
        gamma0_hl = input$gamma0_hl, B0 = input$B0, sigma_epsilon0 = input$sigma_epsilon0, PV = input$PV,
        P_hv = input$P_hv, Pc = input$Pc, Ps = input$Ps, T_val = input$T, theta = input$theta,
        alpha = input$alpha, power = input$power
      )
      # prepend the new record to the existing data frame
      records(bind_rows(new_record, records()))
    })
  })
  
  # remove the selected record from the table
  observeEvent(input$remove_record, {
    selected_rows <- input$records_table_rows_selected
    if (!is.null(selected_rows)) {
      current_records <- records()
      records(current_records[-selected_rows, , drop = FALSE])
    }
  })
  
  # --- 2. core computing logic (data processing & model parameters) ---
  
  # reactive expression to derive all model parameters from user inputs.
  # debounced to prevent recalculation on every single keystroke.
  params <- reactive({
    req(input$lambda_max, input$lambda_min, input$gamma0_hl > 0) # Ensure key inputs are valid
    lmax <- max(input$lambda_max, input$lambda_min)
    lmin <- min(input$lambda_max, input$lambda_min)
    
    lognorm_params <- find_lognorm_params(input$mu_c, input$mu_c * input$prop_c)
    
    list(
      lambda0 = (lmax + lmin) / 2,
      lambda_amp = (lmax - lmin) / 2,
      omega = 2 * pi / 365,
      gamma = log(2) / input$gamma0_hl,
      mu_c = input$mu_c,
      var_c = input$mu_c * input$prop_c,
      E_C2 = (input$mu_c * input$prop_c) + input$mu_c^2,
      mu_c_ln = lognorm_params$mu,
      sigma_c_sq_ln = lognorm_params$sigma_sq,
      B0 = input$B0,
      sigma_epsilon0 = input$sigma_epsilon0,
      PV = input$PV,
      P_hv = input$P_hv / 100,
      Pc = input$Pc,
      Ps = input$Ps / 100,
      T_val = input$T,
      theta = input$theta,
      alpha = input$alpha,
      power = input$power
    )
  }) %>% debounce(800) # wait 800ms after user stops typing before re-calculating
  
  # --- 3. event-driven computations for outPlots ---
  
  # plots are only re-calculated when the user clicks 'Add Record'.
  # use the latest values from the debounced `params` reactive.
  plot_data <- eventReactive(input$add_record, {
    calculate_main_plot_data(params())
  }, ignoreNULL = FALSE)
  
  indiv_traj_data <- eventReactive(input$add_record, {
    calculate_individual_trajectory(params())
  }, ignoreNULL = FALSE)
  
  
  # --- 4. return all reactive objects for the output module to use ---
  
  return(
    list(
      records = records,
      params = params,
      plot_data = plot_data,
      indiv_traj_data = indiv_traj_data
    )
  )
}