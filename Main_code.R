library(DEoptim)
library(mvtnorm)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# 1. OBJECTIVE FUNCTION THAT RETURNS LSE

varmax_sse <- function(parameters, 
                       y, 
                       X, 
                       restricted = FALSE,
                       residuals = FALSE) {
  
  # X should have two columns: Won and TotalAmount
  x <- X  # Keep both predictors
  
  # Extract parameters according to supervisor's specification
  delta <- parameters[1:2]          # intercepts
  
  # Psi matrix (2x2) - flatten by row
  Psi <- matrix(parameters[3:6], nrow = 2, ncol = 2, byrow = TRUE)
  
  # Theta matrix (diagonal 2x2): autoregressive effects
  Theta <- diag(parameters[7:8])
  
  if(restricted) {
    Phi <- -Theta
    # For SSE, we don't need sigma parameters, but we'll keep them for consistency
    # with the parameter count
    if (length(parameters) > 8) {
      Sigma_params <- parameters[9:11]
    }
  } else {
    Phi <- diag(parameters[9:10])
    if (length(parameters) > 10) {
      Sigma_params <- parameters[11:13]
    }
  }
  
  # Compute residuals
  epsilon <- matrix(0, nrow = nrow(y), ncol = 2)
  
  for(i in 1:nrow(y)) {
    if(i == 1) {
      # For first observation
      y_hat <- solve(diag(2) - Theta) %*% delta + Psi %*% x[i, ]
      epsilon[i, ] <- y[i, ] - y_hat
    } else {
      # For subsequent observations
      y_hat <- delta + Psi %*% x[i, ] + Theta %*% y[i-1, ] + Phi %*% epsilon[i-1, ]
      epsilon[i, ] <- y[i, ] - y_hat
    }
  }
  
  # Return residuals or sum of squared errors
  if(residuals) {
    return(epsilon)
  } else {
    sse <- sum(epsilon^2)
    return(sse)
  }
}

# 2. SET UP PARAMETER BOUNDS

# For unrestricted model (varmax): 11 parameters
lower_unrestricted <- c(
  0, 0,                    # delta1, delta2
  -5, -5, -5, -5,         # Psi11, Psi12, Psi21, Psi22
  -(1 - 1e-5), -(1 - 1e-5), # Theta1, Theta2
  -(1 - 1e-5), -(1 - 1e-5), # Phi1, Phi2
  1e-5, -1 + 1e-5, 1e-5    # sigma1, rho, sigma2
)

upper_unrestricted <- c(
  1, 1,                    # delta1, delta2
  5, 5, 5, 5,             # Psi11, Psi12, Psi21, Psi22
  1 - 1e-5, 1 - 1e-5,     # Theta1, Theta2
  1 - 1e-5, 1 - 1e-5,     # Phi1, Phi2
  0.25, 1 - 1e-5, 0.25     # sigma1, rho, sigma2
)

# For restricted model (varmax_restricted): 9 parameters
lower_restricted <- c(
  0, 0,                    # delta1, delta2
  -5, -5, -5, -5,         # Psi11, Psi12, Psi21, Psi22
  -(1 - 1e-5), -(1 - 1e-5), # Theta1, Theta2
  1e-5, -1 + 1e-5, 1e-5    # sigma1, rho, sigma2
)

upper_restricted <- c(
  1, 1,                    # delta1, delta2
  5, 5, 5, 5,             # Psi11, Psi12, Psi21, Psi22
  1 - 1e-5, 1 - 1e-5,     # Theta1, Theta2
  0.25, 1 - 1e-5, 0.25     # sigma1, rho, sigma2
)

# 3. FUNCTION TO CALCULATE AIC AND BIC FOR LEAST SQUARES

calculate_aic_bic_sse <- function(sse, n_params, n_obs, sigma2_est = NULL) {
  # Calculate AIC and BIC for least squares models
  
  # If sigma2_est is not provided, estimate it from SSE
  if (is.null(sigma2_est)) {
    # Estimate sigma^2 = SSE / (n - p) where p = number of parameters
    sigma2_est <- sse / (n_obs - n_params)
  }
  
  # AIC for least squares: n * ln(SSE/n) + 2k
  # But better: AIC = n * ln(SSE/n) + 2k + constant (n * (1 + ln(2π)))
  aic <- n_obs * log(sse/n_obs) + 2 * n_params
  
  # BIC for least squares: n * ln(SSE/n) + k * ln(n)
  bic <- n_obs * log(sse/n_obs) + n_params * log(n_obs)
  
  return(list(AIC = aic, BIC = bic, sigma2_est = sigma2_est))
}

# 4. LOAD AND PREPARE DATA

cat("=== LOADING AND PREPARING DATA ===\n")

# Load the data
data <- read.csv("~/Documents/Invatamant/Uni/Year_4/Thesis/Data/Data_CIAC.csv", stringsAsFactors = FALSE)

# Check the structure of the data
cat("Data structure:\n")
str(data)
cat("\nColumn names:\n")
print(names(data))

# Check for required columns
required_cols <- c("Ppn", "TrialNumber", "PAscore", "NAscore", "Won", "TotalAmount")
missing_cols <- setdiff(required_cols, names(data))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
}

# Convert to appropriate types if needed
data$Ppn <- as.factor(data$Ppn)
data$TrialNumber <- as.integer(data$TrialNumber)

# Get unique participants
participants <- unique(data$Ppn)
n_participants <- length(participants)
cat(sprintf("\nNumber of participants: %d\n", n_participants))

# 5. FUNCTION TO FIT MODELS FOR A SINGLE PARTICIPANT (USING LSE)

fit_participant_models_sse <- function(participant_data, ppn_id) {
  cat(sprintf("\nFitting models for participant %s...\n", ppn_id))
  
  # Prepare data
  y <- as.matrix(participant_data[, c("PAscore", "NAscore")])
  X <- as.matrix(participant_data[, c("Won", "TotalAmount")])
  n_obs <- nrow(y)
  
  # Initialize results storage
  results <- list(
    participant = ppn_id,
    n_obs = n_obs,
    restricted = list(),
    unrestricted = list()
  )
  
  # Fit restricted model
  cat("  Fitting restricted model...\n")
  tryCatch({
    de_restricted <- DEoptim(
      fn = function(params) {
        varmax_sse(
          parameters = params,
          y = y,
          X = X,
          restricted = TRUE,
          residuals = FALSE
        )
      },
      lower = lower_restricted,
      upper = upper_restricted,
      control = DEoptim.control(
        itermax = 500,  
        trace = FALSE,
        strategy = 2,
        NP = 200,       
        parallelType = 0,
        reltol = 1e-8,
        steptol = 100
      )
    )
    
    # Extract results for restricted model
    best_params_restricted <- de_restricted$optim$bestmem
    sse_restricted <- de_restricted$optim$bestval
    
    # Calculate AIC and BIC for SSE
    criteria_restricted <- calculate_aic_bic_sse(
      sse = sse_restricted,
      n_params = length(lower_restricted),
      n_obs = n_obs
    )
    
    # Get residuals
    residuals_restricted <- varmax_sse(
      parameters = best_params_restricted,
      y = y,
      X = X,
      restricted = TRUE,
      residuals = TRUE
    )
    
    # Extract parameter estimates with names (only keep first 8 for restricted)
    param_names_restricted <- c(
      "delta1", "delta2",
      "Psi11", "Psi12", "Psi21", "Psi22",
      "Theta1", "Theta2",
      "sigma1", "rho", "sigma2"  # Keep for consistency but won't be interpreted
    )
    
    # Only store the parameters we actually use (first 8)
    stored_params_restricted <- best_params_restricted[1:8]
    names(stored_params_restricted) <- param_names_restricted[1:8]
    
    results$restricted <- list(
      parameters = stored_params_restricted,
      sse = sse_restricted,
      AIC = criteria_restricted$AIC,
      BIC = criteria_restricted$BIC,
      sigma2_est = criteria_restricted$sigma2_est,
      residuals = residuals_restricted,
      convergence = TRUE
    )
    
  }, error = function(e) {
    cat(sprintf("    ERROR in restricted model: %s\n", e$message))
    results$restricted$convergence <- FALSE
  })
  
  # Fit unrestricted model
  cat("  Fitting unrestricted model...\n")
  tryCatch({
    de_unrestricted <- DEoptim(
      fn = function(params) {
        varmax_sse(
          parameters = params,
          y = y,
          X = X,
          restricted = FALSE,
          residuals = FALSE
        )
      },
      lower = lower_unrestricted,
      upper = upper_unrestricted,
      control = DEoptim.control(
        itermax = 500,
        trace = FALSE,
        strategy = 2,
        NP = 200,
        parallelType = 0,
        reltol = 1e-8,
        steptol = 100
      )
    )
    
    # Extract results for unrestricted model
    best_params_unrestricted <- de_unrestricted$optim$bestmem
    sse_unrestricted <- de_unrestricted$optim$bestval
    
    # Calculate AIC and BIC for SSE
    criteria_unrestricted <- calculate_aic_bic_sse(
      sse = sse_unrestricted,
      n_params = length(lower_unrestricted),
      n_obs = n_obs
    )
    
    # Get residuals
    residuals_unrestricted <- varmax_sse(
      parameters = best_params_unrestricted,
      y = y,
      X = X,
      restricted = FALSE,
      residuals = TRUE
    )
    
    # Extract parameter estimates with names (only keep first 10 for unrestricted)
    param_names_unrestricted <- c(
      "delta1", "delta2",
      "Psi11", "Psi12", "Psi21", "Psi22",
      "Theta1", "Theta2",
      "Phi1", "Phi2",
      "sigma1", "rho", "sigma2"  # Keep for consistency but won't be interpreted
    )
    
    # Only store the parameters we actually use (first 10)
    stored_params_unrestricted <- best_params_unrestricted[1:10]
    names(stored_params_unrestricted) <- param_names_unrestricted[1:10]
    
    results$unrestricted <- list(
      parameters = stored_params_unrestricted,
      sse = sse_unrestricted,
      AIC = criteria_unrestricted$AIC,
      BIC = criteria_unrestricted$BIC,
      sigma2_est = criteria_unrestricted$sigma2_est,
      residuals = residuals_unrestricted,
      convergence = TRUE
    )
    
  }, error = function(e) {
    cat(sprintf("    ERROR in unrestricted model: %s\n", e$message))
    results$unrestricted$convergence <- FALSE
  })
  
  return(results)
}

# 6. RUN ANALYSIS FOR ALL PARTICIPANTS

cat("\n=== FITTING MODELS FOR ALL PARTICIPANTS (USING SSE) ===\n")

all_results <- list()
model_comparison <- data.frame()

# Progress bar
pb <- txtProgressBar(min = 0, max = n_participants, style = 3)

for (i in 1:n_participants) {
  ppn <- participants[i]
  
  # Extract participant data
  participant_data <- data[data$Ppn == ppn, ]
  participant_data <- participant_data[order(participant_data$TrialNumber), ]
  
  # Fit models
  results <- fit_participant_models_sse(participant_data, as.character(ppn))
  
  # Store results
  all_results[[i]] <- results
  
  # Add to comparison table if both models converged
  if (!is.null(results$restricted$convergence) && results$restricted$convergence &&
      !is.null(results$unrestricted$convergence) && results$unrestricted$convergence) {
    
    comparison_row <- data.frame(
      Participant = as.character(ppn),
      n_obs = results$n_obs,
      SSE_restricted = results$restricted$sse,
      SSE_unrestricted = results$unrestricted$sse,
      AIC_restricted = results$restricted$AIC,
      AIC_unrestricted = results$unrestricted$AIC,
      BIC_restricted = results$restricted$BIC,
      BIC_unrestricted = results$unrestricted$BIC,
      Delta_AIC = results$unrestricted$AIC - results$restricted$AIC,
      Delta_BIC = results$unrestricted$BIC - results$restricted$BIC,
      Preferred_AIC = ifelse(results$restricted$AIC < results$unrestricted$AIC, 
                             "Restricted", "Unrestricted"),
      Preferred_BIC = ifelse(results$restricted$BIC < results$unrestricted$BIC, 
                             "Restricted", "Unrestricted"),
      sigma2_restricted = results$restricted$sigma2_est,
      sigma2_unrestricted = results$unrestricted$sigma2_est
    )
    model_comparison <- rbind(model_comparison, comparison_row)
  }
  
  setTxtProgressBar(pb, i)
}

close(pb)

# 7. CREATE DESCRIPTIVE STATISTICS TABLES

cat("\n=== CREATING DESCRIPTIVE STATISTICS ===\n")

# Extract parameter estimates for all participants
extract_parameters_table <- function(all_results, model_type = "restricted") {
  param_table <- data.frame()
  
  for (i in seq_along(all_results)) {
    results <- all_results[[i]]
    
    if (!is.null(results[[model_type]]$convergence) && results[[model_type]]$convergence) {
      params <- results[[model_type]]$parameters
      
      # Create data frame row
      row_data <- as.list(params)
      row_data$Participant <- results$participant
      row_data$SSE <- results[[model_type]]$sse
      row_data$AIC <- results[[model_type]]$AIC
      row_data$BIC <- results[[model_type]]$BIC
      row_data$sigma2_est <- results[[model_type]]$sigma2_est
      
      param_table <- rbind(param_table, data.frame(row_data))
    }
  }
  
  return(param_table)
}

# Create parameter tables
restricted_params <- extract_parameters_table(all_results, "restricted")
unrestricted_params <- extract_parameters_table(all_results, "unrestricted")

# Create descriptive statistics
create_descriptive_stats <- function(param_table, model_name) {
  # Exclude non-numeric columns for statistics
  numeric_cols <- sapply(param_table, is.numeric)
  numeric_data <- param_table[, numeric_cols]
  
  # Remove Participant column if present
  if ("Participant" %in% names(numeric_data)) {
    numeric_data <- numeric_data[, !names(numeric_data) %in% "Participant"]
  }
  
  # Calculate statistics
  stats <- data.frame(
    Parameter = names(numeric_data),
    Mean = apply(numeric_data, 2, mean, na.rm = TRUE),
    SD = apply(numeric_data, 2, sd, na.rm = TRUE),
    Median = apply(numeric_data, 2, median, na.rm = TRUE),
    Min = apply(numeric_data, 2, min, na.rm = TRUE),
    Max = apply(numeric_data, 2, max, na.rm = TRUE),
    Q1 = apply(numeric_data, 2, quantile, probs = 0.25, na.rm = TRUE),
    Q3 = apply(numeric_data, 2, quantile, probs = 0.75, na.rm = TRUE)
  )
  
  stats$Model <- model_name
  return(stats)
}

restricted_stats <- create_descriptive_stats(restricted_params, "Restricted")
unrestricted_stats <- create_descriptive_stats(unrestricted_params, "Unrestricted")

# Combine statistics
all_stats <- rbind(restricted_stats, unrestricted_stats)

# 8. CREATE SUMMARY TABLES

cat("\n=== MODEL COMPARISON SUMMARY ===\n")

# Model selection summary
model_selection_summary <- data.frame(
  n_participants = nrow(model_comparison),
  Restricted_preferred_AIC = sum(model_comparison$Preferred_AIC == "Restricted", na.rm = TRUE),
  Unrestricted_preferred_AIC = sum(model_comparison$Preferred_AIC == "Unrestricted", na.rm = TRUE),
  Restricted_preferred_BIC = sum(model_comparison$Preferred_BIC == "Restricted", na.rm = TRUE),
  Unrestricted_preferred_BIC = sum(model_comparison$Preferred_BIC == "Unrestricted", na.rm = TRUE),
  Mean_Delta_AIC = mean(model_comparison$Delta_AIC, na.rm = TRUE),
  SD_Delta_AIC = sd(model_comparison$Delta_AIC, na.rm = TRUE),
  Mean_Delta_BIC = mean(model_comparison$Delta_BIC, na.rm = TRUE),
  SD_Delta_BIC = sd(model_comparison$Delta_BIC, na.rm = TRUE),
  Mean_SSE_restricted = mean(model_comparison$SSE_restricted, na.rm = TRUE),
  Mean_SSE_unrestricted = mean(model_comparison$SSE_unrestricted, na.rm = TRUE),
  Mean_sigma2_restricted = mean(model_comparison$sigma2_restricted, na.rm = TRUE),
  Mean_sigma2_unrestricted = mean(model_comparison$sigma2_unrestricted, na.rm = TRUE)
)

print(model_selection_summary)

cat("\n=== DESCRIPTIVE STATISTICS FOR PARAMETERS ===\n")

# Print parameter statistics
print(all_stats, row.names = FALSE)

# 9. SAVE RESULTS

cat("\n=== SAVING RESULTS ===\n")

# Save all results
save(all_results, model_comparison, restricted_params, unrestricted_params, all_stats,
     file = "varmax_sse_analysis_results.RData")

# Write tables to CSV files
write.csv(model_comparison, "model_comparison_sse_table.csv", row.names = FALSE)
write.csv(restricted_params, "restricted_model_sse_parameters.csv", row.names = FALSE)
write.csv(unrestricted_params, "unrestricted_model_sse_parameters.csv", row.names = FALSE)
write.csv(all_stats, "parameter_descriptive_statistics_sse.csv", row.names = FALSE)

cat("Results saved to:\n")
cat("- varmax_sse_analysis_results.RData\n")
cat("- model_comparison_sse_table.csv\n")
cat("- restricted_model_sse_parameters.csv\n")
cat("- unrestricted_model_sse_parameters.csv\n")
cat("- parameter_descriptive_statistics_sse.csv\n")

# 10. CREATE VISUALIZATION OF MODEL COMPARISON

cat("\n=== CREATING VISUALIZATIONS ===\n")

# Create AIC/BIC comparison plot
if (nrow(model_comparison) > 0) {
  # AIC comparison plot
  p_aic <- ggplot(model_comparison, aes(x = reorder(Participant, Delta_AIC), y = Delta_AIC)) +
    geom_bar(stat = "identity", aes(fill = Delta_AIC > 0), width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = c(-2, 2), linetype = "dotted", color = "gray50", alpha = 0.5) +
    scale_fill_manual(values = c("#4E79A7", "#E15759"), 
                      labels = c("Restricted better", "Unrestricted better"),
                      name = "Model Preference") +
    labs(
      x = "Participant",
      y = "ΔAIC (Unrestricted - Restricted)",
      title = "Model Comparison by ΔAIC (SSE-based)",
      subtitle = paste("Negative values favor Restricted model, positive values favor Unrestricted",
                       "\nDotted lines at ΔAIC = ±2")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  # SSE comparison plot
  p_sse <- ggplot(model_comparison) +
    geom_density(aes(x = SSE_restricted, fill = "Restricted"), alpha = 0.5) +
    geom_density(aes(x = SSE_unrestricted, fill = "Unrestricted"), alpha = 0.5) +
    scale_fill_manual(values = c("#4E79A7", "#E15759")) +
    labs(
      x = "Sum of Squared Errors (SSE)",
      y = "Density",
      title = "Distribution of SSE Values",
      subtitle = "Lower SSE indicates better model fit",
      fill = "Model"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 11),
      legend.position = "bottom"
    )
  
  # Parameter distributions for key parameters
  if (nrow(restricted_params) > 0 && nrow(unrestricted_params) > 0) {
    # Combine Theta parameter data
    restricted_long <- restricted_params %>%
      select(Participant, Theta1, Theta2) %>%
      mutate(Model = "Restricted") %>%
      pivot_longer(cols = c(Theta1, Theta2), names_to = "Parameter", values_to = "Value")
    
    unrestricted_long <- unrestricted_params %>%
      select(Participant, Theta1, Theta2, Phi1, Phi2) %>%
      mutate(Model = "Unrestricted") %>%
      pivot_longer(cols = c(Theta1, Theta2, Phi1, Phi2), names_to = "Parameter", values_to = "Value")
    
    combined_params <- bind_rows(restricted_long, unrestricted_long)
    
    # Plot parameter distributions
    p_theta <- ggplot(combined_params, aes(x = Value, fill = Model)) +
      geom_density(alpha = 0.5) +
      facet_wrap(~ Parameter, ncol = 2, scales = "free") +
      scale_fill_manual(values = c("#4E79A7", "#E15759")) +
      labs(
        x = "Parameter Value",
        y = "Density",
        title = "Distribution of Theta and Phi Parameters",
        subtitle = "Comparison between Restricted and Unrestricted models"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 11),
        legend.position = "bottom"
      )
    
    # Save plots
    ggsave("Delta_AIC_comparison_sse.png", p_aic, width = 10, height = 6, dpi = 300)
    ggsave("SSE_distribution_comparison.png", p_sse, width = 10, height = 6, dpi = 300)
    ggsave("Parameter_distributions_sse.png", p_theta, width = 10, height = 8, dpi = 300)
    
    cat("Visualizations saved:\n")
    cat("- Delta_AIC_comparison_sse.png\n")
    cat("- SSE_distribution_comparison.png\n")
    cat("- Parameter_distributions_sse.png\n")
  }
  
  # Create a summary plot of model performance
  summary_plot_data <- model_comparison %>%
    summarise(
      AIC_Restricted_better = sum(Preferred_AIC == "Restricted"),
      AIC_Unrestricted_better = sum(Preferred_AIC == "Unrestricted"),
      BIC_Restricted_better = sum(Preferred_BIC == "Restricted"),
      BIC_Unrestricted_better = sum(Preferred_BIC == "Unrestricted")
    ) %>%
    pivot_longer(cols = everything(), names_to = "Criterion", values_to = "Count") %>%
    separate(Criterion, into = c("Criterion", "Model", "better"), sep = "_") %>%
    mutate(Criterion = paste0(Criterion, " better"))
  
  p_summary <- ggplot(summary_plot_data, aes(x = Criterion, y = Count, fill = Model)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    scale_fill_manual(values = c("#4E79A7", "#E15759")) +
    geom_text(aes(label = Count), position = position_dodge(width = 0.7), 
              vjust = -0.5, size = 4, fontface = "bold") +
    labs(
      x = "Information Criterion",
      y = "Number of Participants",
      title = "Model Preference Summary",
      subtitle = sprintf("Total participants: %d", nrow(model_comparison)),
      fill = "Preferred Model"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = 11),
      legend.position = "bottom"
    )
  
  ggsave("model_preference_summary.png", p_summary, width = 8, height = 6, dpi = 300)
  cat("- model_preference_summary.png\n")
}

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("\nKey differences from previous version:\n")
cat("1. Uses Sum of Squared Errors (SSE) instead of negative log-likelihood\n")
cat("2. AIC/BIC calculated using SSE-based formulas: AIC = n*ln(SSE/n) + 2k\n")
cat("3. Sigma parameters are estimated from SSE: sigma² = SSE/(n-p)\n")
cat("4. Only meaningful parameters (delta, Psi, Theta, Phi) are stored\n")
cat("5. Same parameter bounds as supervisor's specification\n")