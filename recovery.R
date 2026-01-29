################################################################################
# RECOVERY
#   - Show that we can estimate back the parameters of each of the models
#   - Show that we can simulate VARMAX data and estimate the exponential 
#     discounting model and vice-versa
################################################################################

# Load in the objective functions of the different models, namely `exponential`, 
# `quasi_hyperbolic`, `double_exponential`, and `varmax`
source(file.path( "~/Documents/Invatamant/Uni/Year_4/Thesis/Data/Niels code/objective_functions.R"))

# Load in the simulation functions for the models
source(file.path("~/Documents/Invatamant/Uni/Year_4/Thesis/Data/Niels code/simulate.R"))

# Combine simulation functions, objective functions, and the bounds of the 
# parameters of the models to each other
models <- list(
    "varmax" = list(
        \(parameters, X) simulate_varmax_ls(parameters, X, restricted = FALSE),
        \(parameters, y, X, residuals = FALSE) varmax_ls(parameters, y, X, restricted = FALSE, residuals = residuals), 
        cbind(
            c(0, 0, -5, -5, -(1 - 1e-5), -(1 - 1e-5), -(1 - 1e-5), -(1 - 1e-5), 1e-5, -1 + 1e-5, 1e-5),
            c(1, 1, 5, 5, 1 - 1e-5, 1 - 1e-5, 1 - 1e-5, 1 - 1e-5, 0.25, 1 - 1e-5, 0.25)
        ),
        function(x) x
    ),
    "varmax_restricted" = list(
        \(parameters, X) simulate_varmax_ls(parameters, X, restricted = TRUE),
        \(parameters, y, X, residuals = FALSE) varmax_ls(parameters, y, X, restricted = TRUE, residuals = residuals), 
        cbind(
            c(0, 0, -5, -5, -(1 - 1e-5), -(1 - 1e-5), 1e-5, -1 + 1e-5, 1e-5),
            c(1, 1, 5, 5, 1 - 1e-5, 1 - 1e-5, 0.25, 1 - 1e-5, 0.25)
        ),
        function(x) x
    )
)

# Perform the actual recovery, save the results, and create a figure displaying
# them. Some important choices:
#   - Simulating 100 trajectories per model
#   - Each trajectory consists of 152 datapoints, just like the VANHASBROECK_2022
#     data
#   - Stimuli are generated according to a uniform distribution with ranges 
#     similar to those used in the paper
#
# Loop over models
set.seed(1)
N <- 100
for(i in seq_along(models)) {
    cat("\n")

    # Simulate the parameters that will be used to simulate this model
    sim_params <- runif(
        nrow(models[[i]][[3]]) * N,
        min = models[[i]][[3]][, 1],
        max = models[[i]][[3]][, 2]
    ) |>
        matrix(
            nrow = N,
            ncol = nrow(models[[i]][[3]]),
            byrow = TRUE
        )

    # Split up the dataset according to row
    splitted <- split(
        as.data.frame(sim_params),
        sort(seq_len(nrow(sim_params)))
    )

    # Loop over each simulation
    est_params <- parallel::mclapply(
        seq_along(splitted),
        function(j) {
            # Give us an update of where they are
            cat(paste0("", names(models)[i], ": ", j, "\n"))

            # Simulate the values of the stimuli. First do a general draw, then 
            # scramble them
            x <- c(
                runif(76, 0.5, 4),
                -runif(76, 0.5, 4)
            )
            x <- rep(
                sample(x, 152, replace = FALSE),
                times = 10
            )

            # Creat the matrix necessary to run the models
            X <- stimulus_matrix(x)

            # Simulate the model itself
            y <- splitted[[j]] |>
                as.numeric() |>
                models[[i]][[1]](X)

            # Estimate the model through DEoptim
            n <- nrow(models[[i]][[3]])
            result <- DEoptim::DEoptim(
                \(x) models[[i]][[2]](x, y, X),
                lower = models[[i]][[3]][1:(n - 3), 1],
                upper = models[[i]][[3]][1:(n - 3), 2],
                control = DEoptim::DEoptim.control(
                    NP = 200,
                    CR = 0.5,
                    F = 0.8,
                    itermax = 1000, 
                    strategy = 6,
                    p = 0.5,
                    trace = FALSE
                )
            )

            # Get the residuals and compute the covariances: Only needed when 
            # doing least-squares estimation, which is only true for the VARMAX 
            # (for now, at least)
            if(grepl("varmax", names(models)[i], fixed = TRUE)) {
                params <- models[[i]][[4]](result$optim$bestmem)
                residuals <- models[[i]][[2]](
                    params, 
                    y, 
                    X,
                    residuals = TRUE
                )
                G <- residuals |>
                    cov() |>
                    chol() |>
                    t()
                G <- G[lower.tri(G, diag = TRUE)] |>
                    as.numeric()
                
                return(c(params, G))

            } else {
                return(models[[i]][[4]](result$optim$bestmem))
            }
        },
        mc.cores = ifelse(
            Sys.info()["sysname"] == "Windows",
            1,
            5
        )
    )
    est_params <- do.call("rbind", est_params)

    # Create a dataframe containing the results
    results <- data.frame(
        iteration = rep(1:nrow(sim_params), each = ncol(sim_params)),
        parameters = rep(1:ncol(sim_params), times = nrow(sim_params)),
        simulated = as.numeric(t(sim_params)),
        estimated = as.numeric(t(est_params))
    )

    
    # Save the results
    output_file <- file.path("results", paste0(names(models)[i], ".csv"))
    cat(paste0("Saving results to: ", output_file, "\n"))
    
    write.csv(
      results,
      file = output_file,
      row.names = FALSE
    )
    
    # Create a plot displaying the results
    plt <- lapply(
        1:ncol(est_params), 
        function(j) {
            plt_data <- data.frame(
                x = sim_params[, j],
                y = est_params[, j]
            )

            return(
                ggplot2::ggplot(
                    data = plt_data, 
                    ggplot2::aes(
                        x = x, 
                        y = y
                    )
                ) +
                    ggplot2::geom_abline(
                        intercept = 0, 
                        slope = 1, 
                        linewidth = 2, 
                        color = "black"
                    ) +
                    ggplot2::geom_point(
                        shape = 21,
                        alpha = 0.25, 
                        color = "black",
                        fill = "cornflowerblue",
                        size = 5
                    ) +
                    ggplot2::annotate(
                        "text",
                        x = models[[i]][[3]][j, 1] + 0.05 * diff(models[[i]][[3]][j, ]),
                        y = models[[i]][[3]][j, 1] + 0.95 * diff(models[[i]][[3]][j, ]),
                        label = format(
                            round(
                                cor(
                                    plt_data$x, 
                                    plt_data$y
                                ),
                                digits = 2
                            ),
                            nsmall = 2
                        ),
                        size = 5,
                        hjust = 0
                    ) +
                    ggplot2::lims(
                        x = models[[i]][[3]][j, ],
                        y = models[[i]][[3]][j, ]
                    ) +
                    ggplot2::labs(
                        title = paste("Parameter", j),
                        x = "Simulated",
                        y = "Estimated"
                    ) +
                    ggplot2::theme(
                        panel.background = ggplot2::element_rect(
                            fill = "white"
                        ),
                        panel.border = ggplot2::element_rect(
                            fill = NA, 
                            color = "black",
                            linewidth = 1
                        ),
                        axis.title = ggplot2::element_text(size = 15),
                        plot.title = ggplot2::element_text(size = 20)
                    )
            )
        }
    )

    ggplot2::ggsave(
        filename = figure_file,,
        ggpubr::ggarrange(
            plotlist = plt, 
            nrow = 1
        ),
        width = length(plt) * 1500, 
        height = 1650,
        unit = "px",
        limitsize = FALSE
    )
}

# VARMAX - exponential seems to well estimated after transformation of the 
# mean/intercept estimated by each model. 