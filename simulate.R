#' Objective function for the VARMAX
#' 
#' @param parameters Numeric vector containing the parameters of the model
#' @param X Numeric matrix denoting the stimuli presented in the experiment
#' @param restricted Logical denoting whether to fit a restricted version of the
#' VARMAX (i.e., one assuming white noise residuals)
#' 
#' @return Numeric matrix containing the simulated observations according to 
#' the model
#' 
#' @export 
simulate_varmax_ls <- function(parameters, 
                               X, 
                               restricted = FALSE) {

    # Adjust the matrix X so that you only use the first column. This should 
    # provide the input needed to evaluate y. [Can be deleted later, providing 
    # a vector of stimuli instead of a matrix]
    x <- X[, 1]

    # Differentiate between a restricted version of the VARMAX or the non-restricted
    # version thereof and extract parameters to give them substantive meaning 
    # (as specified in the paper). Specifically:
    #   - \delta: The intercept of the process
    #   - \psi: The slope of the stimuli
    #   - \theta: The autoregressive effect
    #   - \phi: The slope of the moving average
    #   - \sigma: The residual standard deviation
    delta <- parameters[1:2]
    Psi <- parameters[3:4] |>
        matrix(ncol = 1)

    Theta <- parameters[5:6] |>
        diag()

    if(restricted) {
        Phi <- -Theta
        G <- matrix(0, nrow = 2, ncol = 2)
        G[lower.tri(G, diag = TRUE)] <- parameters[7:9]

    } else {
        Phi <- parameters[7:8] |>
            diag()

        G <- matrix(0, nrow = 2, ncol = 2)
        G[lower.tri(G, diag = TRUE)] <- parameters[9:11]
    }
    Sigma <- G %*% t(G)

    # Compute the expected values of the dependent variable y for each trial. 
    # These will serve as the means of the conditional normal process for the 
    # likelihood function.
    epsilon <- MASS::mvrnorm(
        nrow(X), 
        mu = c(0, 0),
        Sigma = Sigma
    )
    y_hat <- matrix(
        0,
        nrow = nrow(X),
        ncol = 2
    )
    for(i in 1:nrow(X)) {
        # Differentiate between the first lag and the rest. Needed because
        # there is no lag in y nor in the residual yet
        if(i == 1) {
            y_hat[i, ] <- solve(diag(2) - Theta) %*% delta + Psi %*% x[i] + epsilon[i, ]

        } else {
            y_hat[i, ] <- delta + Psi %*% x[i] + Theta %*% y_hat[i - 1, ] + Phi %*% epsilon[i - 1, ] + epsilon[i, ]
        }
    }

    return(y_hat)
}
