#' Create a stimulus matrix based on the input \code{x}
#' 
#' @param x Numeric vector containing the stimuli presented in the experiment
#' 
#' @return Numeric matrix containing the history of all stimuli in the rows, where
#' the first column indicates the most recent stimulus and the last column the 
#' most distant one
#' 
#' @export
stimulus_matrix <- function(x) {
    X <- matrix(
        0, 
        nrow = length(x),
        ncol = length(x)
    )

    for(i in 1:length(x)) {
        X[i, 1:i] <- rev(x[1:i])
    }

    return(X)
}

#' Objective function for the VARMAX, least-squares
#' 
#' @param parameters Numeric vector containing the parameters of the model
#' @param y Numeric matrix of size N x d denoting the dependent variables
#' @param X Numeric matrix denoting the stimuli presented in the experiment, 
#' containing only rows that match the observations in y
#' @param restricted Logical denoting whether to fit a restricted version of the
#' VARMAX (i.e., one assuming white noise residuals)
#' 
#' @return Min-log-likelihood for the observations under the model specified by
#' \code{parameters}
#' 
#' @export 
varmax_ls <- function(parameters, 
                      y, 
                      X, 
                      restricted = FALSE,
                      residuals = FALSE) {

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

    } else {
        Phi <- parameters[7:8] |>
            diag()
    }

    # Compute the expected values of the dependent variable y for each trial. 
    # These will serve as the means of the conditional normal process for the 
    # likelihood function.
    epsilon <- matrix(
        0,
        nrow = nrow(y),
        ncol = 2
    )
    for(i in 1:nrow(y)) {
        # Differentiate between the first lag and the rest. Needed because
        # there is no lag in y nor in the residual yet
        if(i == 1) {
            y_hat <- solve(diag(2) - Theta) %*% delta + Psi %*% x[i]
            epsilon[i, ] <- y[i, ] - y_hat

        } else {
            y_hat <- delta + Psi %*% x[i] + Theta %*% y[i - 1, ] + Phi %*% epsilon[i - 1, ]
            epsilon[i, ] <- y[i, ] - y_hat            
        }
    }

    # Compute the sum of squared error and return this value (residuals = FALSE)
    # or return the residuals themselves (residuals = TRUE)
    if(residuals) {
        return(epsilon)
    } else {
        return(sum(epsilon^2))
    }
}
