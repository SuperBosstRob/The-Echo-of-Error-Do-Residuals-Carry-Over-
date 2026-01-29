# TO DO
#   - For illustration purposes, I include the full computation in the objective
#     functions. However, a lot can be replaced with the simulation functions 
#     to get the y_hat (note however that the simulation functions draw from 
#     a normal distribution: This would need to change to deterministic choice
#     if I want to allow for this approach)

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





#' Objective function for the exponential discounting function
#' 
#' @param parameters Numeric vector containing the parameters of the model
#' @param y Numeric vector denoting the dependent variable
#' @param X Numeric matrix denoting the stimuli presented in the experiment, 
#' containing only rows that match the observations in y
#' 
#' @return Min-log-likelihood for the observations under the model specified by
#' \code{parameters}
#' 
#' @export 
exponential <- function(parameters, 
                        y, 
                        X) {
    
    # Extract parameters to give them substantive meaning (as specified in the
    # paper). Specifically:
    #   - \alpha: The mean of the process
    #   - \beta: The slope of the stimuli
    #   - \gamma: The discounting factor
    #   - \sigma: The residual standard deviation
    alpha <- parameters[1]
    beta <- parameters[2]
    gamma <- parameters[3]
    sigma <- parameters[4]

    # Create a decay vector that specifies the value of the discounting at each 
    # trial
    decay <- gamma^seq(0, ncol(X) - 1, 1)

    # Compute the expected values of the dependent variable y for each trial. 
    # These will serve as the means of the conditional normal process for the 
    # likelihood function.
    #
    # Note that I do this in a vectorized way. Could've easily been done with a 
    # for-loop as well. %*% denotes the matrix multiplication of the stimulus 
    # matrix X with the decay vector, resulting in a vector of a weighted sum 
    # of the stimuli as weighted by the discounting factor.
    y_hat <- alpha + beta * X %*% decay

    # Compute the likelihood assuming normality
    L <- dnorm(
        y, 
        mean = y_hat, 
        sd = sigma
    )

    # Log-transform, summate, and reverse the sign to get the min-log-likelihood.
    # Return this value.
    return(-sum(log(L)))
}





#' Objective function for the quasi-hyperbolic discounting function
#' 
#' @param parameters Numeric vector containing the parameters of the model
#' @param y Numeric vector denoting the dependent variable
#' @param X Numeric matrix denoting the stimuli presented in the experiment, 
#' containing only rows that match the observations in y
#' 
#' @return Min-log-likelihood for the observations under the model specified by
#' \code{parameters}
#' 
#' @export 
quasi_hyperbolic <- function(parameters, 
                             y, 
                             X) {
    
    # Extract parameters to give them substantive meaning (as specified in the
    # paper). Specifically:
    #   - \alpha: The mean of the process
    #   - \beta: The slope of the stimuli
    #   - \nu: The discounting factor after recency has taken effect
    #   - \kappa: The recency effect of the stimuli
    #   - \sigma: The residual standard deviation
    alpha <- parameters[1]
    beta <- parameters[2]
    nu <- parameters[3]
    kappa <- parameters[4]
    sigma <- parameters[5]

    # Create a decay vector that specifies the value of the discounting at each 
    # trial
    decay <- nu^seq(0, ncol(X) - 1, 1)
    recency <- kappa^c(0, rep(1, ncol(X) - 1))

    decay <- decay * recency

    # Compute the expected values of the dependent variable y for each trial. 
    # These will serve as the means of the conditional normal process for the 
    # likelihood function.
    #
    # Note that I do this in a vectorized way. Could've easily been done with a 
    # for-loop as well. %*% denotes the matrix multiplication of the stimulus 
    # matrix X with the decay vector, resulting in a vector of a weighted sum 
    # of the stimuli as weighted by the discounting factor.
    y_hat <- alpha + beta * X %*% decay

    # Compute the likelihood assuming normality
    L <- dnorm(
        y, 
        mean = y_hat, 
        sd = sigma
    )

    # Log-transform, summate, and reverse the sign to get the min-log-likelihood.
    # Return this value.
    return(-sum(log(L)))
}





#' Objective function for the double exponential function
#' 
#' @param parameters Numeric vector containing the parameters of the model
#' @param y Numeric vector denoting the dependent variable
#' @param X Numeric matrix denoting the stimuli presented in the experiment, 
#' containing only rows that match the observations in y
#' 
#' @return Min-log-likelihood for the observations under the model specified by
#' \code{parameters}
#' 
#' @export 
double_exponential <- function(parameters, 
                               y, 
                               X) {
    
    # Extract parameters to give them substantive meaning (as specified in the
    # paper). Specifically:
    #   - \alpha: The mean of the process
    #   - \beta: The slope of the stimuli
    #   - \gamma: The first discounting factor
    #   - \nu: The second discounting factor
    #   - \omega: The weighting between the two types of discounting
    #   - \sigma: The residual standard deviation
    alpha <- parameters[1]
    beta <- parameters[2]
    gamma <- parameters[3]
    nu <- parameters[4]
    omega <- parameters[5]
    sigma <- parameters[6]

    # Create a decay vector that specifies the value of the discounting at each 
    # trial
    decay_1 <- gamma^seq(0, ncol(X) - 1, 1)
    decay_2 <- nu^seq(0, ncol(X) - 1, 1)

    # Compute the expected values of the dependent variable y for each trial. 
    # These will serve as the means of the conditional normal process for the 
    # likelihood function.
    #
    # Note that I do this in a vectorized way. Could've easily been done with a 
    # for-loop as well. %*% denotes the matrix multiplication of the stimulus 
    # matrix X with the decay vector, resulting in a vector of a weighted sum 
    # of the stimuli as weighted by the discounting factor.
    y_hat <- alpha + omega * beta * X %*% decay_1 + (1 - omega) * beta * X %*% decay_2

    # Compute the likelihood assuming normality
    L <- dnorm(
        y, 
        mean = y_hat, 
        sd = sigma
    )

    # Log-transform, summate, and reverse the sign to get the min-log-likelihood.
    # Return this value.
    return(-sum(log(L)))
}







#' Objective function for the ARMAX
#' 
#' @param parameters Numeric vector containing the parameters of the model
#' @param y Numeric vector denoting the dependent variable
#' @param X Numeric matrix denoting the stimuli presented in the experiment, 
#' containing only rows that match the observations in y
#' @param restricted Logical denoting whether to fit a restricted version of the
#' VARMAX (i.e., one assuming white noise residuals)
#' 
#' @return Min-log-likelihood for the observations under the model specified by
#' \code{parameters}
#' 
#' @export 
armax <- function(parameters, 
                  y, 
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
    if(restricted) {
        delta <- parameters[1]
        psi <- parameters[2]
        theta <- parameters[3]
        phi <- -theta
        sigma <- parameters[4]

    } else {
        delta <- parameters[1]
        psi <- parameters[2]
        theta <- parameters[3]
        phi <- parameters[4]
        sigma <- parameters[5]
    }

    # Compute the expected values of the dependent variable y for each trial. 
    # These will serve as the means of the conditional normal process for the 
    # likelihood function.
    epsilon <- numeric(length(y))
    for(i in 1:length(y)) {
        # Differentiate between the first lag and the rest. Needed because
        # there is no lag in y nor in the residual yet
        if(i == 1) {
            y_hat <- delta + psi * x[i]
            epsilon[i] <- y[i] - y_hat

        } else {
            y_hat <- delta + psi * x[i] + theta * y[i - 1] + phi * epsilon[i - 1]
            epsilon[i] <- y[i] - y_hat            
        }
    }

    # Compute the likelihood assuming normality
    L <- dnorm(
        epsilon, 
        mean = 0, 
        sd = sigma
    )

    # Log-transform, summate, and reverse the sign to get the min-log-likelihood.
    # Return this value.
    return(-sum(log(L)))
}







#' Objective function for the VARMAX
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
varmax <- function(parameters, 
                   y, 
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

    Lambda <- parameters[5:6] |>
        diag()
    # angles <- cumsum(parameters[7:8])
    # Q <- cbind(
    #     c(cos(angles[1]), sin(angles[1])),
    #     c(-sin(angles[2]), cos(angles[2]))
    # )
    # Theta <- Q %*% Lambda %*% solve(Q)
    Theta <- Lambda

    if(restricted) {
        Phi <- -Theta
        G <- matrix(0, nrow = 2, ncol = 2)
        G[lower.tri(G, diag = TRUE)] <- parameters[9:11]

    } else {
        Lambda <- parameters[9:10] |>
            diag()
        # angles <- cumsum(parameters[11:12])
        # Q <- cbind(
        #     c(cos(angles[1]), sin(angles[1])),
        #     c(-sin(angles[2]), cos(angles[2]))
        # )
        # Phi <- Q %*% Lambda %*% solve(Q)
        Phi <- Lambda

        G <- matrix(0, nrow = 2, ncol = 2)
        G[lower.tri(G, diag = TRUE)] <- parameters[13:15]
    }
    Sigma <- G %*% t(G)

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
            y_hat <- delta + Psi %*% x[i]
            epsilon[i, ] <- y[i, ] - y_hat

        } else {
            y_hat <- delta + Psi %*% x[i] + Theta %*% y[i - 1, ] + Phi %*% epsilon[i - 1, ]
            epsilon[i, ] <- y[i, ] - y_hat            
        }
    }

    # Compute the likelihood assuming normality
    L <- mvtnorm::dmvnorm(
        epsilon, 
        mean = c(0, 0), 
        sigma = Sigma
    )

    # Log-transform, summate, and reverse the sign to get the min-log-likelihood.
    # Return this value.
    return(-sum(log(L)))
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
