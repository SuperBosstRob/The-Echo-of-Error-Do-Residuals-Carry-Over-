#' Simulate data for the exponential discounting model
#' 
#' @param parameters Numeric vector of parameters for the exponential discounting
#' model.
#' @param X Numeric matrix containing the stimuli presented to the participants
#' 
#' @return Numeric vector containing the simulated observations according to 
#' the model
#' 
#' @export 
simulate_exponential <- function(parameters, 
                                 x) {
    
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
    decay <- gamma^seq(0, nrow(X) - 1, 1)

    # Compute the expected values of the dependent variable y for each trial. 
    # These will serve as the means of the conditional normal process for the 
    # likelihood function.
    #
    # Note that I do this in a vectorized way. Could've easily been done with a 
    # for-loop as well. %*% denotes the matrix multiplication of the stimulus 
    # matrix X with the decay vector, resulting in a vector of a weighted sum 
    # of the stimuli as weighted by the discounting factor.
    return(
        rnorm(
            nrow(X),
            mean = alpha + beta * X %*% decay,
            sd = sigma
        )
    )
}





#' Objective function for the quasi-hyperbolic discounting function
#' 
#' @param parameters Numeric vector containing the parameters of the model
#' @param X Numeric matrix denoting the stimuli presented in the experiment
#' 
#' @return Numeric vector containing the simulated observations according to 
#' the model
#' 
#' @export 
simulate_quasi_hyperbolic <- function(parameters, 
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
    decay <- nu^seq(0, nrow(X) - 1, 1)
    recency <- kappa^c(0, rep(1, nrow(X) - 1))

    decay <- decay * recency

    # Compute the expected values of the dependent variable y for each trial. 
    # These will serve as the means of the conditional normal process for the 
    # likelihood function.
    #
    # Note that I do this in a vectorized way. Could've easily been done with a 
    # for-loop as well. %*% denotes the matrix multiplication of the stimulus 
    # matrix X with the decay vector, resulting in a vector of a weighted sum 
    # of the stimuli as weighted by the discounting factor.
    return(
        rnorm(
            nrow(X),
            mean = alpha + beta * X %*% decay,
            sd = sigma
        )
    )
}





#' Objective function for the double exponential discounting function
#' 
#' @param parameters Numeric vector containing the parameters of the model
#' @param X Numeric matrix denoting the stimuli presented in the experiment
#' 
#' @return Numeric vector containing the simulated observations according to 
#' the model
#' 
#' @export 
simulate_double_exponential <- function(parameters, 
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
    decay_1 <- gamma^seq(0, nrow(X) - 1, 1)
    decay_2 <- nu^seq(0, nrow(X) - 1, 1)

    # Compute the expected values of the dependent variable y for each trial. 
    # These will serve as the means of the conditional normal process for the 
    # likelihood function.
    #
    # Note that I do this in a vectorized way. Could've easily been done with a 
    # for-loop as well. %*% denotes the matrix multiplication of the stimulus 
    # matrix X with the decay vector, resulting in a vector of a weighted sum 
    # of the stimuli as weighted by the discounting factor.
    return(
        rnorm(
            nrow(X),
            mean = alpha + omega * beta * X %*% decay_1 + (1 - omega) * beta * X %*% decay_2,
            sd = sigma
        )
    )
}







#' Objective function for the ARMAX
#' 
#' @param parameters Numeric vector containing the parameters of the model
#' @param X Numeric matrix denoting the stimuli presented in the experiment
#' @param restricted Logical denoting whether to fit a restricted version of the
#' VARMAX (i.e., one assuming white noise residuals)
#' 
#' @return Numeric vector containing the simulated observations according to 
#' the model
#' 
#' @export 
simulate_armax <- function(parameters, 
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
    epsilon <- rnorm(nrow(X), mean = 0, sd = sigma)
    y_hat <- numeric(nrow(X))
    for(i in 1:nrow(X)) {
        # Differentiate between the first lag and the rest. Needed because
        # there is no lag in y nor in the residual yet
        if(i == 1) {
            y_hat[i] <- delta + psi * x[i] + epsilon[i]

        } else {
            y_hat[i] <- delta + psi * x[i] + theta * y_hat[i - 1] + phi * epsilon[i - 1] + epsilon[i]
        }
    }

    return(y_hat)
}







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
simulate_varmax <- function(parameters, 
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
            y_hat[i, ] <- delta + Psi %*% x[i] + epsilon[i, ]

        } else {
            y_hat[i, ] <- delta + Psi %*% x[i] + Theta %*% y_hat[i - 1, ] + Phi %*% epsilon[i - 1, ] + epsilon[i, ]
        }
    }

    return(y_hat)
}









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
