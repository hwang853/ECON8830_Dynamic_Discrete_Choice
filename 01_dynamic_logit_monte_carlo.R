# ECON 8830: Dynamic discrete-choice Monte Carlo experiment
#
# Combined R translation of the MATLAB code shown in PS1DGP.pdf and
# PS1LOGLIK.pdf. The PDFs reference helper functions and data objects that are
# not included in the exported pages. This script supplies a coherent dynamic
# binary-logit DGP and its matching conditional log-likelihood.

simulate_dgp <- function(n, theta, y0 = 0L) {
  stopifnot(n >= 2L, length(theta) == 3L, y0 %in% c(0L, 1L))

  x <- rnorm(n)
  y <- integer(n)
  y_lag <- c(y0, y[-n])

  for (t in seq_len(n)) {
    lag_choice <- if (t == 1L) y0 else y[t - 1L]
    index <- theta[1L] + theta[2L] * x[t] + theta[3L] * lag_choice
    y[t] <- rbinom(1L, size = 1L, prob = plogis(index))
  }

  y_lag <- c(y0, y[-n])
  data.frame(y = y, x = x, y_lag = y_lag)
}

negative_log_likelihood <- function(theta, y, x, y0 = 0L) {
  n <- length(y)
  if (length(x) != n || n < 2L) {
    stop("y and x must have the same length, with at least two observations")
  }

  y_lag <- c(y0, y[-n])
  index <- theta[1L] + theta[2L] * x + theta[3L] * y_lag
  probability <- plogis(index)

  -sum(dbinom(y, size = 1L, prob = probability, log = TRUE))
}

estimate_model <- function(data, theta_start = c(0, 1, -0.5), y0 = 0L) {
  fit <- optim(
    par = theta_start,
    fn = negative_log_likelihood,
    y = data$y,
    x = data$x,
    y0 = y0,
    method = "BFGS",
    hessian = TRUE,
    control = list(reltol = 1e-10, maxit = 1000L)
  )

  inverse_hessian <- tryCatch(
    solve(fit$hessian),
    error = function(e) matrix(NA_real_, nrow = 3L, ncol = 3L)
  )

  list(
    coefficients = fit$par,
    objective = fit$value,
    convergence = fit$convergence,
    inverse_hessian = inverse_hessian
  )
}

run_monte_carlo <- function(
    n = 200L,
    replications = 101L,
    theta_true = c(0, 1, -0.5),
    theta_start = theta_true,
    seed = 8830L) {
  stopifnot(length(theta_true) == 3L, length(theta_start) == 3L)
  set.seed(seed)

  estimates <- matrix(
    NA_real_,
    nrow = replications,
    ncol = 3L,
    dimnames = list(NULL, c("intercept", "x", "lag_choice"))
  )
  objective <- numeric(replications)
  convergence <- integer(replications)
  inverse_hessians <- array(NA_real_, dim = c(3L, 3L, replications))

  for (j in seq_len(replications)) {
    simulated_data <- simulate_dgp(n = n, theta = theta_true)
    fit <- estimate_model(simulated_data, theta_start = theta_start)

    estimates[j, ] <- fit$coefficients
    objective[j] <- fit$objective
    convergence[j] <- fit$convergence
    inverse_hessians[, , j] <- fit$inverse_hessian
  }

  mean_estimate <- colMeans(estimates)
  bias <- mean_estimate - theta_true
  estimate_covariance <- cov(estimates)
  monte_carlo_se <- apply(estimates, 2L, sd) / sqrt(nrow(estimates))
  average_asymptotic_covariance <- apply(
    inverse_hessians,
    MARGIN = c(1L, 2L),
    FUN = function(z) mean(z, na.rm = TRUE)
  )

  list(
    settings = list(
      n = n,
      replications = replications,
      theta_true = theta_true,
      theta_start = theta_start,
      seed = seed
    ),
    estimates = estimates,
    objective = objective,
    convergence = convergence,
    mean_estimate = mean_estimate,
    bias = bias,
    estimate_covariance = estimate_covariance,
    monte_carlo_se = monte_carlo_se,
    average_asymptotic_covariance = average_asymptotic_covariance
  )
}

print_results <- function(results) {
  cat("Dynamic binary-logit Monte Carlo experiment\n")
  cat("Observations per replication:", results$settings$n, "\n")
  cat("Replications:", results$settings$replications, "\n")
  cat("Successful optimizations:", sum(results$convergence == 0L), "\n\n")

  summary_table <- cbind(
    truth = results$settings$theta_true,
    mean = results$mean_estimate,
    bias = results$bias,
    monte_carlo_se = results$monte_carlo_se
  )
  print(round(summary_table, 6L))

  cat("\nCovariance of estimates:\n")
  print(round(results$estimate_covariance, 6L))

  cat("\nAverage inverse Hessian:\n")
  print(round(results$average_asymptotic_covariance, 6L))
}

if (sys.nframe() == 0L) {
  results <- run_monte_carlo()
  print_results(results)

  write.csv(
    results$estimates,
    file = "monte_carlo_estimates.csv",
    row.names = FALSE
  )
}
