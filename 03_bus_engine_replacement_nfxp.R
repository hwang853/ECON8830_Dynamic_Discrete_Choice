# ECON 8830, Part 3: Rust-style bus engine replacement model
#
# R translation of WHRECON8830PS3.pdf. The original code estimates a dynamic
# discrete-choice model with 90 mileage states, replacement decisions, and a
# three-point mileage transition distribution using a nested fixed-point (NFXP)
# algorithm. The source a530875 data are not included in the PDF, so this script
# accepts raw or prepared CSV data and includes a reproducible simulated panel.

log_sum_exp_two <- function(a, b) {
  maximum <- pmax(a, b)
  maximum + log(exp(a - maximum) + exp(b - maximum))
}

transition_matrix <- function(probabilities, n_states = 90L) {
  stopifnot(length(probabilities) == 3L, n_states >= 3L)
  if (any(probabilities < 0) || abs(sum(probabilities) - 1) > 1e-8) {
    stop("transition probabilities must be nonnegative and sum to one")
  }

  matrix_p <- matrix(0, nrow = n_states, ncol = n_states)
  for (state in seq_len(n_states)) {
    for (increment in 0:2) {
      destination <- min(state + increment, n_states)
      matrix_p[state, destination] <-
        matrix_p[state, destination] + probabilities[increment + 1L]
    }
  }
  matrix_p
}

solve_expected_value <- function(
    maintenance_cost,
    replacement_cost,
    probabilities,
    discount = 0.99,
    n_states = 90L,
    tolerance = 1e-10,
    max_iterations = 20000L) {
  matrix_p <- transition_matrix(probabilities, n_states)
  state <- 0:(n_states - 1L)
  expected_value <- numeric(n_states)

  for (iteration in seq_len(max_iterations)) {
    value_keep <- -maintenance_cost * state + discount * expected_value
    value_replace <- -replacement_cost + discount * expected_value[1L]
    integrated_value <- log_sum_exp_two(value_keep, value_replace)
    updated <- as.vector(matrix_p %*% integrated_value)

    if (max(abs(updated - expected_value)) < tolerance) {
      return(list(
        expected_value = updated,
        transition_matrix = matrix_p,
        iterations = iteration,
        converged = TRUE
      ))
    }
    expected_value <- updated
  }

  list(
    expected_value = expected_value,
    transition_matrix = matrix_p,
    iterations = max_iterations,
    converged = FALSE
  )
}

unpack_parameters <- function(parameters) {
  maintenance_cost <- exp(parameters[1L])
  replacement_cost <- exp(parameters[2L])

  probability_scale <- c(parameters[3L], parameters[4L], 0)
  probability_scale <- probability_scale - max(probability_scale)
  probabilities <- exp(probability_scale) / sum(exp(probability_scale))

  list(
    maintenance_cost = maintenance_cost,
    replacement_cost = replacement_cost,
    probabilities = probabilities
  )
}

pack_parameters <- function(
    maintenance_cost,
    replacement_cost,
    probabilities) {
  c(
    log(maintenance_cost),
    log(replacement_cost),
    log(probabilities[1L] / probabilities[3L]),
    log(probabilities[2L] / probabilities[3L])
  )
}

nfxp_negative_log_likelihood <- function(
    parameters,
    panel,
    discount = 0.99,
    n_states = 90L) {
  theta <- unpack_parameters(parameters)
  fixed_point <- solve_expected_value(
    maintenance_cost = theta$maintenance_cost,
    replacement_cost = theta$replacement_cost,
    probabilities = theta$probabilities,
    discount = discount,
    n_states = n_states
  )
  if (!fixed_point$converged) return(1e12)

  state_index <- pmin(pmax(as.integer(panel$state), 0L), n_states - 1L) + 1L
  value_keep <- -theta$maintenance_cost * (state_index - 1L) +
    discount * fixed_point$expected_value[state_index]
  value_replace <- -theta$replacement_cost +
    discount * fixed_point$expected_value[1L]
  probability_replace <- plogis(value_replace - value_keep)
  probability_replace <- pmin(pmax(probability_replace, 1e-12), 1 - 1e-12)

  choice_log_likelihood <- sum(dbinom(
    as.integer(panel$replace),
    size = 1L,
    prob = probability_replace,
    log = TRUE
  ))

  valid_transition <- !is.na(panel$increment)
  increment <- pmin(pmax(as.integer(panel$increment[valid_transition]), 0L), 2L)
  transition_log_likelihood <- sum(log(theta$probabilities[increment + 1L]))

  -(choice_log_likelihood + transition_log_likelihood)
}

estimate_bus_replacement_model <- function(
    panel,
    discount = 0.99,
    n_states = 90L,
    start = list(
      maintenance_cost = 0.05,
      replacement_cost = 8,
      probabilities = c(0.40, 0.55, 0.05))) {
  panel <- validate_bus_panel(panel, n_states)
  start_parameters <- pack_parameters(
    start$maintenance_cost,
    start$replacement_cost,
    start$probabilities
  )

  fit <- optim(
    par = start_parameters,
    fn = nfxp_negative_log_likelihood,
    panel = panel,
    discount = discount,
    n_states = n_states,
    method = "BFGS",
    hessian = TRUE,
    control = list(maxit = 300L, reltol = 1e-9)
  )

  estimate <- unpack_parameters(fit$par)
  covariance_transformed <- tryCatch(
    solve(fit$hessian),
    error = function(e) matrix(NA_real_, 4L, 4L)
  )

  p0 <- estimate$probabilities[1L]
  p1 <- estimate$probabilities[2L]
  jacobian <- matrix(0, nrow = 4L, ncol = 4L)
  jacobian[1L, 1L] <- estimate$maintenance_cost
  jacobian[2L, 2L] <- estimate$replacement_cost
  jacobian[3L, 3L] <- p0 * (1 - p0)
  jacobian[3L, 4L] <- -p0 * p1
  jacobian[4L, 3L] <- -p1 * p0
  jacobian[4L, 4L] <- p1 * (1 - p1)
  covariance <- jacobian %*% covariance_transformed %*% t(jacobian)

  fixed_point <- solve_expected_value(
    estimate$maintenance_cost,
    estimate$replacement_cost,
    estimate$probabilities,
    discount,
    n_states
  )

  list(
    coefficients = c(
      maintenance_cost = estimate$maintenance_cost,
      replacement_cost = estimate$replacement_cost,
      probability_0 = estimate$probabilities[1L],
      probability_1 = estimate$probabilities[2L],
      probability_2 = estimate$probabilities[3L]
    ),
    standard_errors = c(
      maintenance_cost = sqrt(covariance[1L, 1L]),
      replacement_cost = sqrt(covariance[2L, 2L]),
      probability_0 = sqrt(covariance[3L, 3L]),
      probability_1 = sqrt(covariance[4L, 4L]),
      probability_2 = NA_real_
    ),
    objective = fit$value,
    convergence = fit$convergence,
    expected_value = fixed_point$expected_value,
    transition_matrix = fixed_point$transition_matrix,
    fixed_point_iterations = fixed_point$iterations,
    panel = panel
  )
}

validate_bus_panel <- function(panel, n_states = 90L) {
  required <- c("state", "replace", "increment")
  missing_columns <- setdiff(required, names(panel))
  if (length(missing_columns) > 0L) {
    stop("panel is missing columns: ", paste(missing_columns, collapse = ", "))
  }
  panel <- panel[!is.na(panel$state) & !is.na(panel$replace), required]
  panel$state <- pmin(pmax(as.integer(panel$state), 0L), n_states - 1L)
  panel$replace <- as.integer(panel$replace > 0)
  panel$increment <- ifelse(
    is.na(panel$increment),
    NA_integer_,
    pmin(pmax(as.integer(panel$increment), 0L), 2L)
  )
  if (nrow(panel) < 20L) stop("at least 20 observations are required")
  panel
}

prepare_panel_csv <- function(data, n_states = 90L) {
  if (all(c("state", "replace", "increment") %in% names(data))) {
    return(validate_bus_panel(data, n_states))
  }

  required <- c("bus_id", "period", "mileage", "replace")
  missing_columns <- setdiff(required, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      "CSV must contain state/replace/increment or bus_id/period/mileage/replace"
    )
  }

  data <- data[order(data$bus_id, data$period), ]
  data$state <- pmin(floor(data$mileage / 5000), n_states - 1L)
  data$increment <- NA_integer_
  for (bus in unique(data$bus_id)) {
    index <- which(data$bus_id == bus)
    if (length(index) > 1L) {
      next_increment <- pmax(diff(data$state[index]), 0L)
      data$increment[index[-length(index)]] <- pmin(next_increment, 2L)
    }
  }
  validate_bus_panel(data, n_states)
}

prepare_matlab_matrix <- function(raw_matrix, n_states = 90L) {
  if (nrow(raw_matrix) < 128L || ncol(raw_matrix) < 37L) {
    stop("raw MATLAB-style input must have at least 128 rows and 37 columns")
  }
  raw_matrix <- as.matrix(raw_matrix[1:128, 1:37])
  storage.mode(raw_matrix) <- "double"

  mileage <- raw_matrix[12:128, ]
  first_replacement_mileage <- raw_matrix[6, ]
  replacement <- matrix(0L, nrow = 117L, ncol = 37L)

  for (bus in seq_len(37L)) {
    replacement_time <- which(mileage[, bus] > first_replacement_mileage[bus])[1L]
    if (!is.na(replacement_time)) {
      mileage[replacement_time:117L, bus] <-
        mileage[replacement_time:117L, bus] - first_replacement_mileage[bus]
      if (replacement_time > 1L) replacement[replacement_time - 1L, bus] <- 1L
    }
  }
  mileage[mileage[, 20L] > 172100, 20L] <-
    mileage[mileage[, 20L] > 172100, 20L] - 172100
  replacement[80L, 20L] <- 1L

  state <- pmin(pmax(floor(mileage / 5000), 0L), n_states - 1L)
  increment <- apply(state, 2L, function(z) c(pmin(pmax(diff(z), 0L), 2L), NA))
  validate_bus_panel(data.frame(
    state = as.vector(state),
    replace = as.vector(replacement),
    increment = as.vector(increment)
  ), n_states)
}

simulate_bus_panel <- function(
    n_buses = 37L,
    periods = 117L,
    maintenance_cost = 0.0613,
    replacement_cost = 9.6435,
    probabilities = c(0.3996, 0.5876, 0.0138),
    discount = 0.99,
    n_states = 90L,
    seed = 8830L) {
  probabilities <- probabilities / sum(probabilities)
  fixed_point <- solve_expected_value(
    maintenance_cost,
    replacement_cost,
    probabilities,
    discount,
    n_states
  )
  if (!fixed_point$converged) stop("simulation fixed point did not converge")
  set.seed(seed)

  output <- vector("list", n_buses)
  for (bus in seq_len(n_buses)) {
    state <- 0L
    records <- data.frame(
      state = integer(periods),
      replace = integer(periods),
      increment = integer(periods)
    )
    for (period in seq_len(periods)) {
      value_keep <- -maintenance_cost * state +
        discount * fixed_point$expected_value[state + 1L]
      value_replace <- -replacement_cost + discount * fixed_point$expected_value[1L]
      replace <- rbinom(1L, 1L, plogis(value_replace - value_keep))
      increment <- sample(0:2, 1L, prob = probabilities)

      records[period, ] <- c(state, replace, increment)
      state <- if (replace == 1L) {
        increment
      } else {
        min(state + increment, n_states - 1L)
      }
    }
    output[[bus]] <- records
  }
  do.call(rbind, output)
}

save_bus_results <- function(fit, results_dir = "results") {
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  table <- data.frame(
    parameter = names(fit$coefficients),
    estimate = unname(fit$coefficients),
    standard_error = unname(fit$standard_errors)
  )
  write.csv(table, file.path(results_dir, "03_nfxp_estimates.csv"), row.names = FALSE)

  state <- 0:(length(fit$expected_value) - 1L)
  c0 <- fit$coefficients["maintenance_cost"]
  c1 <- fit$coefficients["replacement_cost"]
  value_keep <- -c0 * state + 0.99 * fit$expected_value
  value_replace <- rep(-c1 + 0.99 * fit$expected_value[1L], length(state))
  probability_replace <- plogis(value_replace - value_keep)

  png(file.path(results_dir, "03_value_functions.png"), 900, 650)
  plot(state, value_keep, type = "l", lwd = 2, xlab = "Mileage state",
       ylab = "Choice-specific value", col = "black")
  lines(state, value_replace, lwd = 2, col = "red", lty = 2)
  legend("bottomleft", c("Keep engine", "Replace engine"),
         col = c("black", "red"), lty = c(1, 2), lwd = 2, bty = "n")
  dev.off()

  png(file.path(results_dir, "03_replacement_probability.png"), 900, 650)
  plot(state, probability_replace, type = "l", lwd = 2,
       xlab = "Mileage state", ylab = "Replacement probability", ylim = c(0, 1))
  dev.off()
  invisible(table)
}

print_bus_results <- function(fit) {
  cat("Bus engine replacement NFXP model\n")
  cat("Optimizer convergence code:", fit$convergence, "\n")
  cat("Negative log-likelihood:", fit$objective, "\n")
  cat("Fixed-point iterations:", fit$fixed_point_iterations, "\n\n")
  print(data.frame(
    estimate = fit$coefficients,
    standard_error = fit$standard_errors
  ))
}

if (sys.nframe() == 0L) {
  arguments <- commandArgs(trailingOnly = TRUE)
  if (length(arguments) >= 1L) {
    input <- read.csv(arguments[1L], check.names = FALSE)
    if (nrow(input) >= 128L && ncol(input) >= 37L &&
        !all(c("state", "replace", "increment") %in% names(input))) {
      panel <- prepare_matlab_matrix(input)
      message("Using MATLAB-style 128 x 37 bus matrix: ", arguments[1L])
    } else {
      panel <- prepare_panel_csv(input)
      message("Using prepared bus panel: ", arguments[1L])
    }
  } else {
    panel <- simulate_bus_panel()
    message("No CSV supplied; running the reproducible simulated bus panel.")
  }

  fit <- estimate_bus_replacement_model(panel)
  print_bus_results(fit)
  save_bus_results(fit)
}
