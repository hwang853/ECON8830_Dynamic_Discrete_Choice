# ECON 8830, Part 2: kernel density and local-constant regression
#
# Combined R translation of WHRECON8830HW2Q1.pdf and
# WHRECON8830HW2Q2.pdf. The original PDFs use age, wage, and education vectors
# with 1,000 observations but do not include the source dataset. This script
# accepts a CSV containing those variables and supplies a reproducible example
# dataset when no CSV path is provided.

epanechnikov_kernel <- function(u) {
  0.75 * (1 - u^2) * (abs(u) <= 1)
}

local_constant <- function(y, x, x0, bandwidth, kernel = epanechnikov_kernel) {
  if (length(y) != length(x)) stop("y and x must have the same length")
  if (!is.finite(bandwidth) || bandwidth <= 0) stop("bandwidth must be positive")

  weights <- kernel((x - x0) / bandwidth)
  denominator <- sum(weights)
  if (!is.finite(denominator) || denominator <= .Machine$double.eps) {
    return(NA_real_)
  }
  sum(weights * y) / denominator
}

leave_one_out_cv <- function(
    y,
    x,
    bandwidths,
    kernel = epanechnikov_kernel) {
  n <- length(y)
  if (length(x) != n) stop("y and x must have the same length")

  criterion <- numeric(length(bandwidths))
  for (b in seq_along(bandwidths)) {
    fitted <- numeric(n)
    for (i in seq_len(n)) {
      keep <- seq_len(n) != i
      fitted[i] <- local_constant(
        y = y[keep],
        x = x[keep],
        x0 = x[i],
        bandwidth = bandwidths[b],
        kernel = kernel
      )
    }
    squared_error <- (y - fitted)^2
    criterion[b] <- mean(squared_error, na.rm = TRUE)
  }

  data.frame(bandwidth = bandwidths, cv_mse = criterion)
}

bootstrap_education_effect <- function(
    wage,
    education,
    bandwidth = 4,
    education_high = 16,
    education_low = 12,
    replications = 1000L,
    seed = 8830L) {
  set.seed(seed)
  n <- length(wage)
  effects <- numeric(replications)

  for (d in seq_len(replications)) {
    index <- sample.int(n, size = n, replace = TRUE)
    effects[d] <-
      local_constant(wage[index], education[index], education_high, bandwidth) -
      local_constant(wage[index], education[index], education_low, bandwidth)
  }

  effects <- effects[is.finite(effects)]
  if (length(effects) < 2L) stop("bootstrap produced too few finite estimates")

  # Reproduces the interval formula used in the MATLAB PDF. This is an
  # interval for the simulation mean because it divides by sqrt(B).
  monte_carlo_se <- sd(effects) / sqrt(length(effects))
  matlab_style_ci <- mean(effects) +
    qt(c(0.025, 0.975), df = length(effects) - 1L) * monte_carlo_se

  # Standard percentile interval for the sampling distribution of the effect.
  percentile_ci <- unname(quantile(effects, probs = c(0.025, 0.975)))

  list(
    effects = effects,
    mean_effect = mean(effects),
    bootstrap_sd = sd(effects),
    monte_carlo_se = monte_carlo_se,
    matlab_style_ci = matlab_style_ci,
    percentile_ci = percentile_ci
  )
}

generate_example_data <- function(n = 1000L, seed = 8830L) {
  set.seed(seed)
  age <- pmin(80, pmax(18, rnorm(n, mean = 39.3, sd = 13.4)))
  education <- sample(
    8:20,
    size = n,
    replace = TRUE,
    prob = dnorm(8:20, mean = 13.5, sd = 2.5)
  )
  wage_mean <- 4 + 1.45 * education - 0.018 * (age - 42)^2
  wage <- pmax(0.5, wage_mean + rnorm(n, sd = 7))
  data.frame(age = age, wage = wage, education = education)
}

validate_data <- function(data) {
  required <- c("age", "wage", "education")
  missing_columns <- setdiff(required, names(data))
  if (length(missing_columns) > 0L) {
    stop("CSV is missing required columns: ", paste(missing_columns, collapse = ", "))
  }

  data <- data[complete.cases(data[, required]), required]
  for (name in required) data[[name]] <- as.numeric(data[[name]])
  if (nrow(data) < 10L) stop("at least 10 complete observations are required")
  data
}

run_kernel_analysis <- function(
    data,
    results_dir = "results",
    bootstrap_replications = 1000L,
    seed = 8830L) {
  data <- validate_data(data)
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  age <- data$age
  wage <- data$wage
  education <- data$education
  n <- nrow(data)

  # Q1(a): descriptive statistics and rule-of-thumb bandwidth from the PDF.
  age_sd <- sd(age)
  h_star <- 2.34 * n^(-0.2) * age_sd

  png(file.path(results_dir, "02_age_wage_scatter.png"), 900, 650)
  plot(age, wage, pch = 1, col = "steelblue", xlab = "Age", ylab = "Wage")
  dev.off()

  # Q1(b): Gaussian KDE under three bandwidths.
  age_grid <- seq(20, 60, by = 1)
  kde_bandwidths <- c(h_star / 3, h_star, 3 * h_star)
  kde <- lapply(
    kde_bandwidths,
    function(h) density(age, bw = h, kernel = "gaussian", from = 20, to = 60, n = 512)
  )

  png(file.path(results_dir, "02_age_kde_bandwidths.png"), 900, 650)
  plot(kde[[1]], col = "red", lty = 1, lwd = 2, xlim = range(age_grid),
       xlab = "Age", main = "Age density under alternative bandwidths")
  lines(kde[[2]], col = "black", lty = 3, lwd = 2)
  lines(kde[[3]], col = "blue", lty = 2, lwd = 2)
  legend(
    "topright",
    legend = c("h*/3", "h*", "3h*"),
    col = c("red", "black", "blue"),
    lty = c(1, 3, 2),
    lwd = 2,
    bty = "n"
  )
  dev.off()

  # Q1(c): bandwidth selection for wage on age.
  age_cv <- leave_one_out_cv(wage, age, bandwidths = 1:50)
  best_age_bandwidth <- age_cv$bandwidth[which.min(age_cv$cv_mse)]
  write.csv(age_cv, file.path(results_dir, "02_age_bandwidth_cv.csv"), row.names = FALSE)

  png(file.path(results_dir, "02_age_bandwidth_cv.png"), 900, 650)
  plot(age_cv$bandwidth, age_cv$cv_mse, type = "l", lwd = 2,
       xlab = "Bandwidth", ylab = "Leave-one-out mean squared error")
  abline(v = best_age_bandwidth, col = "red", lty = 2)
  dev.off()

  # Q2(a): bandwidth selection for wage on education and the wage contrast.
  education_cv <- leave_one_out_cv(wage, education, bandwidths = 1:20)
  best_education_bandwidth <- education_cv$bandwidth[which.min(education_cv$cv_mse)]
  write.csv(
    education_cv,
    file.path(results_dir, "02_education_bandwidth_cv.csv"),
    row.names = FALSE
  )

  png(file.path(results_dir, "02_education_bandwidth_cv.png"), 900, 650)
  plot(education_cv$bandwidth, education_cv$cv_mse, type = "l", lwd = 2,
       xlab = "Bandwidth", ylab = "Leave-one-out mean squared error")
  abline(v = best_education_bandwidth, col = "red", lty = 2)
  dev.off()

  # Preserve bandwidth 4 for the contrast, as used in the original PDF.
  education_effect <-
    local_constant(wage, education, 16, bandwidth = 4) -
    local_constant(wage, education, 12, bandwidth = 4)

  # Q2(b): nonparametric bootstrap.
  bootstrap <- bootstrap_education_effect(
    wage = wage,
    education = education,
    bandwidth = 4,
    replications = bootstrap_replications,
    seed = seed
  )

  summary <- data.frame(
    statistic = c(
      "observations",
      "mean_age",
      "sd_age",
      "h_star",
      "best_age_bandwidth",
      "best_education_bandwidth",
      "education_16_minus_12",
      "bootstrap_mean_effect",
      "bootstrap_percentile_lower",
      "bootstrap_percentile_upper",
      "matlab_style_lower",
      "matlab_style_upper"
    ),
    value = c(
      n,
      mean(age),
      age_sd,
      h_star,
      best_age_bandwidth,
      best_education_bandwidth,
      education_effect,
      bootstrap$mean_effect,
      bootstrap$percentile_ci,
      bootstrap$matlab_style_ci
    )
  )
  write.csv(summary, file.path(results_dir, "02_summary.csv"), row.names = FALSE)

  list(
    summary = summary,
    age_cv = age_cv,
    education_cv = education_cv,
    bootstrap = bootstrap
  )
}

print_kernel_results <- function(results) {
  cat("Kernel density and local-constant regression\n\n")
  print(results$summary, row.names = FALSE)
}

if (sys.nframe() == 0L) {
  arguments <- commandArgs(trailingOnly = TRUE)
  if (length(arguments) >= 1L) {
    input_data <- read.csv(arguments[1L])
    message("Using data from: ", arguments[1L])
  } else {
    input_data <- generate_example_data()
    message("No CSV supplied; running the reproducible example dataset.")
  }

  results <- run_kernel_analysis(input_data)
  print_kernel_results(results)
}
