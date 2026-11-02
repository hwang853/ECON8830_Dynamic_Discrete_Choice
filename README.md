# ECON8830 Dynamic Discrete Choice

This repository contains a combined R translation of the MATLAB code shown in
`PS1DGP.pdf` and `PS1LOGLIK.pdf`. It runs a Monte Carlo experiment for a dynamic
binary-logit model with state dependence.

## Model

For period `t`, the binary choice follows

```text
Pr(y_t = 1 | x_t, y_(t-1))
  = logistic(theta_0 + theta_x x_t + theta_lag y_(t-1)).
```

The default parameters preserve the values used in the PDFs:

- sample size: `n = 200`
- true/start value: `theta = (0, 1, -0.5)`
- Monte Carlo replications: `101`

The original DGP PDF used 100 replications, while the likelihood PDF used 101.
Set `replications = 100` in `run_monte_carlo()` to reproduce the former count.

## Running the code

The implementation uses only base R. From the repository directory, run:

```sh
Rscript 01_dynamic_logit_monte_carlo.R
```

The script prints parameter means, bias, Monte Carlo standard errors, the
covariance matrix of the estimates, and the average inverse Hessian. It also
writes `monte_carlo_estimates.csv`.

To use the functions interactively:

```r
source("01_dynamic_logit_monte_carlo.R")

results <- run_monte_carlo(
  n = 200,
  replications = 100,
  theta_true = c(0, 1, -0.5),
  seed = 8830
)

print_results(results)
```

## Translation notes

The PDF exports contain the MATLAB driver code and printed numerical output,
but not the referenced `dgp()` and `loglik()` function definitions or the code
that creates `y` and `x`. To make the R translation complete and runnable, this
repository supplies a standard conditional dynamic-logit DGP and the matching
negative log-likelihood. Here, the three parameters are interpreted as an
intercept, the coefficient on a normally distributed exogenous variable, and
the coefficient on the lagged binary choice.

The R implementation also corrects several apparent issues in the exported
MATLAB snippets:

- it preallocates the full coefficient matrix;
- it generates fresh `y` and `x` data for every replication;
- it retains each optimizer's Hessian instead of only the final Hessian;
- it computes Monte Carlo standard errors using the number of replications,
  rather than the total number of matrix elements;
- it uses stable binomial log-probabilities for the likelihood.

## Files

- `01_dynamic_logit_monte_carlo.R` - Part 1: combined dynamic-logit DGP,
  likelihood, maximum-likelihood estimation, and Monte Carlo analysis.
- `.gitignore` - excludes generated estimates and local R session files.
