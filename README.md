---
title: MSM-VAR
output: github_document
---

# MSM-VAR: Markov-Switching Mean-Adjusted VAR — R Package

An R package implementing the **Markov-Switching Mean-Adjusted VAR (MSM-VAR)** model,
where the mean, autoregressive coefficients, and covariance matrix are interchangeably
state-dependent. Estimation is performed via an exact EM algorithm operating on
augmented states.

---

## Motivation

Standard VAR models assume fixed parameters across time, which is often unrealistic
for economic and financial data that exhibit structural breaks or regime changes
(e.g., expansions vs. recessions, high vs. low volatility periods).
MSM-VAR addresses this by allowing the model's mean, dynamics, and noise
structure to switch across a discrete set of hidden regimes, governed by a
Markov chain.

---

## Features

- **State-dependent means, AR coefficients, and/or covariance matrices**
- **Exact EM algorithm** for parameter estimation with augmented state representation
- Hamilton filtering and smoothing for regime probability inference
- Compatible with both simulated data and real-world multivariate time series
- Designed for reproducibility and extensibility

---

## Installation
```r
# Install from GitHub using devtools
devtools::install_github("ang-fran/MSM-VAR")
```

---

## Usage
```r
library(MSM-VAR)

# Simulate data from a 2-state MSMAH-VAR(1) process
sim_data <- simulate_msmahvar(
  n     = 200,
  k     = 2,      # number of variables
  M     = 2,      # number of regimes
  p     = 1       # lag order
)

# Estimate the model via EM
fit <- estimate_msmahvar(
  data = sim_data$Y,
  M    = 2,
  p    = 1
)

# Inspect results
summary(fit)
plot_smoothed_probs(fit)
```

---

## Methods

Estimation follows the EM algorithm for Markov-switching models:

- **E-step:** Hamilton filter forward pass + smoother backward pass to compute
  regime probabilities
- **M-step:** Closed-form updates for regime means, AR matrices, and covariance matrices,
  weighted by smoothed probabilities
- Iteration continues until log-likelihood convergence

---

## Related Work

This package connects to broader research on regime-switching time series.
See also the [mixed-ms-ar-var](https://github.com/ang-fran/mixed-ms-ar-var) repo
for related mixed-lag Markov-switching AR and VAR implementations.

---

## Author

Angela Ibhade — Graduate Researcher, Time Series & Econometric Modeling
[Portfolio](https://github.com/ang-fran/Angela-data-analysis-portfolio)
