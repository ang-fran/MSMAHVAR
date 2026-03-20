---
title: "MSMAH-VAR"
output: github_document
---

# MSMAH-VAR: EM Estimation

This R package implements **Markov-Switching Mean-Adjusted VAR (MSMAH-VAR)** models:

- State-dependent means, AR coefficients, and covariances
- Exact EM algorithm for augmented states
- Works on simulated or real datasets

```r
library(MSMAHVAR)
# simulate or load data
