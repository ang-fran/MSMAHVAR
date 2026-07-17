#' @importFrom stats rnorm
NULL
# ===============================================
# MSM-VAR R Functions
# ===============================================

#' Log-sum-exp trick for numerical stability
#'
#' @param v Numeric vector
#' @return log(sum(exp(v)))
#' @export
logsumexp <- function(v) {
  m <- max(v)
  m + log(sum(exp(v - m)))
}

#' Log-density of multivariate normal
#'
#' @param x Numeric vector
#' @param mean Mean vector
#' @param Sigma Covariance matrix
#' @param eps Small number to stabilize inversion
#' @return Log-density
#' @export
dmvnorm_log <- function(x, mean, Sigma, eps = 1e-12) {
  x <- as.numeric(x); mean <- as.numeric(mean)
  k <- length(x)

  S <- Sigma
  diag(S) <- diag(S) + eps

  R <- chol(S)
  z <- backsolve(R, x - mean, transpose = TRUE)
  quad <- sum(z^2)
  logdet <- 2 * sum(log(diag(R)))

  -0.5 * (k * log(2*pi) + logdet + quad)
}

#' Generate one multivariate normal sample
#'
#' @param mean Mean vector
#' @param Sigma Covariance matrix
#' @return Numeric vector sample
#' @export
rmvnorm1 <- function(mean, Sigma) {
  k <- length(mean)
  R <- chol(Sigma)
  as.numeric(mean + t(R) %*% rnorm(k))
}

#' Construct a regime-switching specification
#'
#' Declares which parameter blocks are allowed to vary by regime versus
#' being pooled into a single shared estimate across all regimes. Passing
#' the result to \code{ms_var_m_step} (via \code{ms_var_em}) controls how
#' the M-step aggregates weighted sufficient statistics for each block.
#'
#' @param mu Logical; if \code{TRUE} (default), the mean vector varies by
#'   regime. If \code{FALSE}, a single mean is estimated by pooling
#'   responsibility-weighted statistics across all regimes.
#' @param A Logical; if \code{TRUE} (default), the AR coefficient matrices
#'   vary by regime. If \code{FALSE}, a single set of AR matrices is
#'   estimated by pooling weighted least-squares statistics across regimes.
#' @param Sigma Logical; if \code{TRUE} (default), the covariance matrix
#'   varies by regime. If \code{FALSE}, a single covariance matrix is
#'   estimated by pooling weighted residual outer products across regimes.
#' @return A list of class \code{msmvar_switching} with logical elements
#'   \code{mu}, \code{A}, \code{Sigma}.
#' @examples
#' make_switching()                       # fully switching (default)
#' make_switching(A = FALSE)              # AR coefficients pooled/shared
#' make_switching(mu = FALSE, A = FALSE)  # only Sigma switches by regime
#' @export
make_switching <- function(mu = TRUE, A = TRUE, Sigma = TRUE) {
  stopifnot(is.logical(mu), length(mu) == 1,
            is.logical(A), length(A) == 1,
            is.logical(Sigma), length(Sigma) == 1)
  structure(
    list(mu = mu, A = A, Sigma = Sigma),
    class = "msmvar_switching"
  )
}

#' @export
print.msmvar_switching <- function(x, ...) {
  block_status <- function(flag) if (flag) "switching by regime" else "pooled across regimes"
  cat("MSMVAR switching specification:\n")
  cat("  mu    :", block_status(x$mu), "\n")
  cat("  A     :", block_status(x$A), "\n")
  cat("  Sigma :", block_status(x$Sigma), "\n")
  invisible(x)
}

#' Simulate MSM-VAR(p) data
#'
#' @param T Number of time points
#' @param M Number of regimes
#' @param k Number of variables
#' @param p VAR lag order
#' @param mu List of state-dependent means
#' @param A_list List of state-dependent AR matrices
#' @param Sigma_list List of state-dependent covariances
#' @param P Transition matrix
#' @param pi Initial regime probabilities (default uniform)
#' @return List with Y (data) and s (hidden regimes)
#' @export
simulate_ms_var <- function(T, M, k, p, mu, A_list, Sigma_list, P, pi = NULL) {
  stopifnot(length(mu) == M, length(A_list) == M, length(Sigma_list) == M)
  stopifnot(nrow(P) == M, ncol(P) == M)
  for (j in 1:M) {
    stopifnot(length(mu[[j]]) == k)
    stopifnot(length(A_list[[j]]) == p)
    for (lag in 1:p) stopifnot(all(dim(A_list[[j]][[lag]]) == c(k, k)))
    stopifnot(all(dim(Sigma_list[[j]]) == c(k, k)))
  }
  if (is.null(pi)) pi <- rep(1/M, M)
  pi <- pi / sum(pi)

  s <- integer(T)
  Y <- matrix(0, nrow = T, ncol = k)

  # init
  s[1] <- sample(1:M, 1, prob = pi)
  Y[1, ] <- rmvnorm1(mu[[s[1]]], Sigma_list[[s[1]]])

  # early times: unconditional mixture per state
  if (p > 1) {
    for (t in 2:min(p, T)) {
      s[t] <- sample(1:M, 1, prob = P[s[t-1], ])
      Y[t, ] <- rmvnorm1(mu[[s[t]]], Sigma_list[[s[t]]])
    }
  }

  # recursion
  if (T >= p + 1) {
    for (t in (p+1):T) {
      s[t] <- sample(1:M, 1, prob = P[s[t-1], ])
      j <- s[t]

      mean_t <- mu[[j]]
      for (lag in 1:p) {
        mean_t <- mean_t + A_list[[j]][[lag]] %*% (Y[t-lag, ] - mu[[ s[t-lag] ]])
      }

      Y[t, ] <- as.numeric(mean_t) + rmvnorm1(rep(0, k), Sigma_list[[j]])
    }
  }

  list(Y = Y, s = s)
}

#' Decode augmented state index
#' @param idx Integer index
#' @param M Number of regimes
#' @param p Lag order
#' @return Vector (s_t, s_{t-1},...,s_{t-p})
#' @export
decode_z <- function(idx, M, p) {
  x <- idx - 1L
  st <- (x %% M) + 1L
  x <- x %/% M

  base <- M + 1L
  lags <- integer(p)
  for (k in 1:p) {
    lags[k] <- (x %% base)   # 0..M
    x <- x %/% base
  }
  c(st, lags)
}

#' Encode augmented state z_t
#' @param z Vector (s_t, s_{t-1},...,s_{t-p})
#' @param M Number of regimes
#' @param p Lag order
#' @return Integer index
#' @export
encode_z <- function(z, M, p) {
  # z length p+1: z[1] in 1..M, z[2:(p+1)] in 0..M
  st <- z[1]
  lags <- z[-1]
  base <- M + 1L

  x <- (st - 1L)
  mult <- M
  pow <- 1L
  for (k in 1:p) {
    x <- x + mult * (lags[k] * pow)
    pow <- pow * base
  }
  x + 1L
}

#' Compute log-likelihood for all t
#' @param Y Data matrix (T x k)
#' @param mu List of regime means
#' @param A_list List of AR matrices
#' @param Sigma_list List of covariance matrices
#' @param M Number of regimes
#' @param p Lag order
#' @param eps Numerical stability constant
#' @export
compute_log_g_allt <- function(Y, mu, A_list, Sigma_list, M, p, eps = 1e-12) {
  Tn <- nrow(Y)
  k  <- ncol(Y)
  K  <- M * (M + 1)^p

  logg <- matrix(-Inf, nrow = K, ncol = Tn)

  for (t in 1:Tn) {
    q <- min(p, t - 1)
    for (z in 1:K) {
      zt <- decode_z(z, M, p)   # (s_t, s_\{t-1\},...,s_\{t-p\})
      st <- zt[1]

      mean_t <- mu[[st]]
      if (q > 0) {
        for (lag in 1:q) {
          s_lag <- zt[lag + 1]  # 0..M
          if (s_lag == 0) next
          mean_t <- mean_t + A_list[[st]][[lag]] %*% (Y[t-lag, ] - mu[[s_lag]])
        }
      }

      logg[z, t] <- dmvnorm_log(Y[t, ], mean = as.numeric(mean_t),
                                Sigma = Sigma_list[[st]], eps = eps)
    }
  }

  list(logg = logg, K = K)
}

#' Compute fitted values and residuals for an estimated MSMVAR model
#'
#' Derives fitted values and residuals directly from a set of estimated
#' MSMVAR parameters (\code{mu}, \code{A_list}, \code{Sigma_list},
#' \code{P_hat}). Internally re-runs the E-step (\code{ms_var_e_step}) at
#' the estimated parameters to obtain the smoothed augmented-state
#' responsibilities \code{gammaZ}, then constructs, at each t, the
#' responsibility-weighted expectation of the conditional mean over
#' \code{gammaZ[t, ]} — i.e. fitted values marginalize over regime
#' uncertainty rather than conditioning on a single hard-classified regime
#' path. Uses the same augmented-state conditional-mean construction as
#' \code{compute_log_g_allt}.
#'
#' @param Y Data matrix (T x k)
#' @param mu Estimated list of regime means
#' @param A_list Estimated list of regime AR matrices
#' @param Sigma_list Estimated list of regime covariance matrices
#' @param P_hat Estimated transition matrix
#' @param p VAR lag order
#' @param pi Initial regime probabilities used for the E-step (default
#'   uniform)
#' @param eps Small constant for numerical stability, passed to the E-step
#' @return A list with \code{fitted} (T x k matrix of fitted values),
#'   \code{resid} (T x k matrix of residuals, \code{Y - fitted}), and
#'   \code{gammaZ} (the smoothed augmented-state responsibilities used to
#'   construct the fitted values, for reuse/inspection)
#' @export
compute_fitted_msmvar <- function(Y, mu, A_list, Sigma_list, P_hat, p,
                                  pi = NULL, eps = 1e-12) {
  Tn <- nrow(Y)
  k  <- ncol(Y)
  M  <- length(mu)
  stopifnot(nrow(P_hat) == M, ncol(P_hat) == M)

  E <- ms_var_e_step(Y, P_hat, mu, A_list, Sigma_list, p, pi = pi, eps = eps)
  gammaZ <- E$gammaZ
  K <- ncol(gammaZ)

  fitted <- matrix(0, nrow = Tn, ncol = k)

  for (t in 1:Tn) {
    q <- min(p, t - 1)
    mean_t_accum <- rep(0, k)

    for (z in 1:K) {
      w <- gammaZ[t, z]
      if (w <= 0) next
      zt <- decode_z(z, M, p)
      st <- zt[1]

      mean_t <- mu[[st]]
      if (q > 0) {
        for (lag in 1:q) {
          s_lag <- zt[lag + 1]
          if (s_lag == 0) next
          mean_t <- mean_t + A_list[[st]][[lag]] %*% (Y[t-lag, ] - mu[[s_lag]])
        }
      }

      mean_t_accum <- mean_t_accum + w * as.numeric(mean_t)
    }

    fitted[t, ] <- mean_t_accum
  }

  resid <- Y - fitted
  list(fitted = fitted, resid = resid, gammaZ = gammaZ)
}

#' E-step for MSM-VAR EM estimation
#'
#' Runs the forward-backward algorithm (in log space, with scaling) over the
#' augmented state space to produce smoothed regime probabilities and the
#' data log-likelihood, given a fixed set of model parameters. Called
#' internally by \code{ms_var_em} at each EM iteration, but also exported
#' for direct use — e.g. to compute smoothed probabilities for a already-
#' fitted model without re-running the full EM loop.
#'
#' @export
ms_var_e_step <- function(Y, P_hat, mu, A_list, Sigma_list, p, pi = NULL, eps = 1e-12) {
  Tn <- nrow(Y)
  M  <- nrow(P_hat)
  stopifnot(ncol(P_hat) == M)

  if (is.null(pi)) pi <- rep(1/M, M)
  pi <- pi / sum(pi)

  Eg <- compute_log_g_allt(Y, mu, A_list, Sigma_list, M, p, eps = eps)
  logg <- Eg$logg
  K <- Eg$K

  # forward-backward with scaling (in log space)
  logalpha <- matrix(-Inf, nrow = Tn, ncol = K)
  logbeta  <- matrix(-Inf, nrow = Tn, ncol = K)
  csc <- numeric(Tn)

  # initial distribution for z1: only (j,0,...,0) has prob pi[j]
  logpi_z1 <- rep(-Inf, K)
  for (j in 1:M) {
    z1 <- c(j, rep(0L, p))
    logpi_z1[encode_z(z1, M, p)] <- log(pi[j] + eps)
  }

  # t=1
  la1 <- logpi_z1 + logg[, 1]
  csc[1] <- logsumexp(la1)
  logalpha[1, ] <- la1 - csc[1]

  # transitions in augmented space depend only on s_{t-1} -> s_t (shift stack)
  for (t in 2:Tn) {
    for (znew in 1:K) logalpha[t, znew] <- -Inf

    # build from zprev
    for (zprev in 1:K) {
      ap <- logalpha[t-1, zprev]
      if (!is.finite(ap)) next

      zprev_vec <- decode_z(zprev, M, p)
      i <- zprev_vec[1]  # s_{t-1}

      for (j in 1:M) {
        znew_vec <- c(j, zprev_vec[1:p])
        znew <- encode_z(znew_vec, M, p)
        logalpha[t, znew] <- logsumexp(c(logalpha[t, znew], ap + log(P_hat[i, j] + eps)))
      }
    }

    logalpha[t, ] <- logalpha[t, ] + logg[, t]
    csc[t] <- logsumexp(logalpha[t, ])
    logalpha[t, ] <- logalpha[t, ] - csc[t]
  }

  loglik <- sum(csc)

  # backward init
  logbeta[Tn, ] <- 0

  for (t in (Tn-1):1) {
    for (zprev in 1:K) {
      zprev_vec <- decode_z(zprev, M, p)
      i <- zprev_vec[1]

      acc <- -Inf
      for (j in 1:M) {
        znew_vec <- c(j, zprev_vec[1:p])
        znew <- encode_z(znew_vec, M, p)
        acc <- logsumexp(c(acc, log(P_hat[i, j] + eps) + logg[znew, t+1] + logbeta[t+1, znew]))
      }
      logbeta[t, zprev] <- acc - csc[t+1]
    }
  }

  # smoothed gammaZ
  gammaZ <- matrix(0, nrow = Tn, ncol = K)
  for (t in 1:Tn) {
    lg <- logalpha[t, ] + logbeta[t, ]
    denom <- logsumexp(lg)
    gammaZ[t, ] <- exp(lg - denom)
  }

  # smoothed marginals P(s_t=j | Y)
  smoothed_marg <- matrix(0, nrow = Tn, ncol = M)
  for (t in 1:Tn) {
    for (z in 1:K) {
      st <- decode_z(z, M, p)[1]
      smoothed_marg[t, st] <- smoothed_marg[t, st] + gammaZ[t, z]
    }
  }

  # smoothed joint P(s_{t-1}=i, s_t=j | Y) for t>=2
  smoothed_joint <- vector("list", Tn)
  smoothed_joint[[1]] <- NULL

  for (t in 2:Tn) {
    J <- matrix(0, M, M)

    for (zprev in 1:K) {
      zprev_vec <- decode_z(zprev, M, p)
      i <- zprev_vec[1]

      ap <- logalpha[t-1, zprev]
      if (!is.finite(ap)) next

      for (j in 1:M) {
        znew_vec <- c(j, zprev_vec[1:p])
        znew <- encode_z(znew_vec, M, p)

        # proportional to alpha_{t-1}(zprev) * P(i,j) * f(y_t|znew) * beta_t(znew)
        J[i, j] <- J[i, j] + exp(ap) * (P_hat[i, j]) * exp(logg[znew, t]) * exp(logbeta[t, znew])      }
    }

    smoothed_joint[[t]] <- J / (sum(J) + eps)
  }

  list(
    loglik = loglik,
    gammaZ = gammaZ,
    smoothed_marginals = smoothed_marg,
    smoothed_joint = smoothed_joint
  )
}

#' M-step for MSM-VAR EM estimation
#'
#' Updates \code{mu}, \code{A_list}, \code{Sigma_list}, and \code{P_hat}
#' given the smoothed responsibilities from the E-step. The
#' \code{switching} argument (see \code{\link{make_switching}}) controls,
#' block by block, whether weighted sufficient statistics are kept separate
#' per regime or pooled into a single shared estimate. Called internally by
#' \code{ms_var_em} at each EM iteration, but also exported for direct use.
#'
#' @export
ms_var_m_step <- function(Y, gammaZ, smoothed_joint, mu_prev, A_prev, p,
                          switching = make_switching(), eps = 1e-10) {
  Tn <- nrow(Y)
  k  <- ncol(Y)
  M  <- length(mu_prev)
  K  <- ncol(gammaZ)
  stopifnot(K == M * (M + 1)^p)
  stopifnot(inherits(switching, "msmvar_switching"))

  # ---- (A) Accumulate per-regime WLS sufficient statistics for A ----
  XtWX_list <- vector("list", M)
  XtWY_list <- vector("list", M)

  for (j in 1:M) {
    XtWX <- matrix(0, nrow = k*p, ncol = k*p)
    XtWY <- matrix(0, nrow = k*p, ncol = k)

    for (t in 1:Tn) {
      q <- min(p, t - 1)
      for (z in 1:K) {
        w <- gammaZ[t, z]
        if (w <= 0) next
        zt <- decode_z(z, M, p)
        if (zt[1] != j) next

        # build regressor vector X (kp x 1) block by lag
        Xvec <- matrix(0, nrow = k*p, ncol = 1)
        if (q > 0) {
          for (lag in 1:q) {
            s_lag <- zt[lag + 1]
            if (s_lag == 0) next
            block <- (1 + (lag-1)*k):(lag*k)
            Xvec[block, 1] <- (Y[t-lag, ] - mu_prev[[s_lag]])
          }
        }

        ytil <- matrix(Y[t, ] - mu_prev[[j]], ncol = 1)  # k x 1

        XtWX <- XtWX + w * (Xvec %*% t(Xvec))
        XtWY <- XtWY + w * (Xvec %*% t(ytil))            # (kp x k)
      }
    }

    XtWX_list[[j]] <- XtWX
    XtWY_list[[j]] <- XtWY
  }

  solve_A_block <- function(XtWX, XtWY) {
    B <- solve(XtWX + eps * diag(k*p), XtWY)             # (kp x k)
    Aj <- vector("list", p)
    for (lag in 1:p) {
      rows <- (1 + (lag-1)*k):(lag*k)                    # k rows
      Aj[[lag]] <- t(B[rows, , drop = FALSE])            # k x k
    }
    Aj
  }

  A_new <- vector("list", M)
  if (switching$A) {
    for (j in 1:M) A_new[[j]] <- solve_A_block(XtWX_list[[j]], XtWY_list[[j]])
  } else {
    # pool WLS sufficient statistics across regimes -> one shared A for all j
    XtWX_pool <- Reduce(`+`, XtWX_list)
    XtWY_pool <- Reduce(`+`, XtWY_list)
    Aj_shared <- solve_A_block(XtWX_pool, XtWY_pool)
    for (j in 1:M) A_new[[j]] <- Aj_shared
  }

  # ---- (B) Accumulate per-regime sufficient statistics for mu ----
  mu_num <- vector("list", M)
  mu_den <- numeric(M)

  for (j in 1:M) {
    num <- rep(0, k)
    den <- 0

    for (t in 1:Tn) {
      q <- min(p, t - 1)
      for (z in 1:K) {
        w <- gammaZ[t, z]
        if (w <= 0) next
        zt <- decode_z(z, M, p)
        if (zt[1] != j) next

        pred_no_mu <- rep(0, k)
        if (q > 0) {
          for (lag in 1:q) {
            s_lag <- zt[lag + 1]
            if (s_lag == 0) next
            pred_no_mu <- pred_no_mu + A_new[[j]][[lag]] %*% (Y[t-lag, ] - mu_prev[[s_lag]])
          }
        }

        num <- num + w * (Y[t, ] - as.numeric(pred_no_mu))
        den <- den + w
      }
    }

    mu_num[[j]] <- num
    mu_den[j] <- den
  }

  mu_new <- vector("list", M)
  if (switching$mu) {
    for (j in 1:M) mu_new[[j]] <- as.numeric(mu_num[[j]] / (mu_den[j] + eps))
  } else {
    # pool weighted sums across regimes -> one shared mu for all j
    num_pool <- Reduce(`+`, mu_num)
    den_pool <- sum(mu_den)
    shared_mu <- as.numeric(num_pool / (den_pool + eps))
    for (j in 1:M) mu_new[[j]] <- shared_mu
  }

  # ---- (C) Accumulate per-regime sufficient statistics for Sigma ----
  S_list <- vector("list", M)
  S_den  <- numeric(M)

  for (j in 1:M) {
    S <- matrix(0, k, k)
    den <- 0

    for (t in 1:Tn) {
      q <- min(p, t - 1)
      for (z in 1:K) {
        w <- gammaZ[t, z]
        if (w <= 0) next
        zt <- decode_z(z, M, p)
        if (zt[1] != j) next

        mean_t <- mu_new[[j]]
        if (q > 0) {
          for (lag in 1:q) {
            s_lag <- zt[lag + 1]
            if (s_lag == 0) next
            mean_t <- mean_t + A_new[[j]][[lag]] %*% (Y[t-lag, ] - mu_new[[s_lag]])
          }
        }

        e <- as.numeric(Y[t, ] - as.numeric(mean_t))
        S <- S + w * tcrossprod(e)
        den <- den + w
      }
    }

    S_list[[j]] <- S
    S_den[j] <- den
  }

  Sigma_new <- vector("list", M)
  if (switching$Sigma) {
    for (j in 1:M) {
      S <- S_list[[j]] / (S_den[j] + eps)
      diag(S) <- diag(S) + eps
      Sigma_new[[j]] <- S
    }
  } else {
    # pool weighted residual outer products across regimes -> one shared Sigma
    S_pool <- Reduce(`+`, S_list)
    den_pool <- sum(S_den)
    S_shared <- S_pool / (den_pool + eps)
    diag(S_shared) <- diag(S_shared) + eps
    for (j in 1:M) Sigma_new[[j]] <- S_shared
  }

  # ---- (D) Update transition matrix P (always regime-specific) ----
  P_new <- matrix(0, M, M)
  for (i in 1:M) {
    for (j in 1:M) {
      P_new[i, j] <- sum(sapply(2:Tn, function(t) smoothed_joint[[t]][i, j]))
    }
    P_new[i, ] <- P_new[i, ] / (sum(P_new[i, ]) + eps)
  }

  list(mu = mu_new, A_list = A_new, Sigma_list = Sigma_new, P_hat = P_new)
}

#' EM estimation for MSM-VAR model
#'
#' Fits a Markov-Switching Mean-Adjusted VAR(p) model using the EM algorithm.
#' @param Y Data matrix (T x k)
#' @param p VAR lag order
#' @param mu_init Initial mean list
#' @param A_init Initial AR matrices list
#' @param Sigma_init Initial covariance matrices list
#' @param P_init Initial transition matrix
#' @param pi Initial regime probabilities (optional)
#' @param switching A switching specification from \code{\link{make_switching}}
#'   controlling which of \code{mu}, \code{A}, \code{Sigma} vary by regime
#'   versus are pooled across regimes. Defaults to fully switching, i.e.
#'   \code{make_switching()}, which matches the previous (pre-MSMVAR)
#'   behavior of this function.
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @param eps Small constant for stability
#' @param step_size Controls step size between iterations
#' @param verbose Print progress
#' @return List with fitted parameters and log-likelihood
#' @export
ms_var_em <- function(Y, p,
                      mu_init, A_init, Sigma_init, P_init,
                      pi = NULL,
                      switching = make_switching(),
                      max_iter = 300, tol = 1e-6, eps = 1e-10,
                      step_size = 0.5,        # <-- new parameter
                      verbose = TRUE) {

  stopifnot(inherits(switching, "msmvar_switching"))

  Tn <- nrow(Y)
  M <- nrow(P_init)
  stopifnot(ncol(P_init) == M)
  if (is.null(pi)) pi <- rep(1/M, M)

  mu_curr <- mu_init
  A_curr  <- A_init
  S_curr  <- Sigma_init
  P_curr  <- P_init

  ll <- numeric(max_iter)

  for (iter in 1:max_iter) {
    E    <- ms_var_e_step(Y, P_curr, mu_curr, A_curr, S_curr, p, pi = pi, eps = eps)
    ll[iter] <- E$loglik

    Mres <- ms_var_m_step(Y, E$gammaZ, E$smoothed_joint, mu_curr, A_curr, p,
                          switching = switching, eps = eps)

    # ---- Damped update: theta_new = theta_old + step_size * (theta_mstep - theta_old) ----

    # mu: list of k-vectors
    mu_next <- lapply(seq_len(M), function(j)
      mu_curr[[j]] + step_size * (Mres$mu[[j]] - mu_curr[[j]])
    )

    # A_list: list (M) of lists (p) of k x k matrices
    A_next <- lapply(seq_len(M), function(j)
      lapply(seq_len(p), function(lag)
        A_curr[[j]][[lag]] + step_size * (Mres$A_list[[j]][[lag]] - A_curr[[j]][[lag]])
      )
    )

    # Sigma_list: list of k x k matrices
    S_next <- lapply(seq_len(M), function(j)
      S_curr[[j]] + step_size * (Mres$Sigma_list[[j]] - S_curr[[j]])
    )

    # P: M x M matrix
    P_next <- P_curr + step_size * (Mres$P_hat - P_curr)
    # Re-normalise rows so P stays a valid transition matrix
    P_next <- P_next / rowSums(P_next)

    mu_curr <- mu_next
    A_curr  <- A_next
    S_curr  <- S_next
    P_curr  <- P_next

    if (verbose && iter %% 10 == 0) cat("iter:", iter, " loglik:", ll[iter], "\n")

    if (iter > 1 && is.finite(ll[iter]) && abs(ll[iter] - ll[iter-1]) < tol) {
      ll <- ll[1:iter]
      if (verbose) cat("Converged at iteration", iter, "\n")
      break
    }
  }

  list(mu = mu_curr, A_list = A_curr, Sigma_list = S_curr, P_hat = P_curr, loglik = ll,
       switching = switching)
}
#' Compare fitted vs true parameters
#' @param mu_true List of true mean vectors
#' @param A_true List of true AR matrices
#' @param S_true List of true covariance matrices
#' @param P_true True transition matrix
#' @param mu_hat List of estimated mean vectors
#' @param A_hat List of estimated AR matrices
#' @param S_hat List of estimated covariance matrices
#' @param P_hat Estimated transition matrix
#' @param digits Number of decimal places to print
#' @export
compare_ms_var <- function(mu_true, A_true, S_true, P_true,
                           mu_hat,  A_hat,  S_hat,  P_hat,
                           digits = 4) {

  M <- nrow(P_true)
  p <- length(A_true[[1]])

  cost <- function(ord) {
    c_mu <- 0; c_A <- 0; c_S <- 0; c_P <- sum((P_hat[ord, ord] - P_true)^2)
    for (j in 1:M) {
      jj <- ord[j]
      c_mu <- c_mu + sum((mu_hat[[jj]] - mu_true[[j]])^2)
      c_S  <- c_S  + sum((S_hat[[jj]] - S_true[[j]])^2)
      for (lag in 1:p) c_A <- c_A + sum((A_hat[[jj]][[lag]] - A_true[[j]][[lag]])^2)
    }
    c_mu + c_A + c_S + c_P
  }

  ord <- 1:M
  if (M == 2 && cost(c(2,1)) < cost(c(1,2))) ord <- c(2,1)

  mu_hat2 <- mu_hat[ord]
  A_hat2  <- A_hat[ord]
  S_hat2  <- S_hat[ord]
  P_hat2  <- P_hat[ord, ord, drop=FALSE]

  cat("\nP (true vs hat)\n")
  print(round(P_true, digits)); cat("\n"); print(round(P_hat2, digits))

  cat("\nmu (true vs hat)\n")
  for (j in 1:M) {
    cat("\nRegime", j, "\n")
    print(round(cbind(true = mu_true[[j]], hat = mu_hat2[[j]]), digits))
  }

  cat("\nA matrices (true vs hat)\n")
  for (j in 1:M) {
    for (lag in 1:p) {
      cat("\nRegime", j, " Lag", lag, "\n")
      cat("A_true:\n"); print(round(A_true[[j]][[lag]], digits))
      cat("A_hat:\n");  print(round(A_hat2[[j]][[lag]], digits))
    }
  }

  cat("\nSigma (true vs hat)\n")
  for (j in 1:M) {
    cat("\nRegime", j, "\n")
    cat("Sigma_true:\n"); print(round(S_true[[j]], digits))
    cat("Sigma_hat:\n");  print(round(S_hat2[[j]], digits))
  }

  invisible(list(order = ord, mu_hat = mu_hat2, A_hat = A_hat2, S_hat = S_hat2, P_hat = P_hat2))
}
