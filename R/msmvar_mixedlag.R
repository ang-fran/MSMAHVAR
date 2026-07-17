# ================================================================
# MSM-VAR with Regime-Specific (Mixed) Lag Orders
# ================================================================
# This file extends the MSM-VAR package to support per-regime lag
# orders p_vec = c(p_1, ..., p_M).  The conditional mean follows the
# mean-adjusted form from Chapter 4:
#
#   y_t - mu_{s_t} = sum_{l=1}^{p_{s_t}} A_{s_t, l} (y_{t-l} - mu_{s_{t-l}})
#                    + epsilon_t,   epsilon_t ~ N(0, Sigma_{s_t})
#
# The augmented state space uses p_max = max(p_vec) lags for indexing,
# but each regime only populates regressors up to its own p_s.
# The M-step AR update uses pseudoinverse-based weighted OLS (MASS::ginv)
# to handle the rank deficiency that arises when regimes have different
# lag orders.
#
# Public API (mirrors existing package):
#   make_switching()                    -- unchanged, re-exported
#   ms_var_e_step_ml()              -- E-step for mixed-lag model
#   ms_var_m_step_ml()              -- M-step with pseudoinverse OLS
#   ms_var_em_ml()                  -- Full EM loop
#   ms_var_bic_grid_ml()            -- BIC grid search over (p_1,...,p_M)
# ================================================================

# ----------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------

if (!requireNamespace("MASS", quietly = TRUE))
  stop("Package 'MASS' is required for pseudoinverse OLS in ms_var_m_step_ml().")

# Re-use helpers from the installed MSM-VAR package
# (logsumexp, dmvnorm_log, rmvnorm1, decode_z, encode_z, make_switching)


# ----------------------------------------------------------------
# Internal helper: build p_vec from scalar or vector input
# ----------------------------------------------------------------

.resolve_p_vec <- function(p_vec, M) {
  if (length(p_vec) == 1L)
    p_vec <- rep(as.integer(p_vec), M)
  stopifnot(length(p_vec) == M, all(p_vec >= 1L))
  as.integer(p_vec)
}


# ----------------------------------------------------------------
# E-step  (mixed-lag)
# ----------------------------------------------------------------

#' E-step for the mixed-lag MSM-VAR
#'
#' Identical in structure to \code{ms_var_e_step} but accepts a
#' per-regime lag vector \code{p_vec}.  The augmented state space is
#' built with \code{p_max = max(p_vec)} lags; regime \eqn{s} only
#' populates regressors up to its own \code{p_vec[s]} when evaluating
#' the conditional mean.
#'
#' @param Y           Data matrix (\eqn{T \times k}).
#' @param P_hat       Current \eqn{M \times M} transition matrix.
#' @param mu          List of \eqn{M} mean vectors.
#' @param A_list      List of \eqn{M} AR coefficient lists.  Element
#'   \code{A_list[[s]]} has length \code{p_vec[s]}.
#' @param Sigma_list  List of \eqn{M} covariance matrices.
#' @param p_vec       Integer vector of length \eqn{M} giving the lag order
#'   for each regime.  A scalar is broadcast to all regimes (reducing to the
#'   standard equal-lag case).
#' @param switching   Named logical list from \code{make_switching()}.
#' @param pi          Initial regime probability vector.  Default: uniform.
#' @param eps         Numerical stability constant.  Default \code{1e-12}.
#'
#' @return Same structure as \code{ms_var_e_step}: a list with
#'   \code{loglik}, \code{gammaZ}, \code{smoothed_marginals},
#'   \code{smoothed_joint}, and additionally \code{p_max}.
#'
#' @export
ms_var_e_step_ml <- function(Y, P_hat, mu, A_list, Sigma_list,
                                p_vec,
                                switching = make_switching(),
                                pi = NULL, eps = 1e-12) {

  Tn    <- nrow(Y)
  k     <- ncol(Y)
  M     <- nrow(P_hat)
  p_vec <- .resolve_p_vec(p_vec, M)
  p_max <- max(p_vec)

  if (is.null(pi)) pi <- rep(1 / M, M)
  pi <- pi / sum(pi)

  # Total augmented states: M * (M+1)^p_max
  K <- M * (M + 1L)^p_max

  # ---- Observation log-likelihoods for every augmented state ----
  logg <- matrix(-Inf, nrow = K, ncol = Tn)

  for (t in 1:Tn) {
    q_max <- min(p_max, t - 1L)

    for (z in 1:K) {
      zt  <- decode_z(z, M, p_max)   # (s_t, s_{t-1}, ..., s_{t-p_max})
      st  <- zt[1]
      p_s <- p_vec[st]               # lag order for this regime

      mu_idx    <- if (switching$mu)    st else 1L
      A_idx     <- if (switching$A)     st else 1L
      Sigma_idx <- if (switching$Sigma) st else 1L

      mean_t <- mu[[mu_idx]]

      # Only sum lags up to min(p_s, t-1)
      q_s <- min(p_s, t - 1L)
      if (q_s > 0L) {
        for (lag in 1:q_s) {
          s_lag  <- zt[lag + 1L]   # from augmented state history
          if (s_lag == 0L) next
          mu_lag <- if (switching$mu) mu[[s_lag]] else mu[[1L]]
          mean_t <- mean_t +
            A_list[[A_idx]][[lag]] %*% (Y[t - lag, ] - mu_lag)
        }
      }

      logg[z, t] <- dmvnorm_log(Y[t, ],
                                mean  = as.numeric(mean_t),
                                Sigma = Sigma_list[[Sigma_idx]],
                                eps   = eps)
    }
  }

  # ---- Forward pass ----
  logalpha <- matrix(-Inf, nrow = Tn, ncol = K)
  logbeta  <- matrix(-Inf, nrow = Tn, ncol = K)
  csc      <- numeric(Tn)

  # Initial distribution: only augmented states with all-zero lag history
  logpi_z1 <- rep(-Inf, K)
  for (j in 1:M) {
    z1 <- encode_z(c(j, rep(0L, p_max)), M, p_max)
    logpi_z1[z1] <- log(pi[j] + eps)
  }

  la1    <- logpi_z1 + logg[, 1]
  csc[1] <- logsumexp(la1)
  logalpha[1, ] <- la1 - csc[1]

  for (t in 2:Tn) {
    logalpha[t, ] <- -Inf

    for (zprev in 1:K) {
      ap <- logalpha[t - 1L, zprev]
      if (!is.finite(ap)) next

      zprev_vec <- decode_z(zprev, M, p_max)
      i         <- zprev_vec[1]

      for (j in 1:M) {
        znew_vec  <- c(j, zprev_vec[1:p_max])
        znew      <- encode_z(znew_vec, M, p_max)
        logalpha[t, znew] <- logsumexp(
          c(logalpha[t, znew], ap + log(P_hat[i, j] + eps))
        )
      }
    }

    logalpha[t, ] <- logalpha[t, ] + logg[, t]
    csc[t]        <- logsumexp(logalpha[t, ])
    logalpha[t, ] <- logalpha[t, ] - csc[t]
  }

  loglik <- sum(csc)

  # ---- Backward pass ----
  logbeta[Tn, ] <- 0

  for (t in (Tn - 1L):1L) {
    for (zprev in 1:K) {
      zprev_vec <- decode_z(zprev, M, p_max)
      i         <- zprev_vec[1]
      acc       <- -Inf

      for (j in 1:M) {
        znew_vec <- c(j, zprev_vec[1:p_max])
        znew     <- encode_z(znew_vec, M, p_max)
        acc <- logsumexp(c(acc,
                           log(P_hat[i, j] + eps) +
                             logg[znew, t + 1L]   +
                             logbeta[t + 1L, znew]))
      }
      logbeta[t, zprev] <- acc - csc[t + 1L]
    }
  }

  # ---- Smoothed posteriors ----
  gammaZ <- matrix(0, nrow = Tn, ncol = K)
  for (t in 1:Tn) {
    lg          <- logalpha[t, ] + logbeta[t, ]
    denom       <- logsumexp(lg)
    gammaZ[t, ] <- exp(lg - denom)
  }

  smoothed_marg <- matrix(0, nrow = Tn, ncol = M)
  for (t in 1:Tn) {
    for (z in 1:K) {
      st                   <- decode_z(z, M, p_max)[1]
      smoothed_marg[t, st] <- smoothed_marg[t, st] + gammaZ[t, z]
    }
  }

  smoothed_joint      <- vector("list", Tn)
  smoothed_joint[[1]] <- NULL

  for (t in 2:Tn) {
    J <- matrix(0, M, M)
    for (zprev in 1:K) {
      zprev_vec <- decode_z(zprev, M, p_max)
      i         <- zprev_vec[1]
      ap        <- logalpha[t - 1L, zprev]
      if (!is.finite(ap)) next

      for (j in 1:M) {
        znew_vec <- c(j, zprev_vec[1:p_max])
        znew     <- encode_z(znew_vec, M, p_max)
        J[i, j]  <- J[i, j] +
          exp(ap) * P_hat[i, j] *
          exp(logg[znew, t]) *
          exp(logbeta[t, znew])
      }
    }
    smoothed_joint[[t]] <- J / (sum(J) + eps)
  }

  list(
    loglik             = loglik,
    gammaZ             = gammaZ,
    smoothed_marginals = smoothed_marg,
    smoothed_joint     = smoothed_joint,
    p_max              = p_max
  )
}


# ----------------------------------------------------------------
# M-step  (mixed-lag, pseudoinverse OLS)
# ----------------------------------------------------------------

#' M-step for the mixed-lag MSM-VAR
#'
#' Updates all parameter blocks.  The AR update for regime \eqn{s} builds
#' a regressor matrix of width \eqn{k \cdot p_s}, which may differ across
#' regimes.  The weighted normal equations are solved via the Moore-Penrose
#' pseudoinverse (\code{MASS::ginv}) so that the solve remains well-defined
#' even when the effective rank is less than \eqn{k \cdot p_s} (e.g. when
#' regime \eqn{s} has few effective observations or large \eqn{p_s}).
#' An additional ridge term \code{eps * I} is added before inversion for
#' numerical stability.
#'
#' @param Y             Data matrix (\eqn{T \times k}).
#' @param gammaZ        \eqn{T \times K} smoothed augmented-state posteriors.
#' @param smoothed_joint List of smoothed joint transition posteriors.
#' @param mu_prev       Previous mean list (length \eqn{M}).
#' @param A_prev        Previous AR list (length \eqn{M}).
#' @param p_vec         Integer vector of length \eqn{M} (regime lag orders).
#' @param switching     Named logical list from \code{make_switching()}.
#' @param eps           Ridge / stability constant.  Default \code{1e-10}.
#'
#' @return Same structure as \code{ms_var_m_step}: a list with
#'   \code{mu}, \code{A_list}, \code{Sigma_list}, \code{P_hat}.
#'
#' @export
ms_var_m_step_ml <- function(Y, gammaZ, smoothed_joint,
                                mu_prev, A_prev,
                                p_vec,
                                switching = make_switching(),
                                eps = 1e-10) {

  Tn    <- nrow(Y)
  k     <- ncol(Y)
  M     <- length(mu_prev)
  p_vec <- .resolve_p_vec(p_vec, M)
  p_max <- max(p_vec)
  K     <- ncol(gammaZ)
  stopifnot(K == M * (M + 1L)^p_max)

  A_indices     <- if (switching$A)     1:M else 1L
  mu_indices    <- if (switching$mu)    1:M else 1L
  Sigma_indices <- if (switching$Sigma) 1:M else 1L

  target_regimes <- function(j_target, flag)
    if (flag) j_target else 1:M

  # ---- (A) AR matrices via pseudoinverse weighted OLS ----
  # For each target index j_target, build Xvec of width k * p_{j_target}.
  # If non-switching, p_s = p_vec[1] (we use slot 1; see broadcast below).

  A_raw <- vector("list", length(A_indices))
  names(A_raw) <- as.character(A_indices)

  for (jj in seq_along(A_indices)) {
    j_target <- A_indices[jj]
    valid_st <- target_regimes(j_target, switching$A)

    # Lag width for this regime target
    p_j  <- if (switching$A) p_vec[j_target] else max(p_vec[valid_st])
    dim_x <- k * p_j

    XtWX <- matrix(0, nrow = dim_x, ncol = dim_x)
    XtWY <- matrix(0, nrow = dim_x, ncol = k)

    for (t in 1:Tn) {
      q <- min(p_j, t - 1L)

      for (z in 1:K) {
        w  <- gammaZ[t, z]
        if (w <= 0) next
        zt <- decode_z(z, M, p_max)
        if (!(zt[1] %in% valid_st)) next

        # Build regressor vector of length k * p_j
        Xvec <- matrix(0, nrow = dim_x, ncol = 1L)
        if (q > 0L) {
          for (lag in 1:q) {
            s_lag <- zt[lag + 1L]
            if (s_lag == 0L) next
            mu_lag <- if (switching$mu) mu_prev[[s_lag]] else mu_prev[[1L]]
            block  <- (1L + (lag - 1L) * k):(lag * k)
            Xvec[block, 1L] <- Y[t - lag, ] - mu_lag
          }
        }

        mu_st <- if (switching$mu) mu_prev[[zt[1]]] else mu_prev[[1L]]
        ytil  <- matrix(Y[t, ] - mu_st, ncol = 1L)

        XtWX <- XtWX + w * (Xvec %*% t(Xvec))
        XtWY <- XtWY + w * (Xvec %*% t(ytil))
      }
    }

    # Pseudoinverse solve with ridge for stability
    XtWX_reg <- XtWX + eps * diag(dim_x)
    B         <- MASS::ginv(XtWX_reg) %*% XtWY   # dim_x x k

    Aj_list <- vector("list", p_j)
    for (lag in 1:p_j) {
      rows          <- (1L + (lag - 1L) * k):(lag * k)
      Aj_list[[lag]] <- t(B[rows, , drop = FALSE])  # k x k
    }
    A_raw[[jj]] <- Aj_list
  }

  # Broadcast non-switching: replicate single estimate to all M slots,
  # padding shorter lag lists to p_max if needed for state-space consistency
  if (switching$A) {
    A_new <- A_raw
  } else {
    A_new <- rep(list(A_raw[[1]]), M)
  }

  # ---- (B) Means ----
  mu_raw <- vector("list", length(mu_indices))

  for (jj in seq_along(mu_indices)) {
    j_target <- mu_indices[jj]
    valid_st <- target_regimes(j_target, switching$mu)
    p_j      <- if (switching$A) p_vec[j_target] else max(p_vec[valid_st])

    num <- rep(0, k)
    den <- 0

    for (t in 1:Tn) {
      q <- min(p_j, t - 1L)

      for (z in 1:K) {
        w  <- gammaZ[t, z]
        if (w <= 0) next
        zt <- decode_z(z, M, p_max)
        if (!(zt[1] %in% valid_st)) next

        A_use <- A_new[[if (switching$A) zt[1] else 1L]]

        pred_no_mu <- rep(0, k)
        if (q > 0L) {
          for (lag in 1:min(q, length(A_use))) {
            s_lag  <- zt[lag + 1L]
            if (s_lag == 0L) next
            mu_lag     <- if (switching$mu) mu_prev[[s_lag]] else mu_prev[[1L]]
            pred_no_mu <- pred_no_mu +
              A_use[[lag]] %*% (Y[t - lag, ] - mu_lag)
          }
        }

        num <- num + w * (Y[t, ] - as.numeric(pred_no_mu))
        den <- den + w
      }
    }
    mu_raw[[jj]] <- as.numeric(num / (den + eps))
  }

  mu_new <- if (switching$mu) mu_raw else rep(list(mu_raw[[1]]), M)

  # ---- (C) Covariance matrices ----
  Sigma_raw <- vector("list", length(Sigma_indices))

  for (jj in seq_along(Sigma_indices)) {
    j_target <- Sigma_indices[jj]
    valid_st <- target_regimes(j_target, switching$Sigma)
    p_j      <- if (switching$A) p_vec[j_target] else max(p_vec[valid_st])

    S   <- matrix(0, k, k)
    den <- 0

    for (t in 1:Tn) {
      q <- min(p_j, t - 1L)

      for (z in 1:K) {
        w  <- gammaZ[t, z]
        if (w <= 0) next
        zt <- decode_z(z, M, p_max)
        if (!(zt[1] %in% valid_st)) next

        st     <- zt[1]
        mu_use <- mu_new[[if (switching$mu) st else 1L]]
        A_use  <- A_new[[if (switching$A)   st else 1L]]

        mean_t <- mu_use
        if (q > 0L) {
          for (lag in 1:min(q, length(A_use))) {
            s_lag  <- zt[lag + 1L]
            if (s_lag == 0L) next
            mu_lag <- mu_new[[if (switching$mu) s_lag else 1L]]
            mean_t <- mean_t + A_use[[lag]] %*% (Y[t - lag, ] - mu_lag)
          }
        }

        e   <- as.numeric(Y[t, ] - as.numeric(mean_t))
        S   <- S   + w * tcrossprod(e)
        den <- den + w
      }
    }

    S_j       <- S / (den + eps)
    diag(S_j) <- diag(S_j) + eps
    Sigma_raw[[jj]] <- S_j
  }

  Sigma_new <- if (switching$Sigma) Sigma_raw else rep(list(Sigma_raw[[1]]), M)

  # ---- (D) Transition matrix ----
  P_new <- matrix(0, M, M)
  for (i in 1:M) {
    for (j in 1:M)
      P_new[i, j] <- sum(sapply(2:Tn, function(t) smoothed_joint[[t]][i, j]))
    P_new[i, ] <- P_new[i, ] / (sum(P_new[i, ]) + eps)
  }

  list(mu = mu_new, A_list = A_new, Sigma_list = Sigma_new, P_hat = P_new)
}


# ----------------------------------------------------------------
# Full EM loop  (mixed-lag)
# ----------------------------------------------------------------

#' EM estimation for the mixed-lag MSM-VAR
#'
#' Fits a Markov-Switching Mean-Adjusted Heteroskedastic VAR with
#' per-regime lag orders \code{p_vec}.  The conditional mean for regime
#' \eqn{s} is:
#' \deqn{
#'   y_t - \mu_s = \sum_{l=1}^{p_s} A_{s,l}(y_{t-l} - \mu_{s_{t-l}})
#'                + \varepsilon_t, \quad \varepsilon_t \sim N(0, \Sigma_s).
#' }
#' Initialization lists \code{A_init[[s]]} must have length \code{p_vec[s]}.
#'
#' @param Y           Data matrix (\eqn{T \times k}).
#' @param p_vec       Integer vector of length \eqn{M} with per-regime lag
#'   orders.  A scalar is broadcast to all regimes (standard equal-lag case).
#' @param mu_init     List of \eqn{M} initial mean vectors.
#' @param A_init      List of \eqn{M} initial AR lists.  \code{A_init[[s]]}
#'   must have length \code{p_vec[s]}.
#' @param Sigma_init  List of \eqn{M} initial covariance matrices.
#' @param P_init      Initial \eqn{M \times M} transition matrix.
#' @param switching   Named logical list from \code{make_switching()}.
#' @param pi          Initial regime probability vector.  Default: uniform.
#' @param max_iter    Maximum EM iterations.  Default \code{300}.
#' @param tol         Log-likelihood convergence tolerance.  Default \code{1e-6}.
#' @param eps         Ridge / stability constant.  Default \code{1e-10}.
#' @param step_size   EM damping factor in \eqn{(0, 1]}.  Default \code{0.5}.
#' @param verbose     Print progress every 10 iterations.  Default \code{TRUE}.
#'
#' @return A list with \code{mu}, \code{A_list}, \code{Sigma_list},
#'   \code{P_hat}, \code{loglik}, \code{p_vec}, \code{switching}.
#'
#' @examples
#' \dontrun{
#' # Regime 1 uses 1 lag, Regime 2 uses 2 lags
#' A_init <- list(
#'   list(diag(0.3, k)),            # p_1 = 1
#'   list(diag(0.3, k), diag(0.1, k))  # p_2 = 2
#' )
#' fit <- ms_var_em_ml(Y, p_vec = c(1, 2),
#'                        mu_init, A_init, Sigma_init, P_init)
#' }
#'
#' @export
ms_var_em_ml <- function(Y, p_vec,
                            mu_init, A_init, Sigma_init, P_init,
                            switching  = make_switching(),
                            pi         = NULL,
                            max_iter   = 300,
                            tol        = 1e-6,
                            eps        = 1e-10,
                            step_size  = 0.5,
                            verbose    = TRUE) {

  M     <- nrow(P_init)
  p_vec <- .resolve_p_vec(p_vec, M)

  if (is.null(pi)) pi <- rep(1 / M, M)
  switching <- modifyList(make_switching(), switching)

  # Validate A_init lengths
  for (s in 1:M) {
    if (length(A_init[[s]]) != p_vec[s])
      stop(sprintf(
        "A_init[[%d]] has length %d but p_vec[%d] = %d.",
        s, length(A_init[[s]]), s, p_vec[s]
      ))
  }

  mu_curr <- mu_init
  A_curr  <- A_init
  S_curr  <- Sigma_init
  P_curr  <- P_init
  ll      <- numeric(max_iter)

  for (iter in 1:max_iter) {

    E <- ms_var_e_step_ml(Y, P_curr, mu_curr, A_curr, S_curr,
                             p_vec     = p_vec,
                             switching = switching,
                             pi        = pi,
                             eps       = eps)
    ll[iter] <- E$loglik

    Mres <- ms_var_m_step_ml(Y, E$gammaZ, E$smoothed_joint,
                                mu_curr, A_curr,
                                p_vec     = p_vec,
                                switching = switching,
                                eps       = eps)

    # Damped update
    mu_next <- lapply(seq_len(M), function(j)
      mu_curr[[j]] + step_size * (Mres$mu[[j]] - mu_curr[[j]]))

    # AR: only update lags that exist for each regime
    A_next <- lapply(seq_len(M), function(j) {
      p_j <- length(A_curr[[j]])
      lapply(seq_len(p_j), function(lag)
        A_curr[[j]][[lag]] + step_size *
          (Mres$A_list[[j]][[lag]] - A_curr[[j]][[lag]]))
    })

    S_next <- lapply(seq_len(M), function(j)
      S_curr[[j]] + step_size * (Mres$Sigma_list[[j]] - S_curr[[j]]))

    P_next <- P_curr + step_size * (Mres$P_hat - P_curr)
    P_next <- P_next / rowSums(P_next)

    mu_curr <- mu_next
    A_curr  <- A_next
    S_curr  <- S_next
    P_curr  <- P_next

    if (verbose && iter %% 10 == 0)
      cat(sprintf("iter: %d  loglik: %.4f\n", iter, ll[iter]))

    if (iter > 1L && is.finite(ll[iter]) &&
        abs(ll[iter] - ll[iter - 1L]) < tol) {
      ll <- ll[1:iter]
      if (verbose)
        cat(sprintf("Converged at iteration %d  loglik: %.4f\n",
                    iter, ll[iter]))
      break
    }
  }

  list(
    mu         = mu_curr,
    A_list     = A_curr,
    Sigma_list = S_curr,
    P_hat      = P_curr,
    loglik     = ll,
    p_vec      = p_vec,
    switching  = switching
  )
}


# ----------------------------------------------------------------
# BIC grid search over (p_1, ..., p_M)
# ----------------------------------------------------------------

#' BIC grid search for mixed-lag MSM-VAR
#'
#' Fits the mixed-lag model over all combinations of per-regime lag orders
#' drawn from \code{p_grid} and returns a ranked summary table.  The
#' initialization for each grid point uses the parameters from the supplied
#' standard-lag MS-VAR fit (\code{base_fit}), with shorter AR lists
#' truncated and longer ones padded with zero matrices.
#'
#' The BIC penalty is:
#' \deqn{
#'   \mathrm{BIC} = -2\hat{\ell} + \log(T) \cdot \hat{d},
#' }
#' where \eqn{\hat{d}} is the number of freely estimated parameters:
#' \eqn{M(M-1)} transition parameters, \eqn{Mk} means (if switching),
#' \eqn{M \cdot k^2 p_s} AR coefficients (summed over regimes), and
#' \eqn{M \cdot k(k+1)/2} covariance entries (if switching).
#'
#' @param Y         Data matrix (\eqn{T \times k}).
#' @param M         Number of regimes.
#' @param p_grid    Integer vector of candidate lag orders to consider per
#'   regime.  Default \code{1:4}.  All \code{M}-tuples are evaluated.
#' @param base_fit  A fitted object from \code{ms_var_em_ml} or
#'   \code{ms_var_em} (the package's standard EM), used for
#'   warm-start initialization.
#' @param switching Named logical list from \code{make_switching()}.
#' @param pi        Initial regime probability vector.  Default: uniform.
#' @param max_iter  EM iterations per grid point.  Default \code{200}.
#' @param tol       Convergence tolerance.  Default \code{1e-6}.
#' @param eps       Ridge constant.  Default \code{1e-10}.
#' @param step_size EM damping.  Default \code{0.5}.
#' @param verbose   Print progress for each grid point.  Default \code{FALSE}.
#'
#' @return A list with:
#' \describe{
#'   \item{results}{Data frame of grid results sorted by BIC, with columns
#'     \code{p_vec} (as string), \code{loglik}, \code{n_params}, \code{BIC},
#'     \code{converged}, \code{iters}.}
#'   \item{best_fit}{The full fitted object for the BIC-optimal grid point.}
#'   \item{best_p_vec}{Integer vector of the BIC-optimal lag orders.}
#'   \item{all_fits}{Named list of all fitted objects (key = p_vec string).}
#' }
#'
#' @export
ms_var_bic_grid_ml <- function(Y, M, p_grid = 1:4,
                                  base_fit,
                                  switching  = make_switching(),
                                  pi         = NULL,
                                  max_iter   = 200,
                                  tol        = 1e-6,
                                  eps        = 1e-10,
                                  step_size  = 0.5,
                                  verbose    = FALSE) {

  Tn <- nrow(Y)
  k  <- ncol(Y)

  # All M-tuples from p_grid
  grid_list  <- do.call(expand.grid, rep(list(p_grid), M))
  n_grid     <- nrow(grid_list)

  cat(sprintf("BIC grid search: %d combinations over %d regimes\n",
              n_grid, M))

  # Extract base initialization
  base_mu    <- base_fit$mu
  base_Sigma <- base_fit$Sigma_list
  base_P     <- base_fit$P_hat

  # Helper: pad/truncate AR list to target length p_target
  .resize_A <- function(A_s, p_target, k) {
    p_s <- length(A_s)
    if (p_target == p_s) return(A_s)
    if (p_target < p_s)  return(A_s[1:p_target])
    # Pad with zero matrices
    c(A_s, rep(list(matrix(0, k, k)), p_target - p_s))
  }

  # BIC parameter count
  .n_params <- function(p_vec, M, k, switching) {
    n  <- M * (M - 1L)                                        # transition
    n  <- n + if (switching$mu)    M * k           else k     # means
    n  <- n + sum(k^2 * p_vec)                                # AR (per regime)
    n  <- n + if (switching$Sigma) M * k*(k+1)/2  else k*(k+1)/2  # Sigma
    as.integer(n)
  }

  results  <- vector("list", n_grid)
  all_fits <- vector("list", n_grid)

  for (gi in seq_len(n_grid)) {
    p_vec_i <- as.integer(grid_list[gi, ])
    key     <- paste(p_vec_i, collapse = "-")

    if (verbose)
      cat(sprintf("  [%d/%d] p_vec = (%s) ... ",
                  gi, n_grid, paste(p_vec_i, collapse = ", ")))

    # Build A_init for this p_vec
    A_init_i <- lapply(seq_len(M), function(s)
      .resize_A(base_fit$A_list[[s]], p_vec_i[s], k))

    fit_i <- tryCatch(
      ms_var_em_ml(
        Y          = Y,
        p_vec      = p_vec_i,
        mu_init    = base_mu,
        A_init     = A_init_i,
        Sigma_init = base_Sigma,
        P_init     = base_P,
        switching  = switching,
        pi         = pi,
        max_iter   = max_iter,
        tol        = tol,
        eps        = eps,
        step_size  = step_size,
        verbose    = FALSE
      ),
      error = function(e) {
        warning(sprintf("Grid point (%s) failed: %s", key, conditionMessage(e)))
        NULL
      }
    )

    if (is.null(fit_i)) {
      results[[gi]] <- data.frame(
        p_vec     = key,
        loglik    = NA_real_,
        n_params  = NA_integer_,
        BIC       = NA_real_,
        converged = FALSE,
        iters     = NA_integer_,
        stringsAsFactors = FALSE
      )
    } else {
      ll_i      <- tail(fit_i$loglik, 1)
      np_i      <- .n_params(p_vec_i, M, k, switching)
      bic_i     <- -2 * ll_i + log(Tn) * np_i
      conv_i    <- length(fit_i$loglik) < max_iter
      iters_i   <- length(fit_i$loglik)

      if (verbose)
        cat(sprintf("loglik = %.2f  BIC = %.2f  iters = %d\n",
                    ll_i, bic_i, iters_i))

      results[[gi]] <- data.frame(
        p_vec     = key,
        loglik    = ll_i,
        n_params  = np_i,
        BIC       = bic_i,
        converged = conv_i,
        iters     = iters_i,
        stringsAsFactors = FALSE
      )
      all_fits[[gi]] <- fit_i
    }
    names(all_fits)[gi] <- key
  }

  results_df <- do.call(rbind, results)
  results_df <- results_df[order(results_df$BIC, na.last = TRUE), ]

  best_key   <- results_df$p_vec[1]
  best_idx   <- which(names(all_fits) == best_key)[1]
  best_fit   <- all_fits[[best_idx]]
  best_p_vec <- as.integer(strsplit(best_key, "-")[[1]])

  cat(sprintf("\nBIC-optimal p_vec: (%s)  BIC = %.2f\n",
              paste(best_p_vec, collapse = ", "),
              results_df$BIC[1]))

  list(
    results    = results_df,
    best_fit   = best_fit,
    best_p_vec = best_p_vec,
    all_fits   = all_fits
  )
}
