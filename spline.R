# Spline functions copied from [spline.R] (https://github.com/robjhyndman/forecast/blob/master/R/spline.R)

## Set up Sigma of order (n x n)
make.Sigma <- function(n, n0=0) {
  nn <- n + n0
  Sigma <- matrix(0, nrow = nn, ncol = nn)
  for (i in 1:nn)
    Sigma[i, i:nn] <- Sigma[i:nn, i] <- (i * i * (3 * (i:nn) - i)) / 6
  return(Sigma / (n ^ 3))
}

## Compute spline matrices
spline.matrices <- function(n, beta, cc=1e2, n0=0) {
  nn <- n + n0
  Sigma <- make.Sigma(n, n0)
  s <- cbind(rep(1, nn), (1:nn) / n)
  Omega <- cc * s %*% t(s) + Sigma / beta + diag(nn)
  max.Omega <- max(Omega)
  inv.Omega <- solve(Omega / max.Omega, tol = 1e-10) / max.Omega
  P <- chol(inv.Omega)
  return(list(s = s, Sigma = Sigma, Omega = Omega, inv.Omega = inv.Omega, P = P))
}

## Compute smoothing splines
## Return -loglikelihood
# beta multiplied by 1e6 to avoid numerical difficulties in optimization
spline.loglik <- function(beta, y, cc=1e2) {
  n <- length(y)
  mat <- spline.matrices(n, beta / 1e6, cc = cc)
  y.star <- mat$P %*% matrix(y)
  return(-log(det(mat$P)) + 0.5 * n * log(sum(y.star ^ 2)))
}


#' @importFrom stats sd
train_spline <- function(.data, ...) {

  if (length(measured_vars(.data)) > 1) {
    abort("Only univariate responses are supported by SPLINE.")
  }

  y <- unclass(.data)[[measured_vars(.data)]]

  if (all(is.na(y))) {
    abort("All observations are missing, a model cannot be estimated without data.")
  }

  n <- length(y)
  if (n > 100) { # Use only last 100 observations to get beta
    yy <- y[(n - 99):n]
  } else {
    yy <- y
  }
  
  # Find optimal beta using likelihood approach in Hyndman et al paper.
  beta <- optimize(spline.loglik, interval = c(1e-6, 1e7), y = yy)$minimum / 1e6
  
  # Compute spar which is equivalent to beta
  r <- 256 * smooth.spline(1:n, y, spar = 0)$lambda
  lss <- beta.est * n ^ 3 / (n - 1) ^ 3
  spar <- (log(lss / r) / log(256) + 1) / 3
  splinefit <- smooth.spline(1:n, y, spar = spar)
  sfit <- splinefit$y
  
  # Compute matrices for optimal beta
  mat <- spline.matrices(n, beta)
  
  # Get one-step predictors
  yfit <- e <- rep(NA, n)
  if (n > 1000) {
    warning("Series too long to compute training set fits and residuals")
  } else
  {
    for (i in 1:(n - 1))
    {
      U <- mat$Omega[1:i, i + 1]
      Oinv <- solve(mat$Omega[1:i, 1:i] / 1e6) / 1e6
      yfit[i + 1] <- t(U) %*% Oinv %*% y[1:i]
      sd <- sqrt(mat$Omega[i + 1, i + 1] - t(U) %*% Oinv %*% U)
      e[i + 1] <- (y[i + 1] - yfit[i + 1]) / sd
    }
  }
  
  # Compute sigma^2
  sigma2 <- mean(e ^ 2, na.rm = TRUE)
  
  structure(
    list(
      y = y,
      invOmega = mat$inv.Omega,
      beta = beta,
      fitted = sfit,
      resid = y - yfit,
      sigma2 = sigma2
    ),
    class = "SPLINE"
  )
}


#' Spline model
#'
#' \code{SPLINE()} returns an iid model applied to the formula's response variable.
#'
#' The model does not use any specials, and so everything on the formula's
#' right-hand-side will be ignored.
#'
#' @param formula Model specification.
#' @param ... Not used.
#'
#' @return A model specification.
#'
#' @seealso
#' [Forecasting: Principles and Practices, Some nonlinear regression forecasting methods (section 7.7)](https://otexts.com/fpp3/nonlinear-regression.html)
#'
#' @export
SPLINE <- function(formula, ...) {
  spline_model <- new_model_class("SPLINE",
    train = train_spline
  )
  new_model_definition(spline_model, !!enquo(formula), ...)
}

#' @inherit forecast.ARIMA
#'
#' @export
forecast.SPLINE <- function(object, new_data, specials = NULL, bootstrap = FALSE, times = 5000, ...) {
  h <- NROW(new_data)
  y <- object$y
  n <- length(y)
  invOmega = object$invOmega
  beta = object$beta
  res <- residuals(object)
  sigma2 = object$sigma2
  
  # Compute matrices for optimal beta
  fcmat <- spline.matrices(n, beta, n0 = h) 
  
  # Produce forecasts
  U <- fcmat$Omega[1:n, n + (1:h)]
  Omega0 <- fcmat$Omega[n + (1:h), n + (1:h)]
  fc <- t(U) %*% invOmega %*% y
  se <- sqrt(sigma2 * diag(Omega0 - t(U) %*% invOmega %*% U))
  distributional::dist_normal(fc[,1], se)
  
}

#' @inherit fitted.ARIMA
#'
#' @export
fitted.SPLINE <- function(object, ...) {
  object[["fitted"]]
}

#' @inherit residuals.ARIMA
#'
#' @export
residuals.SPLINE <- function(object, ...) {
  object[["resid"]]
}

#' Glance a spline method model
#'
#' Construct a single row summary of the spline method model.
#'
#' Contains the beta of the fitted spline (`beta`), and the variance of residuals (`sigma2`).
#'
#' @inheritParams generics::glance
#'
#' @return A one row tibble summarising the model's fit.
#'
#' @export
glance.SPLINE <- function(object, ...) {
  tibble(beta = object[["beta"]], sigma2 = object[["sigma2"]])
}

#' Refit a SPLINE model
#'
#' Applies a fitted spline method model to a new dataset.
#'
#' @inheritParams refit.ARIMA
#' @param reestimate If `TRUE`, the spline for the fitted model will be re-estimated 
#' to suit the new data. 
#'
#' @export
refit.SPLINE <- function(object, new_data, specials = NULL, reestimate = TRUE, ...) {
  # Update data for re-evaluation
  return(train_spline(new_data, ...))
}