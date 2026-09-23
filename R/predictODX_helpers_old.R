### Calculating over-diagnosis contributions

#' Integral f(t; shape.H, scale.H) S(Y-t; shape.P, scale.P) from L:U
#' 
#' @noRd
#' @param U A scalar numeric. The upper integration boundary.
#' @param L A scalar numeric. The lower integration boundary.
#' @param last A scalar numeric. The final screening age.
#' @param shape.H A scalar numeric. The shape parameter for the Weibull of the
#'   healthy compartment.
#' @param scale.H A scalar numeric. The scale parameter for the Weibull of the
#'   healthy compartment.
#' @param shape.P A scalar numeric. The shape parameter for the Weibull of the
#'   pre-clinical compartment.
#' @param scale.P A scalar numeric. The scale parameter for the Weibull of the
#'   pre-clinical compartment.
#'
#' @return A scalar numeric.
#'
#' @importFrom stats dweibull integrate pweibull
#' @keywords internal
.integral <- function(U, L, last, shape.H, scale.H, shape.P, scale.P) {

  integrand <- function(t, last, shape.H, shape.P, scale.H, scale.P){
    stats::dweibull(t, shape.H, scale.H) * 
      stats::pweibull(last - t, shape.P, scale.P, lower.tail = FALSE)
  }
  
  stats::integrate(integrand, lower = L, upper = U,
                   last = last, 
                   shape.H = shape.H, shape.P = shape.P, 
                   scale.H = scale.H, scale.P = scale.P, 
                   stop.on.error = TRUE)$value
  
}

#' Contribution to the normalization for non-indolent
#' 
#' @noRd
#' @param U A scalar numeric. The upper integration boundary.
#' @param L A scalar numeric. The lower integration boundary.
#' @param shape.H A scalar numeric. The shape parameter for the Weibull of the
#'   healthy compartment.
#' @param scale.H A scalar numeric. The scale parameter for the Weibull of the
#'   healthy compartment.
#' @param shape.P A scalar numeric. The shape parameter for the Weibull of the
#'   pre-clinical compartment.
#' @param scale.P A scalar numeric. The scale parameter for the Weibull of the
#'   pre-clinical compartment.
#'
#' @return A scalar numeric.
#'
#' @importFrom stats pweibull
#' @keywords internal
.normalization_not_indolent <- function(U, L,  
                                        shape.H, scale.H, 
                                        shape.P, scale.P) {
  
  stats::pweibull(U, shape.H, scale.H, lower.tail = FALSE) + 
    .integral(U = U, L = L, last = U, 
              shape.H = shape.H, scale.H = scale.H, 
              shape.P = shape.P, scale.P = scale.P)
}

#' Normalization
#' 
#' Accounts for the time lapse between risk onset and first screen.
#' 
#' @noRd
#' @param psi A scalar numeric. The probability of being indolent.
#' @param U A scalar numeric. The upper integration boundary.
#' @param L A scalar numeric. The lower integration boundary.
#' @param last A scalar numeric. The final screening age.
#' @param shape.H A scalar numeric. The shape parameter for the Weibull of the
#'   healthy compartment.
#' @param scale.H A scalar numeric. The scale parameter for the Weibull of the
#'   healthy compartment.
#' @param shape.P A scalar numeric. The shape parameter for the Weibull of the
#'   pre-clinical compartment.
#' @param scale.P A scalar numeric. The scale parameter for the Weibull of the
#'   pre-clinical compartment.
#'
#' @return A scalar numeric.
#'
#' @keywords internal
.normalization <- function(psi, U, L,
                           shape.H, scale.H, shape.P, scale.P) {
  
  psi + (1.0 - psi) * 
    .normalization_not_indolent(U = U, L = L, 
                                shape.H = shape.H, scale.H = scale.H, 
                                shape.P = shape.P, scale.P = scale.P)
}

#' ALL SCREEN-DETECTED CASES likelihood contribution
#'
#' @noRd
#' @param psi A scalar numeric. The probability of indolence.
#' @param beta A scalar numeric. The screening sensitivity.
#' @param shape.H A scalar numeric. The shape parameter for the Weibull of the
#'   healthy compartment.
#' @param scale.H A scalar numeric. The scale parameter for the Weibull of the
#'   healthy compartment.
#' @param shape.P A scalar numeric. The shape parameter for the Weibull of the
#'   pre-clinical compartment.
#' @param scale.P A scalar numeric. The scale parameter for the Weibull of the
#'   pre-clinical compartment.
#' @param time.points A numeric vector object. The age at risk onset and
#'   ages at each screen
#'   (age_risk_onset, age_screen_1, ..., age_screen_t)
#' @param n.screens An integer object. The number of screens.
#' 
#' @returns A list. Probability of an indolent screen-detected cancer (indolent);
#'   probability of a progressive but overdiagnosed cancer (progressive);
#'   and the normalization accounting for period between risk onset and
#'   first screen (norm).
#'
#' @importFrom stats pweibull
#' @keywords internal
D_j <- function(psi, beta, shape.H, scale.H, shape.P, scale.P, time.points, 
                n.screens) {
    
  # the number of time points 
  l <- n.screens + 1L
  
  risk_onset <- time.points[1L]
  age_first_screen <- time.points[2L]
  
  # initializations
  A <- 0.0
  B <- 0.0

  ## sum over previous screens 
  for (k in 2L:l) {
    tmp <- beta * {(1.0 - beta)^{l - k}}
    A <- A +  psi * tmp * 
      {stats::pweibull(time.points[k] - risk_onset, shape.H, scale.H) - 
          stats::pweibull(time.points[k - 1L] - risk_onset, shape.H, scale.H)}
    
    B <- B + (1.0 - psi) * tmp * 
      .integral(U = time.points[k] - risk_onset, 
                L = time.points[k - 1L] - risk_onset,
                last = time.points[l] - risk_onset,
                shape.H = shape.H, scale.H = scale.H,
                shape.P = shape.P, scale.P = scale.P)
  }
  
  norm <- .normalization(psi = psi, 
                         U = age_first_screen - risk_onset, L = 0.0,
                         shape.H = shape.H, scale.H = scale.H, 
                         shape.P = shape.P, scale.P = scale.P)
  
  list("indolent" = A, "progressive" = B, "norm" = 1.0 / norm)
}  

#' other cause mortality survival function 
#' 
#' @noRd
#' @param age.min A scalar numeric. The age of first screening. 
#'   Note that we are assuming that it is an integer value.
#' @param t A scalar numeric. The time at which the survival function is estimated.
#' @param lambda.other A numeric vector. The other cause hazard.
#' @param lambda.ages A numeric vector. The ages at which lambda.other is
#'   provided.
#' 
#' @returns A numeric vector.
#' @keywords internal
surv.fun <- function(age.min, t, lambda.other, lambda.ages) {

  if (t < age.min) return(1.0)
  
  bounds <- findInterval(c(age.min, t), lambda.ages)
  lambdas <- lambda.other[bounds[1L]:bounds[2L]]
  deltat <- lambda.ages[bounds[1L]:bounds[2L]]
  deltat <- c(deltat[-1L], t) - deltat
  hh <- drop(crossprod(lambdas, deltat))

  exp(-hh)
}

surv.fun <- Vectorize(surv.fun, vectorize.args = "t")

########### Prog.odx.programmatic

#' @noRd
#' @param age.min A numeric scalar. Age at first screening.
#' @param screen.age A numeric scalar. The screen age.
#' @param P.lambda A numeric scalar. The modeled risk
#' @param lambda.other A numeric vector. The other cause hazard
#' @param lambda.ages A numeric vector. The ages at which lambda.other is
#'   provided.
#'
#' @returns A numeric scalar.
#' 
#' @importFrom stats dweibull
#' @keywords internal
Prog.odx.program <- function(age.min, screen.age, P.lambda, lambda.other, lambda.ages) {
  
  help.f <- function(t) {
    surv.fun(age.min = age.min, t = screen.age + t, lambda.other = lambda.other, 
             lambda.ages = lambda.ages) * 
      stats::dweibull(x = t, shape = P.lambda[2L], scale = P.lambda[1L])
  }
  
  out <- integrate(help.f, lower = 0.0, upper = 500.0)

  surv.fun(age.min, screen.age, lambda.other, lambda.ages) - out$value
}

Prog.odx.program <- Vectorize(Prog.odx.program, vectorize.args = "screen.age")



#' Summarize Probability Results
#' 
#' @noRd
#' @param ODX.ind.stat A matrix object. The probability of an indolent 
#'   screen-detected cancer while still alive for each screen/model.
#'   {n_screens x M}
#' @param ODX.prog.stat A matrix object. The probability to have a progressive
#'   and overdiagnosed screen-detected cancer while still alive for each screen/model.
#'   {n_screens x M}
#' @param SD.stat A matrix object. The probability of a screen-detected cancer 
#'   while still alive. {n_screens x M}
#'   
#' @returns A matrix providing the average probability for a screen-detected
#'   cancer being 1) total: indolent or progressive and overdiagnosed;
#'   2) indolent; and 3) mortality: progressive and overdiagnosed and their
#'   95% credible intervals
#'   
#' @importFrom stats quantile
#' @keywords internal
wrap <- function(ODX.ind.stat, ODX.prog.stat, SD.stat) {
  
  #overall
  a <- c(round(100.0 * mean((ODX.ind.stat + ODX.prog.stat) / SD.stat), 1),
         round(100.0 * stats::quantile((ODX.ind.stat + ODX.prog.stat) / SD.stat, 
                                       probs = c(0.025, 0.5, 0.975)), 1))
  #indolent
  b <- c(round(100.0 * mean(ODX.ind.stat / SD.stat), 1),
         round(100.0 * stats::quantile(ODX.ind.stat / SD.stat, 
                                       probs = c(0.025, 0.5, 0.975)), 1))
  #mortality
  c <- c(round(100.0 * mean(ODX.prog.stat / SD.stat), 1),
         round(100.0 * stats::quantile(ODX.prog.stat / SD.stat, 
                                       probs = c(0.025, 0.5, 0.975)), 1))
  
  
  tab <- rbind(a, b, c)
  colnames(tab)[1L] <- "mean"
  tab <- tab[, c(1L, 2L, 4L, 3L)]
  rownames(tab) <- c("total", "indolent", "mortality")
  
  tab
}