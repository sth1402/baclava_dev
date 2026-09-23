#' A single MH update specification
#' @noRd
#' @keywords internal
Adapt <- R6::R6Class(
  "Adapt",
  public = list(
    #' @field delta ('numeric(1)')\cr
    #' Target acceptance rate
    delta = NULL,
    #' @field t0 ('numeric(1)')\cr
    #' stability; larger values essentially delay adaptation
    t0 = NULL,
    #' @field kappa ('numeric(1)')\cr
    #' smoothing; larger values are less smoothing
    kappa = NULL,
    #' @field warmup ('integer(1)')\cr
    #' Number of warmup iterations
    warmup = NULL,
    #' @field warmup_cov ('integer(1)')\cr
    #' adaptive procedure is "approx D/L"; freeze at warmup_cov; adapt s
    warmup_cov = NULL,
    #' @field learn ('character(1)')\cr
    #' Must be one of 's', 'sD', 'sDL', 'none'
    learn = NULL,
    #' @field batch ('integer(1)')\cr
    #' smoothing (update every batch size)
    batch = NULL,
    #' @field s_floor ('numeric(1)')\cr
    #' minimum allowed s
    s_floor = NULL,
    #' @field s_ceil ('numeric(1)')\cr
    #' maximum allowed s
    s_ceil = NULL,
    #' @field rm_gamma (`numeric(1)`)\cr
    #' Robbins-Monro lambda = c / n^{gamma}
    rm_gamma = NULL,
    #' @field rm_gamma (`numeric(1)`)\cr
    #' Robbins-Monro lambda = c / n^{gamma}
    rm_c = NULL,
    #' @field ridge (`numeric(1)`)\cr
    #' ridge penalty
    ridge = NULL,
    #' @field mode (`character(1)`)\cr
    #' fast/slow adaptation
    mode = NULL,
    #' @field mu (`numeric(1)`)\cr
    #' fast mode anchor point
    mu = NULL,
    #' @field gamma (`numeric(1)`)\cr
    #' fast mode gain
    gamma = NULL,

    continue = function() {
      self$warmup <- 0L
    },
    
    #' @description
    #' Defines how adaptive procedure
    #'
    #' @param delta Numeric scalar in (0, 1). Target acceptance rate.
    #' @param warmup Integer scalar. Number of warmup iterations
    #' @param learn Character. Components of proposal to learn. Must be one of
    #'   's' (step size only), 'sD' (step size and shape), 'sDL' (step size, 
    #'   shape, and correlation), or 'none' (no adaptive procedure)
    #' @param batch Integer scalar. Smooth over batch iterations.
    #' @param t0 Positive numeric scalar. Stability; larger values essentially 
    #'   delay adaptation.
    #' @param kappa Numeric scalar in (0,1). Smoothing; larger values are less 
    #'   smoothing.
    #' @param s_floor Positive numeric scalar. Minimum allowed s.
    #' @param s_ceil Positive numeric scalar. Maximum allowed s.
    #' @param rm_gamma Positive numeric scalar. Robbins-Monro lambda =c / n^gamma.
    #' @param rm_c Positive numeric scalar. Robbins-Monro lambda =c / n^gamma.
    #' @param ridge Positive numeric scalar. Ridge penalty in refresh.
    #' @param mode Character string. "fast" uses dual averaging (Nesterov) for 
    #'   step-size adaptation; "slow" uses standard Robbins-Monro.
    #' @param mu Numeric scalar. The anchor point for the "fast" mode -- by
    #'   default set as log(10*s0)
    #' @param gamma Numeric scalar. The gain term controlling aggressiveness
    #'   of the "fast" mode.
    initialize = function(delta = NA_real_, warmup = 0L, learn = 'none',
                          warmup_cov = 0L,
                          batch = 1L, t0 = 10.0, kappa = 0.75,
                          s_floor = 1e-4, s_ceil = 2.5,
                          rm_gamma = 0.7, rm_c = 1.0,
                          ridge = 1e-6, 
                          mode = "fast", mu = NA_real_, gamma = 0.2) {
      
      stopifnot({is.vector(delta, "numeric") && length(delta) == 1L && delta > 0.0} || 
                  all(is.na(delta)))
      stopifnot(is.vector(warmup, "numeric"), length(warmup) == 1L, warmup >= 0)
      stopifnot(is.vector(batch, "numeric") && length(batch) == 1L && batch > 0.0)
      stopifnot(is.vector(t0, "numeric") && length(t0) == 1L && t0 > 0.0)
      stopifnot(is.vector(kappa, "numeric") && length(kappa) == 1L && kappa > 0.0)
      stopifnot(is.vector(s_floor, "numeric") && length(s_floor) == 1L && s_floor > 0.0)
      stopifnot(is.vector(s_ceil, "numeric") && length(s_ceil) == 1L && s_ceil > 0.0)
      stopifnot(s_ceil >= s_floor)
      stopifnot(is.vector(ridge, "numeric"), length(ridge) == 1L, ridge >= 0)
      stopifnot(is.vector(rm_gamma, "numeric"),
                length(rm_gamma) == 1L,
                rm_gamma > 0.5)
      stopifnot(is.vector(rm_c, "numeric"),
                length(rm_c) == 1L,
                rm_c > 0.0)
      
      self$delta <- delta
      self$warmup <- as.integer(warmup)
      self$warmup_cov <- as.integer(warmup_cov)
      self$learn <- learn
      if (!{learn %in% c('s', 'sD', 'sDL', 'none')}) stop("learn misspecified")
      self$batch <- batch
      self$t0 <- t0
      self$kappa <- kappa
      self$s_floor <- s_floor
      self$s_ceil <- s_ceil
      
      self$rm_gamma <- rm_gamma
      self$rm_c <- rm_c
      self$ridge <- as.numeric(ridge)
      self$mode <- mode
      self$mu <- mu
      self$gamma <- gamma

    },
    
    #' Stable shape for Rcpp
    toList = function() {
      list(
        "meta_for" = "Adapt",
        "delta" = self$delta,
        "warmup" = self$warmup,
        "learn" = self$learn,
        "warmup_cov" = self$warmup_cov,
        "batch" = self$batch,
        "t0" = self$t0,
        "kappa" = self$kappa,
        "s_floor" = self$s_floor,
        "s_ceil" = self$s_ceil,
        "rm_gamma" = self$rm_gamma,
        "rm_c" = self$rm_c,
        "ridge" = self$ridge,
        "mode" = self$mode,
        "mu" = self$mu,
        "gamma" = self$gamma
      )
    },
    
    print = function(...) {
      cat("<Adapt>\n")
      cat("mode=", self$mode,
          " delta=", self$delta,
          " t0=", self$t0,
          " kappa=", self$kappa,
          " mu=", self$mu,
          " gamma=",  self$gamma,
          " warmup=", self$warmup, 
          " learn=", self$learn,
          " covstart=", self$warmup_cov,
          " batch=", self$batch,
          " s_floor=", self$s_floor,
          " s_ceil=", self$s_ceil,
          " rm_gamma:", self$rm_gamma,
          " rm_c:", self$rm_c,
          " ridge:", self$ridge, "\n")
      invisible(self)
    }
  )
)

#' Define a Adaptive Procedure
#'
#' Convenience constructor for a \code{Adapt} (R6) describing adaptation settings.
#'
#' @param delta Numeric scalar in (0, 1). Target acceptance rate.
#' @param warmup Integer scalar. Number of warmup iterations
#' @param warmup_cov Integer scalar. Adaptive procedure is approximate D/L;
#'   freeze D/L at warmup_cov; adapt s (only used for learn = 'sD' or learn = 'sDL')
#' @param learn Character. Components of proposal to learn. Must be one of
#'   's' (step size only), 'sD' (step size and shape), 'sDL' (step size, 
#'   shape, and correlation), or 'none' (no adaptive procedure)
#' @param batch Integer scalar in (0,1). Smooth over batch iterations
#' @param t0 Positive numeric scalar. Stability; larger values essentially 
#'   delay adaptation.
#' @param kappa Numeric scalar in (0,1). Smoothing; larger values are less 
#'   smoothing.
#' @param s_floor Positive numeric scalar. Minimum allowed s.
#' @param s_ceil Positive numeric scalar. Maximum allowed s.
#' @param rm_gamma Positive numeric scalar > 0.5. Robbins-Monro lambda = c/ n^{gamma}
#' @param rm_c Positive numeric scalar. Robbins-Monro lambda = c/ n^{gamma}
#' @param ridge Positive numeric scalar. Ridge penalty in refresh.
#' @param mu Numeric scalar. The anchor point for the "fast" mode -- by
#'   default set as log(10*s0)
#' @param gamma Numeric scalar. The gain term controlling aggressiveness
#'   of the "fast" mode.
#'
#' @return An R6 \code{Adapt} object (class \code{"Adapt"}) with a stable
#'   \code{$toList()} representation suitable for downstream C++/Rcpp code.
#'
#' @examples
#'  adaptive = adapt(delta = 0.23, warmup = 2000L, batch = 5.0, t0 = 50, kappa = 0.75,
#'                  s_floor = 1e-3, s_ceil = 3.0, learn = 's')
#'
#' @export
adapt <- function(delta = NA_real_, warmup = 0L, learn = 'none', 
                  warmup_cov = -1L, 
                  batch = 1L, t0 = 10.0, kappa = 0.75,
                  s_floor = 1e-4, s_ceil = 2.5,
                  rm_gamma = 0.7, rm_c = 1.0,
                  ridge = 1e-6, 
                  mode = "fast", mu = NA_real_, gamma = 0.2) {
  
  return(Adapt$new(delta, warmup, learn, warmup_cov, batch, t0, kappa, 
                   s_floor, s_ceil, rm_gamma, rm_c, ridge, 
                   mode, mu, gamma))
}

#' @describeIn adapt Inheritance test
is.Adapt <- function(x) inherits(x, "Adapt")