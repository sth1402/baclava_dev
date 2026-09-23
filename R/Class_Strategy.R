#' A single MH update specification
#' @noRd
#' @keywords internal
#' @include Class_Adapt.R
MHUpdate <- R6::R6Class(
  "MHUpdate",
  public = list(
    #' @field params (`character()`)\cr
    #' Vector of parameter tags included in update step
    params = NULL, 
    #' @field s0 (`numeric(1)`)\cr
    #' Initial guess for update scaling (step = s * (L * z))
    s0 = NULL,
    #' @field adaptive (`list()`)\cr
    #' Adaptive procedure settings Adapt Class object
    adaptive = NULL, 
    #' @field L (`matrix(length(params), length(params))`)
    #' initial guess for covariance matrix
    L = NULL,
    #' @field D (`numeric(length(params))`)
    #' initial guess for shape matrix
    D = NULL,

    continue = function(s, L, D, ...) {
      self$s0 <- s
      self$L <- L
      self$D <- D
      self$adaptive <- Adapt$new(warmup = 0L)
    },
    
    #' @description
    #' Defines how updates for a Compartment Distributions are handled
    #'
    #' @param params A vector of parameter tags included in the MH update
    #' @param s0 An initial guess for s
    #' @param adaptive Adaptive procedure settings (Adapt object)
    #' @param L Initial guess for L 
    #' @param D Initial guess for D 
    initialize = function(params,
                          s0 = 1.0,
                          adaptive = Adapt$new(),
                          L = diag(nrow = length(params), ncol = length(params)),
                          D = rep(1.0, length(params))) {
      
      stopifnot("`params` must be specified" = !missing(params),
                "`params` must be a character vector" = 
                  is.vector(params, "character") && length(params) >= 1L && all(nzchar(params)),
                "`adaptive` must be an Adapt" = is.Adapt(adaptive))
      
      params <- unique(unname(params))
      
      self$params <- params
      
      stopifnot(
        "`s0` must be a positive scalar numeric" = is.vector(s0, "numeric") &&
          length(s0) == 1L && is.finite(s0) && s0 > 0,
        "`L` must be a numeric matrix of dim length(params) x length(params)" = 
          is.matrix(L) && is.numeric(L) && ncol(L) == length(self$params) &&
          nrow(L) == length(self$params) && all(is.finite(L)),
        "`D` must be a numeric vector of length length(params)" = 
          is.vector(D, "numeric") && length(D) == length(self$params) &&
          all(is.finite(D)) && all(D > 0.0))
      
      self$adaptive <- adaptive

      self$s0 <- as.numeric(s0)
      self$L <- L
      self$D <- D
    },
    
    #' Stable shape for Rcpp
    toList = function() {
      list(
        "meta_for" = "MH",
        "params" = self$params,
        "s0" = self$s0,
        "adaptive" = self$adaptive$toList(),
        "L" = self$L,
        "D" = self$D
      )
    },
    
    print = function(...) {
      cat("<MH>\n")
      cat("  params:", paste(self$params, collapse=", "), "\n")
      cat("  adaptive: \n")
      print(self$adaptive)
      cat("  L: \n")
      print(self$L)
      cat("  D: ", self$D, "\n")
      invisible(self)
    }
  )
)

#' Define a single MH update
#'
#' Convenience constructor for a \code{MH} (R6) describing one MCMC MH update
#' step. Use this to specify which parameters are updated together, plus 
#' optional adaptation settings.
#'
#' @param params Character vector of parameter names included in the MH
#'   (duplicates are removed).
#' @param s0 Numeric scalar; initial step-size scale for MH
#'   (step is \eqn{s \cdot (L z)}).
#' @param adaptive An Adapt object.
#' @param L Numeric matrix (length(params)xlength(params)) giving the initial
#'   settings of the covariance matrix. Defaults to identity matrix.
#' @param D Numeric vector (length(params)) giving the initial
#'   settings of the diagonal shape matrix. Defaults to 1s.
#'
#' @return An R6 \code{MH} object (class \code{"MH"}) with a stable
#'   \code{$toList()} representation suitable for downstream C++/Rcpp code.
#'
#' @examples
#' # MH over two parameters with default adaptation
#' b1 <- MH(params = c("mu", "shape"))
#'
#' # MH with custom step-size and defaultadaptation
#' b2 <- MH(
#'   params = c("mu", "shape"),
#'   s0 = 0.8,
#'   adaptive = adapt())
#'
#' @export
MH <- function(params,
               s0 = 1.0,
               adaptive = adapt(),
               L = diag(nrow = length(params), ncol = length(params)),
               D = rep(1.0, length(params))) {
  
  stopifnot("`params` must be provided" = !missing(params))

  return(MHUpdate$new(params, s0, adaptive, L, D))
}

#' @describeIn MH Inheritance test
is.MH <- function(x) inherits(x, "MHUpdate")

#' A single Gibbs update specification
#' @noRd
#' @keywords internal
GibbsUpdate <- R6::R6Class(
  "GibbsUpdate",
  public = list(
    #' @field params (`character()`)\cr
    #' Vector of parameter tags included in update step
    params = NULL, 

    continue = function(...) {},
    
    #' @description
    #' Defines how updates for a Compartment Distributions are handled
    #'
    #' @param params A vector of parameter tags included in the MH update
    initialize = function(params) {
      
      stopifnot("`params` must be specified" = !missing(params),
                "`params` must be a character vector" = 
                  is.vector(params, "character") && length(params) >= 1L && all(nzchar(params)))
      
      params <- unique(unname(params))
      
      self$params <- params
      
    },
    
    #' Stable shape for Rcpp
    toList = function() {
      list(
        "meta_for" = "Gibbs",
        "params" = self$params
      )
    },
    
    print = function(...) {
      cat("<Gibbs>\n")
      cat("  params:", paste(self$params, collapse=", "), "\n")
      invisible(self)
    }
  )
)

#' Define a single Gibbs update
#'
#' Convenience constructor for a \code{Gibbs} (R6) describing one MCMC Gibbs update
#' step. Use this to specify which parameters are updated together.
#'
#' @param params Character vector of parameter names included in the Gibbs Step.
#'
#' @return An R6 \code{Gibbs} object (class \code{"Gibbs"}) with a stable
#'   \code{$toList()} representation suitable for downstream C++/Rcpp code.
#'
#' @examples
#' # Gibbs over two sensitivity parameters
#' b1 <- Gibbs(params = c("beta0", "beta1"))
#'
#' @export
Gibbs <- function(params) {
  
  stopifnot("`params` must be provided" = !missing(params))
  
  return(GibbsUpdate$new(params))
}

#' @describeIn Gibbs Inheritance test
is.Gibbs <- function(x) inherits(x, "GibbsUpdate")



#' @noRd
#' @keywords internal
UpdateStrategies <- R6::R6Class(
  "UpdateStrategies",
  public = list(
    blocks = NULL,          # list of MH and/or Gibbs
    initialize = function(...) {
      blocks <- list(...)

      if (length(blocks) > 0) {
        stopifnot(all(vapply(blocks, function(b) is.MH(b) || is.Gibbs(b), logical(1L))))
      }
      self$blocks <- unname(blocks)
      names(self$blocks) <- lapply(blocks, function(x) paste(x$params, collapse = "/")) |> unlist()
      invisible(self$validate())
    },
    
    validate = function() {
      if (length(self$blocks) == 0) return(invisible(TRUE))
      
      # no duplicate param names across blocks unless allowed
      all_params <- unlist(lapply(self$blocks, function(b) b$params), use.names = FALSE)
      dups <- unique(all_params[duplicated(all_params)])
      if (length(dups)) {
        stop(sprintf("Strategies: parameters appear in multiple update strategies: %s",
                     paste(dups, collapse = ", ")), call. = FALSE)
      }
      invisible(TRUE)
    },
    
    toList = function() {
      list(
        meta_for = "Strategies",
        n_blocks = length(self$blocks),
        blocks = lapply(self$blocks, function(b) b$toList())
      )
    },
    
    print = function(...) {
      cat("<Strategies>", "\n")
      if (length(self$blocks) == 0L) {
        cat("  (no update strategies)\n")
      } else {
        for (i in seq_along(self$blocks)) {
          cat(" ", sprintf("[%02d] ", i)); self$blocks[[i]]$print()
        }
      }
      invisible(self)
    }
  )
)

#' Define MCMC update strategies
#'
#' Convenience constructor for a \code{Strategies} (R6) object that groups
#'   multiple \code{\link[=MH]{MH}} or \code{\link[=Gibbs]{Gibbs}} definitions 
#'   into an ordered strategy.
#'
#' @noRd
#' @param ... One or more \code{\link[=MH]{MH}} or \code{\link[=Gibbs]{Gibbs}} objects. Each
#'   defines one MH/Gibbs update block.
#'
#' @details
#' Use  \code{$toList()} to generate a stable list form for passing to C++.
#'
#' @return An R6 \code{Strategies} object (class \code{"Strategies"})
#'   containing the specified sequence of blocks.
#'
#' @section Methods:
#' \itemize{
#'   \item \code{$validate()} — Ensure no parameter appears in multiple blocks.
#'   \item \code{$toList()} — Convert the strategy to a serializable list for
#'     downstream C++ or storage.
#'   \item \code{$print()} — Compact human-readable summary.
#' }
#'
#' @examples
#' # Single strategy
#' s1 <- Strategies(list(MH(params = c("mu", "shape"))))
#'
#' # Multi-block strategy
#' s2 <- Strategies(list(
#'   MH(params = c("mu", "shape")),
#'   MH(params = "psi")
#' ))
#'
#' # Validate uniqueness of parameter names
#' s2$validate()
#'
#' # Convert to list for downstream use
#' s2$toList()
#'
#' @keywords internal
Strategies <- function(...) {
  UpdateStrategies$new(...)
}

#' @noRd
#' @keywords internal
is.Strategies <- function(x) inherits(x, "UpdateStrategies")

#' @noRd
#' @keywords internal
.validateModel <- function(H, P, beta, psi, strategies) {
  
  stopifnot(is.Compartment(H), is.Compartment(P), 
            is.Compartment(beta), is.Compartment(psi),
            is.Strategies(strategies))
  
  H_vars <- H$getTags()
  P_vars <- P$getTags()
  beta_vars <- beta$getTags()
  psi_vars <- psi$getTags()
  
  # first and foremost, strategies must contain parameters of model components
  all_strategy_params <- lapply(strategies$blocks, function(b) b$params)

  tst <- all(all_strategy_params %in% c(H_vars, P_vars, beta_vars, psi_vars))
  if (!tst) stop("`strategies` contains parameter(s) not found in model components")
  
  # Currently, H, P, and psi must be updated through MH only
  MH_strategies <- lapply(strategies$blocks, is.MH)
  
  Gibbs_params <- unlist(all_strategy_params[!MH_strategies])
  
  if (length(Gibbs_params) > 0L &&
      any(c(H_vars, P_vars, psi_vars) %in% Gibbs_params)) {
    stop("at this time, only beta parameters can employ Gibbs", call. = FALSE)
  }
  
}