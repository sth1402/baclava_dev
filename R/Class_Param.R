#' @noRd
#' @include Class_Prior.R
#' @keywords internal
Param <- R6Class(
  "Param",
  public = list(
    #' @field tag (`character(1)`)\cr
    #' unique tag
    tag = NULL,
    #' @field init (`numeric(1)`)\cr
    #' initial value
    init = NULL, 
    #' @field minv (`numeric(1)`)\cr
    #' minimum allowed value of parameter
    minv = -Inf,
    #' @field maxv (`numeric(1)`)\cr
    #' maximum allowed value of parameter
    maxv = Inf,
    #' @field prop_scale (`character(1)`)\cr
    #' proposal scale `"id", "log", "logit", "logcv"`
    prop_scale = NULL,
    #' @field prior (`Prior`)\cr
    prior = NULL,
    #' @field fixed (`logical(1)`)\cr
    #' true if parameter is fixed (no prior specified)
    fixed = NULL,
    prop_func = NULL,
    t_df = NULL,
    
    #' @description
    #' Defines how a model or distribution parameter. 
    #'
    #' @param init The starting value for the parameter.
    #' @param tag A unique identifier character string.
    #' @param minv The minimum allowed value on the natural scale.
    #' @param maxv The maximum allowed value on the natural scale.
    #' @param prop_scale Scale on which proposals are made one of 'id', 'log', 'logit', 'logcv'
    #' @param prior A `Prior` object specifying the prior distribution for the parameter.
    initialize = function(init, tag, ..., 
                          minv = -Inf, maxv = Inf, logcv = FALSE, 
                          prior = NULL, force_natural = FALSE, 
                          prop_func = "rnorm", t_df = 7.0) {
      
      stopifnot(
        "`init` must be a scalar numeric" = !missing(init) &&
          is.vector(init, "numeric") && length(init) == 1L,
        "`tag` must be a character string" = !missing(tag) &&
          is.vector(tag) && length(tag) == 1L,
        "`minv` must be a scalar numeric" = is.vector(minv, "numeric") && length(minv) == 1L,
        "`maxv` must be a scalar numeric" = is.vector(maxv, "numeric") && length(maxv) == 1L,
        "`minv` must be < `maxv`" = minv < maxv,
        "`logcv` must be a logical" = is.vector(logcv, "logical") && length(logcv) == 1L,
        "`prior` must be a Prior object" = is.null(prior) || is.Prior(prior),
        "`prop_func` must be one of 'rnorm' or 'runif'" = prop_func %in% c("rnorm", "runif", "rt"),
        "`t_df` must be a scalar numeric" = is.vector(t_df, "numeric") && length(t_df) == 1L
      )
      
      self$tag <- as.character(tag)
      self$init <- init
      
      if (minv >= maxv) stop("minv must be < maxv", call. = FALSE)
      self$minv <- minv
      self$maxv <- maxv
      
      if (logcv) {
        if (!isTRUE(all.equal(minv, 0.0)) || is.finite(maxv)) {
          stop("`logcv=TRUE` requires minv = 0; maxv = Inf")
        }
        self$prop_scale <-  "logcv"
      } else {
        if (is.finite(minv) && is.finite(maxv)) {
          self$prop_scale <-  "logit"
        } else if (is.finite(minv) && is.infinite(maxv)) {
          if (isTRUE(all.equal(minv, 0.0))) {
            self$prop_scale <- "log_pos"
          } else self$prop_scale <-  "log_shifted" 
        } else if (is.infinite(minv) && is.finite(maxv)) {
          self$prop_scale <-  "log_flipped"
        } else {
          self$prop_scale <-  "id"
        }
      }
      
      if (force_natural) self$prop_scale <- "id"
      self$prop_func <- prop_func
      self$t_df <- t_df
      
      self$prior <- prior
      self$fixed <- FALSE
      # if no prior given, parameter is assumed fixed
      if (is.null(self$prior)) {
        self$fixed <- TRUE
        self$minv <- init
        self$maxv <- init
        self$prop_scale <- "id"
        self$prior <- none_p()
      }
      
      self$validate()
    },
    
    validate = function() {
      
      if (self$minv > self$maxv) stop("minv > maxv", call. = FALSE)

      if (self$init < self$minv || self$init > self$maxv) {
        stop("initial value must be within specified bounds", call. = FALSE)
      }

    },
    
    #' @description
    #' Return all parameter metadata as a plain list.
    toList = function() {
      list("meta_for" = "Param",
           "tag" = self$tag, 
           "init" = self$init,
           "minv" = self$minv,
           "maxv" = self$maxv,
           "prop_scale" = self$prop_scale, 
           "prop_func" = self$prop_func,
           "t_df" = self$t_df,
           "prior" = self$prior$toList(),
           "fixed" = self$fixed)
    },
    
    setValue = function(value) {
      self$init = value
      invisible(self)
    },
    
    setPropScale = function(value) {
      self$prop_scale = value
      invisible(self)
    },
    
    #' @description
    #' Display parameter details.
    print = function() {
      if (!self$fixed) {
        cat(sprintf(
          "Parameter %s:\n  Initial value : %g\n  Bounds        : (%g, %g)\n  Prior         : %s\n  Update Scale  : %s\n",
          self$tag, self$init, self$minv, self$maxv, self$prior$asString(), self$prop_scale
        ))
      } else {
        cat("Parameter ", self$tag, ": fixed @ ", self$init)
      }
      invisible(self)
    },
    
    getPriorDist = function() { self$prior$getDist() },
    getPriorParams = function() { self$prior$getParams() },
    getTag = function() { self$tag }
  )
)

#' @noRd
#' @keywords internal
Param$fromList <- function(lst) {
  if (lst$meta_for != "Param") stop("must provide a Param list")
  prior <- Prior$fromList(lst$prior)
  param <- Param$new(init = lst$init, tag = lst$tag, minv = lst$minv, maxv = lst$maxv,
                     prop_scale = lst$prop_scale, prior = prior, 
                     prop_func = lst$prop_func, t_df = lst$t_df)
  invisible(param)
}

#' Define a single model or distribution parameter
#'
#' Convenience constructor for a \code{Param} (R6) object representing a single
#'   parameter within a model or distribution. Each parameter has a unique name,
#'   an initial value, and an associated prior distribution, optional bounds,
#'   and scale on which proposals are made if it is estimated rather than fixed.
#'
#' @param tag Character scalar giving the unique parameter tag.
#' @param init Numeric scalar giving the initial (or fixed) value.
#' @param minv Optional numeric scalar giving the minimum allowed value.
#' @param maxv Optional numeric scalar giving the maximum allowed value.
#' @param prop_scale Optional character scalar specifying the transformation used
#'   during MCMC updates. Must be one of:
#'   \itemize{
#'     \item \code{"id"} — identity (no transformation)
#'     \item \code{"log"} — log-transformed updates (for strictly positive parameters)
#'     \item \code{"logit"} — logit-transformed updates (for parameters bounded in (0,1))
#'     \item \code{"logcv"} — log(CV) transformation for shape parameters
#'   }
#' @param prior Optional \code{\link[=Prior]{Prior}} object specifying the prior
#'   distribution. If omitted, the parameter is treated as fixed (i.e.,
#'   not updated during MCMC).
#'
#' @details
#' The \code{Param} class provides a consistent representation for all scalar
#'   model and distribution parameters. Each object stores its initialization,
#'   bounds, prior, and working scale. Parameters with no prior are automatically
#'   marked as fixed and excluded from Metropolis–Hastings or Gibbs updates.
#'
#' @return An R6 \code{Param} object with fields:
#' \itemize{
#'   \item \code{tag} — parameter tag
#'   \item \code{init} — initial value
#'   \item \code{minv} — minimum allowed value
#'   \item \code{maxv} — maximum allowed value
#'   \item \code{prop_scale} — transformation scale
#'   \item \code{prior} — \code{Prior} object
#'   \item \code{fixed} — logical, \code{TRUE} if fixed
#' }
#'
#' @section Methods:
#' \itemize{
#'   \item \code{$fromList(lst)} — Return a Param object from a list generated by $toList().
#'   \item \code{$getPriorDist()} - Retrieve prior distribution name.
#'   \item \code{$getPriorParms()} - Retrieve prior distribution parameters.
#'   \item \code{$getTag()} - Retrieve parameter tag/name.
#'   \item \code{$setValue(value)} - Reset initial value.
#'   \item \code{$setPropScale(value)} - Reset prop_scale.
#'   \item \code{$print()} — Display parameter details.
#'   \item \code{$toList()} — Return all parameter metadata as a plain list.
#' }
#'
#' @examples
#' # Fixed parameter (no prior)
#' p1 <- param(init = 1.0, tag = "p1")
#'
#' # Free parameter with log-scale updates and bounds
#' p2 <- param(init = 0.5,
#'             tag = "sigma", 
#'             minv = 1e-3, maxv = 10,
#'             prop_scale = "log",
#'             prior = lnorm_p(meanlog = 0, sdlog = 1))
#'
#' # View parameter details
#' print(p2)
#'
#' @export
param <- function(init, tag, ...,  minv = -Inf, maxv = Inf, 
                  prop_scale = c("id", "log", "logit", "logcv"), 
                  prior = NULL, force_natural = FALSE, 
                  prop_func = "rnorm", t_df = 7.0) {
  stopifnot(
    "`init` must be provided" = !missing(init),
    "`tag` must be provided" = !missing(tag)
  )
  
  Param$new(init = init, tag = tag, minv = minv, maxv = maxv,
            prop_scale = prop_scale, prior = prior, force_natural = force_natural,
            prop_func = prop_func, t_df = t_df)
}

#' @describeIn param Inheritance test
is.Param <- function(object) inherits(object, "Param")

#' @describeIn param Print
#' @export 
print.Param <- function(x, ...) x$print()