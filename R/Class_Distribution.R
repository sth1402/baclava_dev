#' @noRd
#' @include Class_Model.R Class_Param.R
#' @keywords internal
Distribution <- R6Class(
  "Distribution",
  public = list(
    #' @field type (`character(1)`)\cr
    #' character representation of probability density function
    type = NULL,
    #' @field params (`list()`)\cr
    #' list of distribution parameter specifications
    parameters = NULL,
    
    #' @description
    #' Defines how a model or distribution parameter. 
    #'
    #' @param type The probability density function specified as a string
    #' @param parameters A list of distribution parameter specifications
    initialize = function(type, parameters) {
      
      stopifnot(
        "`type` must be a character" = !missing(type) && 
          is.vector(type, "character") && length(type) == 1,
        "`parameters` must be a list" = !missing(parameters) &&
          is.list(parameters) && length(parameters) > 0L
      )
      
      self$type <- type

      # Shift parameters that are not modeled to be constant Models
      self$parameters <- lapply(parameters, 
                                function(x) { 
                                  if (is.Param(x)) return(const_m(value = x))
                                  if (!is.Model(x))
                                    stop("all distribution parameters must be Param or Model objects")
                                  x
                                })
      names(self$parameters) <- names(parameters)

    },
    
    #' @description
    #' Return all parameter metadata as a plain list.
    toList = function() {
      params <- lapply(self$parameters, function(x) x$toList())
      names(params) <- names(self$parameters)
      list("meta_for" = "Distribution", 
           "class" = class(self)[1L],
           "type" = self$type,
           "params" = params)
    },
    
    setValue = function(values) {
      for (i in seq_along(values)) {
        for (j in seq_along(self$parameters)) {
          self$parameters[[j]]$setValue(value = values[i], tag = names(values)[i])
        }
      }
      invisible(self)
    },

    getPriorDist = function() { lapply(self$parameters, function(x) x$getPriorDist())},
    getPriorParams = function() { lapply(self$parameters, function(x) x$getPriorParams())},
    getTags = function() { lapply(self$parameters, function(x) x$getTags())}
  )
)

#' @noRd
#' @keywords internal
Distribution$fromList <- function(lst) {
  if (lst$meta_for != "Distribution") stop("must provide a Distribution list")
  klass <- getFromNamespace(lst$class, "baclava")
  params <- lapply(lst$params, function(x) {
    switch(x$meta_for,
           "Param"  = Param$fromList(x),
           "Model"  = Model$fromList(x),
           stop("Unknown meta_for in parameter list"))
  })
  klass$new(lst$type, params)
}


#' Specification of Model Components
#'
#' The healthy and preclinical compartments of the model can be specified as
#'   Weibull or Gamma distributions. 
#' The screening sensitivity and indolent fraction are modeled as
#'   Bernoulli distributions. 
#'
#' @param ... Ignored. Included to require named inputs.
#' @param shape A 'Param' object specifying the shape parameters. Use only
#'   if Distribution is Weibull or Gamma. 
#' @param rate If fitting the rate parameter of the distribution, a 'Param' or 
#'   'Model' object. If fitting the mean or scale, do not specify. Use only
#'   if Distribution is Weibull or Gamma. 
#' @param mu If fitting the mean parameter of the distribution, a 'Param' or 
#'   'Model' object. If fitting the rate or scale, do not specify. Use only
#'   if Distribution is Weibull or Gamma. 
#' @param scale If fitting the scale parameter of the distribution, a 'Param' or 
#'   'Model' object. If fitting the rate or mu, do not specify. Use only
#'   if Distribution is Weibull or Gamma. 
#' @param p If fitting the probability parameter of the distribution, a 'Param' or 
#'   'Model' object. Use only if Distribution is Bernoulli. 
#'   
#' @returns An R6 object of class "Distribution". 
#' 
#' @include Class_Prior.R
#' @name Distribution
#' @rdname Distribution
NULL

#' @noRd
#' @keywords internal
WeibullDistribution_mu <- R6Class(
  "WeibullDistribution_mu",
  inherit = Distribution,
  public = list(
    scale = function(tau, theta) {
      mu <- self$parameters$mu$evaluate(tau, theta)
      shape <- theta[length(theta)]
      mu / gamma(1.0 + 1.0 / shape)
    },
    ddist = function(x, tau, theta, ...) {
      shape <- theta[length(theta)]
      scale <- self$scale(tau, theta)
      dweibull(x, shape = shape, scale = scale, ...)      
    },
    pdist = function(x, tau, theta, ...) {
      shape <- theta[length(theta)]
      scale <- self$scale(tau, theta)
      pweibull(x, shape = shape, scale = scale, ...)      
    },
    rdist = function(n, tau, theta, ...) {
      shape <- theta[length(theta)]
      scale <- self$scale(tau, theta)
      rweibull(n, shape = shape, scale = scale)      
    },
    qdist = function(p, tau, theta, ...) {
      shape <- theta[length(theta)]
      scale <- self$scale(tau, theta)
      qweibull(p, shape = shape, scale = scale, ...)      
    }
  )
)

#' @noRd
#' @keywords internal
WeibullDistribution_rate <- R6Class(
  "WeibullDistribution_rate",
  inherit = Distribution,
  public = list(
    scale = function(tau, theta) {
      rate <- self$parameters$rate$evaluate(tau, theta)
      shape <- theta[length(theta)]
      rate^{- 1.0 / shape}
    },
    ddist = function(x, tau, theta, ...) {
      shape <- theta[length(theta)]
      scale <- self$scale(tau, theta)
      dweibull(x, shape = shape, scale = scale, ...)      
    },
    pdist = function(x, tau, theta, ...) {
      shape <- theta[length(theta)]
      scale <- self$scale(tau, theta)
      pweibull(x, shape = shape, scale = scale, ...)      
    },
    rdist = function(n, tau, theta, ...) {
      shape <- theta[length(theta)]
      scale <- self$scale(tau, theta)
      rweibull(n, shape = shape, scale = scale)      
    },
    qdist = function(p, tau, theta, ...) {
      shape <- theta[length(theta)]
      scale <- self$scale(tau, theta)
      qweibull(p, shape = shape, scale = scale, ...)      
    }
  )
)

#' @noRd
#' @keywords internal
WeibullDistribution_scale <- R6Class(
  "WeibullDistribution_scale",
  inherit = Distribution,
  public = list(
    scale = function(tau, theta) {
      self$parameters$scale$evaluate(tau, theta)
    },
    ddist = function(x, tau, theta, ...) {
      shape <- theta[length(theta)]
      scale <- self$scale(tau, theta)
      dweibull(x, shape = shape, scale = scale, ...)      
    },
    pdist = function(x, tau, theta, ...) {
      shape <- theta[length(theta)]
      scale <- self$scale(tau, theta)
      pweibull(x, shape = shape, scale = scale, ...)      
    },
    rdist = function(n, tau, theta, ...) {
      shape <- theta[length(theta)]
      scale <- self$scale(tau, theta)
      rweibull(n, shape = shape, scale = scale)      
    },
    qdist = function(p, tau, theta, ...) {
      shape <- theta[length(theta)]
      scale <- self$scale(tau, theta)
      qweibull(p, shape = shape, scale = scale, ...)      
    }
  )
)

#' @noRd
#' @keywords internal
WeibullDistribution_median <- R6Class(
  "WeibullDistribution_median",
  inherit = Distribution,
  public = list(
    scale = function(tau, theta) {
      median <- self$parameters$median$evaluate(tau, theta)
      shape <- theta[length(theta)]
      median * log(2.0)^{- 1.0 / shape}
    },
    ddist = function(x, tau, theta, ...) {
      shape <- theta[length(theta)]
      scale <- self$scale(tau, theta)
      dweibull(x, shape = shape, scale = scale, ...)      
    },
    pdist = function(x, tau, theta, ...) {
      shape <- theta[length(theta)]
      scale <- self$scale(tau, theta)
      pweibull(x, shape = shape, scale = scale, ...)      
    },
    rdist = function(n, tau, theta, ...) {
      shape <- theta[length(theta)]
      scale <- self$scale(tau, theta)
      rweibull(n, shape = shape, scale = scale)      
    },
    qdist = function(p, tau, theta, ...) {
      shape <- theta[length(theta)]
      scale <- self$scale(tau, theta)
      qweibull(p, shape = shape, scale = scale, ...)      
    }
  )
)

#' @describeIn Distribution Weibull Distribution
#' @export
weibull_d <- function(..., shape, rate, mu, scale, median) {
  
  if (!missing(shape)) {
    if (shape$prop_scale == "logcv") shape$setPropScale("logcv_weibull")
    d_param = list("shape" = shape)
  } else {
    stop("shape must be provided")
  }
  
  args_given <- c(mu = !missing(mu), rate = !missing(rate), 
                  scale = !missing(scale), median = !missing(median))
  if (sum(args_given) != 1L)
    stop("exactly one of mu, rate, scale, or median must be provided")

  if (missing(rate) && missing(scale) && missing(median) && !missing(mu)) {
    WeibullDistribution_mu$new("dweibull", append(d_param, list("mu" = mu)))
  } else if (missing(mu) && missing(scale) && missing(median) && !missing(rate)) {
    WeibullDistribution_rate$new("dweibull", append(d_param, list("rate" = rate)))
  } else if (missing(mu) && missing(rate) && missing(median) && !missing(scale)) {
    WeibullDistribution_scale$new("dweibull", append(d_param, list("scale" = scale)))
  } else if (missing(mu) && missing(rate) && missing(scale) && !missing(median)) {
    WeibullDistribution_median$new("dweibull", append(d_param, list("median" = median)))
  } else {
    stop("one of scale/rate/mu/median must be provided")
  }
  
}

#' @noRd
#' @keywords internal
GammaDistribution_mu <- R6Class(
  "GammaDistribution_mu",
  inherit = Distribution,
  public = list(
    scale = function(tau, theta) {
      mu <- self$parameters$mu$evaluate(tau, theta)
      shape <- theta[length(theta)]
      mu / shape
    },
    ddist = function(x, tau, theta, ...) {
      shape <- theta[length(theta)]
      scale <- self$scale(tau, theta)
      dgamma(x, shape = shape, scale = scale, ...)      
    },
    pdist = function(x, tau, theta, ...) {
      shape <- theta[length(theta)]
      scale <- self$scale(tau, theta)
      pgamma(x, shape = shape, scale = scale, ...)      
    },
    rdist = function(n, tau, theta, ...) {
      shape <- theta[length(theta)]
      scale <- self$scale(tau, theta)
      rgamma(n, shape = shape, scale = scale)      
    },
    qdist = function(p, tau, theta, ...) {
      shape <- theta[length(theta)]
      scale <- self$scale(tau, theta)
      qgamma(p, shape = shape, scale = scale)      
    }
  )
)

#' @noRd
#' @keywords internal
GammaDistribution_rate <- R6Class(
  "GammaDistribution_rate",
  inherit = Distribution,
  public = list(
    scale = function(tau, theta) {
      rate <- self$parameters$rate$evaluate(tau, theta)
      1.0 / rate
    },
    ddist = function(x, tau, theta, ...) {
      shape <- theta[length(theta)]
      scale <- self$scale(tau, theta)
      dgamma(x, shape = shape, scale = scale, ...)      
    },
    pdist = function(x, tau, theta, ...) {
      shape <- theta[length(theta)]
      scale <- self$scale(tau, theta)
      pgamma(x, shape = shape, scale = scale, ...)      
    },
    rdist = function(n, tau, theta, ...) {
      shape <- theta[length(theta)]
      scale <- self$scale(tau, theta)
      rgamma(n, shape = shape, scale = scale)      
    },
    qdist = function(p, tau, theta, ...) {
      shape <- theta[length(theta)]
      scale <- self$scale(tau, theta)
      qgamma(p, shape = shape, scale = scale)      
    }
  )
)

#' @noRd
#' @keywords internal
GammaDistribution_scale <- R6Class(
  "GammaDistribution_scale",
  inherit = Distribution,
  public = list(
    scale = function(tau, theta) {
      self$parameters$scale$evaluate(tau, theta)
    },
    ddist = function(x, tau, theta, ...) {
      shape <- theta[length(theta)]
      scale <- self$scale(tau, theta)
      dgamma(x, shape = shape, scale = scale, ...)      
    },
    pdist = function(x, tau, theta, ...) {
      shape <- theta[length(theta)]
      scale <- self$scale(tau, theta)
      pgamma(x, shape = shape, scale = scale, ...)      
    },
    rdist = function(n, tau, theta, ...) {
      shape <- theta[length(theta)]
      scale <- self$scale(tau, theta)
      rgamma(n, shape = shape, scale = scale)      
    },
    qdist = function(p, tau, theta, ...) {
      shape <- theta[length(theta)]
      scale <- self$scale(tau, theta)
      qgamma(p, shape = shape, scale = scale)      
    }
  )
)


#' @describeIn Distribution Gamma Distribution
#' @export
gamma_d <- function(..., shape, rate, mu, scale) {
  
  if (!missing(shape)) {
    if (shape$prop_scale == "logcv") shape$setPropScale("logcv_gamma")
    d_param = list("shape" = shape)
  } else {
    stop("shape must be provided")
  }
  
  args_given <- c(mu = !missing(mu), rate = !missing(rate), 
                  scale = !missing(scale))
  if (sum(args_given) != 1L)
    stop("exactly one of mu, rate, or scale must be provided")
  
  if (missing(rate) && missing(scale) && !missing(mu)) {
    GammaDistribution_mu$new("dgamma", append(d_param, list("mu" = mu)))
  } else if (missing(mu) && missing(scale) && !missing(rate)) {
    GammaDistribution_rate$new("dgamma", append(d_param, list("rate" = rate)))
  } else if (missing(mu) && missing(rate) && !missing(scale)) {
    GammaDistribution_scale$new("dgamma", append(d_param, list("scale" = scale)))
  } else {
    stop("one of scale/rate/mu must be provided")
  }
}

#' @noRd
#' @keywords internal
BernDistribution <- R6Class(
  "BernDistribution",
  inherit = Distribution,
  public = list(
    ddist = function(x, tau, theta, ...) {
      p <- self$parameters$p$evaluate(tau, theta)
      dbinom(x, size = 1L, prob = p, ...)      
    },
    pdist = function(x, tau, theta, ...) {
      p <- self$parameters$p$evaluate(tau, theta)
      pbinom(x, size = 1L, prob = p, ...)      
    },
    rdist = function(n, tau, theta, ...) {
      p <- self$parameters$p$evaluate(tau, theta)
      rbinom(n, size = 1L, prob = p)      
    },
    qdist = function(p, tau, theta, ...) {
      pr <- self$parameters$p$evaluate(tau, theta)
      qbinom(p, size = 1L, prob = pr)      
    },
    prob = function(tau, theta) {
      self$parameters$p$evaluate(tau, theta)
    }
  )
)

#' @describeIn Distribution Bernoulli Distribution
#' @export
bern_d <- function(..., p) {
  if (missing(p)) stop("p must be provided")
  BernDistribution$new("dbinom", list("p" = p))
}

#' @describeIn Distribution Inheritance test
is.Distribution <- function(x) inherits(x, "Distribution")

#' @describeIn Distribution Print
print.Distribution <- function(x, ...) {
  cat("Distribution:", x$type, "\n")
  cat("Parameters:\n")
  for (nm in names(x$parameters)) {
    cat("  ", nm, ":\n", sep = "")
    print(x$parameters[[nm]])
  }
  invisible(x)
}
