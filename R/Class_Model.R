#' @noRd
#' @include Class_Param.R
#' @keywords internal
Model <- R6Class(
  "Model", 
  public = list(
    #' @field type (`character(1)`)\cr
    #' type of model
    type = NULL,
    #' @field parameters (`list()`)\cr
    #' list of Param objects, one for each model parameter
    parameters = list(),
    
    #' @description
    #' Defines a model; defined in derived classess
    initialize = function(...) {},
    
    #' @description
    #' Return all parameter metadata as a plain list.
    toList = function() {
      # create a list containing all model parameters
      param_list <- lapply(self$parameters, function(x) x$toList())
      param_tags <- lapply(self$parameters, "[[", "tag")
      names(param_list) <- param_tags
      
      list("meta_for" = "Model",
           "type" = self$type, 
           "params" = param_list)
    },
    
    setValue = function(value, tag) {
      idx <- which(vapply(self$parameters, function(x) x$tag == tag, logical(1)))
      if (length(idx) == 0L) return()
      self$parameters[[idx]]$setValue(value)
      invisible(self)
    },
    
    getPriorDist = function() {
      lapply(self$parameters, function(x) x$getPriorDist())
    },
    getPriorParams = function() {
      lapply(self$parameters, function(x) x$getPriorParams())
    },
    getTags = function() { lapply(self$parameters, function(x) x$getTag() )}
  ))

Model$fromList <- function(lst) {
  if (lst$meta_for != "Model") stop("must provide a Model list")
  param_objs <- lapply(lst$params, Param$fromList)
  mdl <- Model$new()
  mdl$type <- lst$type
  mdl$parameters <- param_objs
  invisible(mdl)
}


#' @noRd
#' @keywords internal
SigmoidModel <- R6Class(
  "SigmoidModel", 
  inherit = Model,
  public = list(
    initialize = function(...) {
      self$type = "sigmoid"
      
      args <- list(...)
      stopifnot(all(sapply(args, function(x) is.Param(x))))
      
      expected <- c("r0", "delta", "m", "w")
      
      actual <- names(args)
      stopifnot("provided parameters do not match expectations" = 
                  identical(sort(expected), sort(actual)))
      
      self$parameters <- args[match(expected, names(args))]
    },
    
    evaluate = function(x, theta) {
      # ordering is r0, delta, m, w
      theta[1L] + 
        theta[2L] / (1.0 + exp(-{x - theta[3L]} / theta[4L]))
    }
  ))

#' A Sigmoid Model
#' 
#' The parameter is modeled as
#' \eqn{= r_0 + \delta / (1 + exp(-(\tau_{HP} - m) / w))}{= r0 + delta / (1 + exp(-(tau_HP - m) / w},
#' where \eqn{\tau_{HP}}{tau_HP} is the age at which a participant transitions
#' into the pre-clinical compartment.
#'
#' @param ... Ignored. Used to require named inputs.
#' @param r0 An 'UpdateType' object with type = "MH" or "Fixed".
#' @param delta An 'UpdateType' object with type = "MH" or "Fixed".
#' @param m An 'UpdateType' object with type = "MH" or "Fixed".
#' @param w An 'UpdateType' object with type = "MH" or "Fixed".
#' 
#' @returns An R6 object inheriting from class Model.
#' 
#' @examples
#' 
#' sigmoid_m(r0 = param(1.0, "r0", minv = 0, maxv = Inf, 
#'                      prop_scale = "id", prior = lnormp(log(2), 0.2)),
#'           delta = param(1.0, "delta", prior = normp(2, 1)),
#'           m = param(52, "m"),
#'           w = param(1, "w", minv = 0, maxv = Inf, prior = lnormp(log(1), 0.2)))
#'           
#' @export
sigmoid_m = function(..., r0, delta, m, w) {
  
  stopifnot(
    "all inputs must be specified" = !missing(r0) && !missing(delta) &&
      !missing(m) && !missing(w)
  )
  
  if (length(list(...)) > 0L) warning("extra inputs provided to sigmoidModel")
  
  SigmoidModel$new(r0 = r0, delta = delta, m = m, w = w)
  
}

#' @noRd
#' @keywords internal
LinearModel <- R6Class(
  "LinearModel", 
  inherit = Model,
  public = list(
    initialize = function(...) {
      self$type = "linear"
      
      args <- list(...)
      stopifnot(all(sapply(args, function(x) is.Param(x))))
      
      expected <- c("r0", "delta", "m")
      
      actual <- names(args)
      stopifnot("provided parameters do not match expectations" = 
                  identical(sort(expected), sort(actual)))
      
      self$parameters <- args[expected]
    },
    evaluate = function(x, theta) {
      # ordering is r0, delta, m
      theta[1L] + theta[2L] * {x - theta[3L]}
    }
  ))

#' A Linear Model
#' 
#' The parameter is modeled as 
#' \eqn{= r_0 + \delta \times (\tau_{HP} - m) }{= r0 + delta * (tau_HP - m)},
#' where \eqn{\tau_{HP}}{tau_HP} is the age at which a participant transitions
#' into the pre-clinical compartment.
#' 
#' @param ... Ignored. Used to require named inputs.
#' @param r0 An 'UpdateType' object with type = "MH" or "Fixed".
#' @param delta An 'UpdateType' object with type = "MH" or "Fixed".
#' @param m An 'UpdateType' object with type = "MH" or "Fixed".
#' 
#' @returns An R6 object inheriting from class Model.
#' 
#' @examples
#' 
#' # m is fixed
#' linear_m(r0 = param(1.0, "r0", minv = 0, maxv = Inf, 
#'                     prop_scale = "id", prior = lnormp(log(2), 0.2)),
#'          delta = param(1.0, "delta", prior = normp(2, 1)),
#'          m = param(52, "m"))
#' 
#' @export
linear_m = function(..., r0, delta, m) {
  
  stopifnot(
    "all inputs must be specified" = !missing(r0) && !missing(delta) &&
      !missing(m)
  )
  
  if (length(list(...)) > 0L) warning("extra inputs provided to linearModel")
  
  LinearModel$new(r0 = r0, delta = delta, m = m)
}

#' @noRd
#' @keywords internal
ExpModel <- R6Class(
  "ExpModel", 
  inherit = Model,
  public = list(
    initialize = function(...) {
      self$type = "exp"
      
      args <- list(...)
      stopifnot(all(sapply(args, function(x) is.Param(x))))
      
      expected <- c("r0", "delta", "m")
      
      actual <- names(args)
      stopifnot("provided parameters do not match expectations" = 
                  identical(sort(expected), sort(actual)))
      
      self$parameters <- args[expected]
    },
    evaluate = function(x, theta) {
      # ordering is r0, delta, m
      theta[1L] * exp(theta[2L] * {x - theta[3L]})
    }
  ))

#' An Exponential Model
#' 
#' The parameter is modeled as 
#' \eqn{= r_0 \exp(\delta \times (\tau_{HP} - m))}{= r0 \exp(delta * (tau_HP - m))},
#' where \eqn{\tau_{HP}}{tau_HP} is the age at which a participant transitions
#' into the pre-clinical compartment.
#' 
#' @param ... Ignored. Used to require named inputs.
#' @param r0 An 'UpdateType' object with type = "MH" or "Fixed".
#' @param delta An 'UpdateType' object with type = "MH" or "Fixed".
#' @param m An 'UpdateType' object with type = "MH" or "Fixed".
#' 
#' @returns An R6 object inheriting from class Model.
#' 
#' @examples
#' 
#' # m is fixed
#' exp_m(r0 = param(1.0, "r0", minv = 0, maxv = Inf, 
#'                  prop_scale = "id", prior = lnormp(log(2), 0.2)),
#'       delta = param(1.0, "delta", prior = normp(2, 1)),
#'       m = param(52, "m"))
#' 
#' @export
exp_m = function(..., r0, delta, m) {
  
  stopifnot(
    "all inputs must be specified" = !missing(r0) && !missing(delta) &&
      !missing(m)
  )
  
  if (length(list(...)) > 0L) warning("extra inputs provided to expModel")
  
  ExpModel$new(r0 = r0, delta = delta, m = m)
}

#' @noRd
#' @keywords internal
StepModel <- R6Class(
  "StepModel", 
  inherit = Model,
  public = list(
    initialize = function(...) {
      self$type = "step"
      
      args <- list(...)
      stopifnot(all(sapply(args, function(x) is.Param(x))))
      
      expected <- c("r0", "r1", "m")
      
      actual <- names(args)
      stopifnot("provided parameters do not match expectations" = 
                  identical(sort(expected), sort(actual)))
      
      self$parameters <- args[expected]
    },
    evaluate = function(x, theta) {
      # ordering is r0, r1, m
      theta[1L] * {x < theta[3L]} +
        theta[2L] * {x >= theta[3L]}
    }
  ))

#' A Step Model
#' 
#' The parameter is modeled as 
#' \eqn{= r_0 (\tau_{HP} < m) + r_1 (\tau_{HP} >= m)}{= r_0 (\tau_{HP} < m) + r_1 (\tau_{HP} >= m)}
#' where \eqn{\tau_{HP}}{tau_HP} is the age at which a participant transitions
#' into the pre-clinical compartment.
#'
#' @param ... Ignored. Used to require named inputs.
#' @param r0 An 'UpdateType' object with type = "MH" or "Fixed".
#' @param r1 An 'UpdateType' object with type = "MH" or "Fixed".
#' @param m An 'UpdateType' object with type = "MH" or "Fixed".
#' 
#' @returns An R6 object inheriting from class Model.
#' 
#' @examples
#' 
#' # m is fixed
#' step_m(r0 = param(1.0, "r0", minv = 0, maxv = Inf, 
#'                   prop_scale = "id", prior = lnormp(log(2), 0.2)),
#'        r1 = param(1.0, "r0", minv = 0, maxv = Inf, 
#'                   prop_scale = "id", prior = lnormp(log(2), 0.2)),
#'        m = param(52, "m"))
#'           
#' @export
step_m = function(..., r0, r1, m) {
  
  stopifnot(
    "all inputs must be specified" = !missing(r0) && !missing(r1) &&
      !missing(m)
  )
  
  if (length(list(...)) > 0L) warning("extra inputs provided to stepModel")
  
  StepModel$new(r0 = r0, r1 = r1, m = m)
}

#' @noRd
#' @keywords internal
ConstantModel <- R6Class(
  "ConstantModel", 
  inherit = Model,
  public = list(
    initialize = function(...) {
      self$type = "constant"
      
      args <- list(...)
      stopifnot(all(sapply(args, function(x) is.Param(x))))
      
      self$parameters <- args
    },
    evaluate = function(x, theta) { theta[1L] }
  ))

#' A Constant Model
#' 
#' The parameter is modeled as \eqn{= value}{= value}.
#' 
#' @param ... Ignored. Used to require named inputs.
#' @param value An 'UpdateType' object with type = "MH" or "Fixed".
#' 
#' @returns An R6 object inheriting from class Model.
#' 
#' @examples
#' 
#' const_m(param(1.0, "r0", minv = 0, maxv = Inf, 
#'               prop_scale = "id", prior = lnormp(log(2), 0.2)))
#' 
#' @export
const_m = function(...) {
  
  value <- list(...)
  
  stopifnot(
    "only 1 Param object can be provided" = length(value) == 1L && is.Param(value[[1L]]))
  
  ConstantModel$new(...)
}

#' @noRd
#' @keywords internal
is.Model <- function(x) inherits(x, "Model")

#' @describeIn sigmoid_m Print.
#' @export
print.SigmoidModel <- function(x, ...) {
  cat("Sigmoid Model with parameter: \n")
  lapply(x$parameters, print)
}

#' @describeIn linear_m Print.
#' @export
print.LinearModel <- function(x, ...) {
  cat("Linear Model with parameter: \n")
  lapply(x$parameters, print)
}

#' @describeIn exp_m Print.
#' @export
print.ExpModel <- function(x, ...) {
  cat("Exponential Model with parameter: \n")
  lapply(x$parameters, print)
}

#' @describeIn step_m Print.
#' @export
print.StepModel <- function(x, ...) {
  cat("Step Model with parameter: \n")
  lapply(x$parameters, print)
}

#' @describeIn const_m Print.
#' @export
print.ConstantModel <- function(x, ...) {
  cat("Parameter: \n")
  lapply(x$parameters, print)
}