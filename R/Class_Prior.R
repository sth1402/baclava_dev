#' R6 Class for Prior Distribution Specification
#'
#' Represents a prior distribution with its associated density name and 
#' numeric parameters. Includes validation to ensure that supported 
#' distributions are called with appropriate arguments.
#'
#' @noRd
#' 
#' @return An R6 object of class \code{Prior}.
#'
#' @examples
#' # Normal prior with mean 0 and sd 1
#' p1 <- Prior$new("dnorm", c(0, 1))
#' p1$toList()
#'
#' # Gamma prior with shape 2, scale 0.5
#' p2 <- Prior$new("dgamma", c(2, 0.5))
#'
#' # Attempting to specify an unsupported prior triggers validation
#' \dontrun{
#' Prior$new("unsupported", c(1, 2))
#' }
#'
#' @import R6
#' @keywords internal
Prior <- R6Class(
  "Prior",
  public = list(
    #' @field dist (`character(1)`) \cr 
    #' Name of the distribution (e.g., \code{"dnorm"}, \code{"dgamma"}).
    dist = NULL,
    #' @field params (`Numeric vector`)\cr
    #' Parameters for the specified distribution.
    params = NULL,
    
    #' @description
    #' Checks that the distribution name and parameters conform to the 
    #'   expected shape and range. Called automatically on initialization.
    validate = function() {
      
      if (self$dist == "no_prior") return(invisible(TRUE))
      
      if (!self$dist %in% names(private$prior_spec_schema)) {
        msg <- sprintf("Unknown prior distribution '%s'", self$dist)
        stop(msg, call. = FALSE)
      }
      
      spec <- private$prior_spec_schema[[self$dist]]
      
      if (length(self$params) != spec$n) {
        msg <- sprintf("Prior '%s' expects %d parameters, got %d", 
                       self$dist, spec$n, length(self$params))
        stop(msg, call. = FALSE)
      }
      
      if (!spec$validate(self$params)) {
        msg <- sprintf("Validation failed for prior '%s' with params: %s", 
                       self$dist, paste(self$params, collapse = ", "))
        stop(msg, call. = FALSE)
      }
      names(self$params) <- private$prior_spec_schema[[self$dist]][["names"]]
    },
    
    #' @description
    #' Generates a string showing the 95% coverage area.
    ci_range = function() {
      if (!self$dist %in% names(private$prior_spec_schema)) return(invisible(NULL))
      
      spec <- private$prior_spec_schema[[self$dist]]
      interval <- spec$ci_range(self$params)
      if (!is.null(interval)) {
        sprintf("95%% prior interval [%.6f, %.6f]", 
                interval[1L], interval[2L])
      }
    },
    
    #' @description
    #' Constructor. Validates that the distribution is supported and that 
    #' parameter values meet distribution-specific criteria.
    #' @param dist A character. The prior distribution.
    #' @param params A numeric vector. The prior parameters.
    initialize = function(dist, params = NULL) {
      stopifnot(
        "`dist` must be a character" = !missing(dist) && 
          is.vector(dist, "character") && length(dist) == 1L,
        "all parameters must be numeric" = is.null(params) || 
          {is.vector(params, "numeric") && length(params) > 0L &&
              all(is.finite(params))}
      )
      self$dist <- dist
      self$params <- params
      self$validate()
    },
    
    #' @description
    #' Displays a formatted summary of the prior distribution and its parameters.
    #' @param ... ignored.
    print = function(...) {
      if (!is.null(self$ci_range())) {
        fmt <- "%10s(x,%s) [%s] \n"
        cat(sprintf(fmt, self$dist, paste(paste(names(self$params), self$params, sep = "="), collapse= ","), self$ci_range()))
      } else {
        fmt <- "%10s(x,%s) \n"
        cat(sprintf(fmt, self$dist, paste(paste(names(self$params), self$params, sep = "="), collapse= ",")))
      }
    },
    
    #' @description
    #' Returns a formatted character string of the prior distribution and its parameters.
    #' @param ... ignored.
    asString = function(...) {
      fmt <- "%10s(%s) [%s] \n"
      rng <- self$ci_range()
      if (is.null(rng)) 
        sprintf("%10s(%s)\n", self$dist, paste(self$params, collapse = ", "))
      else
        sprintf("%10s(%s) [%s]\n", self$dist, paste(self$params, collapse = ", "), rng)
    },
    
    #' @description
    #' Returns a named list containing the distribution name and parameter vector.
    toList = function() {
      if (self$dist != "no_prior") {
        list("meta_for" = "Prior",
             "dist" = self$dist, 
             "params" = self$params)
      } else {
        list("meta_for" = "Prior",
             "dist" = self$dist)
      }
    },
    
    #' @description
    #' Returns the distribution name.
    getDist = function() self$dist,
    #' @description
    #' Returns the distribution parameters.
    getParams = function() self$params,
    dfunc = function(x, log = FALSE) {
      if (self$dist == "no_prior") return(NA_real_)
      if (!self$dist %in% c("dinvgamma", "dhalfnorm")) {
        args <- as.list(self$params)
        args$x <- x
        args$log <- log
        out <- do.call(self$dist, args)
      } else if (self$dist == "dinvgamma") {
        logdens <- - self$params["shape"] * log(self$params["scale"]) -
            lgamma(self$params["shape"]) -
            (self$params["shape"] + 1.0) * log(x) -
            1.0 / (self$params["scale"] * x);
          
        if (any(x <= 0)) {
          logdens[x <= 0] <- -Inf
        }
        out <- if (log) { logdens } else { exp(logdens) }
      } else if (self$dist == "dhalfnorm") {
        args <- as.list(self$params)
        args$x <- x
        args$log <- log
        dens <- do.call(stats::dnorm, args)
        if (any(x < 0))  {
          dens[x < 0] <- if (log) -Inf else 0.0
        }
        out <- if (log) { dens + log(2.0) } else { 2.0 * dens }
      }
      out
    },
    qfunc = function(p, lower.tail = TRUE, log.p = FALSE) {
      if (self$dist == "no_prior") return(NA_real_)
      if (!self$dist %in% c("dinvgamma", "dhalfnorm")) {
        args <- as.list(self$params)
        args$p <- p
        args$lower.tail <- lower.tail
        args$log.p <- log.p
        func <- self$dist
        substr(func,1,1) <- "q"
        out <- do.call(func, args)
      } else if (self$dist == "dinvgamma") {
        if (log.p) p <- exp(p)
        if (!lower.tail) p <- 1.0 - p
        out <- 1.0 / qgamma(1 - p, shape = self$params["shape"], scale = self$params["scale"])
      } else if (self$dist == "dhalfnorm") {
        if (log.p) p <- exp(p)
        if (!lower.tail) p <- 1 - p
        out <- self$params["sd"] * stats::qnorm({p + 1.0} / 2.0)
      }
      out
    }), 
  
  private = list(
    prior_spec_schema = list(
      dunif   = list(n = 2, names = c("min", "max"), 
                     validate = function(x) x[1] < x[2], 
                     ci_range = function(x) qunif(c(0.025, 0.975), min=x[1], max=x[2])),
      dnorm   = list(n = 2, names = c("mean", "sd"), 
                     validate = function(x) x[2] > 0, 
                     ci_range = function(x) qnorm(c(0.025, 0.975), mean=x[1], sd=x[2])),
      dhalfnorm = list(n = 1, names = "sd", validate = function(x) x > 0,
                       ci_range = function(x) { c(x * qnorm((1.025) / 2), x * qnorm((1.975) / 2))}),
      dlnorm  = list(n = 2, names = c("meanlog", "sdlog"), 
                     validate = function(x) x[2] > 0, 
                     ci_range = function(x) qlnorm(c(0.025, 0.975), meanlog=x[1], sdlog=x[2])),
      dgamma  = list(n = 2, names = c("shape", "scale"), 
                     validate = function(x) x[1] > 0 && x[2] > 0, 
                     ci_range = function(x) qgamma(c(0.025, 0.975), shape=x[1], scale=x[2])),
      dinvgamma  = list(n = 2, names = c("shape", "scale"), 
                     validate = function(x) x[1] > 0 && x[2] > 0, 
                     ci_range = function(x) {
                       1.0 / qgamma(c(0.975, 0.025), shape=x[1], scale=x[2])
                     }),
      dbeta   = list(n = 2, names = c("shape1", "shape2"), 
                     validate = function(x) all(x > 0), 
                     ci_range = function(x) qbeta(c(0.025, 0.975), shape1=x[1], shape2=x[2])),
      dnbeta  = list(n = 3, names = c("shape1", "shape2", "ncp"), 
                     validate = function(x) all(x[1:2] > 0), 
                     ci_range = function(x) NULL),
      dchisq  = list(n = 1, names = c("df"), 
                     validate = function(x) x[1] > 0, 
                     ci_range = function(x) qchisq(c(0.025, 0.975), df =x[1], ncp = 0)),
      dnchisq = list(n = 2, names = c("df", "ncp"), 
                     validate = function(x) x[1] > 0, 
                     ci_range = function(x) NULL),
      dt      = list(n = 1, names = c("df"), 
                     validate = function(x) x[1] > 0, 
                     ci_range = function(x) qt(c(0.025, 0.975), df=x[1])),
      dnt     = list(n = 2, names = c("df", "ncp"), 
                     validate = function(x) x[1] > 0, 
                     ci_range = function(x) NULL),
      df      = list(n = 2, names = c("df1", "df2"), 
                     validate = function(x) all(x > 0), 
                     ci_range = function(x) qf(c(0.025, 0.975), df1=x[1], df2=x[2])),
      dnf     = list(n = 3, names = c("df1", "df2", "ncp"), 
                     validate = function(x) all(x[1:2] > 0), 
                     ci_range = function(x) NULL),
      dcauchy = list(n = 2, names = c("location", "scale"), 
                     validate = function(x) x[2] > 0, 
                     ci_range = function(x) qcauchy(c(0.025, 0.975), location=x[1], scale=x[2])),
      dlogis  = list(n = 2, names = c("location", "scale"), 
                     validate = function(x) x[2] > 0, 
                     ci_range = function(x) qlogis(c(0.025, 0.975), location=x[1], scale=x[2])),
      dexp    = list(n = 1, names = c("rate"), 
                     validate = function(x) x[1] > 0, 
                     ci_range = function(x) qexp(c(0.025, 0.975), rate=x[1])),
      dweibull= list(n = 2, names = c("shape", "scale"), 
                     validate = function(x) all(x > 0), 
                     ci_range = function(x) qweibull(c(0.025, 0.975), shape=x[1], scale=x[2])),
      dbinom  = list(n = 2, names = c("size", "prob"), 
                     validate = function(x) x[1] >= 0 && x[2] >= 0 && x[2] <= 1, 
                     ci_range = function(x) qbinom(c(0.025, 0.975), size=x[1], prob=x[2])),
      no_prior= list(n = 0, names = character(0), 
                     validate = function(x) TRUE, 
                     ci_range = function(x) NULL)
    )
  )
)

#' @noRd
#' @keywords internal
Prior$fromList <- function(lst) {
  if (lst$meta_for != "Prior") stop("must provide a Prior list")
  Prior$new(lst$dist, lst$params)
}

#' Prior Specifications
#' 
#' Convenience functions for specifying the prior distributions and their
#'   parameter settings for each model parameter.
#'   
#' @param min,max The lower and upper limits of the distribution. Must be finite.
#'
#' @return An S3 object of class `Prior`, which extends list. Element `dist`
#'   is the character name of the density function and `params` is a numeric
#'   vector of the distribution parameters.
#'   
#' @examples
#' 
#' # prior is dunif(x, min = 0, max = 10)
#' unif_p(0, 10)
#' @export
#' @rdname priors
#' @name priors
unif_p <- function(min = 0.0, max = 1.0) {
  Prior$new("dunif", c(min, max))
}

#' @describeIn priors Normal Distribution
#' @param mean,sd The mean and standard deviation of the distribution.
#' @examples
#' 
#' # prior is dnorm(x, mean = 0.2, sd = 1.5)
#' norm_p(0.2, 1.5)
#' @export
norm_p <- function(mean = 0.0, sd = 1.0) {
  Prior$new("dnorm", c(mean, sd))
}

#' @describeIn priors Half-Normal Distribution
#' @param sd The standard deviation of the distribution.
#' @examples
#' 
#' # prior is abs(dnorm(x, mean = 0.0, sd = 1.5))
#' halfnorm_p(1.5)
#' @export
halfnorm_p <- function(sd = 1.0) {
  Prior$new("dhalfnorm", sd)
}

#' @describeIn priors Log-Normal Distribution
#' @param meanlog,sdlog The mean and standard deviation of the distribution on the log scale.
#' @examples
#' 
#' # prior is dlnorm(x, meanlog = 1, sdlog = 1.5)
#' lnorm_p(1, 1.5)
#' @export
lnorm_p <- function(meanlog = 0.0, sdlog = 1.0) {
  Prior$new("dlnorm", c(meanlog, sdlog))
}

#' @describeIn priors Gamma Distribution
#' @param shape,scale The shape and scale parameters. Must be positive, scale strictly.
#' @examples
#' 
#' # prior is dgamma(x, shape = 1, scale = 1.5)
#' gamma_p(1, 1.5)
#' @export
gamma_p <- function(shape, scale = 1.0) {
  stopifnot("`shape` must be specified" = !missing(shape))
  Prior$new("dgamma", c(shape, scale))
}

#' @describeIn priors Inverse-Gamma Distribution
#' @param shape,scale The shape and scale parameters. Must be positive, scale strictly.
#' @examples
#' 
#' # prior is dinvgamma(x, shape = 1, scale = 1.5)
#' invgamma_p(1, 1.5)
#' @export
invgamma_p <- function(shape, scale = 1.0) {
  stopifnot("`shape` must be specified" = !missing(shape))
  Prior$new("dinvgamma", c(shape, scale))
}

#' @describeIn priors Beta Distribution
#' @param shape1,shape2 The non-negative parameters.
#' @examples
#' 
#' # prior is dbeta(x, shape1 = 1, shape2 = 2)
#' beta_p(1, 2)
#' @export
beta_p <- function(shape1, shape2) {
  stopifnot("`shape1` and `shape2` must be specified" = !missing(shape1) && !missing(shape2))
  Prior$new("dbeta", c(shape1, shape2))
}

#' @describeIn priors Noncentral Beta Distribution
#' @param ncp The non-centrality parameter.
#' @examples
#' 
#' # prior is dbeta(x, shape1 = 1, shape2 = 1.5, ncp = 0.4)
#' nbeta_p(1, 1.5, 0.4)
#' @export
nbeta_p <- function(shape1, shape2, ncp) {
  stopifnot("`shape1`, `shape2`, and `ncp` must be specified" = !missing(shape1) 
            && !missing(shape2) && !missing(ncp))
  Prior$new("dnbeta", c(shape1, shape2, ncp))
}

#' @describeIn priors Chi-squared Distribution
#' @param df The degrees of freedom; must be non-negative.
#' @examples
#' 
#' # prior is dchisq(x, df = 2)
#' chisq_p(2)
#' @export
chisq_p <- function(df) {
  stopifnot("`df` must be specified" = !missing(df))
  Prior$new("dchisq", c(df))
}

#' @describeIn priors Noncentral Chi-squared Distribution
#' @examples
#' 
#' # prior is dchisq(x, df = 2, ncp = 0.4)
#' nchisq_p(2, 0.4)
#' @export
nchisq_p <- function(df, ncp) {
  stopifnot("`df` and `ncp` must be specified" = !missing(df) && !missing(ncp))
  Prior$new("dnchisq", c(df, ncp))
}

#' @describeIn priors T Distribution
#' @examples
#' 
#' # prior is dt(x, df = 2)
#' t_p(2)
#' @export
t_p <- function(df) {
  stopifnot("`df` must be specified" = !missing(df))
  Prior$new("dt", c(df))
}

#' @describeIn priors Noncentral T Distribution
#' @examples
#' 
#' # prior is dt(x, df = 2, ncp = 0.4)
#' nt_p(2, 0.4)
#' @export
nt_p <- function(df, ncp) {
  stopifnot("`df` and `ncp` must be specified" = !missing(df) && !missing(ncp))
  Prior$new("dnt", c(df, ncp))
}

#' @describeIn priors F Distribution
#' @param df1,df2 The degrees of freedom; must be non-negative.
#' @examples
#' 
#' # prior is df(x, df1 = 2, df2 = 4)
#' f_p(2, 4)
#' @export
f_p <- function(df1, df2) {
  stopifnot("`df1` and `df2` must be specified" = !missing(df1) && !missing(df2))
  Prior$new("df", c(df1, df2))
}

#' @describeIn priors Noncentral F Distribution
#' @examples
#' 
#' # prior is df(x, df1 = 2, df2 = 4, ncp = 0.4)
#' nf_p(2, 4, 0.4)
#' @export
nf_p <- function(df1, df2, ncp) {
  stopifnot("`df1`, `df2`, and `ncp` must be specified" = !missing(df1) && 
              !missing(df2) && !missing(ncp))
  Prior$new("dnf", c(df1, df2, ncp))
}

#' @describeIn priors Cauchy Distribution
#' @param location,scale The location and scale parameters.
#' @examples
#' 
#' # prior is dcauchy(location = 0.4, scale = 1.4)
#' cauchy_p(0.4, 1.4)
#' @export
cauchy_p <- function(location = 0.0, scale = 1.0) {
  Prior$new("dcauchy", c(location, scale))
}

#' @describeIn priors Exponential Distribution
#' @param rate The rate parameter.
#' @examples
#' 
#' # prior is dexp(rate = 1.4)
#' exp_p(1.4)
#' @export
exp_p <- function(rate = 1.0) {
  Prior$new("dexp", c(rate))
}

#' @describeIn priors Logistic Distribution
#' @examples
#' 
#' # prior is dlogis(location = 0.4, scale = 1.4)
#' logis_p(0.4, 1.4)
#' @export
logis_p <- function(location = 0.0, scale = 1.0) {
  Prior$new("dlogis", c(location, scale))
}

#' @describeIn priors Weibull Distribution
#' @examples
#' 
#' # prior is dweibull(x, shape = 0.4, scale = 1.4)
#' weibull_p(0.4, 1.4)
#' @export
weibull_p <- function(shape, scale = 1.0) {
  stopifnot("`shape` must be specified" = !missing(shape))
  Prior$new("dweibull", c(shape, scale))
}

#' @describeIn priors Binomial Distribution
#' @examples
#' 
#' # prior is dbinom(x, size = 1, prob = 0.5)
#' binom_p(1, 0.5)
#' @export
binom_p <- function(size, prob) {
  stopifnot("`size` and `prob` must be specified" = !missing(size) && !missing(prob))
  Prior$new("dbinom", c(size, prob))
}

#' @describeIn priors General Constructor
#' @examples
#' 
#' # prior is dbinom(x, 1, 0.5)
#' prior_p("dbinom", 1, 0.5)
#' 
#' # prior is dweibull(x, shape = 0.4, scale = 1.4)
#' prior_p("dweibull", 0.4, 1.4)
#' @export
prior_p <- function(dist, ...) {
  params <- unlist(list(...), use.names = FALSE)
  if (length(params) == 0L) stop("must supply at least one parameter")
  Prior$new(dist, params)
}

#' @noRd
#' @keywords internal
none_p <- function() {
  Prior$new("no_prior")
}

is.Prior <- function(x) inherits(x, "Prior")

#' @describeIn priors Print
#' @export 
print.Prior <- function(x, ...) x$print()