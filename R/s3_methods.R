.scale2mu_weibull <- function(scale, shape) {
  scale * {gamma(1.0 + 1.0 / shape)}
}

.scale2mu_gamma <- function(scale, shape) {
  shape / scale
}

.mu2rate_weibull <- function(mu, shape) {
  {gamma(1.0 + 1.0 / shape) / mu}^shape
}

.rate2mu_weibull <- function(rate, shape) {
  {gamma(1.0 + 1.0 / shape) * rate^{-1.0 / shape}}
}

.rate2scale_weibull <- function(rate, shape) {
  rate^{-1.0 / shape}
}

.mu2rate_gamma <- function(mu, shape) {
  shape / mu
}

.rate2mu_gamma <- function(rate, shape) {
  shape / rate
}

.rate2scale_gamma <- function(rate, shape) {
  1.0 / rate
}

#'
#' @describeIn fit_baclava Summary statistics of posterior distribution parameters
#' @param object An object of class \code{baclava}.
#' @param ... Ignored.
#' @importFrom coda effectiveSize
#' @importFrom stats quantile
#' @importFrom tibble as_tibble
#' @export
summary.baclava <- function(object, ..., burnin = 0L, thin = 0L) {
  
  .summary <- function(x) {
    c("mean" = mean(x),
      "q_low" = stats::quantile(x, probs = 0.025),
      "q_high" = stats::quantile(x, probs = 0.975),
      "ESS" = coda::effectiveSize(x))
  }
  
  object$theta <- lapply(
    seq_along(object$theta),
    function(x) {
      lapply(
        seq_along(object$theta[[x]]),
        function(y) {
          colnames(object$theta[[x]][[y]]) <- paste(names(object$theta)[x],
                                                    names(object$theta[[x]])[y], 
                                                    colnames(object$theta[[x]][[y]]), sep = ".")
          object$theta[[x]][[y]]
        })
    })
  
  object$theta <- lapply(object$theta, do.call, what = cbind)
  
  object$theta <- do.call(cbind, object$theta)
  
  n <- nrow(object$theta)
  if (burnin >= n) {
    stop("burnin is too large; number of posterior sets = ", n, call. = FALSE)
  }
  burnin <- max(burnin, 0L)
  
  seq_kept <- seq_len(n)
  seq_kept <- seq_kept[seq(burnin + 1L, n, max(1L, thin))]
  
  object$theta <- object$theta[seq_kept, , drop = FALSE]
  
  res <- apply(object$theta, 2L, .summary)
  
  tibble::as_tibble(t(res), rownames = "param")
}

#'
#' @describeIn fit_baclava Print summary statistics of posterior distribution parameters
#' @param x An object of class \code{baclava}.
#' @export
print.baclava <- function(x, ...) {
  print(summary.baclava(object = x, ...))
  invisible(x)
}

#' Plot Posterior Distribution Parameters
#'
#' Convenience function to facilitate exploration of posterior distributions through
#'   trace plots, autocorrelations, and densities.
#'
#' @param x An object of class \code{baclava}.
#' @param y Ignored
#' @param type A character object. One of \{"density", "trace", "acf"\}. The
#'   type of plot to generate
#' @param burnin An integer object. Optional. The number of burn-in samples.
#'   Used only for \code{type = "trace"}. One trace plot is generated for
#'   the burnin iterations; a second for the post-burnin iterations. Note, this
#'   refers to the number of kept (thinned) samples.
#' @param thin An integer object. If non-zero posterior will be thinned to
#'   n, 2n, 3n, ...
#' @param trace_var A character object. The parameter for which trace plots
#'   are to be generated.
#' @param ... Ignored
#' 
#' @returns A gg object
#' 
#' @examples
#' data(screen_data)
#'
#' # H compartment model
#' 
#' # fitted parameter is the mean of the distribution
#' mu_H <- param(100.0, tag = "mu_H", minv = 0.0, maxv = Inf,
#'               prior = lnorm_p(log(100), 1))
#'               
#' # fixed shape parameter
#' sh_H <- param(2.0, tag = "sh_H")
#' 
#' # Weibull distribution
#' H <- compartment(weibull_d(mean = mu_H, shape = sh_H))
#'
#' # P compartment model
#' 
#' # fitted parameter is the mean of the distribution
#' mu_P <- param(2, tag = "mu_P", minv = 0.0, maxv = Inf,
#'               prior = lnorm_p(log(2), 1))
#'               
#' # fixed shape parameter
#' sh_P <- param(2.0, tag = "sh_P")
#' 
#' # Weibull distribution
#' P <- compartment(weibull_d(mean = mu_P, shape = sh_P))
#'
#' # psi 
#' # probability parameter
#' p_psi <- param(0.2, tag = "p_psi", minv = 0.0, maxv = 1.0, 
#'                prior = beta_p(1.0, 1.0))
#'
#' psi <- compartment(bern_d(p = p_psi))
#' 
#' # beta 
#' p_b <- param(0.8, tag = "p_b", minv = 0.0, maxv = 1.0, prior = beta_p(1, 1))
#' beta <- compartment(bern_d(p = p_b))
#' 
#' #Sampler strategies
#' strategy_H <- MH(params = c("mu_H"), s0 = 0.014, 
#'                  adaptive = adapt(learn = "none"))
#' 
#' strategy_P <- MH(params = c("mu_P"), s0 = 0.06, 
#'                  adaptive = adapt(learn = "none"))
#' 
#' strategy_psi <- MH(params = c("p_psi"), s0 = 0.1,
#'                    adaptive = adapt(learn = "none"))
#' blocks <- list(strategy_H, strategy_P, Gibbs("p_b"), strategy_psi)
#' 
#' # This is for illustration only -- the number of iterations should be
#' # significantly larger and the epsilon values should be tuned.
#' example <- fit_baclava(data.assess = data.screen,
#'                        data.clinical = data.clinical,
#'                        t0 = 30.0,
#'                        H = H, P = P, beta = beta, psi = psi,
#'                        blocks = blocks)
#'                        
#' plot(example)
#' plot(example, type = "trace", trace_var = "p_psi", burnin = 0L) 
#' plot(example, type = "trace", trace_var = "mu_H", burnin = 0L) 
#' plot(example, type = "trace", trace_var = "mu_P", burnin = 0L) 
#' plot(example, type = "trace", trace_var = "p_b", burnin = 0L) 
#' plot(example, type = "acf")
#' 
#' @import ggplot2
#' @importFrom stats acf quantile
#' @export
plot.baclava <- function(x, y, ..., 
                         type = c("density", "trace", "acf"), 
                         burnin = 0L, thin = 0L,
                         trace_var = NULL, scale = c("natural", "proposal")) {
  
  scale <- match.arg(scale)
  
  x$model$blocks <- NULL
  
  type <- match.arg(type)

  xtheta <- x$theta
  xtheta <- lapply(
    seq_along(xtheta),
    function(x1) {
      lapply(
        seq_along(xtheta[[x1]]),
        function(y) {
          colnames(xtheta[[x1]][[y]]) <- paste(names(xtheta)[x1],
                                               names(xtheta[[x1]])[y], 
                                               colnames(xtheta[[x1]][[y]]), sep = ".")
          xtheta[[x1]][[y]]
        })
    })
  xtheta <- lapply(xtheta, do.call, what = cbind)
  xtheta <- do.call(cbind, xtheta)

  alt_names <- lapply(
    seq_along(x$theta),
    function(x1) {
      lapply(
        seq_along(x$theta[[x1]]),
        function(y) {
          colnames(x$theta[[x1]][[y]])
        })
    }) |> unlist()
  
  if (scale == "proposal" && type != "density") {
    for (i in seq_along(x$theta)) {
      mods <- x$model[[names(x$theta)[i]]]
      for (j in seq_along(x$theta[[i]])) {
        mod <- mods$distributions[[j]]
        for (k in seq_along(mod$parameters)) {
          param <- mod$parameters[[k]]
          tos <- switch(
            param$parameters$value$prop_scale,
            "log_pos" = function(x) {log(x)},
            "logit" = function(x) {
              minv = param$parameters$value$minv
              maxv = param$parameters$value$maxv
              tmp = (x - minv) / (maxv - minv)
              log(tmp) - log1p(-tmp)
            },
            "log_shifted" = function(x) log(x -  param$parameters$value$minv),
            "log_flipped" = function(x) log(param$parameters$value$maxv - x)
          )
          thet <- x$theta[[i]][[names(mods$distributions)[j]]][, unlist(param$getTags())]
          idx <- match(unlist(param$getTags()), alt_names)
          xtheta[, idx] <- tos(thet)
        }
      }
    }
  }
  
  n <- nrow(xtheta)
  if (burnin >= n) {
    stop("burnin is too large; number of posterior sets = ", n, call. = FALSE)
  }
  burnin <- max(burnin, 0L)

  if (type == "density") {
    
    seq_kept <- seq_len(n)
    seq_kept <- seq_kept[seq(burnin + 1L, n, max(1L, thin))]
    
    xtheta <- xtheta[seq_kept, , drop = FALSE]
    
    dist_func <- lapply(x$model, function(y) y$getPriorDist()) |> unlist()
    dfuncs <- lapply(x$model, 
                     function(y) lapply(y$distributions$all$parameters,
                                        function(z) z$parameters$value$prior$dfunc)) |> unlist()
    qfuncs <- lapply(x$model, 
                     function(y) lapply(y$distributions$all$parameters,
                                        function(z) z$parameters$value$prior$qfunc)) |> unlist()
    tags <- lapply(x$model, function(y) y$getTags()) |> unlist()
    
    tags_1 <- sapply(names(tags), function(x) {
      tmp <- strsplit(x, ".", fixed = TRUE)[[1]]
      tmp <- tmp[seq_len(min(length(tmp), 2))]
      paste(tmp[!{tmp %in% c("shape", "value")}], collapse=".")
    })

    tags <- paste(tags_1, tags, sep = ".")

    # max_values <- apply(xtheta, 2L, function(x) { v <- density(x); max(v$y)} )

    prior_values <- list()
    for (i in seq_along(dfuncs)) {
      d <- dist_func[i]
      if (d == "no_prior") next
      range <- do.call(qfuncs[[i]], list(p = c(0.025, 0.975)))
      xlow <- min(0, min(xtheta[,tags[i]]))*0.9
      xhigh <- max(xtheta[,tags[i]])*1.1

      prior_values[[i]] <- data.frame(x = seq(xlow, xhigh, length.out = 100),
                                      d = do.call(dfuncs[[i]], c(x = list(seq(xlow, xhigh, length.out = 100)))),
                                      param = tags[i])
      # normalize prior to match scale of posterior density
      #prior_values[[i]]$d <- prior_values[[i]]$d * max_values[tags[i]] / max(prior_values[[i]]$d)
    }
    prior_values <- do.call(rbind, prior_values)
    
    # extract all elements as a data.frame
    res <- data.frame("param" = rep(colnames(xtheta), each = nrow(xtheta)),
                      "value" = c(xtheta))
    
    # get summary statistics for each parameter
    df_summary <- data.frame(summary(x, burnin = burnin, thin = thin))
    
    gg <- ggplot() +
      geom_density(data = res, aes(.data$value)) +
      ggforce::facet_wrap_paginate(.~param, scales = "free", ncol = 2L, nrow = 3L) +
      labs(x = "") +
      geom_vline(data = df_summary, aes(xintercept = mean, group = .data$param), colour = "maroon") +
      geom_vline(data = df_summary, aes(xintercept = .data$`q_low.2.5.`, group = .data$param), colour = "maroon") +
      geom_vline(data = df_summary, aes(xintercept = .data$`q_high.97.5.`, group = .data$param), colour = "maroon") +
      geom_line(data = prior_values, aes(x, d), linetype = "dashed", color = "grey") +
      theme_bw()
    
    n <- ggforce::n_pages(gg)
    
    for (i in seq_len(n)) {
      print(gg + 
              ggforce::facet_wrap_paginate(.~param, scales = "free", ncol = 2L, nrow = 3L, page = i))
      
    }
    
  } else if (type == "trace") {
    thin <- max(1L, thin)
    seq_kept <- seq_len(n)
    seq_kept <- seq_kept[c(seq_len(burnin), seq(burnin + thin, n, thin))]
    xtheta <- xtheta[seq_kept, , drop = FALSE]
    
    if (is.null(trace_var)) trace_var <- colnames(xtheta)
    
    idx <- match(trace_var, colnames(xtheta))
    if (any(is.na(idx))) {
      idx <- match(trace_var, alt_names)
      if (any(is.na(idx))) {
        stop("requested parameter is not present in $theta\n", 
             "available options are: ", paste(colnames(xtheta), collapse = ", "),
             call. = FALSE)
      }
    }
    
    df <- data.frame("sample" = seq_len(nrow(xtheta)),
                     "value" = c(xtheta[, idx]),
                     "param" = rep(trace_var, each = nrow(xtheta)))
    
    df$param <- factor(df$param)
    
    df$group <- factor(c("Recurrent", "Transient")[as.integer(df$sample <= burnin) + 1L],
                       levels = c("Transient", "Recurrent"))
    
    gg <- ggplot(df) + 
      geom_line(aes(x = .data$sample, y = .data$value)) +
      labs(x = "Sample", y = "") +
      scale_x_continuous(trans = "identity") +
      ggforce::facet_wrap_paginate(~group + param, scales = "free", ncol = 2L)
    
    n <- ggforce::n_pages(gg)
    
    if (!is.null(n)) {
      for (i in seq_len(n)) {
        print(gg + 
                ggforce::facet_wrap_paginate(~group + param, scales = "free", ncol = 2L, nrow = 3L, page = i))
      }
      
    } else {
      print(gg)
    }
    
    
  } else if (type == "acf") {
    
    seq_kept <- seq_len(n)
    seq_kept <- seq_kept[seq(burnin + 1L, n, max(1L, thin))]
    
    xtheta <- xtheta[seq_kept, , drop = FALSE]
    
    my_acf <- function(x) {
      tmp <- stats::acf(x = x, plot = FALSE)
      cbind("acf" = as.vector(tmp$acf), "lag" = as.vector(tmp$lag))
    }
    
    df <- apply(xtheta, 2L, my_acf, simplify = FALSE)
    df <- do.call(rbind, lapply(names(df), function(p) {
      if (all(is.na(df[[p]][,1L]))) return(NULL)
      cbind(data.frame(df[[p]]), param = p)
    }))
    df$param <- factor(df$param)
    
    gg <- ggplot() + geom_col(data = df, aes(.data$lag, .data$acf)) + 
      scale_y_continuous(breaks = seq(0.0, 1.0, 0.5)) +
      ggforce::facet_wrap_paginate(~param, ncol = 2L, nrow = 3L) +
      labs(x = "Lags", y= "Auto-correlation\n function")
    
    n <- ggforce::n_pages(gg)
    
    for (i in seq_len(n)) {
      print(gg + 
              ggforce::facet_wrap_paginate(~param, ncol = 2L, nrow = 3L, page = i))
      
    }
    
  }
  gg
}
