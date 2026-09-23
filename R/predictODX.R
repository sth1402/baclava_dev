#' Estimate the Overall and Per Screen Overdiagnosis Rates
#' 
#' Using the posterior parameter distributions, calculates the infinite
#'   population estimates of the probability of overdiagnosis at 
#'   each screening episode due to indolence and/or death by other causes.
#' 
#' Provided birth cohort life table is an all cause tables obtained from the 
#'   CDC Life Tables. Vital Statistics of the United States, 1974 Life Tables, 
#'   Vol. II, Section 5. 1976. Estimated "other cause" mortality will thus be 
#'   overestimated when using these tables. It is recommended that user provide 
#'   data that has been corrected to exclude death due to the disease under analysis.
#'   
#' @param object An object of S3 class 'baclava'. The value object returned
#'   by \code{fit_baclava()}.
#' @param screening.schedule A numeric vector object. A vector of ages at
#'   which screenings occur.
#' @param other.cause.rates A data.frame object. Must contain columns "Rate"
#'   and "Age". 
#' @param burnin An integer object. Optional. The number of burn-in samples.
#'   Note, this refers to the samples of the kept posterior. If > 0L, the 
#'   1L:burnin posteriors are removed.
#' @param thin An integer object. Optional. If > 1, posterior samples 
#'   (after removing burnin samples) are thinned to (thin, 2*thin, 3*thin, ...) 
#' @param model A named list. If the original model assumed a subset structure,
#'   mixed preclinical compartment, or multiple beta models,
#'   this list specifies the model component to use to predict ODX.
#'   If any model components are subset specific,
#'   only 1 subset can be specified. For mixture models or multiple beta models 
#'   apply to the subset, provide names of all needed models.  
#'   Names of the provided list elements must be 'H', 'P', 'psi' and 'beta'.
#'   The options available can be seen by printing names(object$model$X$distributions) where
#'   X is one of P, H, psi, or beta. 
#' @param verbose A logical object. If TRUE, progress bars will be displayed.
#' @param n.core A scalar integer. The number of cores for parallelization (=1 is sequential)
#' @param ... Additional optional inputs.
#'   
#' @returns A list object. For each screen in \code{screening.schedule},
#'   a matrix providing the mean total overdiagnosis and the mean overdiagnosis
#'   due to indolent/progressive tumors, as well as their 95\% prediction intervals.
#'   Similarly, element \code{overall} provides these estimates for the full
#'   screening schedule.
#'
#' @examples
#' 
#' data(screen_data)
#'
#' # H compartment model
#' 
#' mu_H <- param(init = runif(1L, 50.0, 500.0), tag = "mu_H", 
#'               minv = 10.0, prop_scale = "log", 
#'               prior = lnorm_p(log(150.00), 0.4))
#' sh_H <- param(init =  runif(1L, 0.5, 5.0), tag = "sh_H", 
#'               minv =  0.01, prop_scale = "log", 
#'               prior = lnorm_p(log( 1.75), 0.4))
#' dist_H <- weibull_d(mu = mu_H, shape = sh_H)
#' strategy_H <- block(params = c("mu_H", "sh_H"), adaptive = list(warmup = 0))
#' 
#' H <- compartment(dist_H, strategy_H)
#' 
#' # beta model
#' 
#' p_b <- param(init = runif(1L, 0.0, 1.0), tag = "p_b", 
#'              minv = 0.0, maxv = 1.0, 
#'              prior = beta_p(1.0, 1.0))
#' dist_b <- bern_d(p = p_b)
#' 
#' beta <- compartment(dist_b, strategies = "Gibbs")
#' 
#' # P compartment model
#' 
#' mu_P <- param(init = runif(1L, 1.0, 10.0), tag = "mu_P", 
#'               minv = 1e-6, prop_scale = "log", 
#'               prior = lnorm_p(log(2), 0.4))
#' sh_P <- param(init = runif(1L, 0.5, 5.0), tag = "sh_P", 
#'               minv =  0.01, prop_scale = "log", 
#'               prior = lnorm_p(log(1.73), 0.28))
#' dist_P <- gamma_d(mu = mu_P, shape = sh_P)
#' strategy_P <- block(params = c("mu_P", "sh_P"), adaptive = list(warmup = 0))
#' 
#' P <- compartment(dist_P, list(strategy_P))
#' 
#' # psi model
#' 
#' p_psi <- param(init = runif(1L, 0.0, 1.0), tag = "p_psi", 
#'                minv = 0.0, maxv = 1.0, 
#'                prior = beta_p(1.0, 1.0))
#' dist_psi <- bern_d(p = p_psi)
#' strategy_psi <- block(params = c("p_psi"), adaptive = list(warmup = 0))
#' 
#' psi <- compartment(dist_psi, blockStrategies(list(strategy_psi)))
#' 
#' # This is for illustration only -- the number of iterations should be
#' # significantly larger and the epsilon values should be tuned.
#' example <- fit_baclava(data.assess = sim_screen,
#'                        data.clinical = sim_clinical,
#'                        t0 = t0,
#'                        H = H, P = P, psi = psi, beta = beta)
#'                 
#' # if rates are not available, an all cause dataset is provided in the package
#' # NOTE: these predictions will be over-estimated
#'            
#' data(all_cause_rates)
#' all_cause_rates <- all_cause_rates[, c("Age", "both")]
#' colnames(all_cause_rates) <- c("Age", "Rate")
#' 
#' # using single screen for example speed
#' predicted_odx <- predictODX(object = example, 
#'                             other.cause.rates = all_cause_rates,
#'                             screening.schedule = 40, 
#'                             burnin = 10)
#'
#' plot(predicted_odx)
#' 
#' @importFrom utils read.csv
#' @import doParallel
#' @import foreach
#' @import parallel
#' @include predictODX_helpers.R
#' @rdname predictODX
#'
#' @export
predictODX <- function(object, screening.schedule, 
                       other.cause.rates, 
                       model = list("H" = "all", "P" = "all", "psi" = "all", "beta" = "all"), 
                       burnin = 0L, thin = 1L, verbose = TRUE, 
                       n.draws = NULL,
                       n.core = max(1L, parallel::detectCores() - 1L), ...) {
  
  stopifnot(
    "`object` must be an object of class 'baclava'" = !missing(object) &&
      inherits(object, "baclava"),
    "`screening.schedule` must be a non-negative numeric vector" = !missing(screening.schedule) &&
      is.vector(screening.schedule, mode = "numeric") && length(screening.schedule) > 0L && 
      all(screening.schedule >= 0.0),
    "`other.cause.rates` must be a data.frame" = 
      {is.data.frame(other.cause.rates) && 
          all(c("Rate", "Age") %in% colnames(other.cause.rates))},
    "`model` must be a named list" = is.list(model) && !is.null(names(model)),
    "`burnin` must be a non-negative scalar" = 
      is.vector(burnin, mode = "numeric") && length(burnin) == 1L && 
      isTRUE(all.equal(burnin, round(burnin))) && burnin >= 0.0,
    "`verbose` must be a logical" = is.vector(verbose, mode = "logical") && 
      length(verbose) == 1L
  )
  
  # ensure that specified models are present in the baclava model specification
  
  if (!model$H %in% names(object$model$H$distributions)) {
    stop('specified H model not found; provided ', model$H, "; object has ",
         paste(names(object$model$H$distributions), collapse = ", "))
  }
  
  if (!all(model$P %in% names(object$model$P$distributions))) {
    stop('specified P model not found; provided ', model$P, "; object has ",
         paste(names(object$model$P$distributions), collapse = ", "))
  }
  
  if (object$setup$indolent && !model$psi %in% names(object$model$psi$distributions)) {
    stop('specified psi model not found; provided ', model$psi, "; object has ",
         paste(names(object$model$psi$distributions), collapse = ", "))
  }
  
  if (!all(model$beta %in% names(object$model$beta$distributions))) {
    stop('specified beta model not found; provided ', model$beta, "; object has ",
         paste(names(object$model$beta$distributions), collapse = ", "))
  }
  
  # retrieve required information from the analysis object.
  risk_onset <- object$setup$t0

  # sort screening schedule and ensure that other causes age range spans the schedule
  screening.schedule <- sort(screening.schedule)
  if (min(screening.schedule) < min(other.cause.rates$Age))
    stop("'other.cause.rate' must include minimum screening age", call. = FALSE)
  
  if (max(screening.schedule) > max(other.cause.rates$Age)) {
    message("Extended other-cause hazard table to age ", max(screening.schedule), 
            " (Gompertz extrapolation).")
    other.cause.rates <- extend_other_causes(other.cause.rates, 
                                             max_age = max(screening.schedule))
  }
  
  # total number of screens in schedule
  max_n_screens <- length(screening.schedule)

  keptM <- seq_len(nrow(object$theta$H[[1L]]))
  if (burnin > 0L) keptM <- keptM[-seq_len(burnin)]
  
  if (length(keptM) == 0L) stop("burnin removed all posterior draws", call. = FALSE)
  if (thin > 1L) keptM <- keptM[seq(thin, length(keptM), thin)]
  
  if (is.null(n.draws)) n.draws <- length(keptM)
  
  # number of posterior draws to sample
  draws <- keptM
  if (n.draws > length(draws)) {
    warning("requested # of posterior draws > # available after burnin + thinning; duplicate sampling")
  }
  if (n.draws != length(draws)) {
    draws <- sample(draws, n.draws, replace = n.draws > length(draws)) |> sort()
  }
  
  posteriors <- do.call(cbind, object$theta$H)
  posteriors <- cbind(posteriors, do.call(cbind, object$theta$P))
  posteriors <- cbind(posteriors, do.call(cbind, object$theta$beta))
  posteriors <- cbind(posteriors, do.call(cbind, object$theta$psi))
  
  posteriors <- posteriors[draws, , drop = FALSE]
  
  result <- predict_ODX(object$model$H$toList(), object$model$P$toList(), 
                        object$model$psi$toList(), object$model$beta$toList(), object$setup$t0, 
                        screening.schedule,
                        other.cause.rates$Age,
                        other.cause.rates$Rate,
                        posteriors,
                        colnames(posteriors))
  
  res <- list()
  for (i in seq_len(nrow(result$ODX_ind_stat))) {
    res[[i]] <- wrap(result$ODX_ind_stat[i, , drop = FALSE], 
                     result$ODX_prog_stat[i, , drop = FALSE], 
                     result$SD_stat[i, , drop = FALSE])
  }
  names(res) <- paste0("screen.age.", screening.schedule)

  # sum across screens to get overall; should be means, but they cancel
  res$overall <- wrap(colMeans(result$ODX_ind_stat), 
                      colMeans(result$ODX_prog_stat), 
                      colMeans(result$SD_stat))
  
  res$raw.data <- list(
    "ind" = result$ODX_ind_stat[, , drop = FALSE], 
    "prog" = result$ODX_prog_stat[, , drop = FALSE], 
    "stat" = result$SD_stat[, , drop = FALSE]
  )
  
  res$screening.schedule <- screening.schedule
  
  class(res) <- c(class(res), "baclava.ODX.pred")
  
  res
}

#' @param x A an object of S3 class 'baclava.ODX.pred' as returned by \code{predictODX()}.
#' @param y Ignored.
#' @param ... Ignored.
#'
#' @describeIn predictODX Generate column plot of predicted overdiagnosis for each screen.
#'
#' @import ggplot2
#' @export
plot.baclava.ODX.pred <- function(x, y, ...) {
  
  ss <- x[["screening.schedule"]]
  x[["screening.schedule"]] <- NULL
  x[["raw.data"]] <- NULL
  
  n_screens <- length(x) - 1L
  x[["overall"]] <- NULL

  
  df <- do.call(rbind, x) |> as.data.frame()
  colnames(df) <- c("mean", "low", "high", "mid")
  df$screen <- rep(seq_len(n_screens), each = 3L)
  df$screen <- ss[df$screen]
  df$category <- rep(rownames(x[[1L]]), times = n_screens)
  
  df_sub <- subset(df, df$category != "total")
  df_sub$category <- factor(df_sub$category, levels = c("mortality", "indolent"))
  
  df_total <- subset(df, df$category == "total")

  gg <- ggplot(df_sub, aes(x = .data$screen, y = .data$mean)) + 
    geom_col(aes(fill = .data$category)) +
#    geom_col(data = df_total, aes(x = .data$screen, y = .data$mean), fill = NA, color = "gray72") +
    geom_errorbar(data = df_total, aes(x = .data$screen, ymin = .data$low, ymax = .data$high), width = 0.2,
                  position=position_dodge(0.9), color = "gray72") +
    ggtitle("Mean Predicted Overdiagnosis for Each Screen") +
    xlab("Age at Screen") + ylab("Overdiagnosis Rate, %") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
  gg
}

#' @describeIn predictODX print method
#' @export
print.baclava.ODX.pred <- function(x) {
  print(round(x$overall, digits = 1))
  invisible(x)
}
