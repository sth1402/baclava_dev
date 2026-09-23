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
#'   Used only for \code{type = "trace"}. One trace plot is generated for
#'   the burnin iterations; a second for the post-burnin iterations. Note, this
#'   refers to the kept (thinned) samples.
#' @param groups.rateP An integer scalar object. If model included groups with
#'   different sojourn parameters, the group for which overdiagnosis is to
#'   be estimated. Must be one of \code{object$setup$groups.rateP}
#' @param screen.type An integer scalar object. If model included screen-type,
#'   specific sensitivity parameters, the screen-type for which
#'   overdiagnosis is to be estimated. Must be one of \code{object$setup$groups.beta}
#' @param verbose A logical object. If TRUE, progress bars will be displayed.
#'   
#' @returns A list object. For each screen in \code{screening.schedule},
#'   a matrix providing the mean total overdiagnosis and the mean overdiagnosis
#'   due to indolent/progressive tumors, as well as their 95% prediction intervals.
#'   Similarly, element \code{overall} provides these estimates for the full
#'   screening schedule.
#'
#' @examples
#' 
#' data(screen_data)
#' 
#' theta_0 <- list("rate_H" = 7e-4, "shape_H" = 2.0,
#'                 "rate_P" = 0.5  , "shape_P" = 1.0,
#'                 "beta" = 0.9, psi = 0.4)
#' prior <- list("rate_H" = 0.01, "shape_H" = 1,
#'               "rate_P" = 0.01, "shape_P" = 1,
#'               "a_psi" = 1/2 , "b_psi" = 1/2,
#'               "a_beta" = 38.5, "b_beta" = 5.8)
#'
#' # This is for illustration only -- the number of MCMC samples should be
#' # significantly larger and the epsilon values should be tuned.
#' example <- fit_baclava(data.assess = data.screen,
#'                        data.clinical = data.clinical,
#'                        t0 = 30.0,
#'                        theta_0 = theta_0,
#'                        prior = prior)
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
#' @include predictODX_helpers.R
#' @rdname predictODX
#'
#' @export
predictODX_old <- function(object, screening.schedule, 
                           other.cause.rates, n.draws = 1000L,
                           burnin = 1000L, verbose = TRUE) {
  
  stopifnot(
    "`object` must be an object of class 'baclava'" = !missing(object) &&
      inherits(object, "baclava"),
    "`screening.schedule` must be a non-negative numeric vector" = !missing(screening.schedule) &&
      is.vector(screening.schedule, mode = "numeric") && length(screening.schedule) > 0L && 
      all(screening.schedule >= 0.0),
    "`other.cause.rates` must be a data.frame" = 
      {is.data.frame(other.cause.rates) && 
          all(c("Rate", "Age") %in% colnames(other.cause.rates))},
    "`burnin` must be a non-negative scalar" = 
      is.vector(burnin, mode = "numeric") && length(burnin) == 1L && 
      isTRUE(all.equal(burnin, round(burnin))) && burnin >= 0.0,
    "`verbose` must be a logical" = is.vector(verbose, mode = "logical") && 
      length(verbose) == 1L
  )
  
  # retrieve required information from the analysis object.
  risk_onset <- object$setup$t0
  thin <- object$setup$thin
  
  screening.schedule <- sort(screening.schedule)
  if (max(screening.schedule) > max(other.cause.rates) ||
      min(screening.schedule) < min(other.cause.rates)) {
    warning("`other.cause.rates` doesn't span requested screening schedule",
            call. = FALSE)
  }
  max_n_screens <- length(screening.schedule)
  
  keptM <- nrow(object$theta$H[[1L]])
  if (burnin >= keptM) stop("burnin removed all posterior draws", call. = FALSE)
  
  start <- burnin + 1L
  
  # number of posterior draws to sample
  draws <- seq_len(nrow(object$theta$H[[1L]]))
  draws <- draws[start:keptM]
  if (n.draws > length(draws)) {
    warning("requested # of posterior sets > # available after burnin; duplicate sampling")
  }
  if (n.draws != length(draws)) {
    draws <- sample(draws, n.draws, replace = n.draws > length(draws)) |> sort()
  }
  
  screen_ind <- matrix(0.0, nrow = max_n_screens, ncol = n.draws)
  screen_prog <- matrix(0.0, nrow = max_n_screens, ncol = n.draws)
  screen_total <- matrix(0.0, nrow = max_n_screens, ncol = n.draws)
  
  lambda_seq <- matrix(0.0, n.draws, 2L)
  
  if (verbose) cat(paste(0:10, collapse = "----"), "\n", sep = "")
  tick_increment = max(floor(n.draws / 50), 1L)
  
  for (m in seq_len(n.draws)) {
    
    z <- draws[m]
    
    if (m %% tick_increment == 0L && verbose) cat("=")
    
    # model parameters
    if (object$setup$indolent) {
      psi <- object$theta$psi$all[z, 1L]
    } else {
      psi <- 0.0
    }
    
    beta <- object$theta$beta$all[z, 1L]
    scale_H <- object$theta$H$all[z, "s_H"]
    scale_P <- object$theta$P$all[z, "s_P"]
    shape_H <- object$theta$H$all[z, "sh_H"]
    shape_P <- object$theta$P$all[z, "sh_P"]
    
    ## The lambda used
    lambda_seq[m, ] <- c(scale_P, shape_P)
    
    time_points <- risk_onset
    
    for (n_screens in seq_along(screening.schedule)) {
      
      time_points <- c(time_points, screening.schedule[n_screens])
      
      Dj <- D_j(psi = psi, 
                beta = beta, 
                shape.H = shape_H, 
                scale.H = scale_H, 
                shape.P = shape_P, 
                scale.P = scale_P, 
                time.points = time_points,  
                n.screens = n_screens)
      
      # P(SD)
      screen_total[n_screens, m] <- {Dj$indolent + Dj$progressive} * Dj$norm
      
      # P(SD, I=1)
      screen_ind[n_screens, m] <- Dj$indolent * Dj$norm
      
      # P(SD, I=0)
      screen_prog[n_screens, m] <- Dj$progressive * Dj$norm
      
    }
  } 
  if (verbose) cat('\n')
  
  #-- P(screen-detected cancer)
  SD_stat <- matrix(0.0, nrow = max_n_screens, ncol = n.draws)
  #-- P(indolent screen-detected)
  ODX_ind_stat <- matrix(0.0, nrow = max_n_screens, ncol = n.draws)
  #-- P(progressive overdiagnosed cancer)
  ODX_prog_stat <- matrix(0.0, nrow = max_n_screens, ncol = n.draws)
  
  if (verbose) cat(paste(0:10, collapse = "----"), "\n", sep = "")
  
  # probability to survive past each respective screening round
  hh1 <- sapply(screening.schedule, surv.fun, 
                age.min = screening.schedule[1L], 
                lambda.other = other.cause.rates$Rate, 
                lambda.ages = other.cause.rates$Age)
  
  for (z in seq_len(n.draws)) {
    
    if (z %% tick_increment == 0L && verbose) cat("=")
    
    #--Probability to have a screen-detected cancer
    
    # probability of a screen-detected cancer while still alive
    SD_stat[, z] <- (screen_total[, z] * hh1)
    
    ## Probability to have an indolent screen-detected cancer
    
    # probability of an indolent screen-detected cancer while still alive
    ODX_ind_stat[, z] <- (screen_ind[, z] * hh1)
    
    ## Probability to have a progressive and overdiagnosed screen-detected cancer
    hh <- sapply(screening.schedule, 
                 Prog.odx.program, 
                 age.min = screening.schedule[1L], 
                 P.lambda = lambda_seq[z,], 
                 lambda.other = other.cause.rates$Rate, 
                 lambda.ages = other.cause.rates$Age)
    ODX_prog_stat[, z] <- (screen_prog[, z] * hh)
    
  }
  if (verbose) cat('\n')
  
  res <- list()
  for (i in seq_len(nrow(ODX_ind_stat))) {
    res[[i]] <- wrap(ODX_ind_stat[i, , drop = FALSE], 
                     ODX_prog_stat[i, , drop = FALSE], 
                     SD_stat[i, , drop = FALSE])
  }
  names(res) <- paste0("screen.age.", screening.schedule)
  
  # sum across screens to get overall; should be means, but they cancel
  res$overall <- wrap(colMeans(ODX_ind_stat), 
                      colMeans(ODX_prog_stat), 
                      colMeans(SD_stat))
  
  class(res) <- c(class(res), "baclava.ODX.pred.old")
  
  res
}

#' @param x A an object of S3 class 'baclava.PDX.pred' as returned by \code{predictODX()}.
#' @param y Ignored.
#' @param ... Ignored.
#'
#' @describeIn predictODX Generate column plot of predicted overdiagnosis for each screen.
#'
#' @import ggplot2
#' @export
plot.baclava.ODX.pred.old <- function(x, y, ...) {
  
  n_screens <- length(x) - 1L
  
  x[["overall"]] <- NULL
  
  df <- do.call(rbind, x) |> as.data.frame()
  colnames(df) <- c("mean", "low", "high", "mid")
  df$screen <- rep(seq_len(n_screens), each = 3L)
  df$category <- rep(rownames(x[[1L]]), times = n_screens)
  
  df_sub <- subset(df, df$category != "total")
  df_sub$category <- factor(df_sub$category, levels = c("mortality", "indolent"))
  
  df_total <- subset(df, df$category == "total")
  
  gg <- ggplot(df_sub, aes(x = .data$screen, y = .data$mean)) + 
    geom_col(aes(fill = .data$category)) +
    geom_col(data = df_total, aes(x = .data$screen, y = .data$mean), fill = NA, color = "gray72") +
    geom_errorbar(data = df_total, aes(x = .data$screen, ymin = .data$low, ymax = .data$high), width = 0.2,
                  position=position_dodge(0.9), color = "gray72") +
    ggtitle("Mean Predicted Overdiagnosis for Each Screen") +
    xlab("Screen") + ylab("Overdiagnosis Rate, %") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
  gg
}
