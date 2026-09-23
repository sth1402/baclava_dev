#' internal procedure to determine how latent data was provided
#' 
#' @noRd
#' @param object A baclava object.
#' @param preclinical.data.ids A vector of participant ids for preclinical cases
#' @keywords internal
.resolveLatent <- function(object, preclinical.data.ids, verbose = FALSE) {
  
  if (!is.null(object$latent$participants)) {
    col_idx <- match(preclinical.data.ids, object$latent$participants)
  } else {
    if (verbose) 
      message("No `participants` recorded in latent data; ",
              "assuming columns correspond to original data order")
    col_idx <- preclinical.data.ids
  }
  
  if (any(is.na(col_idx)))
    warning("some screen-detected individuals are not present in latent data; ",
            "check `participants` argument used in extractLatent()", call. = FALSE)
  col_idx[!is.na(col_idx)]
}


#' internal procedure to sample from a truncated distribution
#' 
#' @noRd
#' @keywords internal
.rfunc_trunc <- function(a, age.hp, P, params) {
  
  lower_tail <- TRUE
  log_p <- FALSE
  
  n = length(a)
  if (length(age.hp) != n)
    stop("Length mismatch in rfunc_trunc: a, age.hp", call. = FALSE)
  
  out <- numeric(n)
  eps <- 1e-12
  
  for (i in seq_len(n)) {
    ai <- a[i]
    
    # CDF at lower bound
    if (is.finite(ai)) {
      lo <- P$pdist(x = ai, tau = age.hp[i], theta = params, lower.tail = lower_tail, log.p = log_p)
    } else {
      # a = -Inf  -> F(a) = 0
      lo <- 0.0
    }
    
    # degenerate case: truncation beyond support
    if (!is.finite(lo) || lo >= 1 - eps) {
      out[i] <- ai
      next
    }
    
    u  <- runif(1L, min = lo, max = 1.0 - eps)
    
    out[i] <- P$qdist(p = u, tau = age.hp[i], theta = params,
                      lower.tail = lower_tail, log.p = log_p)
    
  }
  out
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
  
  idx_lo <- findInterval(age.min, lambda.ages)
  idx_hi <- findInterval(t, lambda.ages)
  
  idx <- idx_lo:idx_hi
  lambdas <- lambda.other[idx]
  ages <- c(lambda.ages[idx], t)
  ages[1L] <- age.min  # start from age.min not interval start
  deltat <- diff(ages)
  
  exp(-drop(crossprod(lambdas, deltat)))
}

#' @noRd
#' @keywords internal
surv.fun <- Vectorize(surv.fun, vectorize.args = "t")

#' Probability of death in the interval [t0, t1]
#'
#' @noRd
#' @param t0 A numeric scalar. The lower boundary of the time interval in units
#'   equivalent to those used to define \code{prob}.
#' @param t1 A numeric scalar. The upper boundary of the time interval in units
#'   equivalent to those used to define \code{prob}.
#' @param otherCauses A data.frame object. Must contain 2 columns with headers
#'  Age, the age at which the corresponding probability pertains and
#'  Rate, the probability of death from other causes.
#'
#' @returns A scalar numeric. The probability of death from other causes in the
#'   specified interval.
#'   
#' @keywords internal
.computeProbabilityDeath <- function(S0, t1, other.causes) {
  S1 <- baclava:::surv.fun(age.min = 0, t = t1,
                           lambda.other = other.causes$Rate, 
                           lambda.ages = other.causes$Age)
  tt <- 1.0 - (S1 / pmax(S0, 1e-12))
  tt
}



#' Cohort-specific Overdiagnosis
#'
#' Estimates the overall and screening specific overdiagnosis probability
#'   for the cohort of the original analysis.
#'
#' @param object A 'baclava' object. The value object returned by \link{fit_baclava}()
#'   or \link{addLatent()}.
#' @param data.assess A data.frame object. Disease status assessments recorded
#'   during healthy or preclinical compartment, e.g., screenings for disease.
#'   The data must be structured as
#'   \itemize{
#'     \item \code{id}: A character, numeric, or integer object. The unique participant
#'     id to which the record pertains. Multiple records for each id are allowed.
#'     \item \code{age_assess}: A numeric object. The participant's age at time of
#'     assessment.
#'     \item \code{disease_detected}: An integer object. Must be binary 0/1, where
#'     1 indicates that disease was detected at the assessment; 0 otherwise.
#'  }
#'  If the sensitivity parameter (beta) is screen-specific, an additional
#'   column \code{screen_type} is required indicating the type of each
#'   screen.
#'  This input should be identical to that provided to obtain \code{object}.
#' @param data.clinical A data.frame object. The clinical data. The data must
#'   be structured as
#'   \itemize{
#'     \item \code{id}: A character, numeric, or integer object. The unique participant
#'       id to which the record pertains. Note these must include those provided in
#'       \code{data.assess}. Must be only 1 record for each participant.
#'     \item \code{age_entry}: A numeric object. The age at time of entry into the study.
#'       Note that this data is used to calculate a normalization; to expedite
#'       numerical integration, it is recommended that the ages be rounded.
#'       Optional input \code{round.age.entry} can be
#'       set to FALSE if this approximation is not desired; however, the
#'       computation time will significantly increase.
#'     \item \code{endpoint_type}: A character object. Must be one of \{"clinical",
#'       "censored", "preclinical"\}. Type "clinical" indicates that disease
#'       was diagnosed in the clinical compartment (i.e., symptomatic). Type
#'       "preclinical" indicates that disease was diagnosed in the preclinical
#'       compartment (i.e., during an assessment). Type
#'       "censored" indicates disease was not diagnosed prior to end of study.
#'     \item \code{age_endpoint}: A numeric object. The participant's age at the
#'       time the endpoint was evaluated.
#'  }
#'  If the sensitivity parameter (beta) is arm-specific, an additional
#'   column \code{arm} is required indicating the study arm to which each
#'   participant is assigned. Similarly, if the preclinical Weibull distribution is
#'   group-specific, an additional column \code{grp.rateP} is required. See Details
#'   for further information.
#'  This input should be identical to that provided to obtain \code{object}.
#' @param other.cause.rates A data.frame object. Age specific incidence rates
#'   that do not include the disease of interest. Must contain columns "Rate"
#'   and "Age".
#' @param plot A logical object. If TRUE, generates a boxplot of the overdiagnosis
#'   probability for each individual as a function of the screen at which
#'   disease was detected. Includes only the consecutive screens for which more
#'   than 1\% of the screen detected cases were detected.
#'
#' @returns A list object. 
#' \itemize{
#'   \item \code{all} An n x S matrix containing the estimated overdiagnosis
#'     probability for each individual (n) and each posterior parameter set (S).
#'   \item \code{mean.individual} A vector containing the mean across S of
#'     the estimated overdiagnosis for each individual, i.e., \code{rowMeans(all)}.
#'   \item \code{mean.overall} A numeric, the mean overdiagnosis
#'     probability across all posterior parameter sets and screen-detected cases, 
#'     i.e., \code{mean(all)}. 
#'   \item \code{summary.by.screen} A matrix containing the summary statistics of 
#'     \code{mean.individual} for the individuals detected positive at
#'     each screen, i.e., \code{summary(mean.individual[diagnosis_screen_id == i])}. 
#' }
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
#' # This is for illustration only -- the number of Gibbs samples should be
#' # significantly larger and the epsilon values should be tuned.
#' example <- fit_baclava(data.assess = data.screen,
#'                        data.clinical = data.clinical,
#'                        t0 = 30.0,
#'                        theta_0 = theta_0,
#'                        prior = prior,
#'                        save.latent = TRUE)
#'
#' # if rates are not available, an all cause dataset is provided in the package
#' # NOTE: these predictions will be over-estimated
#'
#' data(all_cause_rates)
#' all_cause_rates <- all_cause_rates[, c("Age", "both")]
#' colnames(all_cause_rates) <- c("Age", "Rate")
#'
#' cohort_odx <- cohortODX(object = example,
#'                         data.clinical = data.clinical,
#'                         data.assess = data.screen,
#'                         other.cause.rates = all_cause_rates,
#'                         plot = FALSE)
#'
#' @importFrom stats rweibull
#' @include utilities.R
#' @export
cohortODX <- function(object, data.clinical, data.assess,
                      other.cause.rates, 
                      plot = TRUE, n.samples = 1000L, K = 1L, verbose = TRUE) {
  
  stopifnot(
    "`object` must be an object of class 'baclava'" = !missing(object) &&
      inherits(object, "baclava"),
    "`data.clinical` must be a data.frame with columns id, age_entry, endpoint_type, and age_endpoint" =
      !missing(data.clinical) && is.data.frame(data.clinical) &&
      all(c("id", "age_entry", "age_endpoint", "endpoint_type") %in% colnames(data.clinical)),
    "`data.assess` must be a data.frame with columns id, age_assess, and disease_detected" =
      !missing(data.assess) && is.data.frame(data.assess) &&
      all(c("id", "age_assess", "disease_detected") %in% colnames(data.assess)),
    "`other.cause.rates` must be a data.frame w/ numeric Rate & Age" = 
      is.data.frame(other.cause.rates) && 
      all(c("Rate", "Age") %in% colnames(other.cause.rates)) &&
      is.numeric(other.cause.rates$Rate) && is.numeric(other.cause.rates$Age),
    "`plot` must be logical" = is.vector(plot, mode = "logical") && length(plot) == 1L,
    "`n.samples` must be in an integer > 1" = is.vector(n.samples, "numeric") &&
      length(n.samples) == 1L && isTRUE(all.equal(n.samples, round(n.samples, 0))) &&
      n.samples > 1,
    "`K` must be a positive integer" = is.vector(K, "numeric") && length(K) == 1L &&
      K >= 1L && isTRUE(all.equal(K, as.integer(K)))
  )
  
  # method is not appropriate if there is no indolence indicator
  if (!object$setup$indolent) stop("cohortODX() requires indolence")

  if (is.null(object$latent) || is.vector(object$latent$tau))
    stop("no latent data in object; call addLatent() first if latent data stored externally", call. = FALSE)
  
  # verify inputs
  n.samples <- as.integer(n.samples)
  if (n.samples < 100L)
    warning("`n.samples` = ", n.samples, " is very small; ",
            "posterior summaries may be unreliable. ",
            "Consider using n.samples >= 100", call. = FALSE)
  
  K <- as.integer(K)
  if (K < 10L)
    warning("`K` = ", K, " is very small; ",
            "individual level ODX may be unreliable. ",
            "Consider using K >= 10", call. = FALSE)
  
  # verifying properly formatted data provided  
  data.assess <- .testDataAssess(data.assess)
  data.clinical <- .testDataClinical(data.clinical)
  
  # extracting fully processed models from fitted object
  beta <- object$model$beta
  H <- object$model$H
  P <- object$model$P
  psi <- object$model$psi

  # extend other causes to 120 years
  other.cause.rates <- extend_other_causes(other.cause.rates)
  
  ### Cohort ODX is only calculated using the screen detected cases

  # group data, allowing for only preclinical to be present in provided
  # data (all.types)
  # do not need m_point, defaulting to Inf
  data.objs <- .makeDataObjectsFull(data.clinical = data.clinical,
                                    data.assess = data.assess,
                                    H, P, psi, beta,
                                    t0 = object$setup$t0, 
                                    m_point = Inf, all.types = FALSE)
  
  data.objs <- data.objs$data.objects
  
  # exact all preclinical ids
  preids_in_data_obj <- NULL
  for (i in seq_along(data.objs)) {
    if (data.objs[[i]]$endpoint_type != "preclinical") next
    preids_in_data_obj <- c(preids_in_data_obj, data.objs[[i]]$ids)
  }
  
  # resolve latent
  latent_cols <- .resolveLatent(object, preids_in_data_obj, verbose)
  
  ### Calculate ODX
  
  prob_odx_all <- list()
  n_screens <- NULL
  afs <- NULL
  aad <- NULL
  
  n_iters <- nrow(object$theta$H[[1L]])
  n.samples <- as.integer(min(n.samples, n_iters))
  if (n.samples > n_iters)
    warning("n.samples > number of posterior samples; sampled with replacement")
  
  samples <- sample(n_iters, n.samples, replace = n.samples > n_iters) |> sort()

  # need a pointer just behind to the first column pertaining to the
  # data object under consideration
  track_ids <- 0L
  for (do in seq_along(data.objs)) {
    
    if (data.objs[[do]]$endpoint_type != "preclinical") next
    
    n_screens <- c(n_screens, data.objs[[do]]$ages_screen$lengths)
    afs <- c(afs, data.objs[[do]]$age_entry)
    aad <- c(aad, data.objs[[do]]$endpoint_time)
    
    # clinical ids present in this data grouping
    ids_in_data_obj <- data.objs[[do]]$ids
    
    # preclinical model for this data grouping
    iP <- data.objs[[do]]$iP + 1L
    
    prob_odx <- matrix(NA_real_, nrow = length(ids_in_data_obj), ncol = n.samples,
                       dimnames = list(ids_in_data_obj, NULL))
    
    # This does not change across iterations, calculate once for all individuals
    S_at_SD <- surv.fun(age.min = 0.0, 
                        t = data.objs[[do]]$endpoint_time,
                        lambda.other = other.cause.rates$Rate, 
                        lambda.ages = other.cause.rates$Age)
    
    for (i in seq_along(ids_in_data_obj)) {
      
      # retrieve latent variables
      # need to include a pointer to where we are in the latent variables
      # latent variable constructed in in the order of data objects
      Z_tau_i <- object$latent$tau[, latent_cols[i + track_ids]]
      Z_indolent_i <- object$latent$indolent[, latent_cols[i + track_ids]]
      
      # subset to just the samples under consideration
      Z_tau_i <- Z_tau_i[samples]
      Z_indolent_i <- Z_indolent_i[samples]
      
      # cycle through posterior samples to get get probability
      for (s in seq_along(samples)) {
        
        # if indolent, probability of ODX is 1
        if (Z_indolent_i[s] == 1L) {
          prob_odx[i, s] <- 1.0
          next
        }
        
        # generate a transition time to clinical guaranteeing that it is
        # after the screen-detected age
        ztau <- Z_tau_i[s]
        age_SD <- data.objs[[do]]$endpoint_time[i]
        LB <- max(0, age_SD - ztau)
        
        # preclinical model at this iteration
        useP <- iP
        parsP <- object$theta$P[[useP]][samples[s], ]
        Pd <- object$model$P$distributions[[useP]]
        
        # generate sojourn time through preclinical compartment
        wait <- .rfunc_trunc(a = rep(LB, K),
                             age.hp = rep(ztau, K), 
                             P = Pd, 
                             params = parsP)
        tau_PC_k <- pmax(ztau + wait, age_SD)
        
        prob_k <- .computeProbabilityDeath(S0 = S_at_SD[i],
                                           t1 = tau_PC_k,
                                           other.causes = other.cause.rates)
        prob_odx[i, s] <- mean(prob_k)
        
      }
    }
    
    track_ids <- track_ids + length(ids_in_data_obj)
    
    prob_odx_all <- append(prob_odx_all, list(prob_odx))
  }
  prob_odx_all <- do.call(rbind, prob_odx_all)

  # overall ODX rate
  prob_odx_all <- prob_odx_all * 100.0
  overall <- sum(prob_odx_all, na.rm = TRUE) / {n.samples * nrow(prob_odx_all)}
  
  prob_by_individual <- rowMeans(prob_odx_all, na.rm = TRUE)
  
  res <- list("all" = prob_odx_all, 
              "mean.individual" = prob_by_individual,
              "mean.overall" = overall)
  class(res) <- c("baclava.ODX.cohort", class(res))
  attr(res, "screen_count") <- n_screens
  attr(res, "afs") <- afs
  attr(res, "aad") <- aad
  
  if (plot) print(plot(res))
  
  res
}

#' @param x A an object of S3 class 'baclava.ODX.cohort' as returned by \code{cohortODX()}.
#' @param y Ignored.
#' @param type A character object indicating the x-variable. Must be one of
#'  'screen' (no. of screens having occurred by screen of detection);
#'  'afs' (age of first screen); or 'aad' (age at screen detection).
#' @param ... Ignored.
#'
#' @describeIn cohortODX Generate column plot of cohort overdiagnosis for each screen.
#'
#' @import ggplot2
#' @export
plot.baclava.ODX.cohort <- function(x, y, type = c("screen", "afs", "aad"), ...) {
  
  type <- match.arg(type)
  
  if (type == "screen") {
    n_screen_detected <- nrow(x$all)
    n_screens_vec <- attr(x, "screen_count")
  
    one_percent <- ceiling(n_screen_detected * 0.01)
    large_enough <- which(tabulate(n_screens_vec) >= one_percent)
    large_enough <- large_enough[c(0L, diff(large_enough)) <= 1L]
  
    df <- data.frame(prob_odx = x$mean.individual,
                     screen = factor(n_screens_vec))
  
    gg <- ggplot(subset(df, df$screen %in% large_enough),
                 aes(.data$screen, .data$prob_odx)) +
      geom_boxplot() +
      stat_summary(fun = mean, fun.min = mean, fun.max = mean,
                   geom = "errorbar", color = "red", width = 0.75,
                   show.legend = TRUE) +
      theme(axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 8),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            plot.title = element_text(hjust = 0.5)) +
      ylab("Overdiagnosis Probability") +
      xlab("Screen at which Disease Detected") +
      ggtitle("Overdiagnosis Probability vs Screen at which Disease Detected")
  } else if (type == "afs") {
    n_screen_detected <- nrow(x$all)
    afs <- attr(x, "afs")
    
    df <- data.frame(prob_odx = x$mean.individual,
                     afs = afs)
    
    gg <- ggplot(df,
                 aes(.data$afs, .data$prob_odx)) +
      geom_point() +
      stat_summary(fun = mean, fun.min = mean, fun.max = mean,
                   geom = "errorbar", color = "red", width = 0.75,
                   show.legend = TRUE) +
      theme(axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 8),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            plot.title = element_text(hjust = 0.5)) +
      ylab("Overdiagnosis Probability") +
      xlab("Age of First Screen") +
      ggtitle("Overdiagnosis Probability vs Age of First Screen")    
  } else if (type == "aad") {
    n_screen_detected <- nrow(x$all)
    aad <- attr(x, "aad")
    
    df <- data.frame(prob_odx = x$mean.individual,
                     aad = aad)
    
    gg <- ggplot(df,
                 aes(.data$aad, .data$prob_odx)) +
      geom_point() +
      stat_summary(fun = mean, fun.min = mean, fun.max = mean,
                   geom = "errorbar", color = "red", width = 0.75,
                   show.legend = TRUE) +
      theme(axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 8),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            plot.title = element_text(hjust = 0.5)) +
      ylab("Overdiagnosis Probability") +
      xlab("Age at Screen Detection") +
      ggtitle("Overdiagnosis Probability vs Age at Screen Detection")    
  }
  gg
}