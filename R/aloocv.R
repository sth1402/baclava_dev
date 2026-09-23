#' Approximate Leave-One-Out Cross-Validation
#'
#' Approximate leave-one-out cross-validation computed from the posterior
#'   draws of the Markov chain Monte Carlo sampler as implemented in
#'   \link{fit_baclava}().
#'
#' Computes the predictive fit of a model.
#'   For each individual and each MCMC draw, the function approximates the
#'   marginal likelihood via importance sampling. It samples J.increment values
#'   of the individual's latent variables using the Metropolis-Hastings proposal
#'   distributions and computes the effective sample size (ESS) of the
#'   importance sampling procedure. If the target ESS is not met, J.increment
#'   additional samples are taken, and the ESS is re-evaluated.
#'   This is repeated until either the ESS is
#'   satisfied or J.max samples have been drawn.
#'
#' @param object The value object returned by \link{fit_baclava}().
#' @param data.assess A data.frame. Disease status assessments recorded
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
#' @param data.clinical A data.frame. The clinical data. The data must
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
#' @param J.increment An integer object. The number of replicates of each
#'   participant to generate in each iteration of the importance sampling procedure
#'   to attain desired effective sample size.
#' @param J.max An integer object. The maximum number of samples to be
#'   drawn.
#' @param ess.target An integer object. The target effective sample size in the
#'   importance sampling procedure.
#' @param verbose A logical object. If TRUE, progress information will be
#'   printed. This input will be ignored if n.core > 1.
#' @param lib An optional character vector allowing for library path
#'   to be provided to cluster.
#'
#' @returns The provided baclava object augmented with the ALOOCV results. 
#'   Element aloocv is a list containing \code{summary}: the min, mean, and
#'   the 1% quantile of the ESS; the max, mean, and 99% quantile of J; the
#'   likelihood; and the individual-level and estimated predictive fit; and
#'   \code{result}: the likelihood, ESS, and J for each MCMC sample for each 
#'   participant.
#'
#' @examples
#'
#' data(screen_data)
#'
#' H <- weibull_dis(mu = MH(x = 100,
#'                          walk = reflective(1e-3, min = 0, log = TRUE),
#'                          prior = lnormp(log(100), 1)),
#'                  shape = fixed(2))
#'
#' P <- weibull_dis(mu = MH(x = 2,
#'                          walk = reflective(1e-3, min = 0, log = TRUE),
#'                          prior = lnormp(log(2), 1)),
#'                  shape = fixed(2))
#'
#' psi <- bern_dis(p = MH(x = .1,
#'                        walk = reflective(1e-3, min = 0.0, max = 1.0),
#'                        prior = betap(1, 1)))
#'                 
#' beta <- bern_dis(p = Gibbs(x = .8, prior = betap(1, 1)))
#'                 
#' # This is for illustration only -- the number of iterations should be
#' # significantly larger and the epsilon values should be tuned.
#' example <- fit_baclava(data.assess = data.screen,
#'                        data.clinical = data.clinical,
#'                        t0 = 30.0,
#'                        H = H, P = P, beta = beta, psi = psi)
#'
#' res <- aloocv(example, data.clinical, data.screen)
#'
#' @export
aloocv <- function(object,
                   data.clinical, data.assess,
                   J.increment = 75L,
                   J.max = 225L,
                   ess.target = 50L,
                   verbose = TRUE,
                   lib = NULL) {

  stopifnot(
    "`object` must be of class 'baclava'" = !missing(object) &&
      "baclava" %in% class(object),
    "`data.assess` must be a data.frame with columns id, age_assess, and disease_detected" =
      !missing(data.assess) && is.data.frame(data.assess) &&
      all(c("id", "age_assess", "disease_detected") %in% colnames(data.assess)),
    "`data.clinical` must be a data.frame with columns id, age_entry, endpoint_type, and age_endpoint" =
      !missing(data.clinical) && is.data.frame(data.clinical) &&
      all(c("id", "age_entry", "age_endpoint", "endpoint_type") %in% colnames(data.clinical)),
    "`J.increment` must be a positive integer" = is.vector(J.increment, mode = "numeric") &&
      length(J.increment) == 1L && isTRUE(all.equal(J.increment, round(J.increment))) &&
      J.increment > 0,
    "`J.max` must be a positive integer >= J.increment" = is.vector(J.max, mode = "numeric") &&
      length(J.max) == 1L && isTRUE(all.equal(J.max, round(J.max))) &&
      J.max >= J.increment,
    "`ess.target` must be a positive integer" = is.vector(ess.target, mode = "numeric") &&
      length(ess.target) == 1L && isTRUE(all.equal(ess.target, round(ess.target))) &&
      ess.target > 0,
    "`verbose` must be a logical" = is.vector(verbose, mode = "logical") &&
      length(verbose) == 1L,
    "`lib` must be a character vector" = is.null(lib) || is.character(lib)
  )

  t0 <- object$setup$t0
  
  adaptive = list("warmup" = 0L,
                  "delta" = 0.0,
                  "gamma" = 0.0,
                  "kappa" = 0.0,
                  "m0" = 0L,
                  "smooth" = 5L)
  
  n_mcmc_samples <- nrow(object$tau_hp)
  if (is.null(n_mcmc_samples) || n_mcmc_samples == 0L) 
    stop("cannot perform aloocv without latent data")
  n_participants <- nrow(data.clinical)

  data.assess <- .testDataAssess(data.assess)
  data.clinical <- .testDataClinical(data.clinical)

  # all ids in data.assess must be in data.clinical
  # note this is not true in the other direction
  if (!all(data.assess$id %in% data.clinical$id)) {
    stop("participant ids in `data.assess` not found in `data.clinical`",
         call. = FALSE)
  }
  
  beta <- object$model$beta
  if (is.StratifiedDistribution(beta) && length(beta$dists) == 1L)
    beta <- beta$dists[[1L]]
  H <- object$model$H
  if (is.StratifiedDistribution(H) && length(H$dists) == 1L)
    H <- H$dists[[1L]]
  P <- object$model$P
  if (is.StratifiedDistribution(P) && length(P$dists) == 1L)
    P <- P$dists[[1L]]
  psi <- object$model$psi
  if (is.StratifiedDistribution(psi) && length(psi$dists) == 1L)
    psi <- psi$dists[[1L]]
  
  # this procedure creates a screen_type vector to be added to
  # data.assess and ensures that beta is aligned.
  # Upon completion, beta will ALWAYS be a StratifiedDistribution
  tmp_list <- .processSensitivity(data.assess, data.clinical, beta)
  beta <- tmp_list$beta
  data.assess[[beta$subset.var]] <- tmp_list$assess.beta
  
  ### Healthy Distribution
  
  # this procedure creates a healthy grouping vector to be added to
  # data.clinical and ensures that H is aligned.
  # Upon completion, H will ALWAYS be a StratifiedDistribution
  tmp_list <- .processCohortGrp(data.assess, data.clinical, H, "H")
  H <- tmp_list$dist
  data.clinical[[H$subset.var]] <- tmp_list$clinical.grp
  
  ### Preclinical Distribution

  # this procedure creates a pre-clinical grouping vector to be added to
  # data.clinical and ensures that P is aligned.
  # Upon completion, P will ALWAYS be a StratifiedDistribution
  tmp_list <- .processCohortGrp(data.assess, data.clinical, P, "P")
  P <- tmp_list$dist
  data.clinical[[P$subset.var]] <- tmp_list$clinical.grp
  
  ### Indolence Fraction

  # this procedure creates a pre-clinical grouping vector to be added to
  # data.clinical and ensures that psi is aligned.
  # Upon completion, psi will ALWAYS be a StratifiedDistribution
  tmp_list <- .processCohortGrp(data.assess, data.clinical, psi, "psi")
  psi <- tmp_list$dist
  data.clinical[[psi$subset.var]] <- tmp_list$clinical.grp

  if (object$setup$round.age.entry) {
    data.clinical$age_entry <- as.integer(round(data.clinical$age_entry, 0L))
  }

  data.objs <- .makeDataObjectsFull(data.clinical = data.clinical,
                                    data.assess = data.assess,
                                    H, P, psi, beta,
                                    t0 = t0)
  
  H_list <- lapply(H$dists, function(x) x$toList())
  P_list <- lapply(P$dists, function(x) x$toList())
  psi_list <- lapply(psi$dists, function(x) x$toList())
  beta_list <- lapply(beta$dists, function(x) x$toList())
  
  result_list <- list()
  for (i in seq_along(data.clinical$id)) {
    
    if (verbose && {{i %% floor(nrow(data.clinical)*0.05)} < 0.5}) {
      print(c(i, nrow(data.clinical)))
    }

    # data of i-th subject
    data.clinical_i = data.clinical[i, ]
    data.assess_i   = data.assess[data.assess$id %in% data.clinical_i$id, ]

    # replicate clinical row J.increment times
    data.clinical_i_rep <- data.clinical_i[rep(1L, J.increment), ]
    data.clinical_i_rep$id <- seq_len(J.increment)  # update ID column
    
    # replicate assess rows for each ID in 1:J.increment
    data.assess_i_rep <- do.call(rbind, lapply(seq_len(J.increment), function(j) {
      tmp <- data.assess_i
      tmp$id <- j
      tmp
    }))
    
    # convert to expected data structure
    data.objs <- suppressMessages(
      .makeDataObjectsFull(data.clinical = data.clinical_i_rep,
                           data.assess = data.assess_i_rep,
                           H, P, psi, beta,
                           t0 = t0, all.types = FALSE)[[1L]])
    
    # calculate lik - prop_latent
    result <- model_comparison(data.objs, 
                               H_list, P_list, beta_list, psi_list, t0, 
                               adaptive, 
                               object$theta, 
                               ess.target,
                               J.increment, J.max,
                               n_mcmc_samples)

    result$i <- i
    result_list[[i]] <- result
  }

  i_ordering <- lapply(result_list, "[[", "i") |> unlist()
  lik_all <- lapply(result_list, "[[", "lik") |> do.call(what = rbind)
  ess_all <- lapply(result_list, "[[", "ess") |> do.call(what = rbind)
  J_all <- lapply(result_list, "[[", "J") |> do.call(what = rbind)

  if (any(J_all == J.max & ess_all < ess.target)) {
    warning("Target ESS not met for ",
            sum(J_all == J.max & ess_all < ess.target),
            " MCMC sample(s)/participants",
            call. = FALSE)
  }
  aloocv_results <- list()
  aloocv_results$all <- list(lik = lik_all[i_ordering, ],
                             ess = ess_all[i_ordering, ],
                             J = J_all[i_ordering, ])

  aloocv_results$summary <- list("ess_min" = min(ess_all),
                                 "ess_mean" = mean(ess_all),
                                 "ess_q001" = quantile(ess_all, probs = 0.01),
                                 "J_max" = max(J_all),
                                 "J_mean" = mean(J_all),
                                 "J_q099" = quantile(J_all, probs = 0.99),
                                 "lpd" = rowMeans(aloocv_results$all$lik) |> log() |> sum(),
                                 "elpd_loo_vec" = log(1.0 / rowMeans(1.0 / aloocv_results$all$lik)))
  
  aloocv_results$summary$elpd_loo <- mean(aloocv_results$summary$elpd_loo_vec)

  object$aloocv <- aloocv_results
  object
}