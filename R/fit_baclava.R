#' Bayesian Analysis of Cancer Latency with Auxiliary Variable Augmentation
#'
#' Markov chain Monte Carlo sampler to fit a three-state mixture compartmental
#'  model of cancer natural history to individual-level screening and cancer
#'  diagnosis histories in a Bayesian framework.
#'
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
#' @param baclava.object NULL or a 'baclava' object. To continue a calculation,
#'   provide the object returned by a previous call.
#' @param H An object of R6 class \code{Compartment} specifying the 
#'   distribution(s) and parameter update strategy for the transition time 
#'   through the Healthy compartment. The allowed distributions are currently
#'   limited to weibull_d() and gamma_d().
#' @param P An object of R6 class \code{Compartment} specifying the 
#'   distribution(s) and parameter update strategy for the transition time 
#'   through the Preclinical compartment. The allowed distributions are currently
#'   limited to weibull_d() and gamma_d().
#' @param psi An object of R6 class \code{Compartment} specifying the 
#'   distribution(s) and parameter update strategy for the probability of 
#'   indolence. The allowed distributions are currently
#'   limited to bern_d().
#' @param beta An object of R6 class \code{Compartment} specifying the 
#'   distribution(s) and parameter update strategy for the probability of 
#'   screen detection. The allowed distributions are currently
#'   limited to bern_d() and only Gibbs updates are allowed.
#' @param indolent A logical object. If \code{FALSE}, disease under analysis
#'   does not have an indolent state, i.e., it is always progressive.
#'   This input is provided for convenience; if FALSE, psi is reset as
#'   \code{compartment(distributions = list(bern_d(p = param(tag = "psi.p", init = 0.0))), strategies = NULL)}, 
#'   which sets the value of psi to be 0 with no updates.
#'   If \code{baclava.object} is a 'baclava' object, this input is ignored.
#' @param t0 A non-negative scalar numeric object. The risk onset age. Must be
#'   less than the earliest assessment age, entry age, and endpoint age.
#'   If \code{baclava.object} is a 'baclava' object, this input is ignored.
#' @param M A positive integer object. The number of Monte Carlo samples. This
#'   is the total, i.e., M = adaptive$warmup + n_MCMC.
#' @param thin A positive integer object. Keep each thin-th step of the sampler
#'   after the warmup period, if any, is complete.
#' @param round.age.entry A logical object. If TRUE, the age at time of entry
#'   will be rounded to the nearest integer prior to performing the MCMC.
#'   This data is used to estimate the probability of experiencing clinical
#'   disease prior to entering the study, which is estimated using a
#'   time consuming numerical integration procedure. It is expected that
#'   rounding the ages at time of entry introduces minimal bias. If FALSE,
#'   and ages cannot be grouped, these integrals significantly increase
#'   computation time.
#'   If \code{baclava.object} is a 'baclava' object, this input is ignored.
#' @param verbose A logical object. If \code{TRUE}, a progress bar will be shown
#'   during the MCMC.
#' @param save.latent A logical object. If \code{TRUE}, latent variable tau
#'   and indolence will be returned. These can be very large matrices. To
#'   estimate the cohort overdiagnosis probability using \link{cohortODX}(),
#'   this must be set to TRUE.
#' @param latent_dir A character string to the directory in which latent data
#'   should be stored (ignored if save.latent = FALSE)
#' @param latent_dump = 500L An integer. The frequency at which latent data
#'   should be stored in an external file (ignored if save.latent = FALSE)
#' @param file_prefix A character string. Appended to file names of external
#'   latent data files (ignored if save.latent = FALSE)
#'
#' @returns An object of S3 class \code{baclava}, which extends a list object.
#'   \itemize{
#'     \item theta: A list of the posterior distribution parameters at the 
#'       thinned samples.
#'     \item latent: A list of latent values and acceptance information
#'       \itemize{
#'         \item tau: If \code{save.latent = TRUE}, a matrix.
#'           The age at time of transition from
#'           healthy to preclinical compartment for each participant at the
#'           thinned samples. If \code{save.latent = FALSE}, their values at
#'           the last iteration.
#'.        \item indolent: If \code{save.latent = TRUE && I_marginalize != "ALL"}, a matrix.
#'           The indolent status for each participant at the
#'           thinned samples. If \code{save.latent = FALSE && I_marginalize != "ALL"}, 
#'           their values at
#'           the last iteration. Will be \code{NA} if disease is always progressive.
#'         \item tau_accept: If \code{save.latent = TRUE}, a matrix.
#'           TRUE/FALSE indicating if tau was accepted. Note that NA is used
#'           to indicate an Inf->Inf transition. 
#'           If \code{save.latent = FALSE}, their values at
#'           the last iteration.  
#'         \item accept_rate: Individual level acceptance rates for tau
#'        }
#'     \item accept: A list of the accept indicator of each MH step at the thinned samples.
#'     \item acceptance_prob A list of the acceptance probabilities in the post-burnin
#'       period of the sampler for each MH step.
#'     \item epsilon: A list. The random walk settings used in the post-warmup period.
#'     \item adaptive: A list. The random walk values at each iteration of the
#'       adaptive period.
#'     \item setup: A list of inputs provided to the call.
#'       \itemize{
#'         \item t0: The input age of risk onset.
#'         \item indolent: TRUE if disease is not progressive.
#'         \item round.age.entry: TRUE if age at entry was rounded to the nearest
#'           whole number.
#'         \item M The total number of iterations.
#'         \item thin: The number of samples dropped between kept MCMC iterations.
#'         \item burnin: The number of iterations dropped.
#'      }
#'     \item model: The provided model specification.
#'     \item call: The matched call.
#'   }
#'
#' @examples
#'
#' data(screen_data)
#'
#' # H compartment model
#' # fit parameter is the mean of the distribution
#' mu_H <- param(100.0, tag = "mu_H", minv = 0.0, maxv = Inf,
#'               prior = lnorm_p(log(100), 1))
#' # shape parameter
#' sh_H <- param(2.0, tag = "sh_H"), minv =  0.0, maxv = Inf
#'                prior = lnorm_p(log(2), 0.1))
#' # Weibull distribution
#' H <- compartment(weibull_d(mean = mu_H, shape = sh_H))
#'
#' # P compartment model
#' # fit parameter is the mean of the distribution
#' mu_P <- param(2, tag = "mu_P", minv = 0.0, maxv = Inf,
#'               prior = lnorm_p(log(2), 1))
#' # shape parameter
#' sh_P <- param(2.0, tag = "sh_P"), minv =  0.0, maxv = Inf
#'                prior = lnorm_p(log(2), 0.1))
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
#' strategy_H <- MH(params = c("sc_H"), s0 = 0.014, 
#'                  adaptive = adapt(learn = "none"))
#' 
#' strategy_P <- MH(params = c("r_P"), s0 = 0.06, 
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
#' summary(example)
#' print(example)
#'
#' # To continue this calculation
#' example_continued <- fit_baclava(data.assess = data.screen,
#'                                  data.clinical = data.clinical,
#'                                  baclava.object = example)
#' @importFrom dplyr count
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom tidyr nest
#' @import Rcpp
#' @import RcppNumerical
#' @include computeEndpoints.R Class_Compartment.R RcppExports.R utilities.R Class_Prior.R
#' @useDynLib baclava
#' @export
fit_baclava <- function(data.assess,
                        data.clinical,
                        H, P, psi, beta, 
                        blocks, ...,
                        baclava.object = NULL,
                        M = 100L,
                        thin = 1L,
                        burnin = 0L,
                        t0 = 0,
                        indolent = TRUE,
                        round.age.entry = TRUE,
                        verbose = TRUE,
                        save.latent = FALSE,
                        debug = 0, m_point = Inf,
#                        prop_choice = 1, # 1 = orig; 2 = intermediate; 3 = current
#                        use_orig_prop = FALSE, I_marginalize = NULL,
                        latent_dir = "./", latent_dump = 500L, file_prefix = "Analysis") {
  
  
  
  stopifnot("`baclava.object` must be NULL or of class 'baclava'" = is.null(baclava.object) ||
              "baclava" %in% class(baclava.object))
  
  if (!is.null(baclava.object)) {
    # pull model from prior object
    H <- baclava.object$model$H
    P <- baclava.object$model$P
    beta <- baclava.object$model$beta
    psi <- baclava.object$model$psi
    
    blocks <- baclava.object$model$blocks
    
    t0 <- baclava.object$setup$t0
    indolent <- baclava.object$setup$indolent
    
    # assign initial values as last kept value
    last_it <- nrow(baclava.object$theta$H[[1L]])
    
    H_param_values <- lapply(baclava.object$theta$H, function(x) x[last_it, ]) |> unname() |> unlist()
    names(H_param_values) <- .strip_theta_names(names(H_param_values))
    H$setValue(H_param_values)
    
    P_param_values <- lapply(baclava.object$theta$P, function(x) x[last_it, ]) |> unname() |> unlist()
    names(P_param_values) <- .strip_theta_names(names(P_param_values))
    P$setValue(P_param_values)
    
    beta_param_values <- lapply(baclava.object$theta$beta, function(x) x[last_it, ]) |> unname() |> unlist()
    names(beta_param_values) <- .strip_theta_names(names(beta_param_values))
    beta$setValue(beta_param_values)
    
    psi_param_values <- lapply(baclava.object$theta$psi, function(x) x[last_it, ]) |> unname() |> unlist()
    names(psi_param_values) <- .strip_theta_names(names(psi_param_values))
    psi$setValue(psi_param_values)
    
    # assign epsilon and turn-off adaptation
    for (i in seq_along(baclava.object$epsilon)) {
      blocks$blocks[[names(baclava.object$epsilon)[i]]]$continue(
        s = baclava.object$epsilon[[i]]$s,
        L = baclava.object$epsilon[[i]]$L,
        D = baclava.object$epsilon[[i]]$D, 0.2
      )
      blocks$blocks[[names(baclava.object$epsilon)[i]]]$adaptive$learn <- "none"
    }
    
    I_marginalize = baclava.object$setup$I_marginalize
  }
  
  stopifnot(
    "'H', 'P', and 'beta' must be provided" = !missing(H) && !missing(P) &&
      !missing(beta),
    "`H` must be a Compartment object" = is.Compartment(H),
    "`P` must be a Compartment object" = is.Compartment(P),
    "`beta` must be a Compartment object" = is.Compartment(beta),
    "`blocks` must be a list of Block objects" = {is.list(blocks) && length(blocks) > 0L &&
        all(unlist(lapply(blocks, function(i) is.MH(i) || is.Gibbs(i))))} || is.Strategies(blocks),
    "`data.assess` must be a data.frame with columns id, age_assess, and disease_detected" =
      !missing(data.assess) && is.data.frame(data.assess) &&
      all(c("id", "age_assess", "disease_detected") %in% colnames(data.assess)),
    "`data.clinical` must be a data.frame with columns id, age_entry, endpoint_type, and age_endpoint" =
      !missing(data.clinical) && is.data.frame(data.clinical) &&
      all(c("id", "age_entry", "age_endpoint", "endpoint_type") %in% colnames(data.clinical)),
    "`indolent` must be logical" = is.vector(indolent, mode = "logical") &&
      length(indolent) == 1L,
    "`psi`is required if indolent = TRUE" = !indolent ||
      {!missing(psi) && indolent && is.Compartment(psi)},
    "`t0` must be a non-negative scalar numeric" =
      is.vector(t0, mode = "numeric") && length(t0) == 1L && t0 >= 0.0,
    "`M` must be a positive integer" = is.vector(M, mode = "numeric") &&
      length(M) == 1L && isTRUE(all.equal(M, as.integer(M))) && M > 0.0,
    "`thin` must be a positive integer" = is.vector(thin, mode = "numeric") &&
      length(thin) == 1L && isTRUE(all.equal(thin, as.integer(thin))) && thin > 0.0,
    "`round.age.entry` must be logical" = is.logical(round.age.entry),
    "`verbose` must be a logical" = is.vector(verbose, mode = "logical") &&
      length(verbose) == 1L
  )
  
  if (!is.Strategies(blocks)) blocks <- do.call(Strategies, blocks)
  
  parameters_fit <- lapply(blocks$blocks, function(i) i$params) |> unlist()
  
  not_fixed <- NULL
  for (i in H$distributions) {
    for (k in i$parameters) {
      if (!k$parameters$value$fixed) not_fixed <- c(not_fixed, k$parameters$value$tag)
    }
  }
  
  for (i in P$distributions) {
    for (k in i$parameters) {
      if (!k$parameters$value$fixed) not_fixed <- c(not_fixed, k$parameters$value$tag)
    }
  }
  
  for (i in beta$distributions) {
    for (k in i$parameters) {
      if (!k$parameters$value$fixed) not_fixed <- c(not_fixed, k$parameters$value$tag)
    }
  }
  
  for (i in psi$distributions) {
    for (k in i$parameters) {
      if (!k$parameters$value$fixed) not_fixed <- c(not_fixed, k$parameters$value$tag)
    }
  }
  
  if (!setequal(not_fixed, parameters_fit)) {
    stop("mismatch in parameter specification and strategies: \n\t", 
         "parameters that are not fixed: ", paste(not_fixed, collapse = ", "), "\n\t",
         "parameters in strategies: ", paste(parameters_fit, collapse = ", "), "\n")
  }
  
  ### Verify data.frame inputs
  
  if (debug > 0) print(".testDataAssess")
  data.assess <- .testDataAssess(data.assess)
  if (debug > 0) print(".testDataClinical")
  data.clinical <- .testDataClinical(data.clinical)
  
  if (min(data.assess$age_assess) <= t0) stop("t0 cannot be >= assessment ages", call. = FALSE)
  
  # all ids in data.assess must be in data.clinical
  # note this is not true in the other direction
  if (!all(data.assess$id %in% data.clinical$id)) {
    stop("participant ids in `data.assess` not found in `data.clinical`",
         call. = FALSE)
  }
  
  if (!indolent) {
    message(" * * * disease is always progressive\n")
    psi <- compartment(distributions = list(bern_d(p = param(tag = "psi.p", init = 0.0))))
  }
  
  ## connecting model compartments to data to identify subset structure
  ## if present
  if (debug > 0) print("finalizing")
  H$finalize(data.assess, data.clinical, where="auto")
  P$finalize(data.assess, data.clinical, where="auto")
  psi$finalize(data.assess, data.clinical, where="auto")
  beta$finalize(data.assess, data.clinical, where="assess")
  
  ### misc data prepping
  
  if (debug > 0) print("counting screens")
  n_participants <- nrow(data.clinical)
  
  if (round.age.entry) {
    message(" * * * age at entry rounded to nearest whole number")
    data.clinical$age_entry <- as.integer(round(data.clinical$age_entry, 0L))
  }
  
  if (debug > 0) print("making data")
  data.objs <- .makeDataObjectsFull(data.clinical = data.clinical,
                                    data.assess = data.assess,
                                    H, P, psi, beta,
                                    t0 = t0, m_point = m_point)
  
  stack_to_sorted <- data.objs$stack_to_sorted
  sorted_to_orig <- data.objs$sorted_to_orig
  data.objs <- data.objs$data.objects
  
  ages_hp <- list()
  indolents <- list()
  if (!is.null(baclava.object)) {
    ages_hp_tmp <- drop(baclava.object$latent$last_iter_tau)
    indolents_tmp <- drop(baclava.object$latent$last_iter_indolent)
    ids_orig <- attr(data.clinical, "orig_ord")
    for (i in seq_along(data.objs)) {
      idx <- match(data.objs[[i]]$ids, ids_orig)
      ages_hp[[i]] <- ages_hp_tmp[idx]
      indolents[[i]] <- indolents_tmp[idx]
    }
  }
  
  if (!indolent) {
    indolents <- lapply(data.objs, function(x) rep(0L, x$n))
  }
  
  if (debug > 0) print("making compartment lists")
  # must convert distributions to lists
  H_list <- H$toList()
  P_list <- P$toList()
  psi_list <- psi$toList()
  beta_list <- beta$toList()
  blocks_list <- blocks$toList()
  
  # Need to change MH to MH_with_Z if psi parameters are in strategies
  psi_params <- unlist(psi$getTags())
  cnt <- 0L
  for (i in seq_along(blocks$blocks)) {
    params <- blocks$blocks[[i]]$params
    if (any(psi_params %in% params) && 
        !all(c("CEN", "CLI", "PRE") %in% I_marginalize)) {
      blocks_list$blocks[[i]]$alt_meta_for <- "MH_with_Z"
      cnt <- cnt + 1L
    } else {
      blocks_list$blocks[[i]]$alt_meta_for <- blocks_list$blocks[[i]]$meta_for
    }
  }
  if (cnt > 1) warning("Z will be updated in each MH step involving a psi parameter!")
  
  skip <- max(burnin, max(unlist(lapply(blocks$blocks, function(i) i$adaptive$warmup))))
  
  # clinical-sorted -> original order
  ids_stack <- unlist(lapply(data.objs, "[[", "ids"), use.names = FALSE)
  ids_orig <- attr(data.clinical, "orig_ord")
  
  stack_to_orig <- match(as.character(ids_orig), ids_stack)
  
  if (is.null(I_marginalize)) I_marginalize <- "NONE"
  
  if (debug > 0) print("calling C++")
  tryCatch(
    out_object <- MCMC_cpp_internal(data_objects = data.objs,
                                    indolents = indolents,
                                    taus = ages_hp,
                                    H = H_list, P = P_list, 
                                    psi = psi_list, beta = beta_list,
                                    strategies = blocks_list,
                                    t0 = t0,
                                    M = M,
                                    thin = thin,
                                    burnin = skip,
                                    verbose = as.integer(verbose),
                                    save_latent = as.integer(save.latent),
                                    debug = debug,
                                    # use_bug = prop_choice, 
                                    # use_orig_prop = use_orig_prop,
                                    # I_marginalize_CEN = "CEN" %in% I_marginalize,
                                    # I_marginalize_PRE = "PRE" %in% I_marginalize,
                                    # I_marginalize_CLI = "CLI" %in% I_marginalize,
                                    latent_dir = latent_dir,
                                    latent_dump = latent_dump,
                                    file_prefix = file_prefix,
                                    latent_order = stack_to_orig - 1L),
    error = function(e) {
      stop(e$message)
    })
  
  tryCatch({
    
    ikept <- out_object$ikept
    
    ## add compartment labels
    
    labs <- list(
      "H" = .comp_labels(H),
      "P" = .comp_labels(P),
      "beta" = .comp_labels(beta),
      "psi" = .comp_labels(psi)
    )
    
    names(H$distributions) <- labs$H
    names(P$distributions) <- labs$P
    names(beta$distributions) <- labs$beta
    names(psi$distributions) <- labs$psi
    
    theta_names <- names(out_object$theta)
    
    # limit the compartment names to only those present in theta
    nms_theta <- names(out_object$theta)
    
    # accept elements are on the dimension of block parameters
    # also include tau_hp, which theta does not
    nms_accept <- names(out_object$accept)
    
    out_object$theta <- lapply(out_object$theta, function(cmp) {
      if (is.matrix(cmp)) {
        t(cmp)
      } else {
        lapply(cmp, function(mat) t(mat))
      }
    })
    
    out_object$accept <- lapply(out_object$accept, function(cmp) {
      if (is.matrix(cmp)) {
        t(cmp)
      } else {
        lapply(cmp, function(mat) t(mat))
      }
    })
    
    # the matrices returned are likely to have 1 too many elements
    for (comp in nms_accept) {
      out_object$accept[[comp]] <- unlist(out_object$accept[[comp]])
      if (length(out_object$accept[[comp]]) == ikept) next
      out_object$accept[[comp]] <- utils::head(out_object$accept[[comp]], ikept)
      
      for (i in seq_len(length(out_object$theta[[comp]]))) {
        if (nrow(out_object$theta[[comp]][[i]]) == ikept) next
        out_object$theta[[comp]][[i]] <- out_object$theta[[comp]][[i]][seq_len(ikept), , drop = FALSE]
      }
    }
    
    if (anyNA(stack_to_orig)) {
      warning("ID mismatch between C++ stack and clinical original IDs; returned latent object will be out of order")
    } else {
      out_object$latent$tau_accept_rate <- out_object$latent$tau_accept_rate[stack_to_orig]
      names(out_object$latent$tau_accept_rate) <- ids_orig
    }
    
    out_object$setup <- list("t0" = t0,
                             "indolent" = indolent,
                             "round.age.entry" = round.age.entry,
                             "M" = M,
                             "thin" = thin,
                             "burnin" = skip,
                             "I_marginalize" = I_marginalize, 
                             "cycle" = ifelse(is.null(baclava.object), 1L, baclava.object$setup$cycle + 1L))
    
    out_object$model <- list("H" = H, "P" = P, "psi" = psi, "beta" = beta,
                             "blocks" = blocks)
    out_object$ikept <- NULL
    
    out_object$call <- match.call()
    
    tmp_theta <- do.call(cbind, out_object$theta)
    for (i in H$distributions) {
      for (k in i$parameters) {
        if (k$parameters$value$fixed) {
          tmp_theta <- cbind(tmp_theta, k$parameters$value$init)
          colnames(tmp_theta)[ncol(tmp_theta)] <- k$parameters$value$tag
        }
      }
    }
    
    for (i in P$distributions) {
      for (k in i$parameters) {
        if (k$parameters$value$fixed) {
          tmp_theta <- cbind(tmp_theta, k$parameters$value$init)
          colnames(tmp_theta)[ncol(tmp_theta)] <- k$parameters$value$tag
        }
      }
    }
    
    for (i in beta$distributions) {
      for (k in i$parameters) {
        if (k$parameters$value$fixed) {
          tmp_theta <- cbind(tmp_theta, k$parameters$value$init)
          colnames(tmp_theta)[ncol(tmp_theta)] <- k$parameters$value$tag
        }
      }
    }
    
    for (i in psi$distributions) {
      for (k in i$parameters) {
        if (k$parameters$value$fixed) {
          tmp_theta <- cbind(tmp_theta, k$parameters$value$init)
          colnames(tmp_theta)[ncol(tmp_theta)] <- k$parameters$value$tag
        }
      }
    }
    
    
    theta2 <- list()
    # psi is not always present
    theta2$H <- list()
    pars <- H$getTags()
    for (i in seq_along(labs$H)) {
      theta2$H[[labs$H[i]]] <- tmp_theta[, unlist(pars[[i]]), drop = FALSE]
    }
    
    theta2$P <- list()
    pars <- P$getTags()
    for (i in seq_along(labs$P)) {
      theta2$P[[labs$P[i]]] <- tmp_theta[, unlist(pars[[i]]), drop = FALSE]
    }
    
    theta2$beta <- list()
    pars <- beta$getTags()
    for (i in seq_along(labs$beta)) {
      theta2$beta[[labs$beta[i]]] <- tmp_theta[, unlist(pars[[i]]), drop = FALSE]
    }
    
    theta2$psi <- list()
    pars <- psi$getTags()
    for (i in seq_along(labs$psi)) {
      theta2$psi[[labs$psi[i]]] <- tmp_theta[, unlist(pars[[i]]), drop = FALSE]
    }
    
    
    
    out_object$theta <- theta2
    
    # if (!is.null(baclava.object)) {
    #   for (i in seq_along(out_object$theta)) {
    #     for (j in seq_along(out_object$theta[[i]])) {
    #       out_object$theta[[i]][[j]] <- rbind(baclava.object$theta[[i]][[j]], out_object$theta[[i]][[j]])
    #     }
    #   }
    # }
    
    
    class(out_object) <- "baclava"
  }, error = function(e) {
    message("there was an issue in manipulating the final object", e$message)
    warning("! ! ! ! ! returned object does not respect documented structure ! ! ! ! !")
  }) 
  out_object
}

#' @noRd
#' @keywords internal
is.baclava <- function(x) inherits(x, "baclava")