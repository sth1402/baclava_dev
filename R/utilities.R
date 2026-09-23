.strip_theta_names <- function(nms) {
  if (any(substr(nms, 1, 3) == "mu.")) {
    tst <- substr(nms, 1, 3) == "mu."
    nms[tst] <- gsub("mu.", "", nms[tst], fixed = TRUE)
  } 
  if (any(substr(nms, 1, 5) == "rate.")) {
    tst <- substr(nms, 1, 5) == "rate."
    nms[tst] <- gsub("rate.", "", nms[tst], fixed = TRUE)
  } 
  if (any(substr(nms, 1, 6) == "scale.")) {
    tst <- substr(nms, 1, 6) == "scale."
    nms[tst] <- gsub("scale.", "", nms[tst], fixed = TRUE)
  } 
  if (any(substr(nms, 1, 7) == "median.")) {
    tst <- substr(nms, 1, 7) == "median."
    nms[tst] <- gsub("median.", "", nms[tst], fixed = TRUE)
  }
  nms
}

.continue <- function(comp, eps, posterior) {
  
  if (!is.null(eps)) {
    for (i in seq_along(eps)) {
      H$strategies$blocks[[i]]$continue(eps[[i]]$s, eps[[i]]$L)
    }
  }
  
  last_params <- list()
  for (i in seq_along(posterior)) {
    last_params[[i]] <- posterior[[i]][nrow(posterior[[i]]), ] |> drop()
  }
  H$setValue(last_params)
  H
}

.continuation <- function(baclava.object) {
  
  if (is.null(baclava.object)) return()
  
  if (verbose) message("\n", rep("* ", 10), "\n")
  if (verbose) message("Continuation of previous run")

  H <- object$model$H$toList()
  H_eps <- object$epsilon$H
  for (i in seq_along(H_eps)) {
    H$strategies$blocks[[i]]$s0 <- H_eps[[i]]$s
    H$strategies$blocks[[i]]$L <- H_eps[[i]]$L
    H$strategies$blocks[[i]]$L_diag <- diag(H_eps[[i]]$L)
    H$strategies$blocks[[i]]$adaptive$warmup <- 0L
  }

  P <- object$model$P$toList()
  P_eps <- object$epsilon$P
  for (i in seq_along(P_eps)) {
    P$strategies$blocks[[i]]$s0 <- P_eps[[i]]$s
    P$strategies$blocks[[i]]$L <- P_eps[[i]]$L
    P$strategies$blocks[[i]]$L_diag <- diag(P_eps[[i]]$L)
    P$strategies$blocks[[i]]$adaptive$warmup <- 0L
  }
  
  psi <- object$model$psi$toList()
  psi_eps <- object$epsilon$psi
  for (i in seq_along(psi_eps)) {
    psi$strategies$blocks[[i]]$s0 <- psi_eps[[i]]$s
    psi$strategies$blocks[[i]]$L <- psi_eps[[i]]$L
    psi$strategies$blocks[[i]]$L_diag <- diag(psi_eps[[i]]$L)
    psi$strategies$blocks[[i]]$adaptive$warmup <- 0L
  }
  
  list(H = H, P = P, psi = psi, pia = object$model$pia$toList(), beta = object$model$beta$toList())
}

.testDataAssess <- function(data.assess) {
  
  data.assess$id <- data.assess$id |> as.character()

  if (!is.numeric(data.assess$disease_detected)) {
    stop("data.assess$disease_detected must be binary 0/1", call. = FALSE)
  }
  
  if (!is.integer(data.assess$disease_detected)) {
    data.assess$disease_detected <- as.integer(round(data.assess$disease_detected))
  }

  if (!all(data.assess$disease_detected %in% c(0L, 1L))) {
    stop("`data.assess$disease_detected` must be binary 0/1", call. = FALSE)
  }
  
  # ages cannot be negative
  if (any(data.assess$age_assess < 0.0)) {
    stop("`data.assess$age_assess` cannot be negative", call. = FALSE)
  }
  
  # order the history and clinical datasets
  data.assess <- data.assess[order(data.assess$id, data.assess$age_assess), ]
  rownames(data.assess) <- NULL
  data.assess
}

.testDataClinical <- function(data.clinical) {

  orig_order <- data.clinical$id
  data.clinical$id <- data.clinical$id |> as.character()
  
  # ensure only 1 record per id is present in clinical data
  if (any(duplicated(data.clinical$id))) {
    stop("at least 1 participant id has multiple records in `data.clinical`",
         call. = FALSE)
  }
  
  # ages cannot be negative
  if (any(data.clinical$age_endpoint < 0.0)) {
    stop("`data.clinical$age_endpoint` cannot be negative", call. = FALSE)
  }
  
  if (any(data.clinical$age_entry < 0.0)) {
    stop("`data.clinical$age_entry` cannot be negative", call. = FALSE)
  }
  
  # ensure all endpoint types are as expected
  if (!all(data.clinical$endpoint_type %in% c("clinical", "preclinical", "censored"))) {
    stop("unrecognized endpoint type in `data.clinical$endpoint_type`", call. = FALSE)
  }
  
  # order the history and clinical datasets
  ord_clin <- order(data.clinical$id)
  data.clinical <- data.clinical[ord_clin, ]
  rownames(data.clinical) <- NULL
  
  # stash both permutations as attributes for later
  attr(data.clinical, "orig_ord") <- orig_order
  attr(data.clinical, "ord_clin") <- ord_clin
  inv_clin <- integer(length(ord_clin)); inv_clin[ord_clin] <- seq_along(ord_clin)
  attr(data.clinical, "inv_clin") <- inv_clin
  
  data.clinical
}

.getPreviousResults <- function(data.objects, baclava.object, indolent) {
  
  cnt <- 1L
  
  age_at_tau_hats <- list()
  indolents <- list()
  for (obj in data.objects) {
    idx_obj <- match(obj$ids, colnames(baclava.object$tau_hp))
    if (any(is.na(idx_obj))) stop("data ids do not match previous result", call. = FALSE)
    
    last_row <- nrow(baclava.object$tau_hp)
    while (TRUE) {
      age_at_tau_hats[[cnt]] <- baclava.object$tau_hp[last_row, idx_obj]
      if (!all(age_at_tau_hats[[cnt]] < 1e-8)) break
      last_row <- last_row - 1L
      if (last_row <= 0L) stop("all prior tau_hp values are zero", call. = FALSE)
    }
    
    if (indolent) indolents[[cnt]] <- baclava.object$indolent[last_row, idx_obj]
    
    cnt <- cnt + 1L
  }
  
  list("age.at.tau.hp.hats" = age_at_tau_hats, "indolents" = indolents)      
  
}

# Extend other-cause hazard table to a desired max age
#
# @param tab A data.frame with columns Age and Rate
# @param max_age A scalar integer. The maximum age needed.
# @param method One of 'gompertz' or 'carry'. 'gompertz' assumes log-hazard 
#   increases linearly with age (i.e., a Gompertz tail). A linear model of 
#   log(Rate) versus Age is fitted to the portion of the table with 
#   Age >= tail_fit_min_age. Predicted hazards for later ages are obtained 
#   from this fit. If fewer than five valid points exist in that tail, the 
#   function falls back to the "carry" method. 'carry' assumes the hazard 
#   remains constant beyond the last observed age. The final observed rate is 
#   carried forward.
# @param tail_fit_min_age A scalar numeric. The starting age of the data to
#   to be fit.
#
# A data frame identical in structure to the input, with added rows extending 
#   Age up to max_age and corresponding extrapolated Rate values.
#
#' @noRd
#' @keywords internal
extend_other_causes <- function(tab, max_age = 120L,
                                method = c("gompertz", "carry"),
                                tail_fit_min_age = 60L) {
  
  stopifnot("`tab` must be a data.frame" = !missing(tab) && is.data.frame(tab), 
            "`tab` must contain columns Age and Rate" = all(c("Age","Rate") %in% names(tab)))
  
  method <- match.arg(method)
  tab <- tab[order(tab$Age), ]
  tab$Rate[tab$Rate < 0] <- 0
  last_age <- max(tab$Age)
  if (last_age >= max_age) return(tab)
  
  ages_new <- seq.int(from = last_age + 1L, to = max_age, by = 1L)
  
  if (method == "carry") {
    rate_last <- tail(tab$Rate, 1L)
    add <- data.frame(Age = ages_new, Rate = rate_last)
  } else {
    # Gompertz: log(rate) linear in age on the tail
    tail <- subset(tab, Age >= tail_fit_min_age & Rate > 0)
    if (nrow(tail) < 5L) {  # fallback if not enough points
      rate_last <- tail(tab$Rate, 1L)
      add <- data.frame(Age = ages_new, Rate = rate_last)
    } else {
      fit <- lm(log(Rate) ~ Age, data = tail)
      pred <- exp(predict(fit, newdata = data.frame(Age = ages_new)))
      add <- data.frame(Age = ages_new, Rate = pmax(pred, 0))
    }
  }
  rbind(tab, add)
}
