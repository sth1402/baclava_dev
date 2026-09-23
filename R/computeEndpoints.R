#' Transform age at time of screening data to intervals.
#' 
#' @noRd
#' @param data.obj A list object containing
#'   `$endpoint_type`, the endpoint type; 
#'   `irateP`, the index of the rate_P parameters;
#'   `$n`, the number of cases;
#'   `$endpoint_time`, the censoring times;
#'   `$age_entry`, the ages at time of entry;
#'   `$ages_screen` a list, each element a vector of ages at time of screening; and 
#'   `$n_screen_positive`, an indicator of a positive screen.
#' @param t0 A scalar numeric object. The risk onset age.
#' 
#' @return A list object containing `$values`, a vector of 
#'   endpoints; `$starts` a vector of length n giving the first index of 
#'   `$values` pertaining to each case; `$ends`, a vector of length n giving 
#'   the last index of `$values` pertaining to each case; and `$lengths`, a 
#'   vector of length n giving the number of values of `$values` that pertain 
#'   to each individual. Note that the `$start` and `$stop` values are prepared 
#'   for use in C++, which indexes from 0.
#' 
#' @importFrom utils head
#' @keywords internal
.computeEndpoints <- function(data.obj, t0, m_point) {
  
  required_names <- c("endpoint_type", "n", "endpoint_time", "ages_screen", 
                      "iH", "iP", "ipsi", "ids", "age_entry")
  
  stopifnot(
    "`data.obj` does not contain required data" = !missing(data.obj) &&
      is.list(data.obj) && all(required_names %in% names(data.obj)),
    "`t0` must be a scalar numeric" = !missing(t0) && is.numeric(t0) &&
      is.vector(t0) && length(t0) == 1L
  )
  
  .generateIndices <- function(data.obj, t0, m_point) {
    times <- lapply(
      seq_len(data.obj$n),
      function(i, data, t0) {
        idx = seq.int(data$ages_screen$starts[i] + 1L, data$ages_screen$ends[i] + 1L)
        if (any(data$ages_screen$values[idx] > data$endpoint_time[i])) {
          stop("assessments after endpoint encountered", call. = FALSE)
        }
        if (data$endpoint_time[i] > m_point) {
          screens <- c(data$ages_screen$values[idx], m_point)
        } else {
          screens <- data$ages_screen$values[idx]
        }
        res <- switch(data$endpoint_type,
                      "preclinical" = c(t0, screens),
                      "censored" = c(t0, screens, data$endpoint_time[i], Inf),
                      "clinical" = c(t0, screens, data$endpoint_time[i]),
                      stop("should never happen; contact developer", call. = FALSE))
        res |> unique() |> sort()
      },
      data = data.obj, t0 = t0)
  }
  
  .generateIndicesNoScreens <- function(data.obj, t0, m_point) {
    times <- lapply(
      seq_len(data.obj$n),
      function(i, data, t0) {
        if (data$endpoint_time[i] > m_point) {
          screens <- c(m_point)
        } else {
          screens <- NULL
        }
        res <- switch(data$endpoint_type,
                      "preclinical" = stop("should never happen; contact developer", call. = FALSE),
                      "clinical" = c(t0, screens, data$endpoint_time[i]),
                      "censored" = c(t0, screens, data$endpoint_time[i], Inf))
        res |> unique() |> sort()
      },
      data = data.obj, t0 = t0)
  }
  
  if (is.null(data.obj$ages_screen) || length(data.obj$ages_screen) == 0L) {
    ages <- .generateIndicesNoScreens(data.obj, t0, m_point)
  } else {
    ages <- .generateIndices(data.obj, t0, m_point)
  }
  
  list("values"  = unlist(ages),
       "starts"  = utils::head(c(1L, cumsum(lengths(ages)) + 1L), -1L) - 1L,
       "ends"    = cumsum(lengths(ages)) - 1L,
       "lengths" = lengths(ages))
}

# internal function to vectorize a jagged list
#
# The returned list contains: `$values`, a numeric vector, the input ages 
#   vectorized; `$starts`, an integer vector of length n, each element, i,
#   is the first index of `$values` pertaining to participant i; `$ends`,
#   an integer vector of length n, each element, i, is the last index of
#   `$values` pertaining to participant i; and `$lengths`, an integer
#   vector of length n, each element, i, is the number of elements of
#   `$value` that pertain to participant i. 
# NOTE: The indices provided in `$start` and `$stop` are generated for use 
#   in C++, which indexes vectors starting at 0. To use in R, must add
#   1L.
.compressScreeningAges_routed <- function(data.assess, z_beta, rows) {
  if (length(rows) == 0L) {
    return(list())
  }
  aa <- data.assess$age_assess[rows]
  ids <- data.assess$id[rows]
  typ <- z_beta[rows]
  
  ages_by_id <- split(aa,  ids)
  types_by_id <- split(typ, ids)
  
  list("values" = unname(unlist(ages_by_id)),
       "starts" = unname(head(c(1L, cumsum(lengths(ages_by_id)) + 1L), -1L) - 1L),
       "ends"  = unname(cumsum(lengths(ages_by_id)) - 1L),
       "lengths" = unname(lengths(ages_by_id)),
       "types" = unname(unlist(types_by_id)) - 1L)
}

.basicDataObject_routed <- function(endpoint.type,
                                    clinical_idx,
                                    assess_idx,
                                    data.clinical,
                                    data.assess,
                                    z_H, z_P, z_psi, z_beta,
                                    H, P, psi, beta,
                                    t0, m_point) {
  
  clinical_subset <- data.clinical[clinical_idx, , drop = FALSE]
  assess_subset   <- if (length(assess_idx)) data.assess[assess_idx, , drop = FALSE] else NULL
  
  # constant (by construction of the cohort)
  iH   <- z_H  [clinical_idx[1L]] - 1L
  iP   <- z_P  [clinical_idx[1L]] - 1L
  ipsi <- z_psi[clinical_idx[1L]] - 1L
  
  # endpoint_time (same semantics as before)
  if (endpoint.type == "preclinical") {
    if (is.null(assess_subset))
      stop("one or more participants marked as screen detected do not have screening data", call. = FALSE)
    last_screen <- !duplicated(assess_subset$id, fromLast = TRUE)
    endpoint_time <- assess_subset$age_assess[last_screen]
  } else {
    endpoint_time <- clinical_subset$age_endpoint
  }
  
  ages_screen <- if (!is.null(assess_subset)) {
    .compressScreeningAges_routed(data.assess, z_beta, assess_idx)
  } else list()
  
  tmp_list <- list(
    "endpoint_type" = endpoint.type,
    "iH" = iH,
    "iP" = iP,
    "ipsi" = ipsi,
    "n" = nrow(clinical_subset),
    "ids" = clinical_subset$id,
    "endpoint_time" = endpoint_time,
    "age_entry" = clinical_subset$age_entry,
    "ages_screen" = ages_screen)
  
  tmp_list$endpoints <- .computeEndpoints(tmp_list, t0, m_point)
  tmp_list
}

#' @noRd
#' @param data.clinical A data.frame object. The clinical dataset.
#' @param data.assess A data.frame object. The assessment dataset.
#' @param H A StratifiedDistribution object. The Healthy Compartment model.
#' @param P A StratifiedDistribution object. The Pre-clinical Compartment model.
#' @param psi A StratifiedDistribution object. The psi Compartment model.
#' @param beta A StratifiedDistribution object. The beta Compartment model.
#' @param t0 A scalar numeric. The risk onset age.
#' @param m_point A scalar numeric. Used to force hard boundary cut in
#'   endpoint times to improve tau_hp proposals
#' @param all.types A logical. If TRUE, it is required that all types of endpoints
#'   must be represented in the data.
#' @keywords internal
.makeDataObjectsFull <- function(data.clinical, data.assess,
                                 H, P, psi, beta,
                                 t0, m_point, all.types = TRUE) {
  
  stopifnot(is.data.frame(data.clinical), is.data.frame(data.assess))
  
  if (any(data.assess$age_assess < t0)) {
    stop("assessment ages younger than t0 encountered", call. = FALSE)
  }
  
  # endpoint type → integers 1..3 (order fixed)
  et_levels <- c("preclinical", "censored", "clinical")
  if (!all(data.clinical$endpoint_type %in% et_levels)) {
    stop("unexpected endpoint_type values in clinical data", call. = FALSE)
  }
  z_et <- match(data.clinical$endpoint_type, et_levels)
  if (all.types && (!any(z_et == 1L) || !any(z_et == 2L) || !any(z_et == 3L))) {
    stop("all compartments must be represented in the data; verify `data.clinical$endpoint_type`", call. = FALSE)
  }
  
  # per-compartment strata codes (1..L) aligned to the right domain
  z_H <- if (identical(H$composition, "stratified")) { H$routing$row_to_stratum 
    } else { rep.int(1L, nrow(data.clinical)) }
  z_P <- if (identical(P$composition, "stratified")) { P$routing$row_to_stratum 
    } else { rep.int(1L, nrow(data.clinical)) }
  z_psi <- if (identical(psi$composition, "stratified")) { psi$routing$row_to_stratum 
    } else { rep.int(1L, nrow(data.clinical)) }
  z_beta <- if (identical(beta$composition, "stratified")) { beta$routing$row_to_stratum 
    } else { rep.int(1L, nrow(data.assess)) }
  
  L_H <- max(z_H)
  L_P <- max(z_P)
  L_psi <- max(z_psi)
  L_et <- 3L
  
  # strides for cartesian product (H × P × psi × endpoint_type)
  s_H <- 1L
  s_P <- L_H
  s_psi <- L_H * L_P
  s_et <- L_H * L_P * L_psi
  
  # group id per clinical row (1..G)
  g_clinical <- 1L + (z_H - 1L) * s_H + (z_P - 1L) * s_P +
    (z_psi - 1L) * s_psi + (z_et - 1L) * s_et
  
  rows_by_group_clinical <- split(seq_len(nrow(data.clinical)), g_clinical)
  
  # map clinical id -> group id, then assess rows to those groups
  id_to_group <- setNames(g_clinical, as.character(data.clinical$id))
  if (!all(as.character(data.assess$id) %in% names(id_to_group))) {
    stop("ids in data.assess not present in data.clinical", call. = FALSE)
  }
  g_assess <- unname(id_to_group[as.character(data.assess$id)])
  rows_by_group_assess <- split(seq_len(nrow(data.assess)), g_assess)
  
  # helper to label H/P/psi/beta names for messages
  lab_H <- if (identical(H$composition, "stratified")) H$routing$strata_labels else "all"
  lab_P <- if (identical(P$composition, "stratified")) P$routing$strata_labels else "all"
  lab_psi <- if (identical(psi$composition, "stratified")) psi$routing$strata_labels else "all"
  lab_beta <- if (identical(beta$composition, "stratified")) beta$routing$strata_labels else "all"
  
  # iterate clinical groups in increasing order
  groups <- sort(as.integer(names(rows_by_group_clinical)))
  data.obj <- vector("list", length(groups))
  names(data.obj) <- as.character(groups)
  
  width_n <- max(2L, ceiling(log10(nrow(data.clinical) + 1L)))
  
  for (ii in seq_along(groups)) {
    g  <- groups[ii]
    ic <- rows_by_group_clinical[[as.character(g)]]
    ia <- rows_by_group_assess[[as.character(g)]]
    if (is.null(ia)) ia <- integer(0L)
    
    # cohort labels for messages (constant within group)
    h_i <- z_H  [ic[1L]]
    p_i <- z_P  [ic[1L]]
    psi_i <- z_psi[ic[1L]]
    et_i <- z_et [ic[1L]]
    et_lab <- et_levels[et_i]
    
    beta_labs <- if (length(ia)) {
      unique(lab_beta[ z_beta[ia] ])
    } else character(0)
    
    # messages (keep the spirit of your original)
    if (length(ia)) {
      message(format(nrow(data.clinical[ic, , drop = FALSE]), width = width_n),
              " participants in subset H=", lab_H[h_i],
              ", P=", lab_P[p_i],
              ", psi=", lab_psi[psi_i],
              " with beta in (", paste(sort(beta_labs), collapse = ", "), ") ",
              et_lab)
    } else {
      if (et_lab == "preclinical") {
        stop("one or more participants marked as screen detected do not have screening data", call. = FALSE)
      }
      message(format(nrow(data.clinical[ic, , drop = FALSE]), width = width_n),
              " participants in subset H=", lab_H[h_i],
              ", P=", lab_P[p_i],
              ", psi=", lab_psi[psi_i],
              " with no screening data ",
              et_lab)
    }
    
    # build the exact structure C++ expects
    data.obj[[ii]] <- .basicDataObject_routed(
      endpoint.type = et_lab,
      clinical_idx = ic,
      assess_idx = ia,
      data.clinical = data.clinical,
      data.assess = data.assess,
      z_H = z_H, z_P = z_P, z_psi = z_psi, z_beta = z_beta,
      H = H, P = P, psi = psi, beta = beta,
      t0 = t0, m_point = m_point
    )
  }
  
  if (length(data.obj) == 0L) stop("verify subset structure; no data identified", call. = FALSE)
  
  stack_to_clin <- unlist(lapply(as.character(groups), function(g) rows_by_group_clinical[[g]]),
                          use.names = FALSE)
  
  # pick up the inverse permutation to original clinical order
  inv_clin <- attr(data.clinical, "inv_clin")
  
  res <- list(
    "data.objects" = data.obj,
    "stack_to_sorted" = stack_to_clin,
    "sorted_to_orig" = inv_clin
  )
  res
}
