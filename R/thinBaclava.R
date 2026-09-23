#' Thin based on class of matrix
#' @noRd
#' @keywords internal
.thin_matrix <- function(mat, keep_rows) {
  if (inherits(mat, "big.matrix")) {
    new_mat <- bigmemory::big.matrix(
      nrow = length(keep_rows),
      ncol = ncol(mat),
      type = bigmemory::typeof(mat)
    )
    for (i in seq_along(keep_rows)) {
      new_mat[i, ] <- mat[keep_rows[i], ]
    }
    return(new_mat)
  }
  mat[keep_rows, , drop = FALSE]
}

#' Thin a Baclava Object
#'
#' Applies additional thinning to the posteriors of a fitted baclava object,
#'   and optionally to an associated latent variable list to maintain alignment.
#'
#' @param obj A baclava object as returned by \code{fit_baclava()} or
#'   \code{combineBaclava()}.
#' @param by A scalar positive integer. Keep every \code{by}-th row of the
#'   posteriors. The effective total thinning becomes
#'   \code{obj$setup$thin * by}.
#' @param latent An optional list as returned by \code{extractLatent()} or
#'   \code{combineLatent()}. If provided, the same thinning is applied to
#'   the latent variables to maintain alignment with the posteriors. If
#'   \code{NULL} and latent variables are not stored in \code{obj}, a warning
#'   is issued and a misalignment flag is set in \code{obj$setup}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @returns A reduced baclava object with thinned posteriors and updated
#'   setup. If \code{latent} was provided, the thinned latents are returned
#'   in the obj$latent slot. If thinning was applied without latent
#'   alignment, \code{setup$latent.misalignment.risk} is set to \code{TRUE}.
#'
#' @export
thinBaclava <- function(obj,
                        by,
                        latent = NULL,
                        verbose = TRUE) {
  
  # Input validation
  stopifnot(
    "`obj` must be a baclava object" = !missing(obj) && inherits(obj, "baclava"),
    "`obj` must contain required elements" = all(c("setup", "latent", "theta", "model") %in% names(obj)),
    "`by` must be a positive integer" = !missing(by) &&
      is.vector(by, "numeric") && length(by) == 1L &&
      by >= 1L && isTRUE(all.equal(by, as.integer(by))),
    "`latent` must be a list or NULL" = is.null(latent) || is.list(latent),
    "`verbose` must be logical" = is.vector(verbose, "logical") && length(verbose) == 1L
  )
  
  by <- as.integer(by)
  
  if (by == 1L) {
    if (verbose) message("by = 1: no thinning applied")
    return(obj)
  }
  
  # Determine keep rows
  M_current <- nrow(obj$theta$H[[1L]])
  keep_rows <- seq(1L, M_current, by = by)
  M_new <- length(keep_rows)
  
  if (verbose)
    message("thinning: keeping every ", by, "th row",
            " (", M_current, " \u2192 ", M_new, " iterations;",
            " effective thin = ", obj$setup$thin * by, ")")
  
  # Thin posteriors
  thinned_theta <- obj$theta
  
  for (comp in names(thinned_theta)) {
    for (dist in names(thinned_theta[[comp]])) {
      thinned_theta[[comp]][[dist]] <-
        thinned_theta[[comp]][[dist]][keep_rows, , drop = FALSE]
    }
  }
  
  # Handle latent variables
  thinned_latent <- NULL
  misalignment_risk <- FALSE

  if (is.list(obj$latent) && is.null(obj$latent$latent_dir) && is.matrix(obj$latent$tau)) {
    # no latent directory specified and tau is a matrix means that latent
    # was saved within the baclava object
    if (verbose) message("thinning latent variables stored in object")
    
    if (obj$setup$indolent && is.matrix(obj$latent$indolent))
      thinned_latent$indolent <- obj$latent$indolent[keep_rows, , drop = FALSE]
    
    thinned_latent$tau <- obj$latent$tau[keep_rows, , drop = FALSE]
    
    if (!is.null(latent) && verbose)
      message("Note: `latent` argument ignored as latent variables are ",
              "stored in object")
    thinned_latent$participants <- obj$latent$participants
    
  } else if (is.list(obj$latent) && is.null(obj$latent$latent_dir) && is.vector(obj$latent$tau)) {
    # no latent directory specified and tau is a vector means that latent
    # was not saved
    message("Latent variables were not stored during analysis")
  } else if (!is.null(latent)) {
    # latent provided externally - validate dimensions and thin
    
    if (!{"tau" %in% names(latent)})
      stop("`latent` must contain `tau` element", call. = FALSE)
    
    if (obj$setup$indolent && nrow(latent$indolent) != M_current)
      stop("`latent$ZI` has ", nrow(latent$indolent), " rows but obj has ",
           M_current, " posterior iterations; cannot align", call. = FALSE)
    
    if (nrow(latent$tau) != M_current)
      stop("`latent$tau` has ", nrow(latent$tau), " rows but obj has ",
           M_current, " posterior iterations; cannot align", call. = FALSE)
    
    if (verbose) message("thinning externally provided latent variables")
    
    thinned_latent$indolent <- if (obj$setup$indolent) .thin_matrix(latent$indolent, keep_rows) else NULL
    thinned_latent$tau <- .thin_matrix(latent$tau, keep_rows)
    thinned_latent$participants <- obj$latent$participants

  } else if (is.list(obj$latent) && !is.null(obj$latent$latent_dir)) {
    # no latent provided and none in object
    misalignment_risk <- TRUE
    warning("latent data is stored externally in '", obj$latent$latent_dir, 
            "'; call addLatent() to maintain alignment with posteriors", call. = FALSE)
  } else {
    stop("`obj$latent` does not expected structure", call. = FALSE)
  }
  

  # Update setup
  updated_setup <- obj$setup
  updated_setup$thin <- obj$setup$thin * by
  updated_setup$latent.misalignment.risk <- misalignment_risk
  
  # Assemble output
  structure(
    list(
      "theta" = thinned_theta,
      "latent" = thinned_latent,
      "model" = obj$model,
      "setup" = updated_setup
    ),
    class = c("baclava", "REDUCED")
  )
  
}