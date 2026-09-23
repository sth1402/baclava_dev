#' Combine Baclava Objects Across Continuation Runs
#'
#' Combines a list of baclava objects from sequential continuation runs into
#'   a single object with concatenated posteriors and latent variables.
#'   Objects need not be provided in cycle order — sorting is performed
#'   internally using \code{setup$cycle}.
#'
#' @param ... Two or more baclava objects as returned by
#'   \code{fitBaclava()} or \code{addLatent}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @returns A single baclava object with:
#'   \describe{
#'     \item{\code{$theta}}{Posteriors concatenated across all cycles in
#'     cycle order.}
#'     \item{\code{$latent}}{Latent variables concatenated across cycles if
#'     present in all objects; \code{NULL} otherwise.}
#'     \item{\code{$model}}{Carried forward from the last cycle.}
#'     \item{\code{$setup}}{Updated with summed \code{M}, \code{cycle} as a
#'     vector of all combined cycle numbers, \code{burnin} from cycle 1,
#'     and all other fields verified consistent across cycles.}
#'   }
#'
#' @note No validation is performed to confirm that objects originate from
#'   the same chain. Combining objects from different models will produce
#'   a syntactically valid but scientifically meaningless result.
#'
#' @include extractLatent.R
#' @export
combineBaclava <- function(..., verbose = TRUE) {
  
  baclava.list <- list(...)
  
  # Input validation
  if (!is.list(baclava.list) || length(baclava.list) == 0L)
    stop("`baclava.list` must be a non-empty list of baclava objects")
  
  stopifnot(
    "all elements must be baclava objects" = all(sapply(baclava.list, is.baclava)),
    "all elements must contain required slots" =
      all(sapply(baclava.list, function(x)
        all(c("theta", "latent", "setup", "model") %in% names(x)))),
    "`verbose` must be logical" =
      is.vector(verbose, "logical") && length(verbose) == 1L
  )
  
  if (any(sapply(baclava.list, function(x) inherits(x, "combined"))))
    stop("cannot combine already-combined baclava objects; ",
         "pass all objects in a single call to combineBaclava()", call. = FALSE)
  
  # Sort by cycle
  cycles <- sapply(baclava.list, function(x) x$setup$cycle)
  
  if (any(duplicated(cycles)))
    stop("duplicate cycle numbers detected: ",
         paste(cycles[duplicated(cycles)], collapse = ", "))
  
  ord <- order(cycles)
  baclava.list <- baclava.list[ord]
  cycles <- cycles[ord]
  
  if (verbose)
    message("Combining ", length(baclava.list), " cycles: ",
            paste(cycles, collapse = " \u2192 "))
  
  # Validate consistency across cycles
  .check_consistent <- function(field, label) {
    vals <- sapply(baclava.list, function(x) x$setup[[field]])
    if (length(unique(vals)) > 1L)
      stop(label, " is inconsistent across cycles: ",
           paste(unique(vals), collapse = ", "))
    invisible(vals[1L])
  }
  
  t0 <- .check_consistent("t0", "`t0`")
  thin <- .check_consistent("thin", "`thin`")
  indolent <- .check_consistent("indolent", "`indolent`")
  I_marginalize <- .check_consistent("I_marginalize", "`I_marginalize`")
  round_age <- .check_consistent("round.age.entry", "`round.age.entry`")
  
  # validate theta compartment structure
  compartment_names <- sapply(baclava.list, function(x)
    paste(sort(names(x$theta)), collapse = ","))
  if (length(unique(compartment_names)) > 1L)
    stop("theta compartment structure is inconsistent across cycles")
  
  if (verbose) message("Consistency checks passed")
  
  # Combine posteriors
  combined_theta <- list()
  
  for (comp in names(baclava.list[[1L]]$theta)) {
    combined_theta[[comp]] <- list()
    for (dist in names(baclava.list[[1L]]$theta[[comp]])) {
      mats <- lapply(baclava.list, function(x) x$theta[[comp]][[dist]])
      
      n_cols <- sapply(mats, ncol)
      if (length(unique(n_cols)) > 1L)
        stop("posterior column count inconsistent for ", comp, "$", dist,
             ": ", paste(unique(n_cols), collapse = ", "))
      
      combined_theta[[comp]][[dist]] <- do.call(rbind, mats)
    }
  }
  
  if (verbose)
    message("Posteriors combined: ",
            sum(sapply(baclava.list, function(x) x$setup$M)),
            " total iterations across ", length(baclava.list), " cycles")
  
  # Combine latent variables
  combined_latent <- do.call(.combineLatent, lapply(baclava.list, "[[", "latent"))
  
  # Assemble combined setup
  combined_setup <- list(
    "t0" = t0,
    "indolent" = indolent,
    "round.age.entry" = round_age,
    "M" = sum(sapply(baclava.list, function(x) x$setup$M)),
    "thin" = thin,
    "burnin" = baclava.list[[1L]]$setup$burnin,
    "I_marginalize" = I_marginalize,
    "cycle" = cycles
  )
  
  # Assemble output
  structure(
    list(
      "theta" = combined_theta,
      "latent" = combined_latent,
      "model" = baclava.list[[length(baclava.list)]]$model,
      "setup" = combined_setup
    ),
    class = c("baclava", "combined")
  )
}