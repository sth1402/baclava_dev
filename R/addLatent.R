#' Add Latent Variable Draws from Chunk Files
#'
#' Reads baclava latent variable chunk files from a single directory,
#'   extracting a specified subset of iterations and optionally a subset
#'   of participants. Returns the baclava object with posteriors align
#'   and latent element populated.
#'
#' @param object A baclava object.
#' @param path A character string. Path to the directory containing latent
#'   chunk files.
#' @param iterations An integer vector. Iteration indices to extract.
#'   Note that iterations will be sorted before retrieval; returned object
#'   may not be in the order provided here.
#' @param participants An optional logical or integer vector. If provided,
#'   subsets participants (columns) from the latent matrices. Logical vectors
#'   are applied directly; integer vectors are used as column indices.
#'   Default \code{NULL} retains all participants.
#' @param pattern A character string. Regex pattern identifying chunk files.
#'   Default assumes all *rds files in the directory are latent data files.
#' @param as.fbm A logical. If true, latent variables are returned as objects
#'   of class FBM.
#' @param fbm.path A character string. If \code{as.fbm = TRUE}, the path in
#'   which the FBM objects should be stored.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @returns The provided baclava object modified to include the latent data,
#'   and if iterations is specified, the posteriors are subset to align with
#'   the latent data.
#'   
#' @include extractLatent.R
#' @export
addLatent <- function(object,
                      path,
                      iterations = NULL,
                      participants = NULL,
                      pattern = NULL,
                      as.fbm = FALSE,
                      fbm.path = tempdir(),
                      verbose = TRUE) {
  
  stopifnot(
    "`object` must be a baclava object" = !missing(object) && is.baclava(object),
    "`path` must be a character string" = !missing(path) && 
      is.vector(path, "character") && length(path) == 1L && nzchar(path),
    "`iterations` must be an integer vector" = is.null(iterations) ||
      {is.vector(iterations, "numeric") && length(iterations) > 0L &&
          all(iterations > 0)},
    "`participants` must be a logical or integer vector" = is.null(participants) ||
      {{is.vector(participants, "logical") || is.vector(participants, "numeric")} &&
          length(participants) > 0L},
    "`pattern` must be a character string" = is.null(pattern) || 
      {is.vector(pattern, "character") && length(pattern) == 1L && nzchar(pattern)},
    "`as.fbm` must be logical" = is.vector(as.fbm, "logical") && length(as.fbm) == 1L,
    "`fbm.path` must be a character string" = !as.fbm ||
      {is.vector(fbm.path, "character") && length(fbm.path) == 1L && nzchar(fbm.path)},
    "`verbose` must be logical" = is.vector(verbose, "logical") && length(verbose) == 1L
  )
  
  M_current <- nrow(object$theta$H[[1L]])
  if (any(iterations > M_current))
    stop("some `iterations` values exceed the number of posterior rows (", 
         M_current, ")", call. = FALSE)
  
  # extract latent data
  lat <- .extractLatent(
    path = path,
    iterations = iterations,
    participants = participants,
    pattern = pattern,
    as.fbm = as.fbm,
    fbm.path = fbm.path,
    verbose = verbose
  )
  
  # subset posteriors to align with extracted iterations
  if (!is.null(iterations)) {
    for (comp in names(object$theta)) {
      for (dist in names(object$theta[[comp]])) {
        object$theta[[comp]][[dist]] <- object$theta[[comp]][[dist]][iterations, , drop = FALSE]
      }
    }
    class(object) <- c(class(object), "REDUCED")
  }
  
  # embed latent in object
  object$latent <- lat
  
  object
}