#' @noRd
#' @keywords internal
is.fbm <- function(x) inherits(x, "FBM")

#' Extract Latent Variable Draws from Chunk Files
#'
#' Reads baclava latent variable chunk files from a single directory,
#'   extracting a specified subset of iterations and optionally a subset
#'   of participants. Returns a list suitable for use with \code{cohortODX}.
#'
#' @noRd
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
#' @returns A list with elements:
#'   \describe{
#'     \item{\code{indolent}}{Integer matrix or FBM matrix of dimension
#'       \code{num_iterations} x \code{n_participants}. Indolence indicators.}
#'     \item{\code{tau}}{Numeric matrix of dimension
#'       \code{num_iterations} x \code{n_participants}. Preclinical onset ages.}
#'     \item{\code{iterations}}{Integer vector of the global iteration
#'       indices that were successfully extracted. If all iterations kept, NULL}
#'   }
.extractLatent <- function(path,
                           iterations = NULL,
                           participants = NULL,
                           pattern = NULL,
                           as.fbm = FALSE,
                           fbm.path = tempdir(),
                           object = NULL,
                           verbose = TRUE) {
  
  stopifnot(
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
    "`object` must be a baclava object" = is.null(object) || is.baclava(object),
    "`verbose` must be logical" = is.vector(verbose, "logical") && length(verbose) == 1L
  )
  
  # Input validation
  if (!dir.exists(path)) stop("directory not found: ", path)
  
  if (!is.null(iterations)) iterations <- as.integer(iterations) |> sort()

  if (is.null(pattern)) pattern <- "-latent-chunk-\\d+\\.rds$"

  chunk_files <- dir(path, pattern = pattern, full.names = FALSE)
  
  if (length(chunk_files) == 0L)
    stop("no files matching pattern '", pattern, "' found in: ", path)
  
  # sort by chunk number
  chunk_nums <- as.integer(
    regmatches(chunk_files,
               regexpr("(?<=-latent-chunk-)\\d+(?=\\.rds$)",
                       chunk_files, perl = TRUE))
  )
  chunk_files <- chunk_files[order(chunk_nums)]
  chunk_nums  <- sort(chunk_nums)
  n_chunks <- length(chunk_files)
  
  if (verbose)
    message("Found ", length(chunk_files), " chunk files",
            " (chunks ", min(chunk_nums), "\u2013", max(chunk_nums), ")",
            " in: ", path)
  
  if (as.fbm) {
    
    first_chunk <- readRDS(file.path(path, chunk_files[1]))
    
    n_rows <- if (!is.null(iterations)) {
      length(iterations)
    } else {
      if (is.null(object)) stop("`object` cannot be NULL if iterations are not specified", call. = FALSE)
      nrow(object$theta$H[[1L]])
    }
    
    n_cols <- if (!is.null(participants)) {
      length(participants)
    } else {
      ncol(first_chunk$Z_tau)
    }
  
    Z_tau_obj <- big.matrix(n_rows, n_cols, backingfile = "Z_tau.bin",
                            descriptorfile = "Z_tau.desc", backingpath = fbm.path)
    
    if ("Z_indolent" %in% names(first_chunk))
      Z_ind_obj <- big.matrix(n_rows, n_cols, backingfile = "Z_indolent.bin",
                              descriptorfile = "Z_indolent.desc", backingpath = fbm.path)
  } else {
    Z_tau_obj <- vector("list", n_chunks)
    Z_ind_obj <- vector("list", n_chunks)
  }
  
  # Main extraction loop
  
  iter_offset <- 0L  # rows seen so far
  
  .verifyParticipant <- function(participants, n) {
    if (is.numeric(participants) && !all(participants <= n)) 
      stop("`participants` vector contains values that exceed chunk dimensions", call. = FALSE)
    if (is.logical(participants) && length(participants) != n) 
      stop("`participants` vector does not match chunk dimensions", call. = FALSE)
  }
  
  has_indolent <- logical(n_chunks)
  
  row_cursor <- 1
  for (i in seq_len(n_chunks)) {
    
    chunk_start <- iter_offset + 1L
    
    # read chunk data
    chunk_data <- readRDS(file.path(path, chunk_files[i]))
    has_indolent[i] <- "Z_indolent" %in% names(chunk_data) && !is.null(chunk_data$Z_indolent)
    
    n_iter_in_chunk <- nrow(chunk_data$Z_tau)
    chunk_end <- iter_offset + n_iter_in_chunk
    
    if (n_iter_in_chunk == 0L) {
      warning("  chunk ", chunk_nums[i], " contains no data", call. = FALSE)
      next
    }
    
    ind_mat <- NULL
    tau_mat <- NULL
    if (!is.null(iterations)) {
      # find which requested iterations fall in this chunk
      hits <- which(iterations >= chunk_start & iterations <= chunk_end)
   
      if (length(hits) > 0L) {
        local_rows <- iterations[hits] - iter_offset
      
        # apply participant subsetting
        if (is.null(participants)) {
          if (has_indolent[i]) ind_mat <- chunk_data$Z_indolent[local_rows, , drop = FALSE]
          tau_mat <- chunk_data$Z_tau[local_rows, , drop = FALSE]
        } else {
          .verifyParticipant (participants, ncol(chunk_data$Z_tau))
          if (has_indolent[i]) ind_mat <- chunk_data$Z_indolent[local_rows, participants, drop = FALSE]
          tau_mat <- chunk_data$Z_tau[local_rows, participants, drop = FALSE]
        }
      }
      if (verbose)
        message("  chunk ", chunk_nums[i], ": extracted ", length(hits),
                " draw(s) [rows ", chunk_start, "\u2013", chunk_end, "]")
    } else {
      # apply participant subsetting
      if (is.null(participants)) {
        if (has_indolent[i]) ind_mat <- chunk_data$Z_indolent
        tau_mat <- chunk_data$Z_tau
      } else {
        .verifyParticipant (participants, ncol(chunk_data$Z_tau))
        if (has_indolent[i]) ind_mat <- chunk_data$Z_indolent[, participants, drop = FALSE]
        tau_mat <- chunk_data$Z_tau[, participants, drop = FALSE]
      }
      
      if (verbose)
        message("  chunk ", chunk_nums[i], ": extracted ", n_iter_in_chunk,
                " draw(s) [rows ", chunk_start, "\u2013", chunk_end, "]")
    }

    if (as.fbm) {
      if (!is.null(tau_mat)) {
        rows <- row_cursor:(row_cursor +  nrow(tau_mat) - 1)
        Z_tau_obj[rows, ] <- tau_mat
        if (!is.null(ind_mat)) Z_ind_obj[rows, ] <- ind_mat
        row_cursor <- row_cursor + nrow(tau_mat)
      }
    } else {
      Z_ind_obj[[i]] <- ind_mat
      Z_tau_obj[[i]] <- tau_mat
    }
    
    iter_offset <- chunk_end
    rm(chunk_data); gc()
  }
  
  if (any(has_indolent) && !all(has_indolent)) {
    stop("some latent files have indolent status; others do not", call. = FALSE)
  }
  has_indolent <- all(has_indolent)
  
  if (!as.fbm) {
    if (has_indolent) Z_ind_obj <- do.call(rbind, Z_ind_obj)
    Z_tau_obj <- do.call(rbind, Z_tau_obj)
  }
  
  # Assemble output
  if (!is.null(iterations)) {
    if (nrow(Z_tau_obj) == 0L)
      warning("no iterations were extracted; check that `iterations` indices",
              " fall within [", 1L, ", ", iter_offset, "]")
  
    if (nrow(Z_tau_obj) < length(iterations) && verbose)
      message("Note: requested ", length(iterations), 
              " iterations, extracted ", nrow(Z_tau_obj),
              " (some may fall outside this directory's range)")
  }
  
  if (is.logical(participants)) participants <- which(participants)
  
  structure(
    list(
      "indolent" = if (has_indolent) Z_ind_obj else NULL,
      "tau" = Z_tau_obj,
      "iterations" = iterations,
      "participants" = participants),
    class = "Latent")
}

#' @noRd
#' @keywords internal
is.Latent <- function(x) inherits(x, "Latent")


#' Combine Latent Draw Lists Across Directories
#'
#' Combines the output of multiple \code{\link{extractLatent}} calls into a
#' single list suitable for \code{cohortODX}. Rows are ordered by global
#' iteration index.
#'
#' @noRd
#' @param ... Two or more Latent objects.
#'
#' @returns A list with the same structure as \code{extractLatent} output,
#'   with \code{indolent} and \code{tau} matrices row-bound and sorted by iteration
#'   index.
#'   
.combineLatent <- function(...) {
  
  latent.list <- list(...)
  
  if (length(latent.list) <= 1L)
    stop("functions requires at least 2 Latent objects")
  
  stopifnot(
    "All `latent.list` elements must be Latent objects" =
      all(sapply(latent.list, is.Latent))
  )
  
  participant_list <- lapply(latent.list, "[[", "participants")
  tst <- all(sapply(participant_list, identical, participant_list[[1L]]))
  if (!tst) stop("participant vector is not the same across all Latent objects")
  
  has_indolent <- all(sapply(latent.list, 
                             function (x) all(c("indolent") %in% names(x))))
  
  if (has_indolent) {
    # all ZI matrices must have same number of columns
    n_cols_ZI <- sapply(latent.list, function(x) ncol(x$indolent))
    if (length(unique(n_cols_ZI)) > 1L)
      stop("indolent matrices have inconsistent column counts: ",
           paste(unique(n_cols_ZI), collapse = ", "))
  }
  
  n_cols_ZT <- sapply(latent.list, function(x) ncol(x$tau))
  if (length(unique(n_cols_ZT)) > 1L)
    stop("tau matrices have inconsistent column counts: ",
         paste(unique(n_cols_ZT), collapse = ", "))
  
  n_iters <- sapply(latent.list, function(x) is.null(x$iterations))
  if (any(n_iters) && !all(n_iters))
    stop("cannot combine: some extractions specified `iterations` and others did not")
  
  # combine latent data
  if (has_indolent) ZI_combined <- do.call(rbind, lapply(latent.list, `[[`, "indolent"))
  ZT_combined <- do.call(rbind, lapply(latent.list, `[[`, "tau"))
  
  iters_combined <- unlist(lapply(latent.list, `[[`, "iterations"))
  if (!is.null(iters_combined)) {
    ord <- order(iters_combined)
    if (has_indolent) ZI_combined <- ZI_combined[ord, , drop = FALSE]
    ZT_combined <- ZT_combined[ord, , drop = FALSE]
    iters_combined <- iters_combined[ord]
  }  

  structure(
    list("indolent" = if (has_indolent) ZI_combined else NULL, 
       "tau" = ZT_combined, 
       "iterations" = iters_combined, 
       "participants" = latent.list[[1L]]$participants),
    class = c("Latent", "Combined"))
}