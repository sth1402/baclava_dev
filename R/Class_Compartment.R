#' @noRd
#' @include Class_Distribution.R
#' @keywords internal
Compartment <- R6::R6Class(
  "Compartment",
  public = list(
    #' @field distributions (`list()`)\cr
    #' list of Distribution(s); named for stratified compartments
    distributions = NULL,
    #' @field subset_var (`character(1)`)\cr
    #' NA or character(1) (for stratified)
    subset_var = NA_character_,          
    #' @field routing (`list()`)\cr
    #' filled during finalization -- links compartment to data structures
    routing = NULL,
    
    #' @description
    #' Defines how a model or distribution parameter. 
    #'
    #' @param distributions a Distribution or a named list of Distribution objects
    #' @param subset_var NULL or character(1)
    initialize = function(distributions,       # Distribution OR named list of Distribution
                          subset_var = NA_character_) { # NULL or character(1)

      # Normalize distributions into a *list*
      if (is.Distribution(distributions)) {
        self$distributions <- list(distributions)
        names(self$distributions) <- NULL
      } else {
        stopifnot(is.list(distributions), length(distributions) >= 1L,
                  all(vapply(distributions, is.Distribution, logical(1L))))
        
        # distributions must be named if more than 1 provided
        if (length(distributions) > 1L) {
          if (is.null(names(distributions))) 
            stop("when multiple distributions are specified, they must be named", call. = FALSE)
        } else {
          distributions <- unname(distributions)
        }
        self$distributions <- distributions
      }
      
      tags <- NULL
      is_fixed <- NULL
      for (i in seq_len(length(self$distributions))) {
        distribution <- self$distributions[[i]]
        dist_params <- distribution$parameters # each is a Model
        for (j in seq_len(length(dist_params))) {
          model_params <- dist_params[[j]]$parameters
          for (k in seq_len(length(model_params))) {
            param <- model_params[[k]]
            tags <- c(tags, param$tag)
            is_fixed <- c(is_fixed, param$fixed)
          }
        }
      }
      
      if (any(duplicated(tags))) stop("parameter tags must be unique")
      
      self$subset_var <- as.character(subset_var)
      if (length(self$subset_var) == 0L && length(self$distributions) == 1L)
        self$subset_var <- NA_character_
      
      if (length(self$subset_var) != 1L) 
        stop("subset variable should be NA or a single character string")
      
      # distributions must be named for stratified compartments
      if (!is.na(self$subset_var)) {
        if (is.null(names(self$distributions)) ||
            any(!nzchar(names(self$distributions)))) {
          stop("when subset_var is provided, distributions must be a named list", call. = FALSE)
        }
      }
    },
    
    toList = function() {
      base <- list(
        "meta_for" = "Compartment",
        "subset_var" = self$subset_var, 
        "labels" = names(self$distributions),
        "distributions" = lapply(self$distributions, function(x) x$toList())
      )
      names(base$distributions) <- names(self$distributions)
      
      base
    },
    
    print = function(...) {
      cat("<Compartment>\n")
      if (!is.na(self$subset_var)) cat("  subset_var:", self$subset_var, "\n")
      if (!is.null(names(self$distributions))) {
        cat("  labels:", paste(names(self$distributions), collapse=", "), "\n")
      }
      invisible(self)
    },
    
    setValue = function(value) {
      for (i in seq_along(self$distributions)) {
        self$distributions[[i]]$setValue(value)
      }
      invisible(self)
    },
    
    getPriorDist = function() {
      lapply(self$distributions, function(x) x$getPriorDist())
    },
    getPriorParams = function() {
      lapply(self$distributions, function(x) x$getPriorParams())
    },
    getTags = function() {
      lapply(self$distributions, function(x) x$getTags())
    }
    
  )
)

#' Define a model compartment (single or stratified)
#'
#' Convenience constructor for a \code{Compartment} (R6) object describing how a
#'  model compartment is composed (single distribution or stratified by a data
#'  column).
#'
#' @param name Character(1). Compartment name (e.g., \code{"H"}, \code{"P"}, \code{"psi"}).
#' @param distributions Either a single \code{\link[=Distribution]{Distribution}} object
#'   (for a single compartment) or a list of \code{Distribution} objects.
#' @param subset_var Optional character(1). Name of a column in the data that
#'   defines strata for a \emph{stratified} compartment. If provided, the
#'   \code{distributions} list must be named and those names must match the
#'   unique values of \code{subset_var} (after coercion to character).
#'
#' @details
#' \strong{Finalization / routing}
#' \code{$finalize(data_assess, data_clinical, where, id_cols)} resolves routing for
#' stratified compartments by locating \code{subset_var} in either the assessment data
#' or clinical data (or mapping clinical \code{subset_var} to assessment via id
#' columns when \code{where = "clinical_to_assess"}). The method validates that the
#' observed levels of \code{subset_var} match the names of \code{distributions} and
#' stores a routing map usable by downstream C++.
#'
#' \strong{Serialization}
#' \code{$toList()} returns a stable list with labels, distributions,
#' and routing metadata (for stratified compartments).
#'
#' @return An R6 \code{Compartment} object with fields:
#' \itemize{
#'   \item \code{name}, \code{subset_var}, 
#'   \item \code{distributions} (list of \code{Distribution})
#'   \item \code{routing} (for stratified)
#' }
#'
#' @examples
#' 
#' # define several parameteters to be updated for convenience
#' m1 = param(0.5, "m1", prop_scale = "log", prior = dlnorm_p(log(0.5), 0.2))
#' m2 = param(1.5, "m2", prop_scale = "log", prior = dlnorm_p(log(1.5), 0.2))
#' 
#' s1 = param(2, "s1", prop_scale = "log", prior = dlnorm_p(log(2), 0.2))
#' s2 = param(2, "s2", prop_scale = "log", prior = dlnorm_p(log(2), 0.2))
#' 
#' # Single compartment; defaults to a single MH update include m1 and s1
#' c1 <- compartment(
#'   distributions = weibull_d(mu = m1, shape = s1))
#'
#' # Stratified compartment (two strata defined by a data column with values A/B)
#' dlist <- list(A = gamma_d(mu = m1, shape = s1),
#'               B = gamma_d(mu = m2, shape = s2))
#' c2 <- compartment(distributions = dlist, 
#'                   subset_var = "group")
#'
#' @seealso \code{\link[=Distribution]{Distribution}}
#' @export
compartment <- function(distributions,      # Distribution OR named list of Distribution
                        subset_var = NA_character_) {  # NULL or character(1)

  stopifnot("distributions must be provided" = !missing(distributions))
  
  Compartment$new(distributions, subset_var)
}

#' @describeIn compartment Inheritance test
is.Compartment <- function(x) inherits(x, "Compartment")

# helper: crosswalk clinical -> assess by id
.comp_match_clin_to_assess <- function(assess_ids, clinical_ids, clinical_vec, nm_assess, nm_clin) {
  stopifnot(length(assess_ids) >= 1L, length(clinical_ids) >= 1L, length(clinical_vec) == length(clinical_ids))
  m <- match(assess_ids, clinical_ids)
  if (anyNA(m)) {
    miss <- unique(assess_ids[is.na(m)])
    stop("Compartment finalize(): some ", nm_assess, " ids not found in ", nm_clin, " (e.g., ",
         paste(utils::head(miss, 5L), collapse=", "), ")", call. = FALSE)
  }
  clinical_vec[m]
}

# need to align the provided data with the specified model
# this is called within baclava to link models to individuals
Compartment$set("public", "finalize", function(data_assess,
                                               data_clinical,
                                               where = c("auto", "assess", "clinical"),
                                               id_cols = list(assess="id", clinical="id")) {
  
  # where to look for the subset
  where <- match.arg(where)
  
  # nothing to route for single
  if (length(self$distributions) == 1L) {
    self$routing <- list(
      source = NULL,
      id_cols = NULL,
      strata_labels = NULL,
      row_to_stratum = NULL,
      rows_by_stratum = NULL,
      n_rows = NA_integer_
    )
    return(invisible(self))
  }
  
  # stratified: validate names & locate vector
  dnames <- names(self$distributions)
  if (is.null(dnames) || any(!nzchar(dnames))) {
    stop("stratified distributions must be a named list", call. = FALSE)
  }
  subset_var <- self$subset_var
  if (!is.character(subset_var) || length(subset_var) != 1L || !nzchar(subset_var)) {
    stop("invalid 'subset_var'", call. = FALSE)
  }
  
  # decide source
  pick_source <- function() {
    # if told where to look - just accept it
    if (where != "auto") return(where)
    
    # is the subset_var in the column headers of data_assess
    in_assess <- !is.null(data_assess) && subset_var %in% names(data_assess)
    
    # is the subset_var in the column headers of data_clinical
    in_clinical <- !is.null(data_clinical) && subset_var %in% names(data_clinical)
    
    # if it is in both, stop with error -- ambiguous
    if (in_assess && in_clinical) {
      stop("subset_var ", subset_var, " cannot be present in both data sets", call. = FALSE)
    } else if (in_assess)   return("assess")
    else if (in_clinical)   return("clinical")
    else stop(subset_var, "' not found in provided data", call. = FALSE)
  }
  source <- pick_source()
  
  # extract raw labels
  raw <- switch(source,
                assess = data_assess[[subset_var]],
                clinical = data_clinical[[subset_var]],
                stop("unknown routing source"))
  
  if (length(raw) == 0L) stop("empty routing vector", call. = FALSE)
  if (anyNA(raw)) stop("NA values in '", subset_var, "' are not allowed.", call. = FALSE)
  
  # coerce to character then factor with distribution names as required levels
  raw_chr <- as.character(raw)
  if (!setequal(dnames, unique(raw_chr))) {
    extra <- setdiff(unique(raw_chr), dnames)
    miss  <- setdiff(dnames, unique(raw_chr))
    msg <- sprintf("names(distributions) != levels in '%s'. extra: {%s}; missing: {%s}",
                   subset_var,
                   paste(extra, collapse=", "),
                   paste(miss, collapse=", "))
    stop(msg, call. = FALSE)
  }
  f <- factor(raw_chr, levels = dnames)
  
  # reorder self$distributions to match factor levels (if user handed a different order)
  self$distributions <- self$distributions[levels(f)]
  
  # build integer map and row indices by stratum
  z <- as.integer(f)  # 1...K
  rows_by_stratum <- split(seq_along(z), z)
  K <- length(levels(f))
  if (length(rows_by_stratum) < K) {
    tmp <- vector("list", K)
    for (k in seq_len(K)) tmp[[k]] <- rows_by_stratum[[as.character(k)]]
    rows_by_stratum <- tmp
  }
  
  # stash routing for Rcpp
  self$routing <- list(
    "source" = source,
    "id_cols" = id_cols,
    "strata_labels" = levels(f),
    "row_to_stratum" = z,
    "rows_by_stratum" = rows_by_stratum,
    "n_rows" = length(z)
  )
  
  invisible(self)
})


.comp_labels <- function(comp) {
  if (is.null(comp)) return(character(0))
  
  if (length(comp$distributions) == 1L) {
    "all"
  } else {
    if (is.null(comp$routing$strata_labels)) {
      as.character(1L:length(comp$distributions))
    } else {
      comp$routing$strata_labels
    }
  }

}

