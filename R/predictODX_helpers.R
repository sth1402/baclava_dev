#' Summarize Probability Results
#' 
#' @noRd
#' @param ODX.ind.stat A matrix object. The probability of an indolent 
#'   screen-detected cancer while still alive for each screen/model.
#'   {n_screens x M}
#' @param ODX.prog.stat A matrix object. The probability to have a progressive
#'   and overdiagnosed screen-detected cancer while still alive for each screen/model.
#'   {n_screens x M}
#' @param SD.stat A matrix object. The probability of a screen-detected cancer 
#'   while still alive. {n_screens x M}
#'   
#' @returns A matrix providing the average probability for a screen-detected
#'   cancer being 1) total: indolent or progressive and overdiagnosed;
#'   2) indolent; and 3) mortality: progressive and overdiagnosed and their
#'   95% credible intervals
#'   
#' @importFrom stats quantile
#' @keywords internal
wrap <- function(ODX.ind.stat, ODX.prog.stat, SD.stat) {
  
  SD.stat <- pmax(SD.stat, 1e-12)
  
  # overall
  # among all screen detected cancer (SD.stat), the fraction that or ODX
  # from either indolent tumors (ODX.ind.stat) or progressive tumors that would 
  # die of other causes before clinical presentation (ODX.prog.stat)
  a_term <- c(100.0 * mean((ODX.ind.stat + ODX.prog.stat) / SD.stat),
              100.0 * stats::quantile((ODX.ind.stat + ODX.prog.stat) / SD.stat, 
                                      probs = c(0.025, 0.5, 0.975)))
  # indolent
  # among all screen detected cancer (SD.stat), the fraction that or ODX
  # from indolent tumors (ODX.ind.stat)
  b_term <- c(100.0 * mean(ODX.ind.stat / SD.stat),
              100.0 * stats::quantile(ODX.ind.stat / SD.stat, 
                                      probs = c(0.025, 0.5, 0.975)))
  # mortality
  # among all screen detected cancer (SD.stat), the fraction that or ODX
  # from progressive tumors that would die of other causes before clinical 
  # presentation (ODX.prog.stat)
  c_term <- c(100.0 * mean(ODX.prog.stat / SD.stat),
              100.0 * stats::quantile(ODX.prog.stat / SD.stat, 
                                      probs = c(0.025, 0.5, 0.975)))
  
  
  tab <- rbind(a_term, b_term, c_term)
  colnames(tab)[1L] <- "mean"
  tab <- tab[, c(1L, 2L, 4L, 3L)]
  rownames(tab) <- c("total", "indolent", "mortality")
  
  tab
}