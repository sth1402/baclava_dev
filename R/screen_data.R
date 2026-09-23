#' Toy Dataset
#' 
#' This toy dataset is provided to facilitate examples and provide an
#'   example of the required input format. Though the data were simulated
#'   under a scenario similar to a real-world breast cancer screening trial,
#'   they should not be interpreted as representing true trial data.
#'
#' @usage data(screen_data)
#' 
#' @format Two datasets are provided. 
#' 
#' \code{data.screen} is a data.frame containing
#'   the following screening information for 89 participants (287 assessments)
#'   \itemize{
#'     \item{id:} A character. Participant ids.
#'     \item{age_assess:} A numeric. The participant age as time of assessment.
#'     \item{disease_detected:} An integer. 1 = disease detected at assessment; 
#'     0 otherwise
#'   }
#'   
#' \code{data.clinical} is a data.frame containing the following information for 89 
#'   participants.
#'   \itemize{
#'     \item{id:} A character. Participant ids.
#'     \item{age_entry:} A numeric. The participant age as time of study entry.
#'     \item{endpoint_type:} A character. One of \{"clinical", "preclinical", "censored"\},
#'     indicating if participant was diagnosed with the disease in the clinical
#'     compartment, was diagnosed in the pre-clinical compartment, or was
#'     censored.
#'     \item{age_endpoint:} A numeric. The participant's age at time the endpoint
#'     was ascertained.
#'   }
#'   
#' @name screen_data
#' @aliases data.clinical data.screen
#' @keywords datasets
NULL


#' Toy Dataset Generator
#' 
#' Data provided with package used
#' n = 100
#' t0 = 30
#' theta <- list(
#'  rate_H = 7.5e-5, shape_H = 2.0,
#'  rate_P = 1/2  , shape_P = 1.0,
#'  beta   = 38.5/(38.5+5.8), psi = 0.4
#'  )
#'  with seed 1234L
#' 
#' @noRd
#' @importFrom dplyr case_when
#' @importFrom stats rpois runif rweibull
#' @keywords internal
#' 
.generateToyData <- function(theta, n, t0) {
  # Helper functions ####
  .draw_sojourn_H <- function(theta) {
    stats::rweibull(1, shape = theta$shape_H, scale = theta$scale_H)
  }

  .draw_sojourn_P <- function(theta) {
    stats::rweibull(1, shape = theta$shape_P, scale = theta$scale_P)
  }

  .compute_compartment <- function(t, tau_HP) {
    dplyr::case_when(t < tau_HP ~ "H", TRUE       ~ "P")
  }

  .screen_result <- function(compartment, theta) {
    if(compartment == "H")  return(FALSE)
    if(compartment == "P")  return(stats::runif(1) < theta$beta)
  }

  #
  # Setup ####
  d_process <- data.frame(
    "person_id" = numeric(n), "sojourn_H" = numeric(n), "tau_HP" = numeric(n),
    "indolent"  = numeric(n), "sojourn_P" = numeric(n), "tau_PC" = numeric(n)
  )

  d_obs_screen <- data.frame("id"  = numeric(0L), 
                             "age_assess" = numeric(0L), 
                             "disease_detected" = integer(0L))

  d_obs_censor <- data.frame("id" = numeric(0L), 
                             "age_entry" = numeric(0L), 
                             "endpoint_type" = character(0L), 
                             "age_endpoint" = numeric(0L))

  #
  ## Biological process ####
  for (i in 1L:n) { # person i
    sojourn_H <- .draw_sojourn_H(theta)
    tau_HP <- t0 + sojourn_H
    indolent <- stats::runif(1L) < theta$psi
    sojourn_P <- if(!indolent) { .draw_sojourn_P(theta) } else { Inf }
    tau_PC <- tau_HP + sojourn_P
    d_process[i, ] <- c(i, sojourn_H, tau_HP, indolent, sojourn_P, tau_PC)
  }

  #
  ## Observation process ####
  AFS <- sample(40:80, n, prob = exp(-(40:80)/5), replace = TRUE)
  dAFS <- stats::rexp(n, 1/5)
  age_death <- pmin(AFS + dAFS, 100.0)
  for (i in 1L:n) {  
    tau_HP <- d_process$tau_HP[i]
    tau_PC <- d_process$tau_PC[i]
  
    # we do not observe individuals that develop a clinical cancer before 
    # their first screen.
    if(tau_PC < AFS[i])  next 
  
    age_screen <- AFS[i]
    j <- 1L
  
    while (age_screen < 1e4) { 
      compartment <- .compute_compartment(age_screen, tau_HP)
      screen_detected <- .screen_result(compartment, theta)
    
      d_obs_screen  <- rbind(d_obs_screen,
                             data.frame("id" = i,
                                        "age_assess" = age_screen,
                                        "disease_detected" = as.integer(screen_detected)))

      if (screen_detected) {  # positive screen first
        d_obs_censor <- rbind(d_obs_censor,
                              data.frame("id" = i,
                                         "age_entry" = AFS[i],
                                         "endpoint_type" = "preclinical",
                                         "age_endpoint" = age_screen))
        break
      }
    
      age_screen <- age_screen + 1.0 + stats::rpois(1, 0.5)
      j <- j + 1L
    
      if (tau_PC < age_screen){  # clinical cancer first
        d_obs_censor <- rbind(d_obs_censor,
                              data.frame("id" = i,
                                         "age_entry" = AFS[i],
                                         "endpoint_type" = "clinical",
                                         "age_endpoint" = tau_PC))
        break
      }
      if (age_death[i] < age_screen) {  # death first
        d_obs_censor <- rbind(d_obs_censor,
                              data.frame("id" = i,
                                         "age_entry" = AFS[i],
                                         "endpoint_type" = "censored",
                                         "age_endpoint" = age_death[i]))
        break
      }
    }
  }

  data.screen <- d_obs_screen
  data.clinical = d_obs_censor
  list("data.screen" = data.screen, "data.clinical" = data.clinical)
}

#' All Cause Mortality Rates
#' 
#' All-cause mortality rates for both male and female, male only, and female
#'   only for the 1974 birth cohort. Data were taken from the CDC Life Tables.
#'   Vital Statistics of the United States, 1974 Life Tables, Vol. II, Section 5. 1976.
#'
#' @usage data(all_cause_rates)
#' 
#' @format  
#' \code{all_cause_rates} is a data.frame containing
#'   the following
#'   \itemize{
#'     \item{Age:} An integer.
#'     \item{both:} A numeric. All-cause mortality rate for combined male and female.
#'     \item{male:} A numeric. Male only all-cause mortality rate.
#'     \item{female:} A numeric. Female only all-cause mortality rate.
#'   }
#'   
#'   
#' @name all_cause_rates
#' @keywords datasets
NULL