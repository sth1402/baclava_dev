test_that(".computeEndpoints() returns expected errors", {
  
  data_obj <- list()
  valid <- list()
  valid$endpoint_type = 1L
  valid$n = 5L
  valid$irateP = 2L
  valid$endpoint_time = c(50, 55, 60, 65, 70)
  valid$ages_screen = list(values = c(40, 45, 50, 
                                      40, 55, 
                                      60, 
                                      50, 55, 60, 65,
                                      50, 60, 63, 68, 70),
                           starts = c(1L, 4L, 6L, 7L, 11L) - 1L,
                           ends = c(3L, 5L, 6L, 10L, 15L) - 1L,
                           lengths = c(3L, 2L, 1L, 4L, 5L))
#  valid$age_entry = c(40, 40, 60, 50, 50)
  valid$n_screen_positive = c(0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1)
  
  expect_error(.computeEndpoints(),
               "`data.obj` does not contain required data")
  expect_error(.computeEndpoints(data.obj = valid$censor_time))
  expect_error(.computeEndpoints(data.obj = list("a" = valid, "b" = valid, "c" = valid)),
               "`data.obj` does not contain required data")
  expect_error(.computeEndpoints(data.obj = list(valid, "a" = 10)),
               "`data.obj` does not contain required data")
  
  expect_error(.computeEndpoints(data.obj = valid),
               "`t0` must be a scalar numeric")
  expect_error(.computeEndpoints(data.obj = valid, t0 = c(1,2)),
               "`t0` must be a scalar numeric")
  expect_error(.computeEndpoints(data.obj = valid, t0 = "a"),
               "`t0` must be a scalar numeric")
  
  data_obj <- valid
  
  # problems in screen
  expect_error(.computeEndpoints(data.obj = data_obj, t0 = 70),
               "assessment ages younger than t0 encountered")
  
  data_obj$endpoint_time = c(50, 55, 60, 65, 50)
  expect_error(.computeEndpoints(data.obj = data_obj, t0 = 30),
               "assessments after endpoint encountered")
  data_obj <- valid
  
  # problems in clinical
  valid$endpoint_type <- 3L
  data_obj <- valid
  data_obj$ages_screen$values <- data_obj$ages_screen$values - 20

  expect_error(.computeEndpoints(data.obj = data_obj, t0 = 40),
               "assessment ages younger than t0 encountered")
  data_obj <- valid
  
  data_obj$endpoint_time = c(50, 55, 60, 65, 50)
  expect_error(.computeEndpoints(data.obj = data_obj, t0 = 30),
               "assessments after endpoint encountered")
  data_obj <- valid
  
  # problems in censored
  valid$endpoint_type <- 2L
  data_obj <- valid
  data_obj$ages_screen$values <- data_obj$ages_screen$values - 20
  expect_error(.computeEndpoints(data.obj = data_obj, t0 = 40),
               "assessment ages younger than t0 encountered")
  data_obj <- valid
  
  data_obj$endpoint_time = c(50, 55, 60, 65, 50)
  expect_error(.computeEndpoints(data.obj = data_obj, t0 = 30),
               "assessments after endpoint encountered")
  data_obj <- valid
  
})

test_that(".computeEndpoints() returns expected results", {
  
  data_obj <- list()
  valid <- list()
  valid$endpoint_type = 1L
  valid$irateP = 1L
  valid$n = 5L
  valid$endpoint_time = c(50, 55, 60, 65, 70)
  valid$ages_screen = list(values = c(40, 45, 50, 
                                      40, 55, 
                                      60, 
                                      50, 55, 60, 65,
                                      50, 60, 63, 68, 70),
                           starts = c(1L, 4L, 6L, 7L, 11L) - 1L,
                           ends = c(3L, 5L, 6L, 10L, 15L) - 1L,
                           lengths = c(3L, 2L, 1L, 4L, 5L))
#  valid$age_entry = c(40, 40, 60, 50, 50)
  valid$n_screen_positive = c(0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1)
  
  valid_result_screen <- list("values"  = c(30, 40, 45, 50, 
                                            30, 40, 55, 
                                            30, 60, 
                                            30, 50, 55, 60, 65,
                                            30, 50, 60, 63, 68, 70),
                              "starts"  = c(1L, 5L, 8L, 10L, 15L) - 1L,
                              "ends"    = c(4L, 7L, 9L, 14L, 20L) - 1L,
                              "lengths" = c(4L, 3L, 2L, 5L, 6L))

  valid_result_censored <- list("values"  = c(30, 40, 45, 50, Inf,
                                             30, 40, 55, Inf,
                                             30, 60, Inf,
                                             30, 50, 55, 60, 65, Inf,
                                             30, 50, 60, 63, 68, 70, Inf),
                                "starts"  = c(1L, 6L, 10L, 13L, 19L) - 1L,
                                "ends"    = c(5L, 9L, 12L, 18L, 25L) - 1L,
                              "lengths" = c(5L, 4L, 3L, 6L, 7L))
  
  data_obj <- list("screen" = valid, "clinical" = valid, "censored" = valid)
  
  expect_equal(.computeEndpoints(data.obj = valid, t0 = 30),
               valid_result_screen)
  valid$endpoint_type <- 3L
  expect_equal(.computeEndpoints(data.obj = valid, t0 = 30),
               valid_result_screen)
  valid$endpoint_type <- 2L
  expect_equal(.computeEndpoints(data.obj = valid, t0 = 30),
               valid_result_censored)
})
