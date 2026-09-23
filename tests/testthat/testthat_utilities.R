test_that("`processSensitivity()` returns expected results no screen or am specification equal sizes", {
  
  data.assess <- withr::with_seed(
    1234,
    data.frame("id" = 1:100, 
               "age_assess" = 1:100,
               "disease_detected" = rbinom(100, 1, 0.4))
    )
  
  data.clinical <- withr::with_seed(
    1234,
    data.frame("id" = 1:100, 
               "age_entry" = 1:100,
               "age_endpoint" = 2:101,
               "endpoint_type" = rep("clinical", 100))
  )
  
  theta_0 <- list("beta" = 0.8)
  prior <- list("a_beta" = 1, "b_beta" = 2)
  
  expected <- list("screen.type" = rep("screen", 100),
                   "arm" = rep("all", 100),
                   "theta" = list("beta" = c("screen" = 0.8)),
                   "prior" = list("a_beta" = c("screen" = 1),
                                  "b_beta" = c("screen" = 2)))
  expect_equal(.processSensitivity(data.assess, data.clinical, theta_0, prior), 
               expected)
  
})


test_that("`processSensitivity()` returns expected errors no screen or am specification unequal sizes", {
  
  data.assess <- withr::with_seed(
    1234,
    data.frame("id" = 1:200, 
               "age_assess" = 1:200,
               "disease_detected" = rbinom(200, 1, 0.4))
  )
  
  data.clinical <- withr::with_seed(
    1234,
    data.frame("id" = 1:100, 
               "age_entry" = 1:100,
               "age_endpoint" = 2:101,
               "endpoint_type" = rep("clinical", 100))
  )
  
  theta_0 <- list("beta" = 0.8)
  prior <- list("a_beta" = 1, "b_beta" = 2)
  
  expected <- list("screen.type" = rep("screen", 200),
                   "arm" = rep("all", 100),
                   "theta" = list("beta" = c("screen" = 0.8)),
                   "prior" = list("a_beta" = c("screen" = 1),
                                  "b_beta" = c("screen" = 2)))
  expect_equal(.processSensitivity(data.assess, data.clinical, theta_0, prior), 
               expected)
  
})

test_that("`processSensitivity()` returns expected results screen_type given but not used", {
  
  data.assess <- withr::with_seed(
    1234,
    data.frame("id" = 1:200, 
               "age_assess" = 1:200,
               "disease_detected" = rbinom(200, 1, 0.4),
               "screen_type" = c(rep("1", 100), rep("2", 100)))
  )
  
  data.clinical <- withr::with_seed(
    1234,
    data.frame("id" = 1:100, 
               "age_entry" = 1:100,
               "age_endpoint" = 2:101,
               "endpoint_type" = rep("clinical", 100))
  )
  
  theta_0 <- list("beta" = 0.8)
  prior <- list("a_beta" = 1, "b_beta" = 2)
  
  expected <- list("screen.type" = rep("screen", 200),
                   "arm" = rep("all", 100),
                   "theta" = list("beta" = c("screen" = 0.8)),
                   "prior" = list("a_beta" = c("screen" = 1),
                                  "b_beta" = c("screen" = 2)))
  expect_warning(out <- .processSensitivity(data.assess, data.clinical, theta_0, prior), 
                 "*** only 1 sensitivity type specified in theta; screen_type ignored ***",
                 fixed = TRUE)
  expect_equal(out, expected)
  
})

test_that("`processSensitivity()` returns expected results 1 screen_type given but not used", {
  
  data.assess <- withr::with_seed(
    1234,
    data.frame("id" = 1:200, 
               "age_assess" = 1:200,
               "disease_detected" = rbinom(200, 1, 0.4),
               "screen_type" = c(rep("1", 200)))
  )
  
  data.clinical <- withr::with_seed(
    1234,
    data.frame("id" = 1:100, 
               "age_entry" = 1:100,
               "age_endpoint" = 2:101,
               "endpoint_type" = rep("clinical", 100))
  )
  
  theta_0 <- list("beta" = 0.8)
  prior <- list("a_beta" = 1, "b_beta" = 2)
  
  expected <- list("screen.type" = rep("1", 200),
                   "arm" = rep("all", 100),
                   "theta" = list("beta" = c("1" = 0.8)),
                   "prior" = list("a_beta" = c("1" = 1),
                                  "b_beta" = c("1" = 2)))
  
  expect_equal(.processSensitivity(data.assess, data.clinical, theta_0, prior), expected)
  
})

test_that("`processSensitivity()` returns expected results multiple beta no arm or screen", {
  
  data.assess <- withr::with_seed(
    1234,
    data.frame("id" = 1:200, 
               "age_assess" = 1:200,
               "disease_detected" = rbinom(200, 1, 0.4))
  )
  
  data.clinical <- withr::with_seed(
    1234,
    data.frame("id" = 1:100, 
               "age_entry" = 1:100,
               "age_endpoint" = 2:101,
               "endpoint_type" = rep("clinical", 100))
  )
  
  theta_0 <- list("beta" = c("A" = 0.8, "B" = 0.6))
  prior <- list("a_beta" = c("A" = 0.18, "B" = 0.16), 
                "b_beta" = c("A" = 0.28, "B" = 0.26))
  
  expect_error(.processSensitivity(data.assess, data.clinical, theta_0, prior), 
               "could not connect theta$beta to one of data.clinical$arm or data.assess$screen_type; verify input",
               fixed = TRUE)
  
})

test_that("`processSensitivity()` returns expected results multiple beta arm bad screen", {
  
  data.assess <- withr::with_seed(
    1234,
    data.frame("id" = 1:200, 
               "age_assess" = 1:200,
               "disease_detected" = rbinom(200, 1, 0.4))
  )
  
  data.clinical <- withr::with_seed(
    1234,
    data.frame("id" = 1:100, 
               "age_entry" = 1:100,
               "age_endpoint" = 2:101,
               "endpoint_type" = rep("clinical", 100),
               "arm" = c(rep("A", 50), rep("B", 50)))
  )
  
  theta_0 <- list("beta" = c("A" = 0.8, "B" = 0.6))
  prior <- list("a_beta" = c("A" = 0.18, "B" = 0.16), 
                "b_beta" = c("A" = 0.28, "B" = 0.26))
  
  expected <- list("screen_type" = rep("1", 200),
                   "arm" = rep("all", 100),
                   "theta" = list("beta" = c("1" = 0.8)),
                   "prior" = list("a_beta" = c("1" = 1),
                                  "b_beta" = c("1" = 2)))
  
  expect_error(.processSensitivity(data.assess, data.clinical, theta_0, prior), 
               "some screens do not fall into defined arm designations",
               fixed = TRUE)
  
})

test_that("`processSensitivity()` returns expected results multiple beta arm good screen", {
  
  data.assess <- withr::with_seed(
    1234,
    data.frame("id" = c(1:100, 1:100), 
               "age_assess" = 1:200,
               "disease_detected" = rbinom(200, 1, 0.4))
  )
  
  data.clinical <- withr::with_seed(
    1234,
    data.frame("id" = 1:100, 
               "age_entry" = 1:100,
               "age_endpoint" = 2:101,
               "endpoint_type" = rep("clinical", 100),
               "arm" = c(rep("A", 50), rep("B", 50)))
  )
  
  theta_0 <- list("beta" = c("A" = 0.8, "B" = 0.6))
  prior <- list("a_beta" = c("A" = 0.18, "B" = 0.16), 
                "b_beta" = c("A" = 0.28, "B" = 0.26))
  
  expected <- list("screen.type" = c(rep("A", 50), rep("B", 50), rep("A", 50), rep("B", 50)),
                   "arm" = c(rep("A", 50), rep("B", 50)),
                   "theta" = list("beta" = c("A" = 0.8, "B" = 0.6)),
                   "prior" = list("a_beta" = c("A" = 0.18, "B" = 0.16),
                                  "b_beta" = c("A" = 0.28, "B" = 0.26)))
  
  expect_equal(.processSensitivity(data.assess, data.clinical, theta_0, prior), 
               expected)
  
})


test_that("`processSensitivity()` returns expected results multiple beta screen no arm", {
  
  data.assess <- withr::with_seed(
    1234,
    data.frame("id" = c(1:100, 1:100), 
               "age_assess" = 1:200,
               "disease_detected" = rbinom(200, 1, 0.4),
               "screen_type" = c(rep("A", 100), rep("B", 100)))
  )
  
  data.clinical <- withr::with_seed(
    1234,
    data.frame("id" = 1:100, 
               "age_entry" = 1:100,
               "age_endpoint" = 2:101,
               "endpoint_type" = rep("clinical", 100))
  )
  
  theta_0 <- list("beta" = c("A" = 0.8, "B" = 0.6))
  prior <- list("a_beta" = c("A" = 0.18, "B" = 0.16), 
                "b_beta" = c("A" = 0.28, "B" = 0.26))
  
  expected <- list("screen.type" = c(rep("A", 100), rep("B", 100)),
                   "arm" = c(rep("all", 100)),
                   "theta" = list("beta" = c("A" = 0.8, "B" = 0.6)),
                   "prior" = list("a_beta" = c("A" = 0.18, "B" = 0.16),
                                  "b_beta" = c("A" = 0.28, "B" = 0.26)))
  
  expect_equal(.processSensitivity(data.assess, data.clinical, theta_0, prior), 
               expected)
  
})

test_that("`processSensitivity()` returns expected results multiple beta screen and arm", {
  
  data.assess <- withr::with_seed(
    1234,
    data.frame("id" = c(1:100, 1:100), 
               "age_assess" = 1:200,
               "disease_detected" = rbinom(200, 1, 0.4),
               "screen_type" = c(rep("A", 100), rep("B", 100)))
  )
  
  data.clinical <- withr::with_seed(
    1234,
    data.frame("id" = 1:100, 
               "age_entry" = 1:100,
               "age_endpoint" = 2:101,
               "endpoint_type" = rep("clinical", 100),
               "arm" = rep("C", 100))
  )
  
  theta_0 <- list("beta" = c("A" = 0.8, "B" = 0.6))
  prior <- list("a_beta" = c("A" = 0.18, "B" = 0.16), 
                "b_beta" = c("A" = 0.28, "B" = 0.26))
  
  expected <- list("screen.type" = c(rep("A", 100), rep("B", 100)),
                   "arm" = c(rep("C", 100)),
                   "theta" = list("beta" = c("A" = 0.8, "B" = 0.6)),
                   "prior" = list("a_beta" = c("A" = 0.18, "B" = 0.16),
                                  "b_beta" = c("A" = 0.28, "B" = 0.26)))
  
  expect_equal(.processSensitivity(data.assess, data.clinical, theta_0, prior), 
               expected)
  
})

test_that("`processSensitivity()` returns expected results multiple beta screen and arm", {
  
  data.assess <- withr::with_seed(
    1234,
    data.frame("id" = c(1:100, 1:100), 
               "age_assess" = 1:200,
               "disease_detected" = rbinom(200, 1, 0.4),
               "screen_type" = c(rep("C", 100), rep("D", 100)))
  )
  
  data.clinical <- withr::with_seed(
    1234,
    data.frame("id" = 1:100, 
               "age_entry" = 1:100,
               "age_endpoint" = 2:101,
               "endpoint_type" = rep("clinical", 100),
               "arm" = c(rep("A", 50), rep("B", 50)))
  )
  
  theta_0 <- list("beta" = c("A" = 0.8, "B" = 0.6))
  prior <- list("a_beta" = c("A" = 0.18, "B" = 0.16), 
                "b_beta" = c("A" = 0.28, "B" = 0.26))
  
  expected <- list("screen.type" = c(rep("A", 50), rep("B", 50), rep("A", 50), rep("B", 50)),
                   "arm" = c(rep("A", 50), rep("B", 50)),
                   "theta" = list("beta" = c("A" = 0.8, "B" = 0.6)),
                   "prior" = list("a_beta" = c("A" = 0.18, "B" = 0.16),
                                  "b_beta" = c("A" = 0.28, "B" = 0.26)))
  
  expect_warning(out <- .processSensitivity(data.assess, data.clinical, theta_0, prior), 
                 "*** screen_type designations ignored ***",
                 fixed = TRUE)
  expect_equal(out, expected)
  
})


test_that("`.testGroupBetaInput` returns expected errors", {
  
  prior <- list("a_beta" = 1, "b_beta" = c("a" = 1, "b" = 2))
  theta <- list("beta" = c("a" = 10, "b" = 20))
  
  expect_error(.testGroupBetaInput(c("a", "b", "a"), theta, prior),
               "inappropriate parameter specifications for beta")
  
  prior <- list("a_beta" = c("a" = 1), "b_beta" = c("a" = 1, "b" = 2))
  theta <- list("beta" = c("a" = 10, "b" = 20))
  
  expect_error(.testGroupBetaInput(c("a", "b", "a"), theta, prior),
               "inappropriate parameter specifications for beta")
  
  prior <- list("a_beta" = c("a" = 1, "b" = 2), "b_beta" = c("a" = 1))
  theta <- list("beta" = c("a" = 10, "b" = 20))
  
  expect_error(.testGroupBetaInput(c("a", "b", "a"), theta, prior),
               "inappropriate parameter specifications for beta")
  
  prior <- list("a_beta" = c("a" = 1, "b" = 2), "b_beta" = 1)
  theta <- list("beta" = c("a" = 10, "b" = 20))
  
  expect_error(.testGroupBetaInput(c("a", "b", "a"), theta, prior),
               "inappropriate parameter specifications for beta")
  
  prior <- list("a_beta" = c("a" = 1, "b" = 2), "b_beta" = c("a" = 1, "b" = 2))
  theta <- list("beta" = c("a" = 10, "b" = 20))
  
  expect_error(.testGroupBetaInput(c("a", "a", "a", "a"), theta, prior),
               "names of `theta_0$beta`, `prior$a_beta`, and/or `prior$b_beta` do not match values of `screen_type`",
               fixed = TRUE)

  prior <- list("a_beta" = c("a" = 1), "b_beta" = c("a" = 12))
  theta <- list("beta" = c("a" = 10))
  
  expect_error(.testGroupBetaInput(c("a", "b", "a", "a"), theta, prior),
               "names of `theta_0$beta`, `prior$a_beta`, and/or `prior$b_beta` do not match values of `screen_type`",
               fixed = TRUE)
  
  prior <- list("a_beta" = c("a" = 1, "b" = 2), "b_beta" = c("a" = 1, "b" = 2))
  theta <- list("beta" = c("a" = 10, "c" = 20))
  
  expect_error(.testGroupBetaInput(c("a", "b", "a", "a"), theta, prior),
               "names of `theta_0$beta`, `prior$a_beta`, and/or `prior$b_beta` do not match values of `screen_type`",
               fixed = TRUE)
  
  prior <- list("a_beta" = c("a" = 1, "c" = 2), "b_beta" = c("a" = 1, "b" = 2))
  theta <- list("beta" = c("a" = 10, "b" = 20))
  
  expect_error(.testGroupBetaInput(c("a", "b", "a", "a"), theta, prior),
               "names of `theta_0$beta`, `prior$a_beta`, and/or `prior$b_beta` do not match values of `screen_type`",
               fixed = TRUE)
  
  prior <- list("a_beta" = c("a" = 1, "b" = 2), "b_beta" = c("c" = 1, "b" = 2))
  theta <- list("beta" = c("a" = 10, "b" = 20))
  
  expect_error(.testGroupBetaInput(c("a", "b", "a", "a"), theta, prior),
               "names of `theta_0$beta`, `prior$a_beta`, and/or `prior$b_beta` do not match values of `screen_type`",
               fixed = TRUE)
  
  prior <- list("a_beta" = c("a" = 1, "b" = 2), "b_beta" = c("c" = 1, "b" = 2))
  theta <- list("beta" = c("a" = 10, "b" = 20))
  
  expect_error(.testGroupBetaInput(c("a", "b", "a", "c"), theta, prior),
               "names of `theta_0$beta`, `prior$a_beta`, and/or `prior$b_beta` do not match values of `screen_type`",
               fixed = TRUE)

  prior <- list("a_beta" = c(1, 2), "b_beta" = c("c" = 1, "b" = 2))
  theta <- list("beta" = c("a" = 10, "b" = 20))
  
  expect_error(.testGroupBetaInput(c("a", "b", "a", "c"), theta, prior),
               "names of `theta_0$beta`, `prior$a_beta`, and/or `prior$b_beta` do not match values of `screen_type`",
               fixed = TRUE)
  
  prior <- list("a_beta" = c("a" = 1, "b" = 2), "b_beta" = c(1, 2))
  theta <- list("beta" = c("a" = 10, "b" = 20))
  
  expect_error(.testGroupBetaInput(c("a", "b", "a", "c"), theta, prior),
               "names of `theta_0$beta`, `prior$a_beta`, and/or `prior$b_beta` do not match values of `screen_type`",
               fixed = TRUE)
  
  prior <- list("a_beta" = c("a" = 1, "b" = 2), "b_beta" = c("c" = 1, "b" = 2))
  theta <- list("beta" = c(10, 20))
  
  expect_error(.testGroupBetaInput(c("a", "b", "a", "c"), theta, prior),
               "names of `theta_0$beta`, `prior$a_beta`, and/or `prior$b_beta` do not match values of `screen_type`",
               fixed = TRUE)
  
})

test_that("`.testGroupBetaInput` returns expected results", {
  
  screen_type <- withr::with_seed(1234, sample(c("a","b","c"), 100, TRUE))
  prior <- list("a_beta" = c("b" = 1, "c" = 2, "a" = 3), 
                "b_beta" = c("a" = 3, "c" = 1, "b" = 2))
  theta <- list("beta" = c("a" = 10, "c" = 30, "b" = 20))
  
  expected <- list("screen.type" = as.factor(screen_type),
                   "theta" = list("beta" = c("a" = 10, "b" = 20, "c" = 30)), 
                   "prior" = list("a_beta" = c("a" = 3, "b" = 1, "c" = 2), 
                                  "b_beta" = c("a" = 3, "b" = 2, "c" = 1)))
  
  expect_equal(.testGroupBetaInput(screen_type, prior = prior, theta = theta), expected)
  
})

test_that("`.testGroupRatePInput` returns expected error", {
  prior <- list("rate_P" = c("b" = 1, "c" = 2, "a" = 3), 
                "shape_P" = c("a" = 3, "c" = 1, "b" = 2))
  theta <- list("rate_P" = c("a" = 10, "c" = 30))
  epsilon <- c("a" = 0.1, "b" = 0.2, "c" = 0.3)
  
  expect_error(.testGroupRatePInput(c("a", "a", "a"), 10, theta, prior, epsilon),
               "inappropriate parameter specifications for pre-clinical Weibull; ")
  
  prior <- list("rate_P" = c("b" = 1, "c" = 2, "a" = 3), 
                "shape_P" = c("a" = 3, "c" = 1))
  theta <- list("rate_P" = c("a" = 10, "c" = 30, "b" = 20))
  epsilon <- c("a" = 0.1, "b" = 0.2, "c" = 0.3)
  
  expect_error(.testGroupRatePInput(c("a", "a", "a"), 10, theta, prior, epsilon),
               "inappropriate parameter specifications for pre-clinical Weibull; ")
  
  prior <- list("rate_P" = c("c" = 2, "a" = 3), 
                "shape_P" = c("a" = 3, "c" = 1, "b" = 2))
  theta <- list("rate_P" = c("a" = 10, "c" = 30, "b" = 20))
  epsilon <- c("a" = 0.1, "b" = 0.2, "c" = 0.3)
  
  expect_error(.testGroupRatePInput(c("a", "a", "a"), 10, theta, prior, epsilon),
               "inappropriate parameter specifications for pre-clinical Weibull; ")
  
  prior <- list("rate_P" = c("b" = 1, "c" = 2, "a" = 3), 
                "shape_P" = c("a" = 3, "c" = 1, "b" = 2))
  theta <- list("rate_P" = c("a" = 10, "c" = 30, "b" = 20))
  epsilon <- c("a" = 0.1, "c" = 0.3)
  
  expect_error(.testGroupRatePInput(c("a", "a", "a"), 10, theta, prior, epsilon),
               "inappropriate parameter specifications for pre-clinical Weibull; ")
  
  prior <- list("rate_P" = c("b" = 1, "c" = 2, "a" = 3), 
                "shape_P" = c("a" = 3, "c" = 1, "b" = 2))
  theta <- list("rate_P" = c("a" = 10, "c" = 30, "b" = 20))
  epsilon <- c("a" = 0.1, "b" = 0.2, "c" = 0.3)
  
  expect_error(.testGroupRatePInput(c("a", "a", "a"), 10, 
                                    theta = theta, prior = prior, epsilon = epsilon),
               "names of `theta_0$rate_P`, `prior$rate_P`, `prior$shape_P`, and/or `epsilon_rate_P` do not match values of `grp.rateP`",
               fixed = TRUE)
  
  expect_error(.testGroupRatePInput(c("a", "b", "a"), 10, 
                                    theta = theta, prior = prior, epsilon = epsilon),
               "names of `theta_0$rate_P`, `prior$rate_P`, `prior$shape_P`, and/or `epsilon_rate_P` do not match values of `grp.rateP`",
               fixed = TRUE)
  
  expect_error(.testGroupRatePInput(c("a", "b", "C"), 10, 
                                    theta = theta, prior = prior, epsilon = epsilon),
               "names of `theta_0$rate_P`, `prior$rate_P`, `prior$shape_P`, and/or `epsilon_rate_P` do not match values of `grp.rateP`",
               fixed = TRUE)
  
  expect_error(.testGroupRatePInput(c("a", "b", "c", "d"), 10, 
                                    theta = theta, prior = prior, epsilon = epsilon),
               "names of `theta_0$rate_P`, `prior$rate_P`, `prior$shape_P`, and/or `epsilon_rate_P` do not match values of `grp.rateP`",
               fixed = TRUE)
})

test_that("`.testGroupRatePInput` returns expected messages/warnings", {
  prior <- list("rate_P" = 1, 
                "shape_P" = 1)
  theta <- list("rate_P" = 1)
  epsilon <- 1
  
  expect_message(.testGroupRatePInput(NULL, 10, theta, prior, epsilon),
                 "no grouping variable for pre-clinical Weibull distribution specified")
  
  expect_message(suppressWarnings(.testGroupRatePInput(c("a","b","a","b"), 10, theta, prior, epsilon)),
                 "`grp.rate_P` variable of `data.clinical` ignored; ")
  
  expect_warning(.testGroupRatePInput(c("a","b","a","b"), 10, theta, prior, epsilon), 
                 "`grp.rate_P` variable of `data.clinical` ignored; ")
  
})

test_that("`.testGroupRatePInput` returns expected messages/warnings", {
  prior <- list("rate_P" = 1, 
                "shape_P" = 1)
  theta <- list("rate_P" = 1)
  epsilon <- 1
  
  expected <- list("grp.rateP" = as.factor(rep("all", 10)),
                   "theta" = list("rate_P" = c("all" = 1)), 
                   "prior" = list("rate_P" = c("all" = 1),  "shape_P" = c("all" = 1)), 
                   "epsilon" = c("all" = 1))
  
  expect_equal(.testGroupRatePInput(NULL, n = 10, theta = theta, 
                                    prior = prior, epsilon = epsilon),
               expected)
  
  expect_equal(.testGroupRatePInput(NULL, n = 10, theta = theta, 
                                    prior = prior),
               expected)

  expect_equal(.testGroupRatePInput(c(rep("all", 10)), n = 10, theta = theta, 
                                    prior = prior),
               expected)

  expect_equal(suppressWarnings(
    .testGroupRatePInput(c(rep("all", 5), rep("ALL", 5)), n = 10, theta = theta, 
                         prior = prior)),
               expected)
  
  prior <- list("rate_P" = c("b" = 1, "a" = 2), 
                "shape_P" = c("b" = 10, "a" = 20))
  theta <- list("rate_P" = c("b" = 11, "a" = 21))
  epsilon <- c("b" = 12, "a" = 22)
  
  expected <- list("grp.rateP" = as.factor(c(rep("a", 5), rep("b", 5))),
                   "theta" = list("rate_P" = c("a" = 21, "b" = 11)), 
                   "prior" = list("rate_P" = c("a" = 2, "b" = 1),  "shape_P" = c("a" = 20, "b" = 10)), 
                   "epsilon" = c("a" = 22, "b" = 12))
  
  expect_equal(.testGroupRatePInput(c(rep("a", 5), rep("b", 5)), n = 10, theta = theta, 
                                    prior = prior, epsilon = epsilon),
               expected)
  
  prior <- list("rate_P" = 1, 
                "shape_P" = 2)
  theta <- list("rate_P" = 3)
  epsilon <- 4
  
  expected <- list("grp.rateP" = as.factor(c(rep("a", 10))),
                   "theta" = list("rate_P" = c("a" = 3)), 
                   "prior" = list("rate_P" = c("a" = 1),  "shape_P" = c("a" = 2)), 
                   "epsilon" = c("a" = 4))
  
  expect_equal(.testGroupRatePInput(c(rep("a", 10)), n = 10, theta = theta, 
                                    prior = prior, epsilon = epsilon),
               expected)
})

test_that("`.testAdaptiveInput` return expected errors", {
  
  expect_error(.testAdaptiveInput(1),
               "`adaptive` must be null or a list")

  expect_error(.testAdaptiveInput(data.frame("a" = 1)),
               "`adaptive` must be null or a list")
  
  expect_error(.testAdaptiveInput(list("warm" = 1, "delta" = 1, "gamma" = 1, "kappa" = 1, "m0" = 1)),
               "`adaptive` does not contain required information")
  
  expect_error(.testAdaptiveInput(list("delta" = 1, "gamma" = 1, "kappa" = 1, "m0" = 1)),
               "`adaptive` does not contain required information")
  
  expect_error(.testAdaptiveInput(list("warmup" = "1", "delta" = 1, "gamma" = 1, "kappa" = 1, "m0" = 1)),
               "all elements of `adaptive` must be numeric")
  
  expect_error(.testAdaptiveInput(list("warmup" = 1, "delta" = "1", "gamma" = 1, "kappa" = 1, "m0" = 1)),
               "all elements of `adaptive` must be numeric")
  
  expect_error(.testAdaptiveInput(list("warmup" = 1, "delta" = 1, "gamma" = "1", "kappa" = 1, "m0" = 1)),
               "all elements of `adaptive` must be numeric")
  
  expect_error(.testAdaptiveInput(list("warmup" = 1, "delta" = 1, "gamma" = 1, "kappa" = "1", "m0" = 1)),
               "all elements of `adaptive` must be numeric")
  
  expect_error(.testAdaptiveInput(list("warmup" = 1, "delta" = 1, "gamma" = 1, "kappa" = 1, "m0" = "1")),
               "all elements of `adaptive` must be numeric")
  
  expect_error(.testAdaptiveInput(list("warmup" = 1:2, "delta" = 1, "gamma" = 1, "kappa" = 1, "m0" = 1)),
               "all elements of `adaptive` must be scalar")
  
  expect_error(.testAdaptiveInput(list("warmup" = 1, "delta" = 1:2, "gamma" = 1, "kappa" = 1, "m0" = 1)),
               "all elements of `adaptive` must be scalar")
  
  expect_error(.testAdaptiveInput(list("warmup" = 1, "delta" = 1, "gamma" = 1:2, "kappa" = 1, "m0" = 1)),
               "all elements of `adaptive` must be scalar")

  expect_error(.testAdaptiveInput(list("warmup" = 1.5, "delta" = 1, "gamma" = 1, "kappa" = 1, "m0" = 1)),
               "`adaptive$warmup` must be a positive integer", fixed = TRUE)
  
  expect_error(.testAdaptiveInput(list("warmup" = -1, "delta" = 1, "gamma" = 1, "kappa" = 1, "m0" = 1)),
               "`adaptive$warmup` must be a positive integer", fixed = TRUE)
  
  expect_error(.testAdaptiveInput(list("warmup" = 1, "delta" = 1.2, "gamma" = 1, "kappa" = 1, "m0" = 1)),
               "`adaptive$delta` must in (0, 1)", fixed = TRUE)
  
  expect_error(.testAdaptiveInput(list("warmup" = 1, "delta" = -1, "gamma" = 1, "kappa" = 1, "m0" = 1)),
               "`adaptive$delta` must in (0, 1)", fixed = TRUE)
  
  expect_error(.testAdaptiveInput(list("warmup" = 1, "delta" = 0.8, "gamma" = 1.1, "kappa" = 1, "m0" = 1)),
               "`adaptive$gamma` must in (0, 1]", fixed = TRUE)
  
  expect_error(.testAdaptiveInput(list("warmup" = 1, "delta" = 0.8, "gamma" = -1.1, "kappa" = 1, "m0" = 1)),
               "`adaptive$gamma` must in (0, 1]", fixed = TRUE)

  expect_error(.testAdaptiveInput(list("warmup" = 1, "delta" = 0.8, "gamma" = 1, "kappa" = 1.1, "m0" = 1)),
               "`adaptive$kappa` must in (0, 1]", fixed = TRUE)
  
  expect_error(.testAdaptiveInput(list("warmup" = 1, "delta" = 0.8, "gamma" = 1, "kappa" = -1.1, "m0" = 1)),
               "`adaptive$kappa` must in (0, 1]", fixed = TRUE)

  expect_error(.testAdaptiveInput(list("warmup" = 1, "delta" = 0.8, "gamma" = 1, "kappa" = 1, "m0" = 1.1)),
               "`adaptive$m0` must be a positive integer", fixed = TRUE)
  
  expect_error(.testAdaptiveInput(list("warmup" = 1, "delta" = 0.8, "gamma" = 1, "kappa" = 1, "m0" = -1)),
               "`adaptive$m0` must be a positive integer", fixed = TRUE)
  
})

test_that(".testDataAssess returns expected errors", {
  
  data_assess <- data.frame("id" = 1:5,
                            "disease_detected" = c("dead", "alive", "dead", "alive", "dead"))
  expect_error(.testDataAssess(data_assess),
               "data.assess$disease_detected must be binary 0/1", fixed = TRUE)
  
  data_assess <- data.frame("id" = 1:5,
                            "disease_detected" = c(2, 1, 2, 1, 2))
  expect_error(.testDataAssess(data_assess),
               "`data.assess$disease_detected` must be binary 0/1", fixed = TRUE)
  
  data_assess <- data.frame("id" = 1:5,
                            "disease_detected" = c(1, 0, 1, 0, 1),
                            "age_assess" = c(40, 50, 25, 15, -10))
  expect_error(.testDataAssess(data_assess),
               "`data.assess$age_assess` cannot be negative", fixed = TRUE)
  
})

test_that(".testDataAssess returns expected results", {
  data_assess <- data.frame("id" = c(1:5, 1:5),
                            "disease_detected" = c(1.1, 0.1, 1.1, 0.1, 1.1, 1.1, 0.1, 1.1, 0.1, 1.1),
                            "age_assess" = c(23, 34, 45, 65, 23, 45, 21, 88, 54, 50))
  
  expected <- data.frame("id" = as.character(rep(1:5, each = 2)),
                         "disease_detected" = c(1L, 1L, 0L, 0L, 1L, 1L, 0L, 0L, 1L, 1L),
                         "age_assess" = c(23, 45, 21, 34, 45, 88, 54, 65, 23, 50))
  rownames(expected) <- NULL

  expect_equal(.testDataAssess(data_assess), expected)
})

test_that(".testDataClinical returns expected errors", {
  
  data_clinical <- data.frame("id" = c(1:5, 1:5))
  expect_error(.testDataClinical(data_clinical),
               "at least 1 participant id has multiple records in `data.clinical`")
  
  data_clinical <- data.frame("id" = c(1:5),
                              "age_endpoint" = c(1, 2, 3, 4, -5))
  expect_error(.testDataClinical(data_clinical),
               "`data.clinical$age_endpoint` cannot be negative", fixed = TRUE)

  data_clinical <- data.frame("id" = c(1:5),
                              "age_endpoint" = c(1, 2, 3, 4, 5),
                              "age_entry" = c(1, 2, 3, 4, -5))
  expect_error(.testDataClinical(data_clinical),
               "`data.clinical$age_entry` cannot be negative", fixed = TRUE)
  
  data_clinical <- data.frame("id" = c(1:5),
                              "age_endpoint" = c(1, 2, 3, 4, 5),
                              "age_entry" = c(1, 2, 3, 4, 5),
                              "endpoint_type" = c("clinical", "preclinical", "censored",
                                                  "censored", "prclinical"))
  expect_error(.testDataClinical(data_clinical),
               "unrecognized endpoint type in `data.clinical$endpoint_type`", fixed = TRUE)

})

test_that(".testDataClinical returns expected results", {
  data_clinical <- data.frame("id" = c(2, 3, 4, 5, 1),
                              "age_endpoint" = c(1, 2, 3, 4, 5),
                              "age_entry" = c(1, 2, 3, 4, 5),
                              "endpoint_type" = c("clinical", "preclinical", "censored",
                                                  "censored", "preclinical"))
  
  expected <- data.frame("id" = as.character(c(1, 2, 3, 4, 5)),
                         "age_endpoint" = c(5, 1, 2, 3, 4),
                         "age_entry" = c(5, 1, 2, 3, 4),
                         "endpoint_type" = c("preclinical", "clinical", "preclinical", "censored",
                                             "censored"))
  rownames(expected) <- NULL
  
  expect_equal(.testDataClinical(data_clinical), expected)
})
