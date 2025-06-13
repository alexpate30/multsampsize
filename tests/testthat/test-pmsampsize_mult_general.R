test_that("a very basic test", {
  my_sampsize_calc <- pmsampsize_mult_general(type = "m",
                                              parameters = 17,
                                              shrinkage = 0.9,
                                              cstatistic = c(0.85, 0.92, 0.99, 0.95, 0.75, 0.95, 0.87, 0.87, 0.71, 0.82),
                                              ### New parameters
                                              K = 5,
                                              mult_n_events = c(2557, 186, 176, 467, 120),
                                              mult_nagrsquared_overall = 0.15)

  testthat::expect_equal(my_sampsize_calc$parameters, 17)
  testthat::expect_equal(my_sampsize_calc$sample_size, 13136)
  testthat::expect_equal(my_sampsize_calc$type, "multinomial")
  testthat::expect_equal(class(my_sampsize_calc), "list")
})

