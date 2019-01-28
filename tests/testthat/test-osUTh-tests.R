context("test-uthwigl-osUTh-tests")

data("Hobbit_MH2T_for_iDAD")

output_os <- suppressMessages(osUTh(Hobbit_MH2T_for_iDAD,
                   nbit = 10, # to make the tests run fast
                   fsum_target = 0.01,
                   U48_0_min = 1.265,
                   U48_0_max = 1.270,
                   l = 5.35,
                   U_0 = 25,
                   K_min = 1e-13,
                   K_max = 1e-11,
                   T_min = 1e3,
                   T_max = 20e3,
                   print_summary = FALSE,
                   with_plots = TRUE,
                   save_plots = TRUE))

# here are the tests:

test_that("osUTh() returns a list", {
  expect_true(is.list(output_os))
})

test_that("osUTh() returns a plot", {
  expect_true(file.exists("output_os.png"))
})

# these tests below may fail is the tolerance is set too low:

test_that("osUTh() returns sensible values for the age", {
  
  expect_equal(output_os$results$`Age (ka)`, 7.1, tolerance = 1)
})

test_that("osUTh() returns sensible values for T_final", {
  
  expect_equal(output_os$T_final, 7100, tolerance = 1000)
})


# delete the plot after testing
file.remove("output_os.png")
