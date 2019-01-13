context("test-uthwigl-tests")

data("iolite_export")
 # Only solve for sample MK16

suppressMessages(png(filename = "test-csUTh.png")) # capture the plot
output <-
  suppressMessages(csUTh(
    iolite_export[grepl('MK16', iolite_export$X),],
    nbitchoice = 10,
    detcorrectionchoice = TRUE,
    keepfiltereddata = FALSE,
    print_summary = FALSE
  ))

suppressMessages(dev.off()) # finish capturing the plot

# here are the tests:

test_that("csUTh() returns a data frame", {
  expect_true(is.data.frame(output))
})

test_that("csUTh() returns a plot", {
  expect_true(file.exists("test-csUTh.png"))
})

# these tests below may fail is the tolerance is set too low:

test_that("csUTh() returns sensible values for the age", {
  
  expect_equal(mean(output$`Age (ka)`), 120, tolerance = 10)
})

test_that("csUTh() returns sensible values for (234U/238U)i", {
  
  expect_equal(mean(output$`(234U/238U)i`), 1.131, tolerance = 0.1)
})


# delete the plot after testing
file.remove("test-csUTh.png")
