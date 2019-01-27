context("test-uthwigl-csUTh-tests")

data("Pan2018")
 # Only solve for sample MK16

suppressMessages(png(filename = "test-csUTh.png")) # capture the plot
output_cs <-
  suppressMessages(csUTh(
    Pan2018,
    sample_name = 'YP003',
    nbitchoice = 10,
    detcorrectionchoice = TRUE,
    keepfiltereddata = FALSE,
    print_summary = FALSE
  ))
print(output_cs$plots)
suppressMessages(dev.off()) # finish capturing the plot

# here are the tests:

test_that("csUTh() returns a list", {
  expect_true(is.list(output_cs))
})

test_that("csUTh() returns a plot", {
  expect_true(file.exists("test-csUTh.png"))
})

# these tests below may fail is the tolerance is set too low:

test_that("csUTh() returns sensible values for the age", {
  
  expect_equal(mean(output_cs$results$`Age (ka)`), 120, tolerance = 10)
})

test_that("csUTh() returns sensible values for (234U/238U)i", {
  
  expect_equal(mean(output_cs$results$`(234U/238U)i`), 1.131, tolerance = 0.1)
})


# delete the plot after testing
file.remove("test-csUTh.png")
