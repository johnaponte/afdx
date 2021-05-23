
test_that("AFmodel works", {
  fit <- logitexp(malaria_df1$fever, malaria_df1$density)
  expect_equal(class(fit),"afmodel")
})