test_that("sens_spec works",{
  fit <- logitexp(malaria_df1$fever, malaria_df1$density)
  xx <- senspec(fit,  c(1,100,500,1000,2000,4000,8000,16000, 32000,54000,100000))
  expect_type(xx,"double")
})