
test_that("Data is loaded", {
  expect_s3_class(malaria_df1,"data.frame")
  expect_s3_class(malaria_df2,"data.frame")
})