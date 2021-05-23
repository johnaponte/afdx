test_that("get_latent_model workds",{
  xx <- get_latent_model()
  expect_type(xx, "character")
})