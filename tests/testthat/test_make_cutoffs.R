test_that("make_cutoffs works",{
  expect_type(
    make_n_cutoffs(
      malaria_df1$fever, 
      malaria_df1$density,
      mintot = 50,
      add1 = TRUE),
    "double")

  expect_type(  
  make_cutoffs(
    malaria_df1$fever, 
    malaria_df1$density, 
    add1 = TRUE),
  "double")
  
  df <- data.frame(
    fever = c(0,1,0,1,0,1),
    density = c(0,0,100,150,200,500)
  )
  expect_equal(
    make_cutoffs(df$fever,df$density, add1 = TRUE),
    c(0,1,150,200,500)
  )
  expect_equal(
    make_cutoffs(df$fever,df$density, add1 = FALSE),
    c(0,150,200,500)
  )
})
