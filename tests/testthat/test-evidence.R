context("Check Evidence file")

test_that("checkIfFile",{
  # Data.frame to test
  L3 <- LETTERS[1:3]
  fac <- sample(L3, 10, replace = TRUE)
  d <- data.frame(x = 1, y = 1:10, fac = fac)
  
  expect_true(data.table::is.data.table(checkIfFile(d)))
  expect_error(checkIfFile(L3))
})
