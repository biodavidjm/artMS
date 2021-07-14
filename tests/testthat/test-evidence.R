context("Check Evidence file")

test_that("checkIfFile",{
  # Data.frame to test
  L3 <- LETTERS[1:3]
  fac <- sample(L3, 10, replace = TRUE)
  d <- data.frame(x = 1, y = 1:10, fac = fac)
  
  expect_true(is.data.frame(.artms_checkIfFile(d)))
  expect_error(.artms_checkIfFile(L3))
})

test_that("Check artmsQuantification() output", {
  # Get data ready
  artms_data_ph_config$files$evidence <- artms_data_ph_evidence
  artms_data_ph_config$files$keys <- artms_data_ph_keys
  artms_data_ph_config$files$contrasts <- artms_data_ph_contrast
  artms_data_ph_config$output_extras <- 0
  msresults <- artmsQuantification(yaml_config_file = artms_data_ph_config, 
                                   data_object = TRUE,  
                                   display_msstats = FALSE, 
                                   verbose = TRUE, 
                                   printPDF = FALSE, 
                                   printTables = FALSE)
  # Test
  expect_equal( length(msresults), 5 )
  expect_equal(dim(msresults$ComparisonResult)[1], 2)
  expect_equal(dim(msresults$ComparisonResult)[2], 11)
  expect_equal(dim(msresults$ModelQC)[1], 5)
  expect_equal(dim(msresults$ModelQC)[2], 13)
  expect_equal(length(msresults$FittedModel[[1]]), 13)
})


