context("Average intensity and retention time")

test_that("artmsAvgIntensityRT labels each average with its own column", {
  ev <- data.frame(
    Proteins = c("P1", "P1", "P2"),
    Modified.sequence = c("_AAA_", "_AAA_", "_BBB_"),
    Charge = c(2L, 2L, 3L),
    Intensity = c(1000, 3000, 500),
    Retention.time = c(10, 12, 30),
    Calibrated.retention.time = c(11, 13, 31),
    Raw.file = "r1",
    stringsAsFactors = FALSE
  )
  res <- artmsAvgIntensityRT(ev, verbose = FALSE)

  expect_equal(names(res)[4:6],
               c("Avg_Intensity", "Avg_Retention.time",
                 "Avg_Calibrated.retention.time"))
  p1 <- res[res$Proteins == "P1", ]
  expect_equal(p1$Avg_Intensity, 2000)
  expect_equal(p1$Avg_Retention.time, 11)
  expect_equal(p1$Avg_Calibrated.retention.time, 12)
})
