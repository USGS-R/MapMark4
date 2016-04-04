context("Testing the S3 class for the meta data")

test_that("Object is constructed properly",{
  meta <- Metadata("Lucky Strike PT", "PT001",
                       tractArea = "5.1 km2",
                       nKnownDeposits = 2,
                       description = "Tract for special deposits")
  expect_match(meta$tractName, "Lucky Strike PT")
  expect_match(meta$tractId, "PT001")
  expect_match(meta$tractArea, "5.1 km2")
  expect_equal(meta$nKnownDeposits, 2)
  expect_match(meta$description, "Tract for special deposits")

  expect_match(getTractId(meta), "PT001")
})
