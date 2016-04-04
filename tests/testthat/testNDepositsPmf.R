context("Testing the S3 class for the number of undiscovered deposits")

# Root, D.H., Menzie, W.D., and Scott, W.A., 1992, Computer Monte Carlo
# simulation in quantitative resource estimation: Nonrenewable resources,
# v. 1, no. 2, p. 125-138.
test_that("Pmf A in Tables 3 and 4 of Root, Menzie, and Scott is reproduced",{
  pmf <- getNDepositPmf(NDepositsPmf("Mark3", list(thresholds=c(1,2,4)), ""))
   expect_equal( nrow(pmf), 5 )
  expect_equal( pmf$nDeposits[1], 0 )
  expect_equal( pmf$nDeposits[5], 4 )
  expect_equal( pmf$probs[1], 0.067, tolerance = 0.001, scale = 1 )
  expect_equal( pmf$probs[2], 0.233, tolerance = 0.001, scale = 1 )
  expect_equal( pmf$probs[3], 0.300, tolerance = 0.001, scale = 1 )
  expect_equal( pmf$probs[4], 0.200, tolerance = 0.001, scale = 1 )
  expect_equal( pmf$probs[5], 0.200, tolerance = 0.001, scale = 1 )
})

test_that("Pmf B in Tables 3 and 4 of Root, Menzie, and Scott is reproduced",{
  pmf <- getNDepositPmf(NDepositsPmf( "Mark3", list(thresholds=c(0,2,2)), ""))
  expect_equal( nrow(pmf), 3 )
  expect_equal( pmf$nDeposits[1], 0 )
  expect_equal( pmf$nDeposits[3], 2 )
  expect_equal( pmf$probs[1], 0.2, tolerance = 0.001, scale = 1 )
  expect_equal( pmf$probs[2], 0.2, tolerance = 0.001, scale = 1 )
  expect_equal( pmf$probs[3], 0.6, tolerance = 0.001, scale = 1 )

})

test_that("Mark4 pmf is reproduced",{
  pmf <- getNDepositPmf(NDepositsPmf("Mark4",
                       list(thresholds = c(1, 7, 20),
                            maxNumberOfDeposits = 23), ""))
  expect_equal( nrow(pmf), 24 )
  expect_equal( pmf$nDeposits[1], 0 )
  expect_equal( pmf$nDeposits[24], 23 )
  expect_equal( pmf$probs[1], 0.1, tolerance = 0.001, scale = 1 )
  expect_equal( pmf$probs[24], 0.025, tolerance = 0.001, scale = 1 )
})

test_that("User-specified pmf is reproduced",{
  pmf <- getNDepositPmf(NDepositsPmf("UserSpecified",
                       list(anchorPts = c(1, 4, 10),
                            relProbabilities = c(1, 4, 1)), ""))
  expect_equal( nrow(pmf), 10 )
  expect_equal( pmf$nDeposits[1], 1 )
  expect_equal( pmf$nDeposits[10], 10 )
  expect_equal( pmf$probs[1], 0.0425, tolerance = 0.0001, scale = 1 )
  expect_equal( pmf$probs[10], 0.0425, tolerance = 0.0001, scale = 1 )
})

test_that("Truncated Poisson pmf is reproduced",{
  pmf <- getNDepositPmf(NDepositsPmf("Poisson", list(theMean=5), ""))
  expect_equal( nrow(pmf), 16 )
  expect_equal( pmf$nDeposits[1], 0 )
  expect_equal( pmf$nDeposits[16], 15 )
  expect_equal( pmf$probs[1], 0.00674, tolerance = 0.00001, scale = 1 )
  expect_equal( pmf$probs[16], 0.000157, tolerance = 0.000001, scale = 1 )
})

test_that("Truncated negative binomial pmf is reproduced",{
  pmf <- getNDepositPmf(NDepositsPmf("NegBinomial",
                                     list(theMean = 5, theStdDev = 4), ""))
  expect_equal( nrow(pmf), 32 )
  expect_equal( pmf$nDeposits[1], 0 )
  expect_equal( pmf$nDeposits[32], 31 )
  expect_equal( pmf$probs[1], 7.11e-2, tolerance = 0.01e-2, scale = 1 )
  expect_equal( pmf$probs[32], 4.63e-5, tolerance = 0.01e-5, scale = 1 )
})

test_that("Debug pmf with one non-zero probability is reproduced",{
  pmf <- getNDepositPmf(NDepositsPmf("Debug",
                                     list(nDeposits = 1,
                                         relProbabilities = 1), ""))
  expect_equal( nrow(pmf), 1 )
  expect_equal( pmf$nDeposits[1], 1 )
  expect_equal( pmf$probs[1], 1, tolerance = 0.0001, scale = 1 )
})

test_that("Debug pmf with two non-zero probabilities is reproduced",{
  pmf <- getNDepositPmf(NDepositsPmf("Debug", list(nDeposits = c(2, 3),
                                     relProbabilities = c(0.3, 0.7)), ""))
  expect_equal( nrow(pmf), 2 )
  expect_equal( pmf$nDeposits[1], 2 )
  expect_equal( pmf$nDeposits[2], 3 )
  expect_equal( pmf$probs[1], 0.3, tolerance = 0.0001, scale = 1 )
  expect_equal( pmf$probs[2], 0.7, tolerance = 0.0001, scale = 1 )
})
