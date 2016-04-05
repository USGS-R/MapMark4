context("Testing the S3 class for the pdf for the material tonnages
in a single, undiscovered deposit")

test_that("Object is constructed properly",{
  pdf1 <- TonnagePdf1(ExampleTonnageData, "mt")

  expect_match( getUnits(pdf1), "mt" )
  expect_match( pdf1$pdfType, "empirical" )
  expect_identical( pdf1$isTruncated, TRUE )

  matNames <- colnames(ExampleTonnageData)[-1]
  for(i in 1:length(matNames)){
    expect_match( pdf1$matNames[i], matNames[i] )
  }
})


test_that("Statistics for pdf and observed data match, case 1",{

  # Case 1: One material, empirical distribution, truncation
  pdf1 <- TonnagePdf1(ExampleTonnageData, "mt")

  # number of random samples
  N <- 1000
  rs <- getRandomSamples1(pdf1, 1000, seed = 7, log_rs = TRUE)

  # min
  rsMin <- apply(rs, 2, min)
  pdfMin <- apply(pdf1$logTonnages, 2, min)
  for(i in 1:length(rsMin)){
    expect_true(pdfMin[i] < rsMin[i])
  }

  # max
  rsMax <- apply(rs, 2, max)
  pdfMax <- apply(pdf1$logTonnages, 2, max)
  for(i in 1:length(rsMax)){
    expect_true(pdfMax[i] > rsMax[i])
  }

  # mean
  # The bounds for the mean of the random samples are calculated with the
  # central limit theorem: Five standard deviations is chosen as the bound.
  rsMean <- apply(rs, 2, mean)
  pdfMean <- apply(pdf1$logTonnages, 2, mean)
  pdfSd <- apply(pdf1$logTonnages, 2, sd)
  N <- nrow(pdf1$obsTonnages)
  lb <- pdfMean - 5 * pdfSd  / sqrt(N)
  ub <- pdfMean + 5 * pdfSd  / sqrt(N)
  for(i in 1:length(rsMean)){
    expect_true(lb[i] < rsMean[i])
    expect_true(ub[i] > rsMean[i])
  }

})

test_that("Statistics for pdf and observed data match, case 2",{

  # Case 2: One material, normal distribution, truncation
  pdf1 <- TonnagePdf1(ExampleTonnageData, "mt", pdfType = "normal")

  # number of random samples
  N <- 1000
  rs <- getRandomSamples1(pdf1, 1000, seed = 7, log_rs = TRUE)

  # min
  rsMin <- apply(rs, 2, min)
  pdfMin <- apply(pdf1$logTonnages, 2, min)
  for(i in 1:length(rsMin)){
    expect_true(pdfMin[i] < rsMin[i])
  }

  # max
  rsMax <- apply(rs, 2, max)
  pdfMax <- apply(pdf1$logTonnages, 2, max)
  for(i in 1:length(rsMax)){
    expect_true(pdfMax[i] > rsMax[i])
  }

  # mean
  # The bounds for the mean of the random samples are calculated with the
  # central limit theorem: Five standard deviations is chosen as the bound.
  rsMean <- apply(rs, 2, mean)
  pdfMean <- apply(pdf1$logTonnages, 2, mean)
  pdfSd <- apply(pdf1$logTonnages, 2, sd)
  lb <- pdfMean - 5 * pdfSd / sqrt(N)
  ub <- pdfMean + 5 * pdfSd / sqrt(N)
  for(i in 1:length(rsMean)){
    expect_true(lb[i] < rsMean[i])
    expect_true(ub[i] > rsMean[i])
  }

})

test_that("Statistics for pdf and observed data match, case 3",{

  # Case 3: One material, empirical distribution, no truncation
  pdf1 <- TonnagePdf1(ExampleTonnageData, "mt", isTruncated = FALSE)

  # number of random samples
  N <- 1000
  rs <- getRandomSamples1(pdf1, N, seed = 7, log_rs = TRUE)

  # mean
  # The bounds for the mean of the random samples are calculated with the
  # central limit theorem.
  rsMean <- apply(rs, 2, mean)
  pdfMean <- apply(pdf1$logTonnages, 2, mean)
  pdfSd <- apply(pdf1$logTonnages, 2, sd)
  lb <- pdfMean - 5 * pdfSd / sqrt(N)
  ub <- pdfMean + 5 * pdfSd / sqrt(N)
  for(i in 1:length(rsMean)){
    expect_true(lb[i] < rsMean[i])
    expect_true(ub[i] > rsMean[i])
  }

})

test_that("Statistics for pdf and observed data match, case 4",{

  # Case 4: One material, normal distribution, no truncation
  pdf1 <- TonnagePdf1(ExampleTonnageData, "mt", pdfType = "normal",
                      isTruncated = FALSE)

  # number of random samples
  N <- 1000
  rs <- getRandomSamples1(pdf1, N, seed = 7, log_rs = TRUE)

  # Use the hypothesis test for the multivariate mean. The confidence level
  # is 0.99
  # See p. 210-216 of
  # Johnson, R.A., and Wichern, D.W., 2007, Applied multivariate statistical
  # analysis: Upper Saddle River, New Jersey, Pearson Prentice Hall,
  # 773 p.
  rsMean <- apply(rs, 2, mean)
  pdfMean <- apply(pdf1$logTonnages, 2, mean)
  rsCov <- cov(rs)

  # eqn 5-4, p. 211
  T2 <- N * t(rsMean - pdfMean) %*% solve(rsCov) %*% (rsMean - pdfMean)

  # eqn 505, p. 212
  p <- ncol(rs)
  ub <- (N - 1) * p / (N - p) *qf(0.99, p, N)

  expect_true(T2 < ub)

  # mean
  # The bounds for the mean of the random samples are calculated with the
  # central limit theorem: Three standard deviations is chosen as the bound.
  rsMean <- apply(rs, 2, mean)
  pdfMean <- apply(pdf1$logTonnages, 2, mean)
  pdfSd <- apply(pdf1$logTonnages, 2, sd)
  lb <- pdfMean - 3 * pdfSd / sqrt(N)
  ub <- pdfMean + 3 * pdfSd / sqrt(N)
  for(i in 1:length(rsMean)){
    expect_true(lb[i] < rsMean[i])
    expect_true(ub[i] > rsMean[i])
  }

  # variance
  # The bounds for the variance of the random samples are calculated with the
  # chi-squared distribution. The 0.99 probability interval is chosen as the
  # bounds. See
  # DeGroot, M.H., and Schervish, M.J., 2002, Probability and statistics
  # (3rd ed.): New York, Addison-Wesley, p. 397- 403
  rsVar <- apply(rs, 2, var)
  pdfVar <- apply(pdf1$logTonnages, 2, var)
  testStat <- unname( N * rsVar / pdfVar )
  lb <- qchisq(0.005, N-1)
  ub <- qchisq(0.995, N-1)
  for(i in 1:length(testStat)){
    expect_true(lb < testStat[i])
    expect_true(ub > testStat[i])
  }

})

