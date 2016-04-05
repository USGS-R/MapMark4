context("Testing the S3 class for the pdf for the material tonnages
in all undiscovered deposit in the permissive tract")

test_that("Object is constructed properly",{

  pmf <- NDepositsPmf( "NegBinomial", list(theMean=5,theStdDev=4), "" )
  pdf1 <- TonnagePdf1(ExampleTonnageData, "mt", pdfType = "normal")
  pdfPT <- TonnagePdfPT(pmf, pdf1)

  # The first column is for the names of the deposits.
  matNames <- colnames(ExampleTonnageData)[-1]
  N <- length(matNames)

  expect_match(pdfPT$units, "mt" )
  for(i in 1:N) {
    expect_match(pdfPT$matNames[i], matNames[i])
  }

  for(i in 1:N) {
    expect_equal(pdfPT$theMean[i], pdfPT$predMean[i], tolerance = 0.05)
    expect_equal(pdfPT$theSd[i], pdfPT$predSd[i], tolerance = 0.05)
  }

  if(N > 1) {
    utCor <- pdfPT$theCor[upper.tri(pdfPT$theCor)]
    utPredCor <- pdfPT$predCor[upper.tri(pdfPT$predCor)]
    for(i in 1:length(utCor)) {
      expect_equal(utCor[i], utPredCor[i], tolerance = 0.05)
    }
  }

})


