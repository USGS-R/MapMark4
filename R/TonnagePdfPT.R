#' @title Plot the pdf for the material tonnages in the permissive tract
#'
#' @description Plot the probability density function (pdf)
#' for the material tonnages in the
#' permissive tract.
#'
#' @param object
#' An object of class "TonnagePdf1".
#'
#' @param whichMatTonnage
#' Character vector comprising the names of the materials that will be included
#' in the plot.
#'
#' @param isUsgsStyle
#' Make the plot format similar to, but not identical to, the
#' U.S. Geological Survey style
#'
#' @details
#' An internal call to function reshape2::melt generates the
#' message "No id variables; using all as measure variables", which should
#' be ignored.
#'
#' @examples
#' pmf <- NDepositsPmf( "NegBinomial", list(theMean=5,theStdDev=4), "" )
#' pdf1 <- TonnagePdf1(ExampleTonnageData, "mt", pdfType = "normal")
#' pdfPT <- TonnagePdfPT(pmf, pdf1)
#' plot(pdfPT)
#'
#' @export
#'
plot.TonnagePdfPT <- function(object,
                             whichMatTonnage = object$matNames,
                             isUsgsStyle = TRUE) {

  Need to account for zeros

  for(i in 1:length(whichMatTonnage)) {
    if(!any(whichMatTonnage[i] == object$matNames)) {
      stop( sprintf( "Function plot.TonnagePdf1\n" ),
            sprintf( "Argument whichMatTonnage must match one of the.\n" ),
            sprintf( "materials. The mismatched name is %s\n",
                     whichMatTonnage[i]),
            call. = FALSE )
    }
  }

  df.rs <- reshape2::melt(object$rs[, whichMatTonnage, drop = FALSE])
  # After the melt operation, the first column is row number. It is
  # useless here, so it is removed.
  df.rs <- df.rs[, -1, drop = FALSE]
  colnames(df.rs) <- c("Material", "Tonnage")

  xLabel <- paste("Tonnage (", object$units, ")", sep = "")

  p <- ggplot2::ggplot(df.rs) +
    ggplot2::geom_density(ggplot2::aes(Tonnage, colour = Material)) +
    ggplot2::scale_x_continuous(name = xLabel, trans = "log10")

  if(isUsgsStyle)
    p <- p + ggplot2::xlab(paste("Tonnage, in ", object$units, sep = "")) +
    ggplot2::theme_bw()

  plot(p)

}

#' @title Summarize the pdf for the material tonnages in the permissive tract
#'
#' @description Summarize the probability density function (pdf) for the
#' material tonnages in all undiscovered deposits within the permissive
#' tract.
#'
#' @param object
#' An object of class "TonnagePdfPT"
#'
#' @param nDigits
#' Number of signficant digits.
#'
#' @details
#' The summary statistics for the pdf are calculated from
#' random samples of the pdf.
#'
#' @examples
#' pmf <- NDepositsPmf( "NegBinomial", list(theMean=5,theStdDev=4), "" )
#' pdf1 <- TonnagePdf1(ExampleTonnageData, "mt", pdfType = "normal")
#' pdfPT <- TonnagePdfPT(pmf, pdf1)
#' summary(pdfPT)
#'
#' @export
#'
summary.TonnagePdfPT <- function(object, nDigits = 2) {

  cat( sprintf( "Units for material tonnage: %s\n", object$units ))

  cat( sprintf( "\n\n"))
  cat( sprintf( "Mean vector for the pdf:\n" ))
  tmp <- cbind(object$theMean, object$predMean)
  colnames(tmp) <- c("Actual", "Predicted")
  print(signif(tmp, digits = nDigits))

  cat( sprintf( "\n\n"))
  cat( sprintf( "Standard deviation vector for the pdf:\n" ))
  tmp <- cbind(object$theSd, object$predSd)
  colnames(tmp) <- c("Actual", "Predicted")
  print(signif(tmp, digits = nDigits))


  cat(sprintf( "\n\n"))
  cat(sprintf("Composite correlation matrix for the pdf\n" ))
  cat(sprintf("The upper triangle is the upper triangle of the actual\n"))
  cat(sprintf("correlation matrix. The lower triangle is the lower\n" ))
  cat(sprintf("triangle of the predicted correlation matrix.\n"))
  tmp <- object$theCor
  tmp[lower.tri(tmp)] <- object$predCor[lower.tri(object$predCor)]
  diag(tmp) <- NA
  print(signif(tmp, digits = nDigits))

}

#' @title Construct the pdf for the material tonnages in the
#' permissive tract
#'
#' @description Construct the probability density function (pdf) for the
#' material tonnages in the permissive
#' tract. The pdf is not explicitly specified; instead it is implicitly
#' specified with the random samples.
#'
#' @param oPmf
#' An object of class "NDepositsPmf"
#'
#' @param oPdf
#' An object of class "TonnagePdf1"
#'
#' @param nRandomSamples
#' Number of random samples of the material tonnages in the
#' permissive tract. The actual number may differ slightly because of
#' rounding.
#'
#' @details
#' The random samples are computed using the algorithm in Ellefsen and others
#' (2016). The predicted mean vector and the predicted covariance matrix for
#' the pdf are computed using the analytic formulas in Ellefsen and others
#' (2016). The predicted correlation matrix and the predicted standard
#' deviation vector are calculated from the predicted covariance matrix.
#'
#' @return If the input arguments have an error, the R-value NULL is returned.
#' Otherwise, a list with the following components is returned.
#' @return \item{totalTonnages}{Matrix comprising random samples of
#' the total material tonnages in the permissive tract and the associated
#' number of undiscovered deposits.}
#' @return \item{units}{Input argument units}
#' @return \item{matNames}{Names of the materials.}
#' @return \item{theMean}{Mean vector for the total material tonnages in
#' the permissive tract.}
#' @return \item{theCov}{Covariance matrix for the total material tonnages in
#' the permissive tract.}
#' @return \item{theCor}{Correlation matrix for the total material tonnages in
#' the permissive tract.}
#' @return \item{theSd}{Standard deviation vector for the total material
#' tonnages in the permissive tract.}
#' @return \item{predMean}{Predicted mean vector for the total material
#' tonnages in the permissive tract.}
#' @return \item{predCov}{Predicted covariance matrix for the total material
#' tonnages in the permissive tract.}
#' @return \item{predCor}{Predicted correlation matrix for the total material
#' tonnages in the permissive tract.}
#' @return \item{predSd}{Predicted standard deviation vector for the total
#' material tonnages in the permissive tract.}
#' @return \item{call}{Function call}
#'
#' @references
#' Ellefsen, K.J., Phillips, J.D., Mihalasky, M.J., and Robinson, G.R., Jr.,
#' 2016, Probability calculations for three-part mineral assessments,
#' U.S. Geological Survey Report XXXX
#'
#' @examples
#' pmf <- NDepositsPmf( "NegBinomial", list(theMean=5,theStdDev=4), "" )
#' pdf1 <- TonnagePdf1(ExampleTonnageData, "mt", pdfType = "normal")
#' pdfPT <- TonnagePdfPT(pmf, pdf1)
#'
#' @export
#'
TonnagePdfPT <- function(oPmf, oPdf, nRandomSamples = 20000) {

  if(nRandomSamples < 5000) {
    stop( sprintf( "Function TonnagePdfPT\n" ),
          sprintf( "The number of random samples must be\n" ),
          sprintf( "greater than or equal to 5000.\n" ),
          sprintf( "But the specified number is %d\n", nRandomSamples ),
          call. = FALSE )
  }

  # Number of random samples (in pdfPT) for each non-zero probability in the pmf
  nRs <- round( oPmf$probs * nRandomSamples )

  # Number of random samples (from pdf1) needed for each non-zero probability
  # in the pmf
  nRs1 <- nRs * oPmf$nDeposits

  # All of the random samples that are needed for the computations
  # It is done this way because repeated calls to getRandomSamples can
  # require a lot of time.
  rs1 <- getRandomSamples(oPdf, sum(nRs1))

  iStart <- 1
  internalGetRandomSamples <- function(n) {
    i1 <- iStart
    i2 <- i1 + n - 1
    iStart <<- i2 + 1  # write to variable iStart in the next higher environment
    return(rs1[i1:i2, ])
  }

  totalTonnages <- NULL
  for( i in seq_along( oPmf$nDeposits ) ) {

    # random samples of the tonnage in the permissive tract, for
    # current number of deposits in the pmf
    rs <- matrix(0, nrow = nRs[i],
                 ncol = length(oPdf$matNames),
                 dimnames = list( NULL, oPdf$matNames))

    # If nDeposits[i] == 0, then the tonnage in the permissive tract
    # must be zero. To handle this situation, matrix rs is initialized to zero.

    if( oPmf$nDeposits[i] > 0 ) {
      for( j in 1:oPmf$nDeposits[i] ) {
        rs <- rs + internalGetRandomSamples(nRs[i])
      }
    }

    tmp <- cbind( nDeposits = rep.int( oPmf$nDeposits[i], nRs[i] ), rs )
    totalTonnages <- rbind(totalTonnages, tmp)
  }

  theCov <- cov(totalTonnages[, -1])
  theCor <- cov2cor(theCov)
  theSd <- sqrt(diag(theCov))

  predCov <- oPmf$theMean * oPdf$theCov +
    oPmf$theVar * oPdf$theMean %o% oPdf$theMean
  predCor <- cov2cor(predCov)
  predSd <- sqrt(diag(predCov))

  rval <- list( rs = totalTonnages,
                units = getUnits(oPdf),
                matNames = colnames(totalTonnages)[-1],
                theMean = colMeans(totalTonnages[, -1]),
                theCov = theCov,
                theCor = theCor,
                theSd = theSd,
                predMean = oPmf$theMean * oPdf$theMean,
                predCov = predCov,
                predCor = predCor,
                predSd = predSd,
                call=sys.call() )

  class(rval) <- "TonnagePdfPT"
  return(rval)
}