#' @title Get random samples from the pdf for the
#' material tonnages in
#' all undiscovered deposits within the permissive tract
#'
#' @description Get random samples from the probability density function (pdf)
#' for the material tonnages in all undiscovered deposits within
#' the permissive tract.
#'
#' @param object
#' An object of class "TonnagePdfPT"
#'
#' @return
#' Data frame comprising random samples of
#' the material tonnages in all undiscovered deposits within
#' the permissive tract and the associated
#' number of undiscovered deposits.
#'
#' @examples
#' pmf <- NDepositsPmf( "NegBinomial", list(theMean=5,theStdDev=4), "" )
#' pdf1 <- TonnagePdf1(ExampleTonnageData, "mt", pdfType = "normal")
#' pdfPT <- TonnagePdfPT(pmf, pdf1)
#' getRandomSamples(pdfPT)
#'
#' @export
#'
getRandomSamples.TonnagePdfPT <- function(object) {
  return(object$rs)
}

#' @title Plot the univariate, marginal pdfs for the material tonnages in
#' all undiscovered deposits within the permissive tract
#'
#' @description Plot the unvariate, marginal probability density
#' functions (pdfs)
#' for the material tonnages in all undiscovered deposits within
#' the permissive tract.
#'
#' @param object
#' An object of class "TonnagePdf1".
#'
#' @param whichMatTonnage
#' Character vector comprising the names of the materials that will be included
#' in the plot.
#'
#' @param isUsgsStyle
#' Make the plot format similar to the U.S. Geological Survey style
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

  for(i in 1:length(whichMatTonnage)) {
    if(!any(whichMatTonnage[i] == object$matNames)) {
      stop( sprintf( "Function plot.TonnagePdfPT\n" ),
            sprintf( "Argument whichMatTonnage must match one of the.\n" ),
            sprintf( "materials. The mismatched name is %s\n",
                     whichMatTonnage[i]),
            call. = FALSE )
    }
  }

  rsTonnage <- object$rs[, whichMatTonnage, drop = FALSE]
  nPdfs <- ncol( rsTonnage )
  nSamples <- nrow( rsTonnage )
  countZeros <- unname(colSums(rsTonnage == 0))

  if(nPdfs > 1){
    for( j1 in 1:(nPdfs-1) ) {
      for( j2 in (j1+1):nPdfs ){
        if(!identical(countZeros[j1], countZeros[j2]) ) {
          stop( sprintf( "Function plot.TonnagePdf1\n" ),
                sprintf( "Number of zeros must be equal for all materials\n" ),
                sprintf( "Material: %s  Number of zeros: %d\n",
                         whichMatTonnage[j1], countZeros[j1] ),
                sprintf( "Material: %s  Number of zeros: %d\n",
                         whichMatTonnage[j2], countZeros[j2] ),
                call. = FALSE )
        }
      }
    }
  }

  rsTonnage[rsTonnage == 0] <- NA
  rsLogTonnage <- log10(rsTonnage)

  probZero <- countZeros[1] / nSamples

  theRange <- range(rsLogTonnage, na.rm = TRUE)
  margin <- 0.05 * diff(theRange)
  theRange[1] <- theRange[1] - margin
  theRange[2] <- theRange[2] + margin

  df <- NULL
  for( j in 1:nPdfs ) {

    tmp1 <- density( rsLogTonnage[, j],
                    from = theRange[1], to = theRange[2], na.rm = TRUE)

    tmp2 <- data.frame( Material = rep.int(whichMatTonnage[j], length(tmp1$y)),
                        Tonnage = 10^tmp1$x,
                        Density = tmp1$y * ( 1 - probZero ))
    df <- rbind(df, tmp2)
  }

  if(isUsgsStyle) {
    xLabel <- paste("Tonnage, in ", object$units, sep = "")
  } else {
    xLabel <- paste("Tonnage (", object$units, ")", sep = "")
  }

  caption <- paste( "Prob of zero tonnage = ",
                     round(probZero, digits=3), sep="" )

  p <- ggplot2::ggplot(df) +
    ggplot2::geom_line(ggplot2::aes(x = Tonnage, y = Density, colour = Material)) +
    ggplot2::scale_x_continuous(name = xLabel, trans = "log10") +
    ggplot2::geom_text(ggplot2::aes(x, y, label = caption),
                       data = data.frame(x = min(df$Tonnage), y = 0),
                       hjust = 0, vjust = 1)

  if(isUsgsStyle)
    p <- p + ggplot2::theme_bw()

  plot(p)

}

#' @title Plot the marginal distributions, as a matrix, for the
#' material tonnages
#' in all undiscovered deposits within the permissive tract
#'
#' @description Plot the univariate and bivariate marginal distributions,
#' as a matrix, for the material tonnages
#' in all undiscovered deposit within the permssive tract.
#'
#' @param object
#' An object of class "TonnagePdfPT".
#'
#' @param nPlotSamples
#' Number of samples that represent the pdf in the plots.
#'
#' @param isUsgsStyle
#' Make the plot format similar to the U.S. Geological Survey style
#'
#' @details
#' Recall that the pdf is implicitly specified by random samples. Consequently,
#' the plots within the matrix are generated from \code{nPlotSamples} random
#' samples. A suitable value for parameter \code{nPlotSamples} is
#' approximately 1000.
#'
#' Along the diagonal of the matrix are histograms showing the marginal,
#' univariate distributions. In the upper triangle of the matrix are
#' scatter plots showing the marginal, bivariate distributions.
#' In the lower triangle of the matrix are contour plots
#' showing the marginal, bivariate distributions.
#'
#' Above the matrix of plots is text specifying the probability that
#' the tonnages could be zero.
#'
#' @examples
#' pmf <- NDepositsPmf( "NegBinomial", list(theMean=5,theStdDev=4), "" )
#' pdf1 <- TonnagePdf1(ExampleTonnageData, "mt", pdfType = "normal")
#' pdfPT <- TonnagePdfPT(pmf, pdf1)
#' plotMatrix(pdfPT)
#'
#' @export
#'
plotMatrix.TonnagePdfPT <- function(object, nPlotSamples = 1000,
                                    isUsgsStyle = TRUE) {

  areZero <- object$rs[, "nDeposits"] == 0
  probZero <- sum(areZero) / nrow(object$rs)

  # the first column contains the number of deposits
  tmp <- object$rs[!areZero, -1, drop = FALSE]
  d <- tmp[1:nPlotSamples, , drop = FALSE]
  N <- ncol(d)

  if(isUsgsStyle) {
    tLabel <- paste("tonnage, in ", object$units, sep = "")
  } else {
    tLabel <- paste("tonnage (", object$units, ")", sep = "")
  }

  caption <- paste( "Probability of zero tonnage = ",
                    round(probZero, digits=3), sep="" )

  grid::grid.newpage()
  grid::pushViewport(
    grid::viewport(
      layout =
        grid::grid.layout(N+1, N,
                          heights =
                            grid::unit(rep.int(1,N+1),
                                       c("lines", rep.int("null",N))))))

  grid::grid.text(caption, gp = grid::gpar(fontsize = 15),
                  vp = grid::viewport(layout.pos.row = 1,
                                      layout.pos.col = NULL))

  for (i in 1:N) {
    for(j in 1:N) {
      if(i == j) {
        # diagonal
        matName <- colnames(d)[i]
        p <- ggplot2::ggplot() +
          ggplot2::geom_histogram(ggplot2::aes_string(x = matName),
                                  data = d, bins = 15,
                                  alpha = 0.3, colour = "white",
                                  fill = "blue") +
          ggplot2::scale_x_continuous(name = paste( matName, tLabel, sep = " "),
                                      trans = "log10") +
          ggplot2::scale_y_continuous(name = "", breaks = NULL)
      } else if (j < i) {
        # lower triangle
        xMatName <- colnames(d)[j]
        yMatName <- colnames(d)[i]
        p <- ggplot2::ggplot() +
          ggplot2::geom_density_2d(
            ggplot2::aes_string(x = xMatName, y = yMatName), data = d) +
          ggplot2::scale_x_continuous(
            name = paste( xMatName, tLabel, sep = " "), trans = "log10") +
          ggplot2::scale_y_continuous(
            name = paste( yMatName, tLabel, sep = " "), trans = "log10")

      } else {
        # upper triangle
        xMatName <- colnames(d)[j]
        yMatName <- colnames(d)[i]
        p <- ggplot2::ggplot() +
          ggplot2::geom_point(ggplot2::aes_string(x = xMatName, y = yMatName),
                              data = d, alpha = 0.2, colour = "blue",
                              pch = 16) +
          ggplot2::scale_x_continuous(
            name = paste( xMatName, tLabel, sep = " "), trans = "log10") +
          ggplot2::scale_y_continuous(
            name = paste( yMatName, tLabel, sep = " "), trans = "log10")

      }

      if(isUsgsStyle) {
        p <- p + ggplot2::theme_bw()
      }

      if(N > 1){
        index <- (i - 1) * N + j

        if(isUsgsStyle) {
          figLabel <- LETTERS[index]
          p <- p + ggplot2::ggtitle(figLabel) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0,
                                                              face = "italic"))
        } else {
          figLabel <- paste( "(", letters[index], ")", sep = "")
          p <- p + ggplot2::ggtitle(figLabel) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0))
        }
      }

      print(p, vp=grid::viewport(layout.pos.row=i+1, layout.pos.col=j))
    }
  }
}


#' @title Summarize the univariate, marginal pdfs for the material
#' tonnages in all undiscovered deposits within the permissive tract
#'
#' @description Summarize the univariate, marginal probability density
#' functions (pdfs) for the
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
#' The summary statistics include the
#' 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, and 0.95 quantiles
#' the arithmetic mean, the probability of zero tonnage,
#' and the probability of exceeding the arithmetic mean.
#'
#' @examples
#' pmf <- NDepositsPmf( "NegBinomial", list(theMean=5,theStdDev=4), "" )
#' pdf1 <- TonnagePdf1(ExampleTonnageData, "mt", pdfType = "normal")
#' pdfPT <- TonnagePdfPT(pmf, pdf1)
#' summary(pdfPT)
#'
#' @export
#'
summary.TonnagePdfPT <- function(object, nDigits = 3) {

  CalcTonnageStats <- function( rs )
  {
    probs <- c(0.05,0.10,0.25,0.50,0.75,0.90,0.95)
    labels <- c( paste( "Q_", probs, sep = ""), "Mean", "P(0)", "P(>Mean)" )

    if ( missing( rs ) )
      return( labels )

    quantiles <- quantile( rs, probs, na.rm = TRUE, names = FALSE )
    theMean <- mean( rs )
    p0 <- sum( rs == 0 ) / length( rs )
    p.exceed.mean <- sum( rs > theMean ) / length( rs )

    statistics <- c( quantiles, theMean, p0, p.exceed.mean )
    names( statistics ) <- labels

    return( statistics )

  }

  cat(sprintf("Summary of the pdf for the material tonnages in\n"))
  cat(sprintf("all undiscovered deposits within the permissive tract.\n"))
  cat(sprintf("------------------------------------------------------------\n"))
  cat( sprintf( "Units for material tonnage: %s\n", object$units ))

  theNames <- colnames( object$rs[, -1, drop = FALSE] )
  nRowsInTable <- length( theNames )

  statTable <- matrix( NA, nrow = nRowsInTable,
                       ncol = length( CalcTonnageStats() ) )

  dimnames(statTable) <- list( theNames, CalcTonnageStats() )

  for ( name in theNames )
    statTable[name,] <- signif( CalcTonnageStats( object$rs[ ,name] ),
                                digits = nDigits )
  cat(sprintf("\n\n"))
  print(statTable)

  cat(sprintf("\n\n"))
  cat(sprintf("Explanation\n"))
  cat(sprintf("\"Q_0.05\" is the 0.05 quantile, \"Q_0.1\" is the 0.1 quantile, and so on.\n"))
  cat(sprintf("\"Mean\" is the arithmetic mean. \"P(0)\" is probability of zero tonnage.\n"))
  cat(sprintf("\"P(>Mean)\" is probability that the tonnage exceeds the arithmetic mean.\n"))
  cat(sprintf("\n\n\n\n"))

}

#' @title Print the checks of the pdf
#' for the material tonnages in all undiscovered deposits within the permissive
#' tract
#'
#' @description Print the checks of the probability density function (pdf)
#' for the material tonnages in all undiscovered deposits within the permissive
#' tract.
#'
#' @param object
#' An object of class "TonnagePdfPT"
#'
#' @param nDigits
#' Number of signficant digits.
#'
#' @details
#' The checks were computed in the constructor function
#' TonnagePdfPT. This function prints those checks in a easy-to-read
#' format.
#'
#' The actual and predicted quantities will differ somewhat, usually in the
#' third signficant digit.
#'
#' @examples
#' pmf <- NDepositsPmf( "NegBinomial", list(theMean=5,theStdDev=4), "" )
#' pdf1 <- TonnagePdf1(ExampleTonnageData, "mt", pdfType = "normal")
#' pdfPT <- TonnagePdfPT(pmf, pdf1)
#' printChecks(pdfPT)
#'
#' @export
#'
printChecks.TonnagePdfPT <- function(object, nDigits = 3) {

  cat( sprintf( "Units for material tonnage: %s\n", object$units ))

  cat( sprintf( "\n\n"))
  cat( sprintf( "Mean vectors:\n" ))
  tmp <- cbind(object$theMean, object$predMean)
  colnames(tmp) <- c("Calculated", "Predicted")
  print(signif(tmp, digits = nDigits))
  cat(sprintf("\nExplanation\n" ))
  cat(sprintf("The calculated mean vector is calculated from random\n"))
  cat(sprintf("samples of the pdf.\n"))

  cat( sprintf( "\n\n"))
  cat( sprintf( "Standard deviation vectors:\n" ))
  tmp <- cbind(object$theSd, object$predSd)
  colnames(tmp) <- c("Calculated", "Predicted")
  print(signif(tmp, digits = nDigits))
  cat(sprintf("\nExplanation\n" ))
  cat(sprintf("The calculated standard deviation vector is calculated\n"))
  cat(sprintf(" from random samples of the pdf.\n"))


  cat(sprintf( "\n\n"))
  cat(sprintf("Composite correlation matrix\n" ))
  tmp <- object$theCor
  tmp[lower.tri(tmp)] <- object$predCor[lower.tri(object$predCor)]
  diag(tmp) <- NA
  print(signif(tmp, digits = nDigits))
  cat(sprintf("\nExplanation\n" ))
  cat(sprintf("1. The upper triangle of the composite correlation matrix\n"))
  cat(sprintf("is the upper triangle of the correlation matrix that is\n" ))
  cat(sprintf("calculated from the random samples of the pdf.\n" ))
  cat(sprintf("2. The lower triangle of the composite correlation matrix\n" ))
  cat(sprintf("is the lower triangle of the predicted correlation matrix.\n"))
  cat(sprintf("3. If the number of materials is 1, then the\n"))
  cat(sprintf("composite correlation matrix is irrelevant.\n"))

}

#' @title Construct the pdf for the material tonnages in
#' all undiscovered deposits within the
#' permissive tract
#'
#' @description Construct the probability density function (pdf) for the
#' material tonnages in all undiscovered deposits within the permissive
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
#' Number of random samples representing the pdf for the  material tonnages
#' in all undiscovered deposits within the
#' permissive tract.
#'
#' @details
#' The random samples are computed using the algorithm in Ellefsen and others
#' (2016). The predicted mean vector and the predicted covariance matrix for
#' the pdf are computed using the analytic formulas in Ellefsen and others
#' (2016). The predicted correlation matrix and the predicted standard
#' deviation vector are calculated from the predicted covariance matrix.
#'
#' @return If the input arguments have an error, the R-value NULL is returned.
#' Otherwise, a list with the following components is returned. Note that
#' these components pertain to the material tonnages in all undiscovered
#' deposits within the permissive tract.
#' @return \item{rs}{Data frame comprising the random samples of
#' the material tonnages and the associated number of undiscovered deposits.}
#' @return \item{units}{Input argument units}
#' @return \item{matNames}{Names of the materials.}
#' @return \item{theMean}{Mean vector of the material tonnages.}
#' @return \item{theCov}{Covariance matrix of the material tonnages.}
#' @return \item{theCor}{Correlation matrix of the material tonnages.}
#' @return \item{theSd}{Standard deviation vector of the material tonnages.}
#' @return \item{predMean}{Predicted mean vector of the material tonnages.}
#' @return \item{predCov}{Predicted covariance matrix of the material tonnages.}
#' @return \item{predCor}{Predicted correlation matrix of the material
#' tonnages.}
#' @return \item{predSd}{Predicted standard deviation vector of the
#' material tonnages.}
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

  if(nRandomSamples < 2500) {
    warning( sprintf( "Function TonnagePdfPT\n" ),
          sprintf( "The number of random samples should be\n" ),
          sprintf( "greater than 2500.\n" ),
          sprintf( "But the specified number is %d\n", nRandomSamples ),
          call. = FALSE )
  }

  # Number of random samples (in pdfPT) for each non-zero probability in the pmf
  # The number of random samples nRs[i] may equal zero when
  # oPmf$probs[i] * nRandomSamples < 0.5. This may occur in the tails of
  # the pmf.
  nRs <- round( oPmf$probs * nRandomSamples )

  # Number of random samples (from pdf1) needed for each non-zero probability
  # in the pmf
  nRs1 <- nRs * oPmf$nDeposits

  # All of the random samples that are needed for the computations
  # It is done this way because repeated calls to getRandomSamples can
  # require a lot of time.
  rs1 <- getRandomSamples(oPdf, sum(nRs1))

  # Get a subset of rs1
  iStart <- 1
  internalGetRandomSamples <- function(n) {
    i1 <- iStart
    i2 <- i1 + n - 1
    iStart <<- i2 + 1  # write to variable iStart in the next higher environment
    return(rs1[i1:i2, , drop = FALSE])
  }

  totalTonnages <- data.frame()
  for( i in seq_along( oPmf$nDeposits ) ) {

    # See the comment associated with variable nRs. For this special case,
    # there is nothing to add.
    if(oPmf$nDeposits[i] > 0 && nRs[i] == 0) next

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

    tmp <- cbind( nDeposits = oPmf$nDeposits[i], rs )
    totalTonnages <- rbind(totalTonnages, tmp)
  }

  # make the ordering of the rows random
  N <- nrow(totalTonnages)
  totalTonnages <- data.frame(totalTonnages[sample.int(N, size = N), ],
                              row.names = 1:N)

  theCov <- cov(totalTonnages[, -1, drop = FALSE])
  theCor <- cov2cor(theCov)
  theSd <- sqrt(diag(theCov))

  predCov <- oPmf$theMean * oPdf$theCov +
    oPmf$theVar * oPdf$theMean %o% oPdf$theMean
  predCor <- cov2cor(predCov)
  predSd <- sqrt(diag(predCov))

  rval <- list( rs = totalTonnages,
                units = getUnits(oPdf),
                matNames = colnames(totalTonnages)[-1],
                theMean = colMeans(totalTonnages[, -1, drop = FALSE]),
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