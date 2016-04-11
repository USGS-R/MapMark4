#' @title Get the units for tonnage
#'
#' @description Get the units for tonnage.
#'
#' @param object
#' An object of class "TonnagePdf1"
#'
#' @return
#' A character string containing the units.
#'
#' @examples
#' pdf1 <- TonnagePdf1(ExampleTonnageData, "mt")
#' cat(sprintf("Units for tonnage: %s\n", getUnits(pdf1)))
#'
#' @export
#'
getUnits <- function(object) {
  return(object$units)
}

#' @title Get random samples from the pdf for material tonnages in a
#' single, undiscovered deposit
#'
#' @description Get random samples from the probability density function (pdf)
#' for the
#' material tonnages in a single, undiscovered deposit within the permissive
#' tract.
#'
#' @param object
#' An object of class "TonnagePdf1"
#'
#' @param nSamples
#' Integer specifying the number of random samples
#'
#' @param seed
#' Integer containing the seed.
#'
#' @param log_rs
#' Logical variable indicating whether the random sample are transformed
#' with the natural logarithm.
#'
#' @details
#' To generate random samples, the known material tonnages are
#' transformed with the natural logarithm.
#' If the \code{type} is \code{empirical}, then
#' a multivariate kernel density estimate
#' of the log-transformed tonnages is used to
#' generate log-transformed random samples
#' (Duong, 2007; Hastie and others, 2009, p. 208-209;
#' Shalizi, 2016, p. 308-330).
#' In contrast, if the \code{type} is \code{normal},
#' then a multivariate normal distribution is used to generate
#' log-transformed random samples.
#'
#' The log-transformed random samples are converted to random samples of
#' material tonnage with exponentiation.
#'
#' @return Matrix containing random samples from the pdf for the
#' material tonnages for a single, undiscovered deposit in the permissive
#' tract.
#'
#' @references
#' Duong, Tarn, 2007, ks - Kernel density estimation and kernel discriminant
#' analysis for multivariate data in R: Journal of Statistical Software,
#' v. 21, issue 7, \url{http://www.jstatsoft.org/}
#'
#' Hastie, Tevor, Tibshirani, Robert, and Friedman, Jerome, 2009,
#' The elements of statistical learning - Data mining, inference, and
#' prediction (2nd ed.): New York, Springer Science + Business Media, LLC, 745 p.
#'
#' Shalizi, C.R., 2016, Advanced data analysis from an elementary point of
#' view: Draft book manuscript publicly available at
#' \url{http://www.stat.cmu.edu/~cshalizi/ADAfaEPoV/}
#'
#' @examples
#' pdf1 <- TonnagePdf1(ExampleTonnageData, "mt")
#' rs <- getRandomSamples(pdf1, 2518)
#'
#' @export
#'
getRandomSamples.TonnagePdf1 <- function(object, nSamples, seed = NULL, log_rs = FALSE) {

  # The number of random samples that are generated is 2 * nSamples because,
  # if the random samples are truncated, then there will be enough remaining
  # so that nSamples random samples can be returned.

  N <- 2 * nSamples

  set.seed(seed)

  if(object$pdfType == "empirical") {
    fhat <- ks::kde(x = object$logTonnages, H = ks::Hpi(x = object$logTonnages))

    indices <- sample( 1:nrow(object$logTonnages), N, replace=TRUE )

    # For each element of indices, a random sample could be drawn. However,
    # this results in very slow code. Instead, determine the unique indices
    # (uniqueIndices) and the number of occurrences for each unique index
    # (counts). Then for uniqueIndices[i], draw counts[i] random samples.
    tmp <- table(indices)
    uniqueIndices <- as.integer(names(tmp))
    counts <- as.vector(tmp, mode = "integer")

    theCov <- as.matrix(fhat$H)
    rsLogTonnage <- NULL
    for(i in seq_along(uniqueIndices)) {
      index <- uniqueIndices[i]
      rs <- mvtnorm::rmvnorm( counts[i],
                              mean = object$logTonnages[index, , drop = FALSE],
                              sigma = theCov)
      rsLogTonnage <- rbind(rsLogTonnage, rs)
    }

    # Make the draws random
    rsLogTonnage <- rsLogTonnage[sample.int(N, size = N), ]

  } else {
    rsLogTonnage <- mvtnorm::rmvnorm(N,
                                     mean = colMeans(object$logTonnages),
                                     sigma = cov(object$logTonnages))
  }

  if(object$isTruncated == TRUE){
    areWithinBnds <- rep.int(TRUE, N)
    for(j in 1:ncol(rsLogTonnage)) {
      theRange <- range(object$logTonnages[, j])
      areWithinBnds <- areWithinBnds &
        theRange[1] <= rsLogTonnage[, j] &
        rsLogTonnage[, j] <= theRange[2]
    }
    rsLogTonnage <- rsLogTonnage[areWithinBnds, ]
  }

  rsLogTonnage <- rsLogTonnage[1:nSamples, , drop = FALSE]
  colnames(rsLogTonnage) <- object$matNames

  if(log_rs == FALSE) {
    return(exp(rsLogTonnage))
  } else {
    return(rsLogTonnage)
  }

}

#' @title Plot the univariate, marginal cdfs for the material tonnages
#' in a single, undiscovered deposit
#'
#' @description Plot the univariate, marginal cumulative distribution
#' functions (cdfs)
#' for the material tonnages in a single, undiscovered deposit within the
#' permissive tract. Overlaid on the plot is the
#' empirical cumulative distribution function (ecdf)
#' for the known material tonnages.
#'
#' @param object
#' An object of class "TonnagePdf1".
#'
#' @param isUsgsStyle
#' Make the plot format similar to the U.S. Geological Survey style
#'
#' @details
#' An internal call to function reshape2::melt generates the
#' message "No id variables; using all as measure variables", which should
#' be ignored.
#'
#' In the plot, the solid line(s) represent the cdf, and the dots represent
#' the ecdf.
#'
#' If there is only one material tonnage, the deviance, which measures the
#' misfit between the known material tonanges and the pdf, is printed at the
#' at the top. If there are two or more material tonnages, then the sum
#' of the deviances for each material tonnage is printed.
#'
#' @examples
#' pdf1 <- TonnagePdf1(ExampleTonnageData, "mt")
#' plot(pdf1)
#'
#' @export
#'
plot.TonnagePdf1 <- function(object,
                             isUsgsStyle = TRUE) {

  df.rs <- reshape2::melt(object$rs)
  # After the melt operation, the first column is row number. It is
  # useless here, so it is removed.
  df.rs <- df.rs[, -1, drop = FALSE]
  colnames(df.rs) <- c("Material", "Tonnage")

  df.obs <- reshape2::melt(object$knownTonnages[, -1, drop = FALSE])
  colnames(df.obs) <- c("Material", "Tonnage")

  xLabel <- paste("Tonnage (", object$units, ")", sep = "")

  if(ncol(object$rs) == 1){
    title <- paste("Deviance = ",
                   signif(object$sumDeviance, digits = 3), sep = "")
  } else {
    title <- paste("Sum of the deviances = ",
                   signif(object$sumDeviance, digits = 3), sep = "")
  }

  p <- ggplot2::ggplot(df.rs) +
    ggplot2::stat_ecdf(ggplot2::aes(Tonnage, colour = Material)) +
    ggplot2::stat_ecdf(ggplot2::aes(Tonnage, colour = Material),
                       data = df.obs, geom = "point") +
    ggplot2::scale_x_continuous(name = xLabel, trans = "log10") +
    ggplot2::ylab("Probability") +
    ggplot2::ggtitle(title)

  if(isUsgsStyle)
    p <- p + ggplot2::xlab(paste("Tonnage, in ", object$units, sep = "")) +
    ggplot2::theme_bw()

  plot(p)

}

#' @title Summarize the pdf for the material tonnages in a single, undiscovered
#' deposit
#'
#' @description Summarize the probability density function (pdf) for the
#' material tonnages in a single, undiscovered deposit within the permissive
#' tract.
#'
#' @param object
#' An object of class "TonnagePdf1"
#'
#' @param nDigits
#' Number of signficant digits.
#'
#' @details
#' It is common that the summary statistics for the known material
#' tonnages differ somewhat from the summary statistics for the pdf, especially
#' when the pdf is truncated. The reason for the difference is that tonnages
#' typically have an enormous range, making the summary statistics somewhat
#' non-robust.
#'
#' @examples
#' pdf1 <- TonnagePdf1(ExampleTonnageData, "mt")
#' summary(pdf1)
#'
#' @export
#'
summary.TonnagePdf1 <- function(object, nDigits = 2) {

  cat(sprintf("Summary of the pdf for the material tonnages in a single,\n"))
  cat(sprintf("undiscovered deposit within the permissive tract.\n"))
  cat(sprintf("------------------------------------------------------------\n"))
  cat( sprintf( "Units for material tonnage: %s\n", object$units ))
  cat( sprintf( "Pdf type: %s\n", object$pdfType ))
  if(object$isTruncated) {
    cat( sprintf("Pdf is truncated at the lowest and the highest values\n"))
    cat( sprintf("of the material tonnages in the discovered deposits.\n"))
  } else {
    cat( sprintf("Pdf is not truncated.\n"))
  }
  cat( sprintf( "Number of discovered deposits: %d\n", nrow(object$knownTonnages) ))
  cat( sprintf( "Number of materials: %d\n", length(object$matNames) ))
  cat( sprintf( "Materials: " ))
  for(i in 1:length(object$matNames)) {
    cat( sprintf( "%s    ", object$matNames[i] ))
  }
  cat( sprintf( "\n" ))

  if(length(object$matNames) == 1) {
    cat(sprintf("Deviance = %g\n", object$sumDeviance))
  } else {
    cat(sprintf("Sum of deviances = %g\n", object$sumDeviance))
  }

  cat( sprintf( "\n\n"))
  cat( sprintf( "Summary statistics for the material tonnages (discovered deposits):\n" ))
  print(summary(object$knownTonnages[, -1]))
  cat( sprintf( "\n"))
  cat( sprintf( "Summary statistics for the pdf:\n" ))
  print(summary(object$rs))

  cat(sprintf("\n\n"))
  cat(sprintf("Explanation\n"))
  cat(sprintf("\"1st Qu.\" refers to the first quartile.\n"))
  cat(sprintf("\"3rd Qu.\" refers to the third quartile.\n"))
  cat(sprintf("\n\n\n\n"))


}

#' @title Print checks of the
#' pdf for the material tonnages in a single, undiscovered
#' deposit
#'
#' @description Print checks of the
#' probability density function (pdf) for the
#' material tonnages in a single, undiscovered deposit within the permissive
#' tract. Summary statistics are calculated for both
#' the known material tonnages and the pdf that represents those tonnages.
#'
#' @param object
#' An object of class "TonnagePdf1"
#'
#' @param nDigits
#' Number of signficant digits.
#'
#' @details
#' The statistics for the pdf are calculated from
#' random samples of that pdf.
#'
#' It is common that the statistics for the known material
#' tonnages differ somewhat from the statistics for the pdf, especially
#' when the pdf is truncated. The reason for the difference is that tonnages
#' typically have an enormous range, making the statistics somewhat
#' non-robust.
#'
#' @examples
#' pdf1 <- TonnagePdf1(ExampleTonnageData, "mt", pdfType = "empirical")
#' printChecks(pdf1)
#'
#' @export
#'
printChecks.TonnagePdf1 <- function(object, nDigits = 3) {

  cat( sprintf( "\n\n"))
  cat( sprintf( "Mean vectors for material tonnages:\n" ))
  tmp <- cbind(object$theKnownMean, object$theMean)
  colnames(tmp) <- c("Discovered", "Pdf")
  print(signif(tmp, digits = nDigits))

  cat( sprintf( "\n\n"))
  cat( sprintf( "Standard deviation vectors for material tonnages:\n" ))
  tmp <- cbind(object$theKnownSd, object$theSd)
  colnames(tmp) <- c("Discovered", "Pdf")
  print(signif(tmp, digits = nDigits))

  cat(sprintf( "\n\n"))
  cat(sprintf("Composite correlation matrix\n" ))
  tmp <- object$theKnownCor
  tmp[lower.tri(tmp)] <- object$theCor[lower.tri(object$theCor)]
  diag(tmp) <- NA
  print(signif(tmp, digits = nDigits))
  cat(sprintf("\nExplanation\n" ))
  cat(sprintf("The upper triangle is the upper triangle of the\n"))
  cat(sprintf("correlation matrix for the material tonnages.\n"))
  cat(sprintf("from the discovered deposits.\n"))
  cat(sprintf("The lower triangle is the lower triangle of the\n" ))
  cat(sprintf("correlation matrix for the material tonnages\n"))
  cat(sprintf("that are represented by the pdf.\n"))


}


#' @title Construct the pdf for the material tonnages in a single,
#' undiscovered deposit
#'
#' @description Construct the probability density function (pdf) for the
#' material tonnages in a single, undiscovered deposit within the permissive
#' tract. The pdf is not explicitly specified; instead it is implicitly
#' specified with the random samples that are generated from it.
#'
#' @param knownTonnages
#' Dataframe containing the known material tonnages. (See details.)
#'
#' @param units
#' Character string containing units for the material tonnages.
#'
#' @param pdfType
#' Character string containing the type of pdf for
#' the log-transformed material tonnages. The choices are either
#' "empirical" or "normal".
#'
#' @param isTruncated
#' Logical variable indicating whether the pdf is
#' truncated at the lowest and highest values of the known material tonnages.
#'
#' @param minNDeposits
#' Minimum number of deposits with known material tonnages.
#'
#' @param nRandomSamples
#' Number of random samples used to compute summary statistics and the
#' marginal cumulative distribution functions.
#'
#' @param seed
#' Seed for the random number generator.
#'
#' @details
#' Data frame knownTonnages comprises the material tonnages for the known
#' deposits. Each row comprises the data for one deposit.
#' The first column lists the names of the deposits.
#' The second column and subsequent columns, if any, list
#' the tonnages for each type of material.
#' All material tonnages in the data frame must be greater than zero and
#' must not be missing.
#' The minimum number of
#' deposits (that is, rows in the data frame) should be greater than or
#' equal to 20.
#'
#' The columns of the data frame have headings. The heading for the
#' first column is "Name" or "Deposit name".
#' The headings for the second column and subsequent columns, if any, are the
#' names of the materials. For example, they might be "Ore" and "Cu".
#' The headings for the second and subsequent columns are important because
#' they are used in plots and tables.
#'
#' The misfit between the known material tonnages and the pdf that represents
#' those material tonnages is quantified with the deviance (McElreath, 2016,
#' p. 177-182). If there are two or more material tonnages, then the deviance is
#' calculated for each material tonnage and the individual deviances are summed.
#' Smaller values of the deviance indicate a better fit.
#'
#' @return If the input arguments have an error, the R-value NULL is returned.
#' Otherwise, a list with the following components is returned.
#' @return \item{knownTonnages}{Input argument knownTonnages.}
#' @return \item{units}{Input argument units}
#' @return \item{pdfType}{Input argument pdfType.}
#' @return \item{isTruncated}{Input argument isTruncated.}
#' @return \item{matNames}{Names of the materials.}
#' @return \item{logTonnages}{Matrix with the natural logarithm
#' transformation of the known material tonnages.}
#' @return \item{theKnownMean}{Mean vector for the known tonnages.}
#' @return \item{theKnownCov}{Covariance matrix for the known tonnages.}
#' @return \item{theKnownSd}{Standard deviation vector for the known tonnages.}
#' @return \item{theKnownCor}{Correlation matrix for the known tonnages.}
#' @return \item{rs}{Random samples of the pdf.}
#' @return \item{theMean}{Mean vector for the pdf.}
#' @return \item{theCov}{Covariance matrix for the pdf.}
#' @return \item{theSd}{Standard deviation vector for the pdf.}
#' @return \item{theCor}{Correlation matrix for the pdf.}
#' @return \item{sumDeviance}{Sum of the deviances measuring the relative fit of
#' the pdf.}
#' @return \item{call}{Function call.}
#'
#' @references
#' McElreath, Richard, 2016, Statistical rethinking - A Bayesian course
#' with examples in R and Stan: New York, CRC Press, 469 p.
#'
#' @examples
#' pdf1 <- TonnagePdf1(ExampleTonnageData, "mt")
#'
#' @export
#'
TonnagePdf1 <- function(knownTonnages,
                        units,
                        pdfType = "empirical",
                        isTruncated = TRUE,
                        minNDeposits = 20,
                        nRandomSamples = 20000, seed = 7) {

  CalcSumDeviance <- function(rsPdf, rsData, nBins = 30) {

    logRsPdf <- log(rsPdf)
    logRsData <- log(rsData)

    N <- ncol(logRsPdf)
    deviance <- 0
    for(i in 1:N) {
      theBreaks <- seq(from = min(logRsData[, i], logRsPdf[, i]),
                       to = max(logRsData[, i], logRsPdf[, i]),
                       length.out = nBins)

      tmpPdf <- hist(logRsPdf[, i], breaks = theBreaks, plot = FALSE)

      tmpData <- hist(logRsData[, i], breaks = theBreaks, plot = FALSE)

      deviance <- deviance - 2 * sum(tmpPdf$density * tmpData$counts)
    }

    return(deviance)
  }

  if(ncol(knownTonnages) < 2) {
    stop( sprintf( "Function TonnagePdf1\n" ),
          sprintf( "The number of columns in dataframe knownTonnage must be.\n" ),
          sprintf( "greater than or equal to 2.\n" ),
          sprintf( "But the actual number is 1.\n" ),
          call. = FALSE )

  }
  if(nrow(knownTonnages) < minNDeposits) {
    stop( sprintf( "Function TonnagePdf1\n" ),
          sprintf( "The number of rows in dataframe knownTonnage must be.\n" ),
          sprintf( "greater than or equal to %d.\n", minNDeposits ),
          sprintf( "But the actual number is %d.\n", nrow(knownTonnages) ),
          call. = FALSE )
  }

  for(j in 2:ncol(knownTonnages)){
    if(any(is.na(knownTonnages[, j]))) {
      stop( sprintf( "Function TonnagePdf1\n" ),
            sprintf( "Column %d of dataframe knownTonnage has\n", j ),
            sprintf( "one or more missing values\n"),
            call. = FALSE )
    }
    if(any(knownTonnages[, j] <= 0.0)) {
      stop( sprintf( "Function TonnagePdf1\n" ),
            sprintf( "Column %d of dataframe knownTonnage has\n", j ),
            sprintf( "one or more zero values\n"),
            call. = FALSE )
    }
  }

  if(!any(pdfType == c("empirical", "normal"))) {
    stop( sprintf( "Function TonnagePdf1\n" ),
          sprintf( "Argument type must be either empirical or normal.\n", j ),
          sprintf( "It is specified as %s\n", pdfType),
          call. = FALSE )
  }

  theKnownMean <- colMeans(knownTonnages[, -1])
  theKnownCov <- cov(knownTonnages[, -1])
  theKnownSd <- sqrt(diag(theKnownCov))
  theKnownCor <- cov2cor(theKnownCov)

  rval <- list( knownTonnages = knownTonnages,
                units = units,
                pdfType = pdfType,
                isTruncated = isTruncated,
                matNames = colnames(knownTonnages)[-1],
                logTonnages = as.matrix(log(knownTonnages[, -1])),
                theKnownMean = theKnownMean,
                theKnownCov = theKnownCov,
                theKnownSd = theKnownSd,
                theKnownCor = theKnownCor)

  class(rval) <- "TonnagePdf1"

  rval$rs <- getRandomSamples(rval, nRandomSamples, seed = seed)
  rval$theMean <- colMeans(rval$rs)
  rval$theCov <- cov(rval$rs)
  rval$theSd <- sqrt(diag(rval$theCov))
  rval$theCor <- cov2cor(rval$theCov)
  rval$sumDeviance <- CalcSumDeviance(rval$rs, rval$knownTonnages[, -1])
  rval$call <- sys.call()

  return(rval)
}