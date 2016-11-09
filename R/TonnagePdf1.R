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

#' @title Get random samples from the pdf for material tonnages in one
#' undiscovered deposit
#'
#' @description Get random samples from the probability density function (pdf)
#' for the
#' material tonnages in one undiscovered deposit within the permissive
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
#' material tonnages for one undiscovered deposit in the permissive
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

    # Ensure that random samples are in a matrix --- there is only one
    # pathologic case
    if(ncol(object$logTonnages) == 1){
      rsLogTonnage <- matrix(rsLogTonnage, ncol = 1)
    }

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
    rsLogTonnage <- rsLogTonnage[areWithinBnds, , drop = FALSE]
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
#' in one undiscovered deposit
#'
#' @description Plot the univariate, marginal cumulative distribution
#' functions (cdfs)
#' for the material tonnages in one undiscovered deposit within the
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
plot.TonnagePdf1 <- function(object, isUsgsStyle = TRUE) {

  # If there are too many random samples, subset them. Then change
  # the data structure
  n <- nrow(object$rs)
  if(n > 5000) {
    indices <- sample.int(n, size = 5000)
    df.rs <- reshape2::melt(object$rs[indices, , drop = FALSE])
  } else {
    df.rs <- reshape2::melt(object$rs)
  }

  # After the melt operation, the first column is row number. It is
  # useless here, so it is removed.
  df.rs <- df.rs[, -1, drop = FALSE]
  colnames(df.rs) <- c("Material", "Tonnage")

  df.obs <- reshape2::melt(object$knownTonnages[, -1, drop = FALSE])
  colnames(df.obs) <- c("Material", "Tonnage")

  if(isUsgsStyle) {
    xLabel <- paste("Tonnage, in ", object$units, sep = "")
  } else {
    xLabel <- paste("Tonnage (", object$units, ")", sep = "")
  }

  if(ncol(object$rs) == 1){
    caption <- paste("Deviance = ",
                   signif(object$sumDeviance, digits = 3), sep = "")
  } else {
    caption <- paste("Sum of the deviances = ",
                   signif(object$sumDeviance, digits = 3), sep = "")
  }

  p <- ggplot2::ggplot(df.rs) +
    ggplot2::stat_ecdf(ggplot2::aes(Tonnage, colour = Material),
                       pad = FALSE, geom = "step") +
    ggplot2::stat_ecdf(ggplot2::aes(Tonnage, colour = Material),
                       data = df.obs, geom = "point") +
    ggplot2::scale_x_continuous(name = xLabel, trans = "log10") +
    ggplot2::ylab("Probability") +
    ggplot2::geom_text(ggplot2::aes(x, y, label = caption),
                       data = data.frame(x = min(df.rs$Tonnage), y = 0),
                       hjust = 0, vjust = 1)

  if(isUsgsStyle) {
    p <- p + ggplot2::theme_bw()
  }

  plot(p)

}

#' @title Plot the marginal distributions, as a matrix, for the
#' material tonnages
#' in one undiscovered deposit
#'
#' @description Plot the univariate and bivariate marginal distributions,
#' as a matrix, for the material tonnages
#' in one undiscovered deposit. This matrix shows how well the
#' probability density function (pdf) represents the known material tonnages.
#'
#' @param object
#' An object of class "TonnagePdf1".
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
#' univariate
#' distributions. Along the horizontal axis of each histgram is a (red)
#' rug plot of the corresponding, known material tonnages.
#'
#' In the upper triangle of the matrix are scatter plots
#' showing the marginal, bivariate distributions. These scatter plots use
#' blue dots. Overlaid on each plot is
#' another scatter plot for the corresponding, known material tonnages.
#' These overlaying plots use red dots.
#'
#' In the lower triangle of the matrix are contour plots
#' showing the marginal, bivariate distributions. Overlaid on each plot is
#' a scatter plot for the corresponding, known material tonnages. These
#' overlaying plots use red dots.
#'
#' @examples
#' pdf1 <- TonnagePdf1(ExampleTonnageData, "mt")
#' plotMatrix(pdf1)
#'
#' @export
#'
plotMatrix.TonnagePdf1 <- function(object, nPlotSamples = 1000,
                                   isUsgsStyle = TRUE) {

  d <- as.data.frame(object$rs[1:nPlotSamples, , drop = FALSE])
  N <- ncol(d)
  f <- object$knownTonnages[, -1, drop = FALSE]

  if(isUsgsStyle) {
    tLabel <- paste("tonnage, in ", object$units, sep = "")
  } else {
    tLabel <- paste("tonnage (", object$units, ")", sep = "")
  }

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(N,N)))

  for (i in 1:N) {
    for(j in 1:N) {
      if(i == j) {
        # diagonal
        matName <- colnames(d)[i]
        p <- ggplot2::ggplot() +
          ggplot2::geom_histogram(ggplot2::aes_string(x = matName,
                                                      y = "..density.."),
                                  data = d, bins = 15,
                                  alpha = 0.3, colour = "white",
                                  fill = "blue") +
          ggplot2::geom_rug(ggplot2::aes_string(x = matName),
                            data = f, alpha = 0.3, colour = "red") +
          ggplot2::scale_x_continuous(name = paste( matName, tLabel, sep = " "),
                                      trans = "log10") +
          ggplot2::scale_y_continuous(name = "Density")
      } else if (j < i) {
        # lower triangle
        xMatName <- colnames(d)[j]
        yMatName <- colnames(d)[i]
        p <- ggplot2::ggplot() +
          ggplot2::geom_density_2d(ggplot2::aes_string(x = xMatName,
                                                       y = yMatName),
                                   data = d) +
          ggplot2::geom_point(ggplot2::aes_string(x = xMatName, y = yMatName),
                              data = f, alpha = 0.5, colour = "red", pch = 16) +
          ggplot2::scale_x_continuous(name =
                                        paste( xMatName, tLabel, sep = " "),
                                      trans = "log10") +
          ggplot2::scale_y_continuous(name =
                                        paste( yMatName, tLabel, sep = " "),
                                      trans = "log10")

      } else {
        # upper triangle
        xMatName <- colnames(d)[j]
        yMatName <- colnames(d)[i]
        p <- ggplot2::ggplot() +
          ggplot2::geom_point(ggplot2::aes_string(x = xMatName, y = yMatName),
                              data = d, alpha = 0.2,
                              colour = "blue", pch = 16) +
          ggplot2::geom_point(ggplot2::aes_string(x = xMatName, y = yMatName),
                              data = f, alpha = 0.5,
                              colour = "red", pch = 16) +
          ggplot2::scale_x_continuous(name =
                                        paste( xMatName, tLabel, sep = " "),
                                      trans = "log10") +
          ggplot2::scale_y_continuous(name =
                                        paste( yMatName, tLabel, sep = " "),
                                      trans = "log10")

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

      print(p, vp=grid::viewport(layout.pos.row=i, layout.pos.col=j))
    }
  }
}


#' @title Summary comparison of the pdf for the material tonnages in one
#' undiscovered deposit and the known material tonnages in the discovered
#' deposits
#'
#' @description Using summary statistics, this function compares
#' the probability density function (pdf)
#' for the material tonnages in one
#' undiscovered deposit and the known material tonnages in the discovered
#' deposits. The pdf is represented by random samples drawn from it.
#' The summary statistics are minimum, 0.25 quantile, median,
#' mean, 0.75 quantile, maximum, standard deviation, and correlation matrix.
#' The comparison is involves both the log-transformed material tonnages
#' and the untransformed material tonnages.
#'
#' @param object
#' An object of class "TonnagePdf1"
#'
#' @param nDigits
#' Number of signficant digits.
#'
#' @details
#' It is common that the corresponding statistics differ somewhat. For
#' statistics calculated with log-transformed material tonnages, the differences
#' should be small. There are two, possible exceptions:
#' The minimum and the maximum might differ a lot, if the pdf is not
#' truncated (see function TonnagePdf1).
#' These statistics are particularly important because the pdf is fit to the
#' log-transformed material tonnages from the discovered deposits.
#'
#' For statistics calculated with untransformed material tonnages, the
#' differences between corresponding statistics may be large. The reason is
#' that the pdf is fit to the log-transformed material tonnages, not the
#' untransformed material tonnages (which would be too difficult).
#'
#' @examples
#' pdf1 <- TonnagePdf1(ExampleTonnageData, "mt")
#' summary(pdf1)
#'
#' @export
#'
summary.TonnagePdf1 <- function(object, nDigits = 3) {

  PrintStat <- function(x, y, nDigits, FUNC, ...){
    tmp <- cbind(apply(x, 2, FUNC, ...), apply(y, 2, FUNC, ...))
    colnames(tmp) <- c("Known", "Pdf")
    print(signif(tmp, digits = nDigits))
  }

  PrintCor <- function(x, y, nDigits){
    tmp <- cor(x)
    tmp1 <- cor(y)
    tmp[lower.tri(tmp)] <- tmp1[lower.tri(tmp1)]
    diag(tmp) <- NA
    print(signif(tmp, digits = nDigits))
  }

  cat(sprintf("Summary comparison of the pdf for the material tonnages in\n"))
  cat(sprintf("one undiscovered deposit and the known material tonnages\n"))
  cat(sprintf("in the discovered deposits.\n"))
  cat(sprintf("------------------------------------------------------------\n"))
  cat( sprintf( "Units for material tonnage: %s\n", object$units ))
  cat( sprintf( "Pdf type: %s\n", object$pdfType ))
  if(object$isTruncated) {
    cat( sprintf("Pdf is truncated at the lowest and the highest values\n"))
    cat( sprintf("of the known material tonnages in the discovered deposits.\n"))
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

  cat( sprintf( "\n###############################################################\n"))
  cat( sprintf( "This section pertains to the log-transformed material tonnages.\n"))

  logRs <- log(object$rs)

  cat( sprintf( "\nMean\n" ))
  PrintStat(object$logTonnages, logRs, nDigits, mean)

  cat( sprintf( "\nStandard deviation\n" ))
  PrintStat(object$logTonnages, logRs, nDigits, sd)

  cat( sprintf( "\nMinimum\n" ))
  PrintStat(object$logTonnages, logRs, nDigits, min)

  cat( sprintf( "\n0.25 quantile\n" ))
  PrintStat(object$logTonnages, logRs, nDigits, quantile, prob = 0.25)

  cat( sprintf( "\nMedian\n" ))
  PrintStat(object$logTonnages, logRs, nDigits, median)

  cat( sprintf( "\n0.75 quantile\n" ))
  PrintStat(object$logTonnages, logRs, nDigits, quantile, prob = 0.75)

  cat( sprintf( "\nMaximum\n" ))
  PrintStat(object$logTonnages, logRs, nDigits, max)

  cat(sprintf("\nComposite correlation matrix\n" ))
  PrintCor(object$logTonnages, logRs, nDigits)

  cat( sprintf( "\n###############################################################\n"))
  cat( sprintf( "This section pertains to the (untransformed) material tonnages.\n"))

  knownTonnages <- object$knownTonnages[, -1, drop = FALSE]

  cat( sprintf( "\nMean\n" ))
  PrintStat(knownTonnages, object$rs, nDigits, mean)

  cat( sprintf( "\nStandard deviation\n" ))
  PrintStat(knownTonnages, object$rs, nDigits, sd)

  cat( sprintf( "\nMinimum\n" ))
  PrintStat(knownTonnages, object$rs, nDigits, min)

  cat( sprintf( "\n0.25 quantile\n" ))
  PrintStat(knownTonnages, object$rs, nDigits, quantile, prob = 0.25)

  cat( sprintf( "\nMedian\n" ))
  PrintStat(knownTonnages, object$rs, nDigits, median)

  cat( sprintf( "\n0.75 quantile\n" ))
  PrintStat(knownTonnages, object$rs, nDigits, quantile, prob = 0.75)

  cat( sprintf( "\nMaximum\n" ))
  PrintStat(knownTonnages, object$rs, nDigits, max)

  cat(sprintf("\nComposite correlation matrix\n" ))
  PrintCor(knownTonnages, object$rs, nDigits)

  cat( sprintf( "\n###############################################################\n"))

  cat(sprintf("\n\n"))
  cat(sprintf("Explanation\n"))
  cat(sprintf("1. The composite correlation matrix has two parts: its upper\n"))
  cat(sprintf("triangle and its lower triangle. The upper triangle pertains\n"))
  cat(sprintf("to the (log-transformed/untransformed) known material\n"))
  cat(sprintf("tonnages in the discovered deposits. The lower triangle\n" ))
  cat(sprintf("pertains to the (log-transformed/untransformed) material\n"))
  cat(sprintf("tonnages that are represented by the pdf.\n"))
  cat(sprintf("2. If the number of materials is 1, then the\n"))
  cat(sprintf("composite correlation matrix is irrelevant.\n"))
  cat(sprintf("\n\n"))

}


#' @title Construct the pdf for the material tonnages in one
#' undiscovered deposit
#'
#' @description Construct the probability density function (pdf) for the
#' material tonnages in one undiscovered deposit within the permissive
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
#' marginal cumulative distribution functions. It should be a large value
#' to ensure that the summary statistics are precise.
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
#' If the \code{pdfType} is \code{empirical}, then
#' a multivariate kernel density estimate
#' of the log-transformed, known tonnages is used to
#' generate log-transformed random samples
#' (Duong, 2007; Hastie and others, 2009, p. 208-209;
#' Shalizi, 2016, p. 308-330).
#' In contrast, if the \code{pdfType} is \code{normal},
#' then a multivariate normal distribution is used to generate
#' log-transformed random samples. For both cases,
#' the log-transformed random samples are converted to random samples of
#' material tonnage with exponentiation.
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
#' Duong, Tarn, 2007, ks - Kernel density estimation and kernel discriminant
#' analysis for multivariate data in R: Journal of Statistical Software,
#' v. 21, issue 7, \url{http://www.jstatsoft.org/}
#'
#' Hastie, Tevor, Tibshirani, Robert, and Friedman, Jerome, 2009,
#' The elements of statistical learning - Data mining, inference, and
#' prediction (2nd ed.): New York, Springer Science + Business Media, LLC, 745 p.
#'
#' McElreath, Richard, 2016, Statistical rethinking - A Bayesian course
#' with examples in R and Stan: New York, CRC Press, 469 p.
#'
#' Shalizi, C.R., 2016, Advanced data analysis from an elementary point of
#' view: Draft book manuscript publicly available at
#' \url{http://www.stat.cmu.edu/~cshalizi/ADAfaEPoV/}
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
                        nRandomSamples = 1000000, seed = 7) {

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

  theKnownMean <- colMeans(knownTonnages[, -1, drop = FALSE])
  theKnownCov <- cov(knownTonnages[, -1, drop = FALSE])
  theKnownSd <- sqrt(diag(theKnownCov))
  theKnownCor <- cov2cor(theKnownCov)

  rval <- list( knownTonnages = knownTonnages,
                units = units,
                pdfType = pdfType,
                isTruncated = isTruncated,
                matNames = colnames(knownTonnages)[-1],
                logTonnages = as.matrix(log(knownTonnages[, -1, drop = FALSE])),
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
  rval$sumDeviance <- CalcSumDeviance(rval$rs, rval$knownTonnages[, -1, drop = FALSE])
  rval$call <- sys.call()

  return(rval)
}