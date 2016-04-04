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
#' exData <- ExampleTonnageData
#' pdf1 <- TonnagePdf1(exData, "mt")
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
#' To generate random samples, the observed material tonnages are
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
#' prediction, 2nd ed., Springer Science + Business Media, LLC, 745 p.
#'
#' Shalizi, C.R., 2016, Advanced data analysis from an elementary point of
#' view: Draft book manuscript publicly available at
#' \url{http://www.stat.cmu.edu/~cshalizi/ADAfaEPoV/}
#'
#' @examples
#' exData <- ExampleTonnageData
#' pdf1 <- TonnagePdf1(exData, "mt")
#' rs <- getRandomSamples(pdf1, 2518)
#'
#' @export
#'
getRandomSamples <- function(object, nSamples, seed = NULL, log_rs = FALSE) {

  # The number of random samples that are generated is 2 * nSamples because,
  # if the random samples are truncated, then there will be enough remaining
  # so that nSamples random samples can be returned.

  N <- 2 * nSamples

  set.seed(seed)

  if(object$pdfType == "empirical") {
    fhat <- ks::kde(x = object$logTonnages, H = ks::Hpi(x = object$logTonnages))

    indices <- sample( 1:nrow(object$logTonnages), N, replace=TRUE )
    theMeans <- object$logTonnages[indices, , drop = FALSE]
    theCov <- as.matrix(fhat$H)

    rsLogTonnage <- matrix(NA_real_, nrow = N, ncol=ncol(fhat$x) )
    for(i in 1:N) {
      rsLogTonnage[i, ] <- mvtnorm::rmvnorm( 1, mean = theMeans[i, ],
                                             sigma = theCov)
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

#' @title Plot the cdf for the material tonnages in a single,
#' undiscovered deposit
#'
#' @description Plot the cumulative distribution function (cdf)
#' for the material tonnages in a single, undiscovered deposit within the
#' permissive tract. Overlaid on the plot is the
#' empirical cumulative distribution function (ecdf)
#' for the observed material tonnages.
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
#' In the plot, the solid line(s) represent the cdf, and the dots represent
#' the ecdf.
#'
#' @examples
#' pdf1 <- TonnagePdf1(ExampleTonnageData, "mt")
#' plot(pdf1)
#'
#' @export
#'
plot.TonnagePdf1 <- function(object,
                             whichMatTonnage = object$matNames,
                             isUsgsStyle = TRUE) {

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

  df.obs <- reshape2::melt(object$obsTonnages[, whichMatTonnage, drop = FALSE])
  colnames(df.obs) <- c("Material", "Tonnage")

  xLabel <- paste("Tonnage (", object$units, ")", sep = "")

  p <- ggplot2::ggplot(df.rs) +
    ggplot2::stat_ecdf(ggplot2::aes(Tonnage, colour = Material)) +
    ggplot2::stat_ecdf(ggplot2::aes(Tonnage, colour = Material),
                       data = df.obs, geom = "point") +
    ggplot2::scale_x_continuous(name = xLabel, trans = "log10") +
    ggplot2::ylab("Probability")

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
#' tract. Summary statistics are calculated for both
#' the observed material tonnages and
#' the pdf that represents those tonnages.
#'
#' @param object
#' An object of class "TonnagePdf1"
#'
#' @param nDigits
#' Number of signficant digits.
#'
#' @details
#' The summary statistics for the pdf are calculated from
#' random samples of that pdf.
#'
#' It is common that the summary statistics for the observed material
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

  cat( sprintf( "Units for material tonnage: %s\n", object$units ))
  cat( sprintf( "Pdf type: %s\n", object$pdfType ))
  if(object$isTruncated) {
    cat( sprintf("Pdf is truncated at the lowest and "))
    cat( sprintf("highest values of the observed material tonnages.\n"))
  } else {
    cat( sprintf("Pdf is not truncated.\n"))
  }
  cat( sprintf( "Number of observed deposits: %d\n", nrow(object$obsTonnages) ))
  cat( sprintf( "Number of materials: %d\n", length(object$matNames) ))
  cat( sprintf( "Materials: \n" ))
  for(i in 1:length(object$matNames)) {
    cat( sprintf( "%s    ", object$matNames[i] ))
  }

  cat( sprintf( "\n\n"))
  cat( sprintf( "Summary statistics for the observed material tonnages:\n" ))
  print(summary(object$obsTonnages[, -1]))
  cat( sprintf( "\n"))
  cat( sprintf( "Summary statistics for the pdf:\n" ))
  print(summary(object$rs))

  cat( sprintf( "\n\n"))
  cat( sprintf( "Mean vector:\n" ))
  tmp <- cbind(object$theObsMean, object$theMean)
  colnames(tmp) <- c("Observed", "Pdf")
  print(signif(tmp, digits = nDigits))

  cat( sprintf( "\n\n"))
  cat( sprintf( "Standard deviation vector:\n" ))
  tmp <- cbind(object$theObsSd, object$theSd)
  colnames(tmp) <- c("Observed", "Pdf")
  print(signif(tmp, digits = nDigits))

  cat(sprintf( "\n\n"))
  cat(sprintf("Composite correlation matrix\n" ))
  cat(sprintf("The upper triangle is the upper triangle of the\n"))
  cat(sprintf("correlation matrix for the observed tonnages.\n"))
  cat(sprintf("The lower triangle is the lower triangle of the\n" ))
  cat(sprintf("correlation matrix for the pdf.\n"))
  tmp <- object$theObsCor
  tmp[lower.tri(tmp)] <- object$theCor[lower.tri(object$theCor)]
  diag(tmp) <- NA
  print(signif(tmp, digits = nDigits))

}

#' @title Construct the pdf for the material tonnages in a single,
#' undiscovered deposit
#'
#' @description Construct the probability density function (pdf) for the
#' material tonnages in a single, undiscovered deposit within the permissive
#' tract. The pdf is not explicitly specified; instead it is implicitly
#' specified with the random samples that are generated from it.
#'
#' @param obsTonnages
#' Dataframe containing the observed material tonnages. (See details.)
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
#' truncated at the lowest and highest values of the observed material tonnages.
#'
#' @param minNDeposits
#' Minimum number of deposits with observed material tonnages.
#'
#' @param nRandomSamples
#' Number of random samples used to compute summary statistics and the
#' cumulative distribution function.
#'
#' @param seed
#' Seed for the random number generator.
#'
#' @details
#' Dataframe obsTonnages comprises the material tonnages for the known
#' deposits. Each row comprise the data for one deposit. Each column comprises
#' a particular aspect of the data. The first column always lists the names
#' of the deposits. The name of the first column is always "Name".
#' The second column and subsequent columns, if any, list
#' the tonnages for each type of material. For example, if there is only
#' one material (which would correspond to only one resource), then its
#' tonnages would be listed in the second column. The name of this column
#' is the name of the resource. For example, it might be "Sand".
#' If there are two materials
#' (which might correspond to ore and a resource), then the ore tonnages would
#' be listed in the second column and the resource tonnages in the third
#' column. Again, the names of the second and third columns are the names of the
#' materials. For example, they might be "Ore" and "Copper".
#' If there are three or more materials, then the format of the
#' dataframe is similar.
#'
#' All material tonnages in the dataframe must be greater than zero and
#' must not be missing.
#'
#' The minimum number of
#' deposits (that is, rows in the dataframe) should be greater than or
#' equal to 20.
#'
#' @return If the input arguments have an error, the R-value NULL is returned.
#' Otherwise, a list with the following components is returned.
#' @return \item{obsTonnages}{Input argument obsTonnages.}
#' @return \item{units}{Input argument units}
#' @return \item{pdfType}{Input argument pdfType.}
#' @return \item{isTruncated}{Input argument isTruncated.}
#' @return \item{matNames}{Names of the materials.}
#' @return \item{logTonnages}{Matrix with the natural logarithm
#' transformation of the observed material tonnages.}
#' @return \item{theObsMean}{Mean vector for the observed tonnages.}
#' @return \item{theObsCov}{Covariance matrix for the observed tonnages.}
#' @return \item{theObsSd}{Standard deviation vector for the observed tonnages.}
#' @return \item{theObsCor}{Correlation matrix for the observed tonnages.}
#' @return \item{rs}{Random samples of the pdf.}
#' @return \item{theMean}{Mean vector for the pdf.}
#' @return \item{theCov}{Covariance matrix for the pdf.}
#' @return \item{theSd}{Standard deviation vector for the pdf.}
#' @return \item{theCor}{Correlation matrix for the pdf.}
#' @return \item{call}{Function call.}
#'
#' @examples
#' pdf1 <- TonnagePdf1(ExampleTonnageData, "mt")
#'
#' @export
#'
TonnagePdf1 <- function(obsTonnages,
                        units,
                        pdfType = "empirical",
                        isTruncated = TRUE,
                        minNDeposits = 20,
                        nRandomSamples = 20000, seed = 7) {

  if(ncol(obsTonnages) < 2) {
    stop( sprintf( "Function TonnagePdf1\n" ),
          sprintf( "The number of columns in dataframe obsTonnage must be.\n" ),
          sprintf( "greater than or equal to 2.\n" ),
          sprintf( "But the actual number is 1.\n" ),
          call. = FALSE )

  }
  if(nrow(obsTonnages) < minNDeposits) {
    stop( sprintf( "Function TonnagePdf1\n" ),
          sprintf( "The number of rows in dataframe obsTonnage must be.\n" ),
          sprintf( "greater than or equal to %d.\n", minNDeposits ),
          sprintf( "But the actual number is %d.\n", nrow(obsTonnages) ),
          call. = FALSE )
  }

  for(j in 2:ncol(obsTonnages)){
    if(any(is.na(obsTonnages[, j]))) {
      stop( sprintf( "Function TonnagePdf1\n" ),
            sprintf( "Column %d of dataframe obsTonnage has\n", j ),
            sprintf( "one or more missing values\n"),
            call. = FALSE )
    }
    if(any(obsTonnages[, j] <= 0.0)) {
      stop( sprintf( "Function TonnagePdf1\n" ),
            sprintf( "Column %d of dataframe obsTonnage has\n", j ),
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

  theObsMean <- colMeans(obsTonnages[, -1])
  theObsCov <- cov(obsTonnages[, -1])
  theObsSd <- sqrt(diag(theObsCov))
  theObsCor <- cov2cor(theObsCov)

  rval <- list( obsTonnages = obsTonnages,
                units = units,
                pdfType = pdfType,
                isTruncated = isTruncated,
                matNames = colnames(obsTonnages)[-1],
                logTonnages = as.matrix(log(obsTonnages[, -1])),
                theObsMean = theObsMean,
                theObsCov = theObsCov,
                theObsSd = theObsSd,
                theObsCor = theObsCor)

  rval$rs <- getRandomSamples(rval, nRandomSamples, seed = seed)
  rval$theMean <- colMeans(rval$rs)
  rval$theCov <- cov(rval$rs)
  rval$theSd <- sqrt(diag(rval$theCov))
  rval$theCor <- cov2cor(rval$theCov)
  rval$call <- sys.call()

  class(rval) <- "TonnagePdf1"
  return(rval)
}