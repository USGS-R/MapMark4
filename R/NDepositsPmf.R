

# The function argument(s) is (are) defined in the documentation for
# function NDepositsPmf.
#
#' @useDynLib ProMinerCore
#'
CalcMark3Pmf <- function( thresholds ){

  # These variables are used in fortran subroutine Mark3Pmf, so they are used here too.
  ND <- as.integer( c( 0, thresholds ) )
  sizeND <- length( ND )
  INOD <- length( thresholds )
  XX <- vector( mode="numeric", length=(ND[INOD+1]+1) )
  sizeXX <- length( XX )
  status <- 0

  if( INOD != 3 && INOD != 5 && INOD != 7 && INOD != 9 ) {
    stop( sprintf( "Function CalcMark3Pmf\n" ),
          sprintf( "The number of thresholds must be 3, 5, 7, or 9.\n" ),
          sprintf( "But the actual number is %d.\n", INOD ),
          call. = FALSE )
  }

  if( any( diff( ND ) < 0 ) ) {
    stop( sprintf( "Function CalcMark3Pmf\n" ),
          sprintf( "The specified thresholds must be in non-decreasing order.\n" ),
          call. = FALSE )
  }

  tmp <- .Fortran( "Mark3Pmf", as.integer(ND), as.integer(sizeND), as.integer(INOD),
                   as.double(XX), as.integer(sizeXX), as.integer(status) )

  if( tmp[[6]] != 0 ) {
    stop( sprintf( "Function CalcMark3Pmf\n" ),
          sprintf( "Fortran subroutine Mark3Pmf failed.\n" ),
          sprintf( "Status = %d\n", tmp[[6]] ),
          call. = FALSE )
  }

  specifiedAccdf <- (c( 0.90, 0.50, 0.10, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001 ))[1:INOD]

  probs <- tmp[[4]]
  nDeposits <- 0:(length(probs)-1)

  return( list( nDeposits=nDeposits,
                probs=probs,
                specifiedAccdf=specifiedAccdf ) )

}

# The function argument(s) is (are) defined in the documentation for
# function NDepositsPmf.
#
#' @useDynLib ProMinerCore
#'
CalcMark4Pmf <- function( thresholds, maxNumberOfDeposits ){

  # These variables are used in fortran subroutine Mark4Pmf, so they are used here too.
  ND <- as.integer( c( 0, thresholds ) )
  sizeND <- length( ND )
  INOD <- length( thresholds )
  maxd <- maxNumberOfDeposits
  XX <- vector( mode="numeric", length=(maxd+1) )
  sizeXX <- length( XX )
  status <- 0

  if( INOD != 3 && INOD != 5 && INOD != 7 && INOD != 9 ) {
    stop( sprintf( "Function CalcMark4Pmf\n" ),
          sprintf( "The number of thresholds must be 3, 5, 7, or 9.\n" ),
          sprintf( "But the actual number is %d.\n", INOD ),
          call. = FALSE )
  }

  if( any( diff( ND ) < 0 ) ) {
    stop( sprintf( "Function CalcMark4Pmf\n" ),
          sprintf( "The specified thresholds must be in non-decreasing order.\n" ),
          call. = FALSE )
  }

  if( maxd < tail( thresholds, n=1 ) ) {
    stop( sprintf( "Function CalcMark4Pmf\n" ),
          sprintf( "The maximum number of deposits must be >= %d\n", tail( thresholds, n=1 ) ),
          sprintf( "But the actual number is %d.\n", maxd ),
          call. = FALSE )
  }


  tmp <- .Fortran( "Mark4Pmf", as.integer(ND), as.integer(sizeND), as.integer(INOD), as.integer(maxd),
                   as.double(XX), as.integer(sizeXX), as.integer(status) )

  if( tmp[[7]] != 0 ) {
    stop( sprintf( "Function CalcMark4Pmf\n" ),
          sprintf( "Fortran subroutine Mark4Pmf failed.\n" ),
          sprintf( "Status = %d\n", tmp[[7]] ),
          call. = FALSE )
  }

  specifiedAccdf <- (c( 0.90, 0.50, 0.10, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001 ))[1:INOD]

  probs <- tmp[[5]]
  nDeposits <- 0:(length(probs)-1)

  return( list( nDeposits=nDeposits,
                probs=probs,
                specifiedAccdf=specifiedAccdf ) )

}


# The function argument(s) is (are) defined in the documentation for
# function NDepositsPmf.
#
CalcPoissonPmf <- function( theMean, relProbabilityThreshold=0.001 ) {

  if( theMean <= 0.0 ) {
    stop( sprintf( "Function CalcPoissonPmf\n" ),
          sprintf( "The expected number of deposits must be > 0.0\n" ),
          sprintf( "The actual value is %g)\n", theMean ),
          call. = FALSE )
  }

  if( relProbabilityThreshold <= 0.0 || 0.1 < relProbabilityThreshold ) {
    stop( sprintf( "Function CalcPoissonPmf\n" ),
          sprintf( "Argument relProbabilityThreshold should be within the interval (0.0,0.1]\n" ),
          sprintf( "The actual value is %g)\n", relProbabilityThreshold ),
          call. = FALSE )
  }

  # probability at the mode
  mode.prob <- dpois( floor(theMean), theMean )

  # Why did I define the tail probability in this manner? That is, why not
  # specify the tail probability directly (e.g., as 0.001). If the mean is
  # large, then all probabilities in the pmf are small. So, if the
  # tail probability is specified as 0.001, then the associated quantile might
  # be very close to the mode. Indeed, if the mean is large enough, then
  # a value of 0.001 might be greater than the probability of the mode. In
  # this case, nDepositRange probably (pun intended) could not be be computed.
  # However, these problems are avoided by using this definition.
  tail.prob <- relProbabilityThreshold * mode.prob

  # nDepositRange comprises two quantiles. The first quantile defines the left
  # tail of the pmf for which the probability is less than or equal to
  # tail.prob. Similarly, the second quantile defines the right tail of the
  # pmf for which the probability is less than or equal to tail.prob
  nDepositRange <- qpois( c( tail.prob, 1-tail.prob ), theMean )
  nDeposits <- nDepositRange[1]:nDepositRange[2]

  tmp <- dpois( nDeposits, theMean )
  probs <- tmp / sum( tmp )

  return( list( nDeposits=nDeposits, probs=probs ) )
}


# The function argument(s) is (are) defined in the documentation for
# function NDepositsPmf.
#
CalcUserSpecifiedPmf <- function( anchorPts, relProbabilities ){

  if(  length( anchorPts ) < 2 ) {
    stop( sprintf( "Function CalcUserSpecifiedPmf\n" ),
          sprintf( "The number of anchorPts must be >= 2.\n" ),
          sprintf( "But the actual number is %d.\n", length( anchorPts ) ),
          call. = FALSE )
  }

  if( length( anchorPts ) != length( relProbabilities ) ) {
    stop( sprintf( "Function CalcUserSpecifiedPmf\n" ),
          sprintf( "The length of anchorPts must equal the length of relProbabilities.\n" ),
          call. = FALSE )
  }

  if( any( anchorPts < 0 ) ) {
    stop( sprintf( "Function CalcUserSpecifiedPmf\n" ),
          sprintf( "All elements of anchorPts must be >= 0.\n" ),
          call. = FALSE )
  }

  if( any( diff( anchorPts ) <= 0 ) ) {
    stop( sprintf( "Function CalcUserSpecifiedPmf\n" ),
          sprintf( "All elements of anchorPts must be in strictly ascending order.\n" ),
          call. = FALSE )
  }

  if( any( relProbabilities < 0.0 ) ) {
    stop( sprintf( "Function CalcUserSpecifiedPmf\n" ),
          sprintf( "All elements of relProbabilities must be >= 0.\n" ),
          call. = FALSE )
  }

  nDeposits <- anchorPts[1]:anchorPts[length(anchorPts)]

  result <- approx( anchorPts, relProbabilities, nDeposits )

  probs <- result$y / sum( result$y )

  return( list( nDeposits=nDeposits, probs=probs ) )
}


# The function argument(s) is (are) defined in the documentation for
# function NDepositsPmf.
#
CalcNegBinomialPmf <- function( theMean, theStdDev,
                                relProbabilityThreshold=0.001 ){

  if( theMean <= 0.0 ) {
    stop( sprintf( "Function CalcNegBinomialPmf\n" ),
          sprintf( "The expected number of deposits must be > 0.0\n" ),
          sprintf( "The actual value is %g)\n", theMean ),
          call. = FALSE )
  }

  if( theMean >= theStdDev^2 ) {
    stop( sprintf( "Function CalcNegBinomialPmf\n" ),
          sprintf( "The expected number of deposits must be < the square of standard deviation.\n" ),
          call. = FALSE )
  }

  if( relProbabilityThreshold <= 0.0 || 0.1 < relProbabilityThreshold ) {
    stop( sprintf( "Function CalcNegBinomialPmf\n" ),
          sprintf( "Argument relProbabilityThreshold should be within the interval (0.0,0.1]\n" ),
          sprintf( "The actual value is %g)\n", relProbabilityThreshold ),
          call. = FALSE )
  }

  theProb = theMean / theStdDev^2
  theSize <- theMean * theProb / ( 1 - theProb )

  # The formula for the mode is taken from
  # https://en.wikipedia.org/wiki/Negative_binomial_distribution. Variable p
  # maps to 1-theProb.
  if( theSize <= 1 ) {
    theMode <- 0
  } else {
    theMode <- floor( (1-theProb)*(theSize-1)/theProb )
  }

  # probability at the mode
  mode.prob <- dnbinom( theMode, theSize, theProb )

  # The rationale for this definition of the tail probability is
  # given in function CalcPoissonPmf.
  tail.prob <- relProbabilityThreshold * mode.prob

  # The meaning of nDepositRange is given in function CalcPoissonPmf
  nDepositRange <- qnbinom( c( tail.prob, 1-tail.prob ), theSize, theProb )
  nDeposits <- nDepositRange[1]:nDepositRange[2]

  tmp <- dnbinom( nDeposits, theSize, theProb )
  probs <- tmp / sum( tmp )

  return( list( nDeposits=nDeposits, probs=probs ) )
}

# The function argument(s) is (are) defined in the documentation for
# function NDepositsPmf.
#
CalcDebugPmf <- function( nDeposits, relProbabilities ) {

  if( length( nDeposits ) != length( relProbabilities ) ) {
    stop( sprintf( "Function CalcDebugPmf\n" ),
          sprintf( "The length of nDeposits must equal the length of relProbabilities.\n" ),
          call. = FALSE )
  }

  if( any( nDeposits < 0 ) ) {
    stop( sprintf( "Function CalcDebugPmf\n" ),
          sprintf( "All elements of nDeposits must be >= 0.\n" ),
          call. = FALSE )
  }

  if( any( diff( nDeposits ) <= 0 ) ) {
    stop( sprintf( "Function CalcDebugPmf\n" ),
          sprintf( "All elements of nDeposits must be in strictly ascending order.\n" ),
          call. = FALSE )
  }

  if( any( relProbabilities < 0.0 ) ) {
    stop( sprintf( "Function CalcDebugPmf\n" ),
          sprintf( "All elements of relProbabilities must be >= 0.\n" ),
          call. = FALSE )
  }

  probs <- relProbabilities / sum( relProbabilities )

  return( list( nDeposits=nDeposits, probs=probs ) )

}

#' @title Plot the pmf for the number of undiscovered deposits
#'
#' @description Plots the probability mass function (pmf) for the number of
#' undiscovered deposits in the permissive tract. If the type is Mark3 or Mark4,
#' then the alternative complementary cumulative distribution function
#' is plotted too.
#'
#' @param object
#' An object of class "NDepositsPmf"
#'
#' @param isMeanPlotted
#' Logical variable indicating whether the mean is plotted on the pmf.
#'
#' @param barWidth
#' Width of the bars, which represent the probabilities for the numbers of
#' undiscovered deposits.
#'
#' @param isUsgsStyle
#' Make the plot format similar to the U.S. Geological Survey style
#'
#' @details
#' The alternative complementary cumulative distribution function is
#' G(x) = Pr( X >= x ) for - infinity < x < infinity. Between
#' the discrete values of X, the function G(x) is defined and would
#' be plotted as a horizontal lines. The problem is that these horizonal
#' lines will be confusing to people who are unfamilar with cumulative
#' distribution functions for discrete random variables. Consequently,
#' they are omitted.
#'
#' @references
#' Grimmett, G., and Stirzaker, D., 2001, Probability and random processes,
#' 3rd ed., Oxford University Press.
#'
#' @examples
#' pmf1 <- NDepositsPmf("Poisson", list(theMean=5), "")
#' plot(pmf1)
#'
#' pmf2 <- NDepositsPmf("Mark3", list(thresholds=c(1,7,20)), "")
#' plot(pmf2)
#'
#' @export
#'
plot.NDepositsPmf <- function( object, isMeanPlotted = TRUE,
                               barwidth = 0.1, isUsgsStyle = TRUE) {

  df <- data.frame(nDeposits = object$nDeposits,
                   probs = object$probs)

  p <- ggplot2::ggplot(df) +
    ggplot2::geom_bar(ggplot2::aes(x = nDeposits, y = probs),
                      stat = "identity", width = barwidth) +
    ggplot2::scale_x_continuous("Number of undiscovered deposits",
                                limits = range(df$nDeposits)+c(-1,1)) +
    ggplot2::ylab("Probability")

  if(isMeanPlotted == TRUE) {
    p <- p + ggplot2::geom_vline(xintercept = object$theMean, colour = "red")
  }

  if(isUsgsStyle)
    p <- p + ggplot2::theme_bw()

  if( object$type == "Mark3" || object$type == "Mark4" ) {
    df2 <- data.frame(nDeposits = object$nDeposits,
                     accdf = object$accdf)
    df3 <- data.frame(thresholds = object$pmf.args$thresholds,
                      specifiedAccdf = object$specifiedAccdf)

    q <- ggplot2::ggplot(df2) +
      ggplot2::geom_point(ggplot2::aes(x = thresholds, y = specifiedAccdf),
                          data = df3, colour = "red", size = 4) +
      ggplot2::geom_point(ggplot2::aes(x = nDeposits, y = accdf)) +
      ggplot2::ylim(0,1) +
      ggplot2::xlab("Number of undiscovered deposits") +
      ggplot2::ylab("Probability")

    if(isUsgsStyle)
      q <- q + ggplot2::theme_bw()

    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout=grid::grid.layout(1,2)))
    plot(p, vp=grid::viewport(layout.pos.row=1, layout.pos.col=1))
    plot(q, vp=grid::viewport(layout.pos.row=1, layout.pos.col=2))

  } else {
    plot(p)
  }

}


#' @title Summarize the pmf for the number of undiscovered deposits
#'
#' @description Summarize the probability mass function (pmf) for
#' the number of undiscovered deposits in the permissive tract.
#'
#' @param object
#' An object of class "NDepositsPmf"
#'
#' @examples
#' pmf <- NDepositsPmf( "Poisson", list(theMean=5), "" )
#' summary( pmf )
#'
#' @export
#'
summary.NDepositsPmf <- function( object) {

  cat(sprintf("Summary of the pmf for the number of undiscovered deposits\n"))
  cat(sprintf("within the permissive tract.\n"))
  cat(sprintf("------------------------------------------------------------\n"))
  cat( sprintf( "Type: %s\n", object$type ))
  cat( sprintf( "Description: %s\n", object$description ))
  cat( sprintf( "Mean: %g\n", object$theMean ))
  cat( sprintf( "Variance: %g\n", object$theVar ))
  cat( sprintf( "Standard deviation: %g\n", sqrt(object$theVar) ))
  cat( sprintf( "Coefficent of variation (\\%): %g\n",
                sqrt(object$theVar)*100/object$themean ))

  cat( sprintf( "Mode: %d\n", object$nDeposits[which.max(object$probs)] ))
  cat( sprintf( "Smallest number of deposits in the pmf: %d\n",
                min(object$nDeposits) ))
  cat( sprintf( "Largest number of deposits in the pmf: %d\n",
                max(object$nDeposits) ))
  cat(sprintf("\n\n\n\n"))

}

#' @title Get the pmf for the number of undiscovered deposits
#'
#' @description Get the probability mass function (pmf)
#' for the number of undiscovered deposits in the permissive tract.
#'
#' @param object
#' An object of class "NDepositsPmf"
#'
#' @return Dataframe containing the pmf. The dataframe has two columns: The
#' first lists the numbers of deposits. The second lists the probabilities
#' associated with those numbers.
#'
#' @examples
#' nDepositsPmf <- NDepositsPmf( "Poisson", list(theMean=5), "" )
#' getNDepositPmf( nDepositsPmf )
#'
#' @export
#'
getNDepositPmf <- function(object) {
  return(data.frame(nDeposits = object$nDeposits,
                    probs = object$probs))
}



#' @title Construct the pmf for the number of undiscovered deposits
#'
#' @description Construct the probability mass function (pmf) for the number of
#' undiscovered deposits within the permissive tract
#'
#' @param type
#' Character string with the type of pmf (See Details).
#'
#' @param pmf.args
#' List with the arguments characterizing the specified type (See Details).
#'
#' @param description
#' Character string with a short description of the pmf.
#'
#' @details
#' The type must be one of "Mark3", "Mark4",
#' "UserSpecified", "Poisson", "NegBinomial", or "Debug". Each type has
#' different arguments, which are specified in pmf.args. The arguments for
#' each type are described below.
#'
#' \emph{Mark3}
#'
#' List pmf.args has one element, which is named "thresholds". Thresholds
#' is a vector of integers. The meaning of each element in the vector
#' is described in this table. To make the table concise,
#' threshold[1] is represented by t[1], and so on. The size of vector thresholds
#' can be only 3, 5, 7 or 9.
#'
#' \tabular{lllll}{
#' Element \tab Size = 3             \tab Size = 5              \tab Size = 7               \tab Size = 9                \cr
#' t[1]    \tab P( N >= t[1] ) = 0.9 \tab P( N >= t[1] ) = 0.9  \tab P( N >= t[1] ) = 0.9   \tab P( N >= t[1] ) = 0.9    \cr
#' t[2]    \tab P( N >= t[2] ) = 0.5 \tab P( N >= t[2] ) = 0.5  \tab P( N >= t[2] ) = 0.5   \tab P( N >= t[2] ) = 0.5    \cr
#' t[3]    \tab P( N >= t[3] ) = 0.1 \tab P( N >= t[3] ) = 0.1  \tab P( N >= t[3] ) = 0.1   \tab P( N >= t[3] ) = 0.1    \cr
#' t[4]    \tab ---                  \tab P( N >= t[4] ) = 0.05 \tab P( N >= t[4] ) = 0.05  \tab P( N >= t[4] ) = 0.05   \cr
#' t[5]    \tab ---                  \tab P( N >= t[5] ) = 0.01 \tab P( N >= t[5] ) = 0.01  \tab P( N >= t[5] ) = 0.01   \cr
#' t[6]    \tab ---                  \tab ---                   \tab P( N >= t[6] ) = 0.005 \tab P( N >= t[6] ) = 0.005  \cr
#' t[7]    \tab ---                  \tab ---                   \tab P( N >= t[7] ) = 0.001 \tab P( N >= t[7] ) = 0.001  \cr
#' t[8]    \tab ---                  \tab ---                   \tab ---                    \tab P( N >= t[8] ) = 0.0005 \cr
#' t[9]    \tab ---                  \tab ---                   \tab ---                    \tab P( N >= t[9] ) = 0.0001
#' }
#'
#' The expression P( N >= t[1] ) = 0.9 means that the probability that the
#' number of the deposits (N)
#' will be at least t[1] is 0.9. That is, there is a 0.9 probability
#' (90 percent chance) of finding t[1] or more
#' deposits. The elements of t (thresholds) must be nondecreasing. For example,
#' if the size is 3, then t[1] <= t[2] <= t[3]. Although it seems that the
#' elements of t (thresholds) should
#' be strictly increasing, this is not required by the algorithm implemented
#' in Fortran subroutine Mark3Pmf
#'
#' If t[1] is 0, then, according to the previous definition, P( N >= 0 ) = 0.9.
#' Of course,
#' this is wrong because P( N >= 0 ) = 1. To address this problem, the
#' algorithm is implemented such that P( N >= 1 ) < 0.9.
#'
#' The p quantile of the cumulative distribution function F is defined as
#' the smallest x such that F(x) >= p. (DeGroot and Schverish, 2002, p. 115).
#' That is, the p quantile is the smallest x such that P( X <= x ) >= p where
#' X is a random variable.
#' This definition obviously differs from that for thresholds. Hence, quantiles
#' are not thresholds.
#'
#' The probabilities in each column of the above table appear to constitute
#' a complementary cumulative
#' distribution function (ccdf). However, a ccdf G is defined as
#' G(x) = P( X > x). (See the wikipedia reference.)
#' Because the table lists P( X >= x ), the probabilities in each column do
#' not constitute a ccdf. To
#' prevent any confusion, these probabilites are called "alternative
#' complementary cumulative distribution function,"
#' which will be abbreviated "accdf".
#'
#' Fortran subroutine Mark3Pmf calculates a pmf such that its aacdf
#' approximates (but rarely equals) the aacdf that is specified with
#' the thresholds.
#' This subroutine was extracted from program Mark3B (Root and others, 1992;
#' Root and others, 1998) and then compiled as a direct link library.
#'
#' \emph{Mark4}
#'
#' List pmf.args has two elements, which are named "thresholds" and
#' "maxNumberOfDeposits". Thresholds was described in the previous section.
#' Variable maxNumberOfDeposits is an integer; it must be greater than or
#' equal to the last element in thresholds.
#'
#' Fortran subroutine Mark4Pmf calculates a pmf such that its aacdf
#' equals exactly the aacdf that is specified with the thresholds.
#' The algorithm was developed and implemented by Jeffery D. Phillips.
#'
#' \emph{UserSpecified}
#'
#' List pmf.args has two elements, which are named "anchorPts" and
#' "relProbabilities". anchorPts is an integer vector, and relProbabilities is
#' a real vector with the same size as anchorPts.
#'
#' A common way to use this function is to have three elements in vector
#' anchorPts and
#' three elements in vector relProbabilities. Regarding vector anchorPts,
#' the first element
#' is the minimum number of deposits in the tract, the second is the mode,
#' and the third is the maximum
#' number of deposits in the tract. The first element must be greater than
#' or equal to 0, and all elements
#' must be in strictly ascending order. That is,
#' 0 <= anchorPts[1] < anchorPts[2] < anchorPts[3].
#' Regarding vector relProbabilities, the first, second, and third elements
#' are the relative probabilities
#' associated with anchorPts[1], anchorPts[2], and anchorPts[3]. The relative
#' probability is any postive,
#' real-valued number. Because anchorPts[2] is the mode, relProbabilities[2]
#' must be
#' greater than relProbabilities[1] and relProbabilities[3]. The sum of the
#' elements in relProbabilities
#' is not required to be 1.
#'
#' The relative probabilites between anchorPts[1] and anchorPts[3] are
#' calculated by linearly interpolating the
#' elements of relProbabilities. Then these interpolated, relative
#' probabilities are scaled, so that their sum is 1.
#'
#' This function is very flexibile: The minimum number of elements in vector
#' anchorPts is 2; there is no maximum.
#' (Vector relProbabilities and vector anchorPts must have the same length.)
#' If this flexibility is exploited, then
#' the user should carefully check the calculated probability mass function.
#'
#' \emph{Poisson}
#'
#' List pmf.args has one element, which is named "theMean". theMean is
#' real-valued and is the mean of the Poisson pmf.
#'
#' For the Poisson pmf, the random variable that represents the number of
#' undiscovered deposits extends from 0 to infinity. Of course, the Monte Carlo
#' simulation cannot be performed for an infinite range. When the number of
#' deposits is far from the mean, the probabilities are so small
#' that they have no practical effect on the simulation results.
#' The range associated with these small probabilites is delimited
#' by a lower bound and an upper bound, which are calulated from
#' variable relProbabilityThreshold. To understand the meaning of this variable,
#' assume that it is 0.001. The lower bound is the quantile for which
#' the probability in the left tail is 0.001 times the probability of the
#' mode of the pmf. Similarly, the right bound is the quantile for which
#' the probability in the right tail is 0.001 times the probability of the
#' mode of the pmf.
#' The pmf is truncated at the lower and the upper bounds, and then rescaled
#' so that the sum of the probabilities is 1. This rescaling changes the
#' mean of the pmf slightly.
#' Usually, the mean of the pmf is small, so that the lower bound is 0.
#' In this (common) case, the lower bound does not affect
#' the truncation of the Poisson pmf.
#'
#' Variable relProbabilityThreshold has a default value, but it can be changed
#' by specifyting the new value in list pmf.args.
#'
#' \emph{NegBinomial}
#'
#' List pmf.args has two elements, which are named "theMean" and "theStdDev".
#' Both theMean and theStdDev are real-valued. theMean and theStdDev are
#' the mean and the standard deviation of the negative binomial pmf. theMean
#' must be less than the square of variable theStdDev.
#'
#' For the negative binomial pmf, the random variable that represents the
#' number of
#' undiscovered deposits extends from 0 to infinity. Of course, the Monte Carlo
#' simulation cannot be performed for an infinity range. The procedure that
#' solves for this problem is described in section for the Poisson pmf.
#'
#' \emph{Debug}
#'
#' List pmf.args has two elements, which are named "nDeposits" and
#' "relProbabilities". nDeposits is an integer vector, and relProbabilities is
#' a real vector with the same size as anchorPts.
#'
#' Regarding vector nDeposits, all elements must be in strictly ascending order.
#' That is,
#' 0 <= nDeposits[1] < nDeposits[2] < nDeposits[3] and so on. Regarding
#' vector relProbabilities, the elements are postive,
#' real-valued numbers. The sum of the elements in relProbabilities is
#' not required to be 1. Vector nDeposits
#' and vector relProbabilites must have the same length.
#'
#' Here are three examples: If nDeposits = c(1) and relProbabilities=c(1),
#' then the probability mass function will be 1
#' when nDeposits is 1. If nDeposits = c(10) and relProbabilities=c(1),
#' then the probability mass function will be 1
#' when nDeposits is 10. If nDeposits = c(0,1) and relProbabilities=c(1.0,0.5),
#' then the probability mass function will
#' be 0.667 and 0.333 when nDeposits is 0 and 1, respectively.
#'
#' @return If the input arguments have an error, the R-value NULL is returned.
#' Otherwise, a list with the following components is returned.
#' @return \item{Type}{Input argument type}
#' @return \item{pmf.args}{Input argument pmf.args}
#' @return \item{description}{Input argument description}
#' @return \item{call}{Function call}
#' @return \item{probs}{Vector containing just the non-zero probabilities
#' in the pmf.}
#' @return \item{nDeposits}{Vector containing the number of undiscovered
#' deposits that are associated with the non-zero probabilities in the pmf.
#' The size of vector nDeposits is equal to the size of vector probs.}
#' @return \item{specifiedAccdf}{If the type is "Mark3" or "Mark4", then
#' this entry is a vector containing probabilities from the alternative
#' complementary cumulative distribution function. These probabilites are
#' associated with the thresholds, which are specified in pmf.arg. Otherwise,
#' this entry is NULL.}
#' @return \item{theMean}{The expected value (mean) of the number of
#' undiscovered deposits within the permissive tract. The number might differ
#' slightly from the input specification.}
#' @return \item{theVar}{The variance of the number of undiscovered deposits
#' within the permissive tract.}
#' @return \item{accdf}{The alternative complementary cumulative distribution
#' function for the pmf. Although it is calculated for all types of pmfs,
#' it is useful only for the "Mark3" and "Mark4" types.}
#' @return \item{entropy}{The information entropy, which is calculated with
#' the natural logarithm.}
#'
#' @references
#' DeGroot, M.H., and Schervish, M.J., 2002, Probability and statistics: Addison-Wesley.
#'
#' Root, D.H., Menzie, W.D., and Scott, W.A., 1992, Computer Monte Carlo simulation in quantitative
#' resource estimation: Nonrenewable resources, v. 1, no. 2, p. 125-138.
#'
#' Root, D.H., Scott, W.A., Jr., Schruben, P.G., 1998, MARK3B Resource assessment program for Macintosh:
#' U.S. Geological Survey, Open-file report 98-356.
#'
#' http://en.wikipedia.org/wiki/Cumulative_distribution_function
#'
#'
#' @examples
#' # Mark3 pmf
#' nDepositsPmf1 <- NDepositsPmf( "Mark3", list(thresholds=c(1,7,20)),
#' "Test Mark3" )
#' plot( nDepositsPmf1 )
#'
#' # Mark4 pmf
#' nDepositsPmf2 <- NDepositsPmf( "Mark4",
#' list(thresholds=c(1,7,20),maxNumberOfDeposits=23), "Test Mark4" )
#' plot( nDepositsPmf2 )
#'
#' # User specified pmf
#' nDepositsPmf3 <- NDepositsPmf( "UserSpecified",
#' list(anchorPts=c(1,4,10),relProbabilities=c(1,4,1)),
#' "Test UserSpecified" )
#' plot( nDepositsPmf3 )
#'
#' # Poisson pmf
#' nDepositsPmf4 <- NDepositsPmf( "Poisson", list(theMean=5), "Test Poisson" )
#' plot( nDepositsPmf4 )
#'
#' # Negative binomial pmf
#' nDepositsPmf5 <- NDepositsPmf( "NegBinomial", list(theMean=5,theStdDev=4),
#' "Test Negative binomial" )
#' plot( nDepositsPmf5 )
#'
#' # Debug pmf
#' nDepositsPmf6 <- NDepositsPmf( "Debug", list(nDeposits=1,relProbabilities=1),
#' "Test Debug 1" )
#' plot( nDepositsPmf6 )
#'
#' # Debug pmf
#' nDepositsPmf7 <- NDepositsPmf( "Debug",
#' list(nDeposits=c(2,3),relProbabilities=c(0.3,0.7)), "Test Debug 2" )
#' plot( nDepositsPmf7 )
#'
#' @export
#'
NDepositsPmf <- function( type, pmf.args, description="" ) {


  # Remove any parts of the pmf for which the number of draws will be 0. These
  # parts are associated with very small probabilities.
  AdjustPmf <- function( probs, nDeposits, nTotalDraws )
  {

    nDraws <- round( probs * nTotalDraws )
    indices <- which( nDraws == 0, arr.ind=TRUE )
    if( length( indices ) > 0 ) {
      cat( sprintf( "WARNING: Function AdjustPmf\n" ) )
      cat( sprintf( "Some probabilities are so small that the number of draws are 0.\n" ) )
      cat( sprintf( "These probabilities are omitted from the simulation.\n" ) )
      cat( sprintf( "The problems are\n" ) )
      cat( sprintf( "Probability    nDraws    nDeposits\n" ) )
      for( i in 1:length(indices) ) {
        cat( sprintf( "%g       %g        %g\n",
                      probs[indices[i]], nDraws[indices[i]],
                      nDeposits[indices[i]] ) )
      }
      nDeposits <- nDeposits[-indices]
      probs <- probs[-indices]
      probs <- probs / sum( probs )
    }

    return( list( probs=probs, nDeposits=nDeposits) )

  }

  CalcSummaryStats <- function( probs, nDeposits, base=exp(1) )
  {

    theMean <- sum( probs * nDeposits )
    theVar <- sum( probs * nDeposits^2 ) - theMean^2
    tmp <- c( 1, 1-cumsum( probs ) )
    theAccdf <- tmp[-length(tmp)]  # exclude the last element, which is 0
    theEntropy <- - sum( probs * log( probs, base=base ) )

    return( list( theMean=theMean, theVar=theVar, theAccdf=theAccdf,
                  theEntropy=theEntropy ) )

  }


  pmf <- switch( type,
                 Mark3 = do.call( CalcMark3Pmf, pmf.args ),
                 Mark4 = do.call( CalcMark4Pmf, pmf.args ),
                 UserSpecified = do.call( CalcUserSpecifiedPmf, pmf.args ),
                 NegBinomial = do.call( CalcNegBinomialPmf, pmf.args ),
                 Poisson = do.call( CalcPoissonPmf, pmf.args ),
                 Debug = do.call( CalcDebugPmf, pmf.args ) )


  # Obviously, it must be greater than 0. Greater than 5000 is preferable.
  nTotalDraws <- 20000

  adjusted.pmf <- AdjustPmf( pmf$probs, pmf$nDeposits, nTotalDraws )

  stats <- CalcSummaryStats( adjusted.pmf$probs, adjusted.pmf$nDeposits )

  # If the type is neither Mark3 nor Mark4, then specifiedAccdf is
  # set to NULL.
  rval <- list( type=type,
                pmf.args=pmf.args,
                description=description,
                call=sys.call(),
                probs=adjusted.pmf$probs,
                nDeposits=adjusted.pmf$nDeposits,
                specifiedAccdf=pmf$specifiedAccdf,
                theMean=stats$theMean,
                theVar=stats$theVar,
                accdf=stats$theAccdf,
                entropy=stats$theEntropy )

  class(rval) <- "NDepositsPmf"
  return(rval)
}
