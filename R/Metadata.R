#' @title Get the tract identifier
#'
#' @description Get the tract identifier, which is a unique character string
#' specifying the permissive tract.
#'
#' @param object
#' An object of class "Metadata"
#'
#' @return
#' A character string containing the tract identifier.
#'
#' @examples
#' meta <- Metadata("Lucky Strike PT", "PT001")
#' cat(sprintf("Tract identifier: %s\n", getTractId(meta)))
#'
#' @export
#'
getTractId <- function(object) {
  return(object$tractId)
}

#' @title Print the meta data for the probability calculations
#'
#' @description Print the meta data for the probability calculations for one
#' permissive tract.
#'
#' @param object
#' An object of class "Metadata"
#'
#' @examples
#' metaData <- Metadata("Lucky Strike PT", "PT001")
#' print(metaData)
#'
#' @export
#'
print.Metadata <- function(object) {

  cat(sprintf("Meta data for the permissive tract\n"))
  cat(sprintf("Name: %s\n", object$tractName))
  cat(sprintf("Identifier: %s\n", object$tractId))
  cat(sprintf("Area: %s\n", object$tractArea))
  cat(sprintf("Number of known deposits: %d\n", object$nKnownDeposits))
  cat(sprintf("Description: %s\n", object$description))

}

#' @title Construct the metadata for the probability calculations
#'
#' @description Construct the metadata for the probability calculations
#' for one permissive tract.
#'
#' @param tractName
#' Character string containing the name of the permissive tract.
#'
#' @param tractId
#' Character string containing the tract identifier, which is a unique
#' identifier for the permissive tract. (See Details).
#'
#' @param tractArea
#' Character string containing the area of the permissive tract and the units
#' for the area.
#'
#' @param nKnownDeposits
#' Integer containing the number of known deposits in the permissive tract.
#'
#' @param description
#' A description of the permissive tract.
#'
#' @details
#' Parameter \code{tractId}, which is a character string, may be used to file
#' names associated with the probability calculations. Consequently, it is
#' recommended that \code{tractId} be suitable for a file name. That is,
#' \code{tractId}
#' should contain neither spaces nor odd characters such as !, $, and &.
#'
#' @return A list with the following components is returned.
#' @return \item{tractName}{Input argument tractName}
#' @return \item{tractId}{Input argument tractId}
#' @return \item{tractArea}{Input argument tractArea}
#' @return \item{nKnownDeposits}{Input argument nKnownDeposits}
#' @return \item{description}{Input argument description}
#' @return \item{call}{Function call}
#'
#' @examples
#' meta <- Metadata("Lucky Strike PT", "PT001", "mt")
#'
#' @export
#'
Metadata <- function(tractName, tractId,
                     tractArea = "",
                     nKnownDeposits = 0,
                     description = "") {

  rval <- list( tractName = tractName,
                tractId = tractId,
                tractArea = tractArea,
                nKnownDeposits = nKnownDeposits,
                description = description,
                call=sys.call() )

  class(rval) <- "Metadata"
  return(rval)
}