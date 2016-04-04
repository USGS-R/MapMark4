#' @title Summarize the probability calculations
#'
#' @description Summarize the probability calculations for the
#' mineral resource assessment of the permissive tract.
#'
#' @param oMetadata
#' An object of class "Metadata"
#'
#' @param oPmf
#' An object of class "NDepositsPmf"
#'
#' @param oPdf1
#' An object of class "TonnagePdf1"
#'
#' @param oPdfPT
#' An object of class "TonnagePdfPT"
#'
#' @details
#' This function calls the summary functions for the respective S3 classes.
#'
#' @examples
#' meta <- Metadata("Lucky Strike PT", "PT001")
#' pmf <- NDepositsPmf( "NegBinomial", list(theMean=5,theStdDev=4), "" )
#' pdf1 <- TonnagePdf1(ExampleTonnageData, "mt", pdfType = "normal")
#' pdfPT <- TonnagePdfPT(pmf, pdf1)
#' summary(meta, pmf, pdf1, pdfPT)
#'
#' @export
#'
SummaryProbCalc <- function(oMetaData, oPmf, oPdf, oPdfPT) {

  summary.Metadata(oMetaData)

  summary.NDepositsPmf(oPmf)

  summary.TonnagePdf1(oPdf)

  summary.TonnagePdfPT(oPdfPT)

}
