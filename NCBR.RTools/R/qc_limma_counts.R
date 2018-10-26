####################
#
#  Frederick National Laboratory
#  Susan Huse, susan.huse@@nih.gov
#  August 13, 2018
#
####################
#' Prints out raw, normalized, and differential expression data comparing two sets of samples for a particular gene
#'
#' For a given gene symbol and two vectors of sample names, 
#' prints out the raw data, the normalized data, the differential expression data, 
#' and test results of the model
#'  
#' @param symbolname the gene symbol serving as rownames in the data objects
#' @param A a vector of sample names used as column names in the data objects, to compare with B
#' @param B a second vector of sample names used as column names in the data objects, to compare with A
#' @param C contrast name to be used in calling topTable (coef="C")
#' @param theRaw the raw counts data matrix
#' @param theEset the eset of filtered and normalized expression data
#' @param theFit the limma eBayes fit object
#' @param testResults the result of decideTests
#' 
#' @return list object with [[1]] = numeric matrix of raw count data, 
#'                          [[2]] = dataframe with sample names and covariates
#'
#' @author Susan Huse \email{susan.huse@@nih.gov}
#' @keywords RNASeq Pipeliner limma
#'
#' @examples
#' eset <- make_filtered_eset(x=raw, y=covarMtx, g=groupsDF$Group, 
#'                            minc = filter_minCPM, mins = filter_minSampleCount, 
#'                            nmethod="quantile")
#' fit <- lmFit(eset, design)
#' cont.matrix.dTGFb <- makeContrasts(
#'                      WT_dTGFb = WNP-WTP,
#'                      KO_dTGFb = KNP-KTP,
#'                      Diff=(KNP-KTP)-(WNP-WTP), 
#'                      levels=design
#'                      )
#' fit.dTGFb <- contrasts.fit(fit, cont.matrix.dTGFb)
#' fit.dTGFb <- eBayes(fit.dTGFb)
#' results.dTGFb <- decideTests(fit.dTGFb, method="nestedF")
#' qc_counts(symbolname=i, 
#'           A = rownames(groupsDF[groupsDF$All == "WTP",]), 
#'           B = rownames(groupsDF[groupsDF$All == "WNP",]), 
#'           theRaw=raw, theEset=eset, theFit=fit.dTGFb, testresults=dTGFb, isSymbol=TRUE)
#' 
#' @export
qc_limma_counts <- function(symbolname, A, B, C, theRaw, theEset, theFit, testresults) {
  print("---------------------------------------------------------------")
  print(paste("Differential Expression Data for ", symbolname))
  print("---------------------------------------------------------------")
  
  print("Raw Expression Data")
  theRaw <- theRaw[symbolname,c(A,B)]
  print(theRaw)
  print(paste("A:", format(mean(theRaw[A]), digits=3), 
              ", B:", format(mean(theRaw[B]), digits=3), 
              "Diff:", format(mean(theRaw[A]) - mean(theRaw[B]), digits=3)))
  
  print("Normalized Expression Data")
  theEset <- theEset$E[symbolname,c(A,B)]
  print(theEset)
  print(paste("A:", format(mean(theEset[A]), digits=3), 
              ", B:", format(mean(theEset[B]), digits=3), 
              "Diff:", format(mean(theEset[A]) - mean(theEset[B]), digits=3)))
  
  if (! is.null(C)) {
    print(paste0("Differential Expression Data: topTable (", C, ")"))
    theTable <- topTable(theFit, coef=C, number=nrow(theFit$coefficients))
    print(theTable[symbolname,])
  }
  
  if(! is.null(testresults)) {
    print("Test Results (Diff)")
    theTable <- topTable(theFit, coef="Diff", number=nrow(theFit$coefficients))
    print(theTable[symbolname,])
    theTest <- testresults[symbolname,]
    print(theTest)
  }
}