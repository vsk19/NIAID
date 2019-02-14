####################
#
#  Frederick National Laboratory
#  Susan Huse, susan.huse@@nih.gov
#  October 16, 2018
#
####################
#' Filter topTable results by fold change, significance, and n
#'
#' Runs limma topTable on an input fit object, 
#' filters by minimum q-value, logFC, and returns up to n rows.
#'
#' Runs the limma topTable function, returning all rows, 
#' then filters out non-hits by minimum q-value (FDR), 
#' and a minimum, absolute fold change (abs(log2FC)).
#' The results are ordered by abs(logFC) and then the first "n" rows are returned.
#' 
#' The addGene option will parse out the gene name from the gene ID reported by the CCBR Pipeliner, 
#' and add the gene name as a new column to the output matrix.
#' 
#' Returns a topTable style dataframe
#'  
#' @param theFit a limma fit (MArrayLM) object
#' @param theCoef the coefficient in the fit object to analyze (default=NULL)
#' @param q minimum FDR / q / adjusted p-value to include (default=0.05)
#' @param n maximum number of rows to include in the final output dataframe. Use "all" to include all significant rows (default = 10)
#' @param lFC minimum log2 fold change to include either up or down (abs(logFC)) (default = 0.58 = FC of 1.5)
#' @param addGene parse out the gene name from the gene ID reported with CCBR Pipeliner (default = TRUE)
#' @param miRNA boolean, specifying if fit object is an miRNA DGELRT object (default=FALSE)
#' 
#' @return dataframe, with the format provided by limma::topTable
#'
#' @author Susan Huse \email{susan.huse@@nih.gov}
#' @keywords RNASeq Pipeliner limma
#'
#' @examples
#' myFit <- limma_model(x=eset, y=covarMtx, cols=covars, interact=FALSE)
#' myHits <- get_filtered_topTable(myFit, theCoef="isKO", q=0.1, n=50, lFC=1)
#' @export
get_filtered_topTable <- function(theFit, theCoef=NULL, q=0.05, n=10, lFC=0.58, addGene=TRUE, miRNA=FALSE) {
  if(miRNA) {
    myDF <- topTags(theFit, n=nrow(theFit$coefficients))$table
    addGene=FALSE
    myDF <- myDF[myDF$FDR <= q,]
  } else {
    myDF <- topTable(theFit, coef=theCoef, number=nrow(theFit$coefficients))
    myDF <- myDF[myDF$adj.P.Val <= q,]
  }
  
  myDF <- myDF[abs(myDF$logFC) >= lFC,]
  if(addGene) { myDF <- insert_EnsemblGene(myDF) }
  myDF <- myDF[order(abs(myDF$logFC), decreasing=TRUE),]
  if (n != "all") {
    if (nrow(myDF) > n) {myDF <- myDF[1:n,]}
  }
  return(myDF)
}
