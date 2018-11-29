####################
#
#  Frederick National Laboratory
#  Susan Huse, susan.huse@@nih.gov
#  October 22, 2018
#
####################
#' Print limma::topTable results to Excel
#'
#' Reads a limma fit (MArrayLM) object, runs topTable with the specified parameters
#' the exports the dataframe to xlsx or csv format.
#'
#' For the given fit object, as well as the coefficient to analyze, 
#' topTable2Excel will run topTable returning all results, before filtering by
#' the maximum FDR (q, adj.P.Val), a minimum logFC, and a maximum number of rows.
#' If more rows match the FDR and the logFC filtering criteria than the number of rows requested
#' the data will be sorted by absolute value of logFC in reverse order before selecting the top n lines.
#' 
#' Data can be exported to either .xlsx or to .csv format
#' Data can be appended to an Excel spreadsheet by creating a new sheet (tab), but data cannot be appended to an existing csv file.
#' 
#' To export all data, set q=1.0
#' 
#' If you run into java out of memory errors:
#' ** Error in .jcall("RJavaTools", "Ljava/lang/Object;", "invokeMethod", cl,  : 
#' java.lang.OutOfMemoryError: GC overhead limit exceeded **
#' you can increase the java memory:
#' options(java.paramters = "-Xmx8000m")
#' 
#' 
#' @param theFit a limma fit (MArrayLM) object
#' @param theCoef the coefficient in the fit object to analyze (default=NULL)
#' @param theFile the name of the output file for the exported data
#' @param theSheet the name of the sheet (tab) in the output Excel file for the data
#' @param q maximum FDR / q / adjusted p-value to include (default=0.05)
#' @param sortby ordering of the rows exported.  If "logFC" then sort by descreasing abs(logFC), otherwise sort by FDR.
#' @param append should the data be added as a new sheet to an existing file, or should the output file be reinitialized with this as the only sheet (default=TRUE)
#' @param excel export to xlsx format, otherwise export to csv format (default = TRUE)
#' @param sep if the file format is excel (excel=FALSE), determines what delimiter to use (default=",").  Use "\t" for tab.
#' 
#' @return Excel (xlsx or csv) file is saved to disk.
#'
#' @author Susan Huse \email{susan.huse@@nih.gov}
#' @keywords RNASeq limma topTable excel xlsx csv
#'
#' @examples
#' options(java.parameters = "-Xmx8000m")
#' top2Excel(theFit=myFit, theCoef=NULL, theFile=Project_results.xlsx, theSheet="MyProject", append=FALSE)
#'
#' @export
topTable2Excel <- function(theFit, theCoef, theFile, theSheet="topTable", q=0.05, sortby="logFC", append=TRUE, excel=TRUE, sep=",") {
  if(append == TRUE & excel == FALSE) {
    return("Error in topTable2Excel: cannot append data to an existing csv file.  Append is only available for Excel (xlsx) files.  Exiting...")
  }
  myDF <- topTable(theFit, coef=theCoef, number=nrow(theFit$coefficients))
  myDF <- myDF[myDF$adj.P.Val <= q,]
  myDF <- insert_EnsemblGene(myDF)
  if(sortby == "logFC") {
    myDF <- myDF[order(abs(myDF$logFC), decreasing=TRUE),]
  } else if (! is.null(sortby)) {
    myDF <- myDF[order(myDF$adj.P.Val, decreasing=TRUE),]
  }
  # print(dim(myDF))
  if(excel){
    write.xlsx2(myDF, file=theFile, sheetName=theSheet, append=append, col.names=TRUE, row.names=TRUE)
  } else {
    write.table(myDF, file=theFile, quote=FALSE, col.names=NA, row.names=TRUE, sep=sep)
  }
}
