####################
#
#  Frederick National Laboratory
#  Susan Huse, susan.huse@@nih.gov
#  November 9, 2018
#
####################
#' Check topTable RNASeq results against the expression set for QC
#'
#' Compares the differential expression values of two sets of samples (treatment and control) 
#' as reported in topTable against the underlying expression set to quality control the 
#' magnitude and direction of the computed fold change.  
#' 
#' QC checks: 
#'     1) does the sample names reported in trtExpr and ctrlExpr match the expected sample names for each set?
#'     2) does Eset_AveExpr equal Fit_AveExpr exactly?  They should be equal to 6 decimal places
#'     3) does the direction (+ or -) of Eset_logFC match the direction of Fit_logFC?
#'     4) does the magnitude of Eset_logFC match the magnitude of Fit_logFC?  
#'        These will be different but certainly less than 0.5.
#'     5) Does the variance level in the treatment and in the control samples warrant Fit_pVal reported?
#' 
#' Returns a list object with three dataframes:
#'     checkTable - the comparison table, with the first n rownames from topTable, 
#'                  Mean_Trt = average expression level of treatment samples
#'                  Mean_Ctrl = average expression level of control samples
#'                  Eset_AveExpr = average expression level of treatment and control samples
#'                  Fit_AveExpr = average expression level read from the topTable results
#'                  Eset_logFC = the predicted log2 fold change calculated directly from the expression data
#'                  Fit_logFC = the calculate log2 fold change for the contrast reported in topTable.
#'                  Fit_pVal = the p-value reported in topTable
#'                  
#'     trtExpr -    a dataframe of the expression levels of each of the treatment samples 
#'                  for the rows reported in checkTable
#'                  
#'     ctrlExpr -   a dataframe of the expression levels of each of the control samples 
#'                  for the rows reported in check Table
#' 
#' 
#' 
#' @param eset an expression set matrix used for creating the fit to be tested
#' @param fit the limma fit to be evaluated
#' @param coef the fit coefficient (column name from fit$coefficients) to evaluate
#' @param trt a vector of column names or indices for subsetting the treatment samples from eset
#' @param ctrl a vector of column names or indices for subsetting the control samples from eset
#' @param n the number of rows to pull from topTable and report on [default=10]
#' 
#' @return a list object of three dataframes: checkTable, trtExpr, ctrlExpr
#'
#' @author Susan Huse \email{susan.huse@@nih.gov}
#' @keywords RNASeq limma topTable check
#'
#' @examples
#' check_topTable(eset=theEset, fit=theFit, coef="Trt_Ctrl", trt=trtSamples, ctrl=ctrlSamples)
#' check_topTable(eset=theEset, fit=theFit, coef="KOvsWT",
#'                trt=theSamples[grepl("KO_", theSamples)], ctrl=theSamples[grepl("WT_", theSamples)])
#'
#' @export
check_topTable <- function(eset, fit, coef, trt, ctrl, n=10) {
  # run topTable
  topDF <- topTable(fit, coef=coef, p.value=1, number=n)
  # head2(topDF)
  
  # Subset the eset
  e_trt <- eset[rownames(topDF), trt]
  e_ctrl <- eset[rownames(topDF), ctrl]
  
  # calculate logFC from eset and create the for return
  df <- data.frame(Mean_Trt=numeric(), 
                   Mean_Ctrl=numeric(), 
                   Eset_AveExpr=numeric(),
                   Fit_AveExpr=numeric(),
                   Eset_logFC= numeric(), 
                   Fit_logFC = numeric(), 
                   Fit_pVal = numeric())
  for(row in rownames(topDF)) {
    # print(topDF[row,])
    mean_t <- (mean(e_trt[row,]))
    mean_c <- (mean(e_ctrl[row,]))
    mean_tc <- mean(eset[row, c(trt, ctrl)])
    lfc <- mean_t - mean_c
    df[row,] <- c(mean_t, mean_c, mean_tc, topDF[row, "AveExpr"], 
                  lfc, topDF[row, "logFC"], topDF[row, "P.Value"])
  }
  return(list(checkTable=df, trtExpr=e_trt, ctrlExpr=e_ctrl))
}
