####################
#
#  Frederick National Laboratory
#  Susan Huse, susan.huse@@nih.gov
#  October 16, 2018
#
####################
#' Create a Venn diagram matrix of significant genes
#'
#' Creates matrices of significant differentially expressed genes from a list of topTable results for drawing Venn diagrams
#'
#' Reads a list of differentially expressed gene comparisons as limma topTable objects,
#' and creates a unified matrix of differentially expressed genes across the list of contrasts.
#' The output is designed for use with the limma::vennDiagram function.
#' NB: the vennDiagram function is limited to 5 sets (contrasts)
#' 
#' The input list of topTable dataframes will be combined into a unified dataframe of fold change values.
#' NB: It is very important to FILTER FOR SIGNIFICANCE BEFORE running get_venn_matrix.
#' get_filtered_topTable has been designed to do this prefiltering.
#' 
#' This function will merge the logFC data from each of the input filtered, topTable dataframes
#' to create a single logFC matrix of all topTable results. 
#' Genes significant in one contrast but not another will have NA values.
#' 
#' The logFC matrix will be converted to a Venn diagram matrix with values of 1,0,-1, where
#' significant logFC > 1 = 1, 
#' significant logFC < 1 = -1, and
#' not significant logFC = 0.
#' This Venn matrix can then be used by limma::vennDiagram
#' 
#' Returns a list of two dataframes: [[1]] = Venn matrix, [[2]] = logFC matrix
#'  
#' @param dflist a list of topTables results, e.g., from get_filtered_topTable
#' @param theNames a vector of names for labeling each of the sets.
#' @param q minimum FDR / q / adjusted p-value to include (default=0.05)
#' @param n maximum number of rows to include in the final output dataframe (default = 10)
#' @param lFC minimum log2 fold change to include either up or down (abs(logFC)) (default = 1.5)
#' @param addGene parse out the gene name from the gene ID reported with CCBR Pipeliner (default = TRUE)
#' 
#' @return a list of two dataframes of the format provided by limma::topTable, the Venn counts and the logFC
#'
#' @author Susan Huse \email{susan.huse@@nih.gov}
#' @keywords RNASeq Pipeliner limma
#'
#' @examples
#' Fit1 <- limma_model(x=eset, y=covarMtx, cols=covars, interact=FALSE)
#' Set1 <- get_filtered_topTable(Fit1, theCoef="isKO", q=0.1, n=50, lFC=1)
#' Set2 <- get_filtered_topTable(Fit2, theCoef="isKO", q=0.1, n=50, lFC=1)
#' Set3 <- get_filtered_topTable(Fit3, theCoef="isKO", q=0.1, n=50, lFC=1)
#' myVenns = get_venn_matrix(list(
#'                get_filtered_topTable(Fit1, theCoef="isKO_TRUE", q=0.1, n=10000, lFC=1), 
#'                get_filtered_topTable(Fit2, theCoef="isKO_TRUE", q=0.1, n=10000, lFC=1), 
#'                get_filtered_topTable(Fit3, theCoef="isKO_TRUE", q=0.1, n=10000, lFC=1)), 
#'             theNames=c("Set1", "Set2", "Set3"))
#' vennDiagram(myVenns[[1]], main="Differential Expression across sets", 
#'             circle.col=c("red", "blue", "green"), include="both")
#' 
#' @export
get_venn_matrix <- function(dflist, theNames) {
  # Get all the sets of significant genes by set
  
  # initialize the merged dataframe using the first 2 dfs in the list
  df <- merge(dflist[[1]], dflist[[2]], all=TRUE, by="row.names")
  colnames(df)[which(colnames(df) == "logFC.x")] <- theNames[1]
  colnames(df)[which(colnames(df) == "logFC.y")] <- theNames[2]
  
  # add any additional dataframes in the dflist to the new dataframe
  if(length(dflist) > 2) {
    for(i in 3:length(dflist)) {
      df <- merge(df, dflist[[i]], all=TRUE, by.x="Row.names", by.y=0)
      colnames(df)[which(colnames(df) == "logFC")] <- theNames[i]
    }
  }
  
  # set to the rownames
  rownames(df) <- df$Row.names
  df <- df[, theNames]
  
  # Create a binary version of the dataframe
  df2 <- df
  df2[is.na(df2)] <- 0
  df2[df2 > 0] <- 1
  df2[df2 < 0] <- -1
  return(list(df2, df))
}
