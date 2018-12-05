####################
#
#  Frederick National Laboratory
#  Susan Huse, susan.huse@@nih.gov
#  October 23, 2018
#
####################
#' Merge limma:topTable results for multiple comparisons
#'
#' Reads a list of limma fit (MArrayLM) objects, runs topTable with the specified parameters
#' and then merges all of the tables together by geneID
#'
#' The list of fits, must have a matching list of names for columns (the comparisons), and 
#' a list of coefficient column names to provide to topTable.
#' The resulting matrix of genes from all comparisons will then be filtered by adj.P.Val (q), 
#' sorted, and then the top n rows returned.
#' The output column names will be designated by the fit name, "_", and the topTable column name, 
#' e.g., myFitName_logFC, myFitName_P.Value
#'
#' Filtering can be limited by selecting q=1.0 to accept all genes, irrespective of significance.
#' The number of rows returned can be set, or to return all rows set n=NULL

#' Before selecting n rows and before returning the dataframe, the data will be sorted.
#' Default sorting is abs(logFC) in descending order.  
#' To sort on a different column, use the standard topTable output column names (e.g., "adj.P.Val" and "P.Value"
#' 
#' @param fits a list of limma fit (MArrayLM) objects
#' @param names a vector of names matching the fit objects to be appended with the topTable column names
#' @param coefs a vector of coefficient names matching the fit objects to be used by topTable
#' @param q maximum FDR / q / adjusted p-value to include (default=0.05)
#' @param logFC minimum logFC to include (default=1)
#' @param n maximum number of rows to include (default=10), use NULL to include all filtered rows.
#' @param sortby ordering of the rows exported. If "logFC" then sort by descreasing abs(logFC), otherwise sort by the specified topTable column.
#' @param miRNA boolean, specifying if fit objects are miRNA DGELRT objects (default=FALSE)
#'
#' @return dataframe merged from all input topTable(fit, coef) results
#'
#' @author Susan Huse \email{susan.huse@@nih.gov}
#' @keywords RNASeq limma topTable merge
#'
#' @examples
#' allfits <- list()
#' for(i in names(DE3_list)){allfits[[i]] <- DE3_list[[i]]$fit}
#' for(i in names(DE4_list)){allfits[[i]] <- DE4_list[[i]]$fit}
#' allDF <- merge_topTables(allfits, 
#'                          names=c(names(DE3_list), names(DE4_list)),
#'                          coefs=c(DE3_coefs, DE4_coefs), 
#'                          q=0.10, n=NULL, sortby="logFC")
#'
#' @export
merge_topTables <- function(fits, names, coefs, q=0.05, logFC=1, n=10, sortby="logFC", miRNA=FALSE) {
  # Check the inputs
  if(length(fits) != length(names)) {stop("fits list and names vector must be the same length")}
  if(length(fits) < 2) {stop("Need at least 2 fits to merge")}
  topHits <- c()
  allDFs <- list()
  
  # Get all of the top tables
  for(i in 1:length(fits)) {
    if(miRNA) {
      df <- topTags(fits[[i]], n=nrow(fits[[i]]$coefficients))$table
      colnames(df)[grep("FDR", colnames(df))] <- "adj.P.Val"
    } else {
      df <- topTable(fits[[i]], coef=coefs[i], n=nrow(fits[[i]]$coefficients))
    }

    colnames(df) <- paste(names[i], colnames(df), sep="_")
    allDFs[[length(allDFs) + 1 ]] <- df
    # create a running list of the q significant results
    topHits <- c(topHits, rownames(df[(df[,paste(names[i], "adj.P.Val", sep="_")] <= q) &
                                        (abs(df[,paste(names[i], "logFC", sep="_")]) >= logFC), ]))
  }
  
  # Get the list of hits in common and create a combined DF
  theHits <- unique(topHits)
  
  newDF <- cbind(allDFs[[1]][topHits,], allDFs[[2]][topHits,])
  for(i in 3:length(fits)) {
    newDF <- cbind(newDF, allDFs[[i]][topHits,])
  }
  
  # Sort the best to the top
  cols <- paste(names, sortby, sep="_")
  if(sortby == "logFC") {
    newDF <- newDF[order(rowSums(abs(newDF[,cols])), decreasing=TRUE),]
  } else {
    newDF <- newDF[order(rowSums(newDF[,cols])),]
  }
  
  # Take top n
  if(! is.null(n)) {
    if(n < nrow(newDF)) {
      newDF <- newDF[1:n,]
    }
  }
  # if((! is.null(n)) & (n < nrow(newDF))) {newDF <- newDF[1:n,]}
  
  return(newDF)
}

