####################
#
#  Frederick National Laboratory
#  Susan Huse, susan.huse@@nih.gov
#  August 9, 2018
#
####################
#' Creates a filtered expression set from raw counts for limma modeling
#'
#' Filters a raw counts matrix using edgeR::filterByExpr, 
#'  normalizes the filtered data using limma::voom, and 
#' 
#' Returns a filtered, normalized, expression set.
#'  
#' @param x raw counts matrix (can be created with load_RawCountsGroup)
#' @param d design - passed directly to edgeR::filterByExpr (use NULL if using group)
#' @param g group - passed directly to edgeR::filterByExpr (default=NULL), if using group set d=NULL
#' @param minc integer value for minimum CPM count filtering threshold (default = 0.5)
#' @param mins integer value for minimum number of samples passing minc to be included in filtered set (default = 3)
#' @param nmethod normalization method to be used by voom (default = "quantile")
#' 
#' @return Elist (limma): expression matrix, weights, design, targets
#'
#' @author Susan Huse \email{susan.huse@@nih.gov}
#' @keywords RNASeq Pipeliner limma
#'
#' @examples
#' dfList <- load_RawCountsGroup(counts_file = file.path(dir_data, file_data), 
#'                               groups_file = file.path(dir_data, file_groups))
#' raw <- dfList[[1]]
#' groupsDF <- dfList[[2]]
#' rm(dfList)
#' names(covar_levels) <- covars
#' covarMtx <- factordf_to_mtx(groupsDF[, covars], cnames=covars, clevels=covar_levels)
#' eset <- make_filtered_eset(x=raw, y=covarMtx, g=groupsDF$Group, 
#'                            minc = filter_minCPM, mins = filter_minSampleCount, 
#'                            nmethod="quantile")
#' myFit <- limma_model(x=eset, y=covarMtx, cols=covars, interact=FALSE)
#' 
#' @export
make_filtered_eset <- function(x, d, g=NULL, minc=0.5, mins=3, nmethod="quantile", plot=TRUE) {
  # Filter low expression genes using edgeR
  dge <- DGEList(counts=x)
  keep <- filterByExpr(dge, design=d, min.count=minc, min.total.count=(minc * mins))
  dge <- dge[keep,,keep.lib.sizes=FALSE]
  
  # Normalize and convert to expressionset
  exprSet <- voom(dge, design=d, normalize.method="quantile", plot=plot)
  return(exprSet)
}
