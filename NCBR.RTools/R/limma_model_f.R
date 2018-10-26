####################
#
#  Frederick National Laboratory
#  Susan Huse, susan.huse@@nih.gov
#  August 9, 2018
#
####################
#' Run a limma model on an expressionset and formula
#'
#' Reads an expression set and formula
#'  and runs the model design, fits the model, and runs eBayes.  
#'  Returns the MArrayLM object ("fit")
#'  
#' @param x ExpressionSet (can be created with make_filtered_eset)
#' @param f formula
#' 
#' @return MArrayLM object ("fit")
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
#' Genotype <- covarMtx[,"Genotype"]
#' Treatment <- covarMtx[, "Treatment"]
#' myF <- as.formula("~ Genotype + Treatment")
#' myFit <- limma_model(x=eset, f=myF)
#' 
#' @export
limma_model_f <- function(x, f) {
  design <- model.matrix(f)
  fit <- lmFit(x, design)
  fit <- eBayes(fit)
  return(fit)
}
