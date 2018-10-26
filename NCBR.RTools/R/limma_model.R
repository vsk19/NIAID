####################
#
#  Frederick National Laboratory
#  Susan Huse, susan.huse@@nih.gov
#  August 9, 2018
#
####################
#' Run a limma model on an expressionset and covariates matrix
#'
#' Reads an expression set, covariate matrix, covariate names
#'  and runs the model design (with and without interaction), 
#'  fits the model, and runs eBayes.  
#'  Returns the MArrayLM object ("fit")
#'
#' @param x ExpressionSet (can be created with make_filtered_eset)
#' @param y covariate matrix (can be created with factordf_to_mtx)
#' @param cols character vector of column names to use from the covariate matrix
#' @param interact boolean value for whether to include interaction terms between covariates (default=FALSE)
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
#' myFit <- limma_model(x=eset, y=covarMtx, cols=covars, interact=FALSE)
#' 
#' @export
limma_model <- function(x, y, cols, interact=FALSE) {
  if(interact) {joiner=" * "} else {joiner=" + "}
  
  # Create the formula with covariates for the model
  for(varname in cols) {
      assign(varname, y[,varname])
    }
  f <- as.formula(paste("~ ", paste(cols, collapse=joiner), collapse=""))

  print(paste("Creating design model from formula: ", paste(f, collapse=" ")))
  
  design <- model.matrix(f)
  fit <- lmFit(x, design)
  fit <- eBayes(fit)
  return(fit)
}
