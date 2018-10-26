####################
#
#  Frederick National Laboratory
#  Susan Huse, susan.huse@@nih.gov
#  August 9, 2018
#
####################
#' Convert a dataframe with factors to a numeric matrix
#'
#' Reads in a dataframe with along with factor column names and levels
#'   and returns a numeric matrix populated with the factor level integers
#'
#' @param x dataframe containing factor data
#' @param cnames character vector of column names to be used as factors
#' @param clevels character vector of factor levels, determining the order of the factor levels
#' 
#' @return numeric matrix
#'
#' @author Susan Huse \email{susan.huse@@nih.gov}
#' @keywords RNASeq Pipeliner limma
#'
#' @examples
#' factordf_to_mtx(groupsDF, cnames=covars, clevels=covar_levels)
#' 
#' @export
factordf_to_mtx <- function(x, cnames=covars, clevels=covar_levels) {
  for(n in cnames) {
    f <- factor(x[,n], levels=clevels[[n]])
    x[,n] <- as.integer(f)
  }
  return(as.matrix(x))
}
