####################
#
#  Frederick National Laboratory
#  Susan Huse, susan.huse@@nih.gov
#  August 9, 2018
#
####################
#' Splits the combined EnsemblID and gene name
#'
#' Splits a character or character vector containing EnsemblID|Gene
#'  and returns either the EnsemblID (without decimal values) or the gene name
#'  
#' @param x string or character vector of "symbol" from RNASeq pipeline EnsemblID|Gene
#' @param s either "E" or "G" to return the cleaned EnsemblID or the Gene name
#' 
#' @return string or character vector of the requested names
#'
#' @author Susan Huse \email{susan.huse@@nih.gov}
#' @keywords RNASeq Pipeliner limma
#'
#' @examples
#' split_symbol("ENSMUSG00000000001.4|Gnai3", "E")
#' 
#' @export
split_symbol <- function(x, s=c("E","G")) {
  if(! s %in% c("E", "G")) {
    stop("output symbol must be either 'E' for EnsemblID or 'G' for Gene")
  }
  
  # Split on pipe - returns a list
  splitvalues <- strsplit(x, "\\|")
  
  # If E, take first element and remove trailing .# 
  # if G, take the second element
  if (s == "E") {
    y <- sapply(splitvalues, function(x) {
      y <- gsub(".[0-9]*$", "", x[1], perl=TRUE)
      return(y)
    })
    
  } else if(s == "G") {
    y <- sapply(splitvalues, function(x) {return(x[2])})
  }
  return(y)
}
