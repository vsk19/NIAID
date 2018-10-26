####################
#
#  Frederick National Laboratory
#  Susan Huse, susan.huse@@nih.gov
#  October 23, 2018
#
####################
#' print the head(obj) in 2 dimensions
#'
#' Prints out the head of an object in multiple dimensions.
#' The head() function prints out the first 6 rows, but all columns of a data.frame or matrix
#' head2 will print out a limited number of rows and columns.  
#' 
#' If the input object is a matrix or data.frame, it will be truncated in both rows (x) and columns (y).
#' If the input object is a list, the list object will be limited to the first x items, 
#' and each item within it will be printed using head2 recursively.
#' If the input object is a vector, the first x items will be printed.
#' 
#' @param obj an R object to print
#' @param x the number of items in the first dimension (length or rows) to print
#' @param y the number of items in the second dimension (columns or list subitems) to print
#' 
#' @return beginning of the input object is printed, returns TRUE if end is reached
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
head2 <- function(obj, x=5, y=5) {
  # Matrix or Dataframe
  if(class(obj) %in% c("matrix", "data.frame")) {
    print(obj[1:min(x, nrow(obj)), 1:min(y, ncol(obj))])

  # Vectors
  } else if (class(obj) %in% c("vector")) {
    print(obj[1:min(x, length(obj))])

  # List objects
  } else if (class(obj) == "list") {
    for(i in 1:min(x, length(obj))) {
      print(head2(obj[[i]]))
    }

  # For now, punt on everything else
  } else {
    print(head(obj))
  }

  # exit
  return(TRUE)
}
