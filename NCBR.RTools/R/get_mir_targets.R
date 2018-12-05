####################
#
#  Frederick National Laboratory
#  Susan Huse, susan.huse@@nih.gov
#  October 12, 2018
#
####################
#' Creates a target gene dataframe from a list of mature miRNA
#'
#' Reads a vector of mature miRNA names, 
#' compares them to a database of filtered targets, and
#' returns a dataframe with two columns: miRNA_ID, target_ID.
#' Download the miRDB database from http://mirdb.org/index.html, and gunzip the file.
#' The miRDB download (e.g., miRDB_v5.0_prediction_result.txt.gz) is a three-column, tab-delimited file:
#'   miRNA_ID, gene_ID, Target_Prediction_Score
#'   optionally, to initial load performance, pre-filter the file on the first three letters of each line that specifies the species
#'   e.g., hsa = Homo sapiens miRNA, or mmu = Mus musculus
#'   
#' The miRNA IDs in the miRDB file will be preceded by xxx- denoting the species name.  
#'   If the input vector of miRNA IDs does not have the species prefix, 
#'   you can specify one for the entire vector with the sp parameter.
#'   If sp is used to specify a prefix, do not include the "-", e.g., use "hsa" or "mmu".
#'
#' miRDB includes a target prediction score [http://mirdb.org/faq.html#How_to_interpret_the_target_prediction_score]
#'   A prediction score > 80 is recommended.  Use minscore to select a different threshold.
#'   
#' keepsp is used to remove a species prefix that has been added.  If TRUE, keep the species prefix from miRDB, 
#'   if FALSE, remove the prefix, returning the miRNA_IDs to the original vector input.
#' 
#' 
#' Returns a dataframe with three columns miRNA_ID, Target_ID, and Tgt_Pred_Score
#'  
#' @param x vector of mature miRNA IDs
#' @param f target_file a text file download of miRDB targets and prediction scores (see Details)
#' @param sp species prefix for miRNA IDs if not present in the input mis vector (see Details)
#' @param minscore integer value for minimum target prediction score to include (default = 80)
#' @param keepsp keep a species prefix that has been added to miRNA data to match miRDB (default=TRUE)
#' 
#' @return dataframe: a three-column dataframe: miRNA_ID, Target_ID, Tgt_Pred_Score
#'
#' @author Susan Huse \email{susan.huse@@nih.gov}
#' @keywords miRNA miRDB
#'
#' @examples
#' refdir <- "/Users/myself/miRDB"
#' topdf <- topTags(exactTest(known_mir))
#' get_mir_targets(rownames(topdf$table), file.path(refdir, "miRDB_v5.0_prediction_result.txt"), sp="mmu", minscore=80)
#' 
#' @export
get_mir_targets <- function(x, f, minscore=80, sp=NULL, keepsp=TRUE) {
  # known_targets = get_mir_targets(rownames(topdf$table), "/Users/husesm/FNL/Projects/Colitis/NCBR-5/miRDB/mmu.miRDB.pred.txt", sp="mmu", minscore=80)
  # Load the miRDB file of targets
  mirdb = read.table(f, header=FALSE, sep="\t")
  colnames(mirdb) <-  c("miRNA_ID", "Target_ID", "Tgt_Pred_Score")
  
  # update the vector if sp is not null
  if(! is.null(sp)) {
    x <- paste(sp, x, sep="-")
  }
  
  df <- mirdb[mirdb[,"miRNA_ID"] %in% x, ]
  df <- df[df[, "Tgt_Pred_Score"] >= minscore, ]
  
  if(! keepsp) {
    df$miRNA_ID <- sub(paste0("^", sp, "-"), "", df$miRNA_ID, perl=TRUE)
  }
  
  return(df)
}
  