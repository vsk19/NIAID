####################
#
#  Frederick National Laboratory
#  Susan Huse, susan.huse@@nih.gov
#  August 9, 2018
#
####################
#' Creates an ExpressionSet and covariates dataframe from RNASeq pipeline
#'
#' Imports a raw counts file from the RNASeq pipeliner output (RawCountFile_RSEM_genes.txt)
#'  and a groups or targets file with covariates, 
#'  and creates an numeric matrix of raw counts and covariates dataframe for linear modeling with limma.
#'  Returns a list containing 1) raw counts matrix and 2) groups/covariates dataframe
#'  
#' @param counts_file raw counts file from the RNASeq pipeline (RawCountFile_RSEM_genes.txt)
#' @param groups_file groups file from RNASeq pipeline, tab-delimited
#'                    Col1 = SampleID, Col2 = group, Col3 = SampleName, Col4..Col# covariates
#' 
#' @return list object with [[1]] = numeric matrix of raw count data, 
#'                          [[2]] = dataframe with sample names and covariates
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
load_RawCountsGroup <- function(counts_file, groups_file) {
  # Load Raw Counts file and the groups file
  if(! file_test("-f", counts_file)) {stop(paste("Unable to locate counts file:", counts_file))}
  if(! file_test("-f", groups_file)) {stop(paste("Unable to locate groups file:", groups_file))}

  df <- read.csv(file=counts_file, header=TRUE, sep="\t", as.is=TRUE)
  grpDF <- read.csv(file=groups_file, header=TRUE, sep="\t", as.is=TRUE)
  colnames(grpDF)[1] <- "SampleID"
  colnames(grpDF)[3] <- "SampleName"
  
  # Clean up the sample IDs
  IDs <- sapply(colnames(df), 
                function(x) {gsub(".RSEM.genes.*$", "", gsub("^.*DEG_RSEM_genes.", "", x))})
  unname(IDs)
  
  # Be sure the groups dataframe is in the same order as the expression dataframe
  rownames(grpDF) <- grpDF$SampleID
  names <- grpDF[IDs[2:length(IDs)],3]
  
  # reset rownames to SampleName
  rownames(grpDF) <- grpDF$SampleName
  
  # Set the df column names to SampleName rather than SampleID
  colnames(df) <- c("Symbol", names)
  
  # Convert to matrix and reset the rownames to Symbol
  mtx <- as.matrix(df[,2:ncol(df)], dimnames=list(df$Symbol, colnames(df[2:ncol(df)])))
  rownames(mtx) <- df$Symbol
  
  # Factorize the grpDF
  if(ncol(grpDF) > 3) {colset <- c(2, 4:ncol(grpDF))} else {colset <- 2}
  for(i in colset) {
    # print(i)
    if(is.character(grpDF[,i])) 
    {
      # print("is character")
      # cvname <- colnames(grpDF)[i]
      # # print(paste("CVName: ",cvname))
      # if(cvname %in% names(cvlevels)){
      #   grpDF[,cvname] <- factor(grpDF[,cvname], levels=cvlevels[[cvname]])
      #   print(grpDF[,cvname])
      # } else {
        grpDF[,i] <- factor(grpDF[,i])
        # print(grpDF[,i])
      # }
    }
  }
  
  # Return the count matrix and the groups Dataframe
  return(list(mtx, grpDF))
}
