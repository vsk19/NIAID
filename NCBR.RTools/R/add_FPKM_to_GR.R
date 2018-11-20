##################
#
# Frederick National Laboratory (NCBR)
# Tovah Markowitz, tovah.markowitz@nih.gov
# November 19, 2018
#
##################
#' Add FPKM values to a genomic ranges structure
#'
#' Combine a data.frame with FPKM values and a genomic ranges object
#'
#' Requires: GenomicRanges
#'
#' Returns a genomic ranges object with FPKM values (mean of replicates, if given)
#' Gene symbols or ensemblIDs can also be added to to the object if requested and
#' and available in the FPKM data.frame
#' NOTE: current version assumes that length of the genomic ranges object and the 
#' number of rows in the data.frame are equal and that the gene sets are identical, just
#' in a different order.
#'
#'
#' @param GR a genomic ranges object containing gene names in metadata column 'name'
#' @param FPKM a data.frame containing FPKM data and gene names of the same format as those in GR$name
#' @param FPKMgene the name of the column of FPKM to be compared to GR$name for merging
#' @param columns a vector of column names from FPKM with FPKM values to be analyzed
#' @param FPKMsymbol the name of the column with gene symbols/names to be added to GR [default:""]
#' @param FPKMensembl the name of the column with ensembl gene IDs to be added to GR [default:""]
#' @param outFile the name of the file to be created with the results of this function, optional.
#'
#' @return a genomic ranges object with FPKM data
#'
#' @author Tovah Markowitz \email{tovah.markowitz@nih.gov}
#' @keywords RNAseq GenomicRanges FPKM
#'
#' @examples
#' bed <- rtracklayer::import.bed("geneinfo.bed")
#' RSEM <- read.table("RSEM.genes.FPKM.all_samples.txt", header=T, stringsAsFactors=F)
#' out <- add_FPKM_to_GR(GR=bed, FPKM=RSEM, FPKMgene="gene_id", 
#'    columns=c("S35_FPKM","S36_FPKM"), FPKMsymbol="GeneName",outFile= "geneinfo_FPKM.txt")
#'
#' @export
add_FPKM_to_GR <- function(GR, FPKM, FPKMgene, columns, FPKMsymbol=FALSE, 
		FPKMensembl=FALSE, outFile) {

	# reorder FPKM to match GR
	FPKMsorted <- FPKM[order(match(FPKM[[FPKMgene]], GR$name)), ]

	# add other gene names if requested/available
	if ( FPKMsymbol != "" ) {
		GR$symbol <- FPKMsorted[[FPKMsymbol]]
	}
	
	if ( FPKMensembl != "" ) {
		GR$ensembl <- FPKMsorted[[FPKMensembl]]
	}

	# get data from each column of FPKM values requested
	# if multiple columns, assume replicates and average
	if (length(columns) > 1) {
   		tmp <- FPKMsorted[[columns[1]]]
   		for (i in 2:length(columns)) {
       		tmp <- cbind(tmp, FPKMsorted[[columns[i]]])
   		}
   		GR$meanFPKM <- apply(tmp, 1, mean)
	} else {
 		GR$FPKM <- FPKMsorted[[columns]]
	}

	# if an output file name is given, write output
	if (!(missing(outFile))) {
		write.table(GR, outFile, quote=F, sep="\t", row.names=F)
	}
	
	return(GR)
}