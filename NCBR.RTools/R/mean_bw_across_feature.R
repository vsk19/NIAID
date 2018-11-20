##################
#
# Frederick National Laboratory (NCBR)
# Tovah Markowitz, tovah.markowitz@nih.gov
# November 19, 2018
#
##################
#' 
#' Calculate the mean bigwig score across a series of features
#' 
#' Requires: rtracklayer, GenomicRanges
#'
#' This function takes a set of regions(features) given in the form of a genomic ranges
#' object. For each region, the bigwig values are extracted using rtracklayer. Using this
#' information, the following data is saved: 1) the mean of the extracted values, 2) the
#' sum of the extracted values, 3) the number of non-zero values, and 4) the "score" that
#' is calculated as the sum divided by the non-zero length.
#'
#'
#' @param GR a genomic ranges object
#' @param bigwigFile name/location of the bigwig file
#' @param progressBar whether or not to run the progress bar, suggested when running interactive [default: TRUE]
#' @param outFile the name of the file to be created with the results of this function, optional.
#'
#' @return a genomic ranges object with new columns: bwMean, bwTotal, nNotZero, score
#'
#' @author Tovah Markowitz \email{tovah.markowitz@nih.gov}
#' @keywords bigwig bw GenomicRanges rtracklayer
#'
#' @examples
#' bed <- rtracklayer::import.bed("geneinfo.bed")
#' out <- mean_bw_across_feature(GR=bed, bigwigFile="mH2A_D.inputnorm.bigwig", 
#'     outFile= "mH2A_D_at_genes.txt") 
#' 
#'
#' @export
mean_bw_across_feature <- function(GR, bigwigFile, progressBar=TRUE, outFile) {

	library(rtracklayer)

	# prepare to add bw data to GR file
	GR$bwTotal <- NA
	GR$nNotZero <- NA
	GR$bwMean <- NA

	total <- length(GR)
	# to create a progress bar
	if (progressBar == TRUE) {
		pb <- txtProgressBar(min=0, max=total, style=3)
	}

	# for each feature: read in bw region, calculate mean bw score (bwMean),
	# calculate total score (bwTotal), calculate number of nonzero bases (nNotZero),
	# calculate score based upon total score and number of nonzero bases (score) 
	for (i in 1:total) {
		bwregion <- import(BigWigFile(bigwigFile), as="NumericList",
			selection=BigWigSelection(GR[i]))
		GR$bwMean[i] <- mean(bwregion[[1]])
		GR$bwTotal[i] <- sum(bwregion[[1]])
		GR$nNotZero[i] <- sum(bwregion[[1]] != 0)
    	GR$score[i] <- GR$bwTotal[i] / GR$nNotZero[i]
    	if (progressBar == TRUE) {
    		setTxtProgressBar(pb,i)
    	}
	}

	if (progressBar == TRUE) {
		close(pb)
	}

	if (!(missing(outFile))) {
		write.table(GR, outFile, quote=F, sep="\t", row.names=F)
	}

	return(GR)
}