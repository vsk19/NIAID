####################
#
#  Frederick National Laboratory
#  Susan Huse, susan.huse@@nih.gov
#  October 16, 2018
#
####################
#' Create a volcano plot from a limma fit object
#'
#' Creates a volcano plot from a limma fit (MArrayLM) object, using FDR as the y-axis, and coloring significant genes.
#'
#' The limma fit is plotted using ggplot, 
#' with a horizontal line for the FDR threshold, 
#' and two vertical lines for the logFC threshold.  
#' The plot will have a main title and a subtitle with counts of differentially expressed genes.
#'  
#' @param f input limma fit object
#' @param coef name of the column in theFit$coefficients to plot
#' @param title text to use for the plot main title
#' @param maxq maximum FDR value of genes considered statistically significant (default=0.1)
#' @param minLFC minimum log2 foldchange of genes considered biologically significant (default=1)
#' @param sigcol color name for the lines and points highlighting significant differential expression (default="red")
#' @param showsub boolean, if TRUE show subtitle with counts, if FALSE no subtitling
#' @param miRNA boolean, specifying if fit object is an miRNA DGELRT object (default=FALSE)
#' 
#' @return ggplot object for printing
#'
#' @author Susan Huse \email{susan.huse@@nih.gov}
#' @keywords RNASeq limma volcano plot
#'
#' @examples
#' print(plot_de_volcano(f=myFit, coef="TreatmntTRUE", title="Experiment 1"))
#' 
#' @export
plot_de_volcano <- function(f, coef, title="Volcano Plot", maxq=0.1, minLFC=1, sigcol="red", showsub=TRUE, miRNA=FALSE) {
  if(miRNA) {
    x <- topTags(f, n=nrow(f$coefficients))$table
    colnames(x)[grep("FDR", colnames(x))] <- "adj.P.Val"
  } else {
    x <- topTable(f, coef=coef, n=nrow(f$coefficients))
  }
  
  x$logq <- -log10(x$adj.P.Val)
  minFDR <- -log10(maxq)
  sigs <- x$logq >= minFDR & abs(x$logFC) >= minLFC
  
  cnt_plus <- sum(x$logq >= minFDR & x$logFC >= minLFC)
  cnt_neg <- sum(x$logq >= minFDR & x$logFC <= -minLFC)
  if(showsub) {
    subtitle <- paste0(sum(sigs), " DE genes: ", cnt_plus, " Pos and ", cnt_neg, " Neg",
                       " (FDR<", maxq, ", log2FC>", minLFC, ")")
  } else {
    subtitle <- NULL
  }
  
  p <- ggplot(x) +
    geom_point(aes(x=logFC, y=logq), size=0.5) + 
    geom_point(data=x[sigs,], aes(x=logFC, y=logq), size=0.5, col=sigcol) +
    ggtitle(title, subtitle=subtitle) + xlab("log2 Fold Change") + ylab("log10 FDR") +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) + 
    geom_vline(xintercept=minLFC, col=sigcol, linetype="dashed", size=0.5) + 
    geom_vline(xintercept=-minLFC, col=sigcol, linetype="dashed", size=0.5) + 
    geom_hline(yintercept=minFDR, col=sigcol, linetype="dashed", size=0.5)
  
  return(p)
}
