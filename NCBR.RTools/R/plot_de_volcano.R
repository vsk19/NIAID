####################
#
#  Frederick National Laboratory
#  Susan Huse, susan.huse@@nih.gov
#  October 16, 2018
#
#  updated: 1/7/19 by Tovah Markowitz to improve plot coloring
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
#' @param interactive boolean, if TRUE use plot_ly to create an interactive plot, otherwise use ggplot (default=TRUE)
#' @param sigVal column used for determining significance, "adj.P.Val" or "P.Value" (default="adj.P.Val")
#' @param maxp maximum FDR or P value of genes considered statistically significant (default=0.1)
#' @param minLFC minimum log2 foldchange of genes considered biologically significant (default=1)
#' @param colSigHigh color name for the lines and points highlighting significant differential expression (default="red")
#' @param colSigLow color name for the points that are significant but have a low FC (default="green")
#' @param colNonLow color name for the points that are not significant and have a low FC (default="grey")
#' @param colNonHigh color name for the points that are not significant but have a high FC (default="grey")
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
plot_de_volcano <- function(f, coef, title="Volcano Plot", interactive=TRUE,
                            sigVal="adj.P.Val", maxp=0.1, minLFC=1.5, 
                            colSigHigh="orangered4", colSigLow="steelblue", 
                            colNonLow = "grey50", colNonHigh="palegreen4",
                            showsub=TRUE, miRNA=FALSE) {
  if(miRNA) {
    x <- topTags(f, n=nrow(f$coefficients))$table
    colnames(x)[grep("FDR", colnames(x))] <- "adj.P.Val"
  } else {
    x <- topTable(f, coef=coef, n=nrow(f$coefficients))
  }
  # print(colnames(x))
  # test p or FDR
  if (!sigVal %in% c("adj.P.Val", "P.Value")) {
    stop('significance value column, sigVal, must be either "adj.P.Val" or "P.Value"')
  }
  x$Gene <- rownames(x)
  x$logq <- -log10(x[, sigVal])
  # x$logFC <- x[, coef]
  minFDR <- -log10(maxp)
  sigHigh <- x$logq >= minFDR & abs(x$logFC) >= minLFC
  sigLow <- x$logq >= minFDR & abs(x$logFC) < minLFC
  nonHigh <- x$logq < minFDR & abs(x$logFC) >= minLFC
  nonLow <- x$logq < minFDR & abs(x$logFC) < minLFC

  x$Significant <- "Not Significant, LowFC"
  x$Significant[sigHigh] <- "Significant, HighFC"
  x$Significant[sigLow] <- "Significant, LowFC"
  x$Significant[nonHigh] <- "Not Significant, HighFC"
  sigLevels=c("Significant, HighFC", "Significant, LowFC", 
              "Not Significant, HighFC", "Not Significant, LowFC")
  # want to set the factor levels, but if one is missing it can crash
  x$Significant <- factor(x$Significant, levels=sigLevels[sigLevels %in% unique(x$Significant)]) 
  
  cnt_plus <- sum(x$logq >= minFDR & x$logFC >= minLFC)
  cnt_neg <- sum(x$logq >= minFDR & x$logFC <= -minLFC)
  if(showsub) {
    subtitle <- paste0(sum(x$logq >= minFDR), " DE genes: ", cnt_plus, " Pos and ", cnt_neg, " Neg",
                       " (FDR<", maxp, ", log2FC>", round(minLFC,2), ")")
  } else {
    subtitle <- NULL
  }
  # xaxislist <- list(title="Fold Change", range=c(-5,5), 
  #                   tickvals=seq(-5,5,1), ticktext=c("-32","-16","-8","-4","-2","1","2","4","8","16","32"))
  xaxislist <- list(title="log2 Fold Change")
  
  if(sigVal=="adj.P.Val") {
    yaxislist <- list(title="-Log10 FDR") #, range=c(0,10))
  } else {
    yaxislist <- list(title="-Log10 p-value") #, range=c(0,10))
  }
  
  if(interactive) {
    colvect <- c(colSigHigh, colSigLow, colNonHigh, colNonLow)
    p <- plot_ly(data=x, x=~logFC, y=~logq, text=~Gene, color=~Significant, colors=colvect,
                 type="scatter", mode="markers", marker=list(size=3)) %>%
      add_lines(x = -minLFC, showlegend=FALSE, line = list(color = I(colSigHigh), width = 2)) %>%
      add_lines(x = minLFC, showlegend=FALSE, line = list(color = I(colSigHigh), width = 2)) %>%
      add_lines(y = minFDR, showlegend=FALSE, line = list(color = I(colSigHigh), width = 2)) %>%
      layout(title = title, xaxis=xaxislist, yaxis=yaxislist)
  } else {
    p <- ggplot(x, aes(x=logFC, y=logq)) +
      geom_point(aes(x=logFC, y=logq), size=0.5) +
      geom_point(data=x[sigHigh,], size=0.5, col=colSigHigh) +
      geom_point(data=x[sigLow,], size=0.5, col=colSigLow) +
      geom_point(data=x[nonHigh,], size=0.5, col=colNonHigh) +
      geom_point(data=x[nonLow,], size=0.5, col=colNonLow) +
      ggtitle(title, subtitle=subtitle) + xlab("log2 Fold Change") + ylab("-log10 FDR") +
      theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
      geom_vline(xintercept=minLFC, col=colSigHigh, linetype="dashed", size=0.5) +
      geom_vline(xintercept=-minLFC, col=colSigHigh, linetype="dashed", size=0.5) +
      geom_hline(yintercept=minFDR, col=colSigHigh, linetype="dashed", size=0.5)
  }
  return(p)
}
