####################
#
#  Frederick National Laboratory
#  Justin Lack, justin.lack@@nih.gov
#  January 7, 2019
#
####################
#' Annotate genes using bioDBnet
#'
#' Return annotations from multiple databases included in bioDBnet for a list of genes
#'
#' Takes a list of gene symbols, constructs a curl command that accesses 
#' the FNLCR Biological Database Network, bioDBnet, and returns a dataframe
#' containing annotations from each database for each gene.
#' Optionally writes out the annotations to a text file.
#'  
#' @param glist input character vector of gene symbols (default=NULL)
#' @param gfile input text file, one column, each line includes one gene symbol (default=NULL)
#' @param outfile optional output text file to write annotation dataframe to (default=NULL)
#' @param delim field delimiter used in output text file (default=tab)
#' @param taxonID taxonomic ID for the gene list (default=9606 for Homo sapiens, use 10090 for Mus musculus)
#' @param maxN maximum number of genes to send in each curl command (default=1000)
#' 
#' @return annotation dataframe with gene symbols are row names and annotation databases are column names
#' additionally, will output to file as comma or tab delimited fields
#'
#' @author Justin Lack \email{justin.lack@@nih.gov}
#' @keywords bioDBnet annotation
#'
#' @examples
#' biodbnet_annot(gfile="mygenes.txt", outfile="myannots.txt")
#' myannots <- (genes=mygenes, taxonID=10090)
#' 
#' @importFrom jsonlite flatten fromJSON
#' @importFrom RCurl getURL
#' @import dplyr
#' @export
biodbnet_annot <- function(genes=NULL, gfile=NULL, outfile=NULL, taxonID=9606, maxN=1000, delim="\t") {
  # Load up the vector of gene names
  # first check how to read the genes
  if(! is.null(genes) & ! is.null(gfile) ) {
    stop("Provide input genes using either a vector object or a filename, but not both")
  }
  if(is.null(genes) & is.null(gfile)) {
    stop("Provide input genes using either a vector object or a filename")
  }
  
  # read in the genes from a file if necessary
  if(! is.null(gfile)) {
    genes <- read.table(gfile, stringsAsFactors = FALSE)[,1]
  }

  # # genes <- genes[1:10]
  # genes <- c("SPON2","LOC100130872","CTBP1-AS","CTBP1","CTBP1-AS2","MAEA","UVSSA",
  #            "CRIPAK","NKX1-1","PPP2R2C","CPZ","USP17L10","USP17L10","USP17L11",
  #            "TP53","BRAF","BRCA2","KRAS")
  
  # set up the curl text
  curlText1 <- "https://biodbnet-abcc.ncifcrf.gov/webServices/rest.php/biodbnetRestApi.json?method=db2db&format=row&input=genesymbol&inputValues="
  # SPON2,LOC100130872,CTBP1-AS,CTBP1,CTBP1-AS2,MAEA,UVSSA,CRIPAK,NKX1-1,PPP2R2C,CPZ,USP17L10,USP17L10,USP17L11,TP53,BRAF,BRCA2,KRAS
  curlText2 <- "&outputs=biocartapathwayname,go-biologicalprocess,go-cellularcomponent,go-molecularfunction,hprdproteincomplex,keggdiseaseid,keggpathwayid,keggpathwayinfo,keggpathwaytitle,ncipidpathwayname,ncipidproteincomplex,pantherid,pfamid,reactomepathwayname,tigerfamsid"
  curlSpecies <- paste0("&taxonId=", taxonID)

  # Loop for maxN genes at a time
  loopN <- length(genes)%/%maxN
  if(length(genes)%%maxN > 0) {loopN <- loopN + 1} # add the last one if there is a remainder
  
  print(paste("Sending ", loopN, " curl commands to bioDBnet, this may take a while..."))
  annotsList <- list()
  for(i in 1:loopN) {
    # get the gene vector parameters for this iteration
    first <- (i-1)*maxN + 1
    last <- i*maxN
    if(last > length(genes)) {last <- length(genes)}
    # print(paste(i, first, last))
  
    # Create the curl text for the gene subset
    curlGenes <- paste(genes[first:last], collapse=",")
    curlText <- paste0(curlText1, curlGenes, curlText2, curlSpecies)
    
    # send the curl command and output the dataframe
    curlResponse <- RCurl::getURL(curlText)
    annotsList[[i]] <- jsonlite::flatten(jsonlite::fromJSON(curlResponse))
  }
  
  # combine all of the annotation tables and export them
  annots <- dplyr::bind_rows(annotsList)
  colnames(annots)[1] <- "Genes"
  
  if(! is.null(outfile)) {
    write.table(annots, file=outfile, quote = FALSE,row.names = FALSE, sep=delim)
  }
  return(annots)
}