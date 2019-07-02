#' @title Plot Gene Expression Data
#' @description Vissualize gene expression data.
#'
#' @details
#' This is a major to do.
#'
#' @param x R data object; Most typically this is a seqExpressionSet there is support for other datatypes as well.
#' @param gene character; Gene or vector of gene names. These can either be rownames from the gene expression data or looked up in the feature data.
#' @param by positive integer; sets the cex multiplier for point size.
#' @param plotType character; method to be used for ploting dots. Can be set to "jitter", "linear", "beeswarm" or "distribution".
#' @param main character; The main plot title. Defaults to true for automated generation.
#' @param symbol numeric; scaling factor controlling the width of the boxes.
#' @param na.rm positive integer; sets pty for plotting data points. Can be a vector to support additional graphical customization.
#' @param shiny logical; if a p-value can be easily calculated for your data, it will be displayed using the \code{sub} annotation setting.
#' @param ... additional options for S3 method variants
#'
#' @examples
#' ToDo<-1
#' @import dplyr
#' @importFrom purrr map
#' @importFrom Biobase exprs pData fData
#' @export
#' @seealso \code{\link{boxplot}}
genePlot <- function(x, gene, by=NULL, plotType=c("box","dot","bar","violin","density"), symbol="GeneSymbol", main=TRUE, na.rm=FALSE, shiny=FALSE, ...) {UseMethod("genePlot",x)}

#' @import dplyr
#' @importFrom purrr map
#' @importFrom Biobase exprs pData fData
#' @export
genePlot.default <- function(x, gene, by=NULL, plotType=c("box","dot","bar","violin","density"), symbol="GeneSymbol", main=TRUE, na.rm=FALSE, shiny=FALSE, ...) {
  npOptions<-list(...)
  data<-1
  if(sum(gene %in% fData(x)[,symbol])>0){
    if(length(gene)>1) {
      data<-map(gene, function(g) exprs(x)[fData(x)[,symbol]==g,]) %>% as.data.frame()
      colnames(data)<-gene
    } else {
      data<-exprs(x)[which(fData(x)[,symbol]==gene),]
    }
  } else if(sum(gene %in% rownames(exprs(x)))>0) {
    if(length(gene)>1) {
      data<-t(exprs(x)[gene,])
    } else {
      data<-exprs(x)[gene,]
    }
  }else {
    stop("unable to identify gene listed in annotation or rownames of the data provided")
  }
  if(!is.null(by) & is.character(by)) {
    print("hi!")
    by<-pData(x)[,by]
  } else{
    stopifnot((is.data.frame(by) | is.null(by) | is.factor(by)))
  }
  if(main==TRUE) {
    if(length(gene)>1) {
      main<-paste0(c("Gene Expression:",paste0(gene,collapse=", ")),collapse=" ")
    } else {
      main<-paste0(gene, " Expression")
    }
  }
  npOptions<-append(list(x=data,by=by,na.rm=na.rm,main=main),npOptions)
  if(plotType[1]=="box"){
    print(npOptions)
    do.call("niceBox",npOptions)
  } else if (plotType[1]=="dot") {
    do.call("niceDots",npOptions)
  } else if (plotType[1]=="violin") {
    do.call("niceVio",npOptions)
  } else if (plotType[1]=="bar") {
    do.call("niceBar",npOptions)
  } else if (plotType[1]=="density") {
    do.call("niceDensity",npOptions)
  } else {
    stop("invalid plot type")
  }
}
