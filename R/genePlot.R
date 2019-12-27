#' @title Plot Gene Expression Data
#' @description Vissualize gene expression data.
#'
#' @details
#' The \code{genePlot} is designed to make vissualizatin of gene expression data simple and easy for R novices and bioinformaticians alike.
#' The function is an S3 generic that accept various R and Bioconductor datasets as input and exracts the expression, factor and annotation data from them according to type.
#' The factors allow for spliting expression data from one or more genes into groups and for plot types with data point overlays, points can be colored by factors levels as well.
#' If the input data is a Bioconductor data set such as an \code{\link[Biobase]{ExpressionSet}} and the \code{gene} option is used, \code{genePlot} will attempt to look up the genes in the associated feature annotation data (e.g. \code{\link[Biobase]{fData}}) according to the data input type and look for the gene symbol column indicated by the \code{symbol} option (defaults to 'GeneSymbol').
#' If no matches are found the row names of are checked of the expression data are check for matches as well.
#' If charactar values are given for factor input, \code{genePlot} will attempt to look up assocation phenotype data (e.g. \code{\link[Biobase]{pData}}).
#' One can also pass raw data vectors/data frames and/or factors to \code{genePlots} to bypass this feature, which is critical for data sets and data formats where integrated phenotype and feature data is not available.
#' The \code{genePlot} uses the \code{NicePlots} graphics library and any \code{NicePlots} option and/or theme can be used in conjuction with options detailed below.
#' The \code{plotType} options supported correspond to \code{NicePlots} functions and include box plots (\code{\link[NicePlots]{niceBox}}), dot plots (\code{\link[NicePlots]{niceDots}}), violin plots (\code{\link[NicePlots]{niceVio}}), bar plots (\code{\link[NicePlots]{niceBar}}) as well as both one/two dimentional kernal density plots (\code{\link[NicePlots]{niceDensity}}).
#' Supported data input types include: \code{\link[Biobase]{ExpressionSet}}, \code{\link[EDASeq]{SeqExpressionSet-class}}, \code{\link[limma]{EList-class}}, \code{\link[DESeq2]{DESeqTransform}}, as well as standard R data types such as \code{\link[base]{vector}}, \code{\link[base]{matrix}}, \code{\link[base]{data.frame}}, and \code{\link[tibble]{tibble}}.
#' \code{genePlot} silently returns a list of class \code{npData} that conatains a summarized findings, p-values (if indicated), extracted plotting data, and plotting options.
#' All npData objects can be replotted using  the \code{\link[graphics]{plot}} function, \code{genePlot} or any of the \code{NicePlots} functions.
#' Options passed to any of these, including \code{plotType} will override the options for the \code{npData} object.
#'
#' @param x R data object; Most typically this is an \code{ExpressionSet} there is support for other datatypes as well.
#' @param gene character; Gene or vector of gene names. These can either be rownames from the gene expression data or looked up in the feature data.
#' @param group factor or name of factor to be exracted from \code{x} (e.g. \code{\link[Biobase]{pData}}). Used as the primary grouping factor.
#' @param subGroup factor or name of factor to be exracted from \code{x} (e.g. \code{\link[Biobase]{pData}}). Used to subgroup data unless multiple genes are selected in which case \code{subGroup} is ignored.
#' @param highlight factor or name of factor to be exracted from \code{x} (e.g. \code{\link[Biobase]{pData}}). Used to color data points by factor levels. Only valid for graphs with point overlays.
#' @param facet factor or name of factor to be exracted from \code{x} (e.g. \code{\link[Biobase]{pData}}). Split the data into multiple smaller graphs.
#' @param plotType character; Can be set to "box", "violin, "dot", "bar" or "denisity" for boxplots, violin plots, dot plots, bar plots, and kernal desity plots, respectively.
#' @param main character; The main plot title. Defaults to true for automated generation.
#' @param symbol character; Colname of of gene symbols in the feature data of \code{x} (\code{fData}).
#' @param na.rm logical; Removes \code{\link{NA}} values prior to ploting.
#' @param shiny logical; Use \code{\link[shiny]{shiny}} interfaces if available.
#' @param groupByGene logical; If more then one gene is listed and \code{grouByGene} is \code{TRUE}
#' @param useNormCounts logical; By default \code{genePlot} will try to use normCounts instead of counts in \code{SeqExpressionSets}. Set to FALSE to use raw counts instead, though this will generate a warning about useing non-normalized data.
#' @param ... Any paramenter recognized by \code{NicePlots} functions.
#'
#' @examples
#' ToDo<-1
#'
#' @importFrom purrr map
#' @importFrom Biobase exprs pData fData
#' @export
#' @seealso \code{\link[NicePlots]{niceBox}} \code{\link[NicePlots]{niceVio}} \code{\link[NicePlots]{niceBar}} \code{\link[NicePlots]{niceDots}} \code{\link[NicePlots]{niceDensity}}
genePlot <- function(x, gene, plotType=c("box","dot","bar","violin","density","suface"), symbol="GeneSymbol", main=TRUE, na.rm=TRUE, group=NULL, subGroup=NULL, highlight=NULL, facet=NULL, shiny=FALSE, groupByGene=TRUE, useNormCounts=TRUE, ...) {UseMethod("genePlot",x)}

#' @importFrom purrr map
#' @importFrom NicePlots niceBox niceVio niceBar niceDensity
#' @importFrom Biobase exprs pData fData
#' @export
genePlot.default <- function(x, gene=NULL, plotType=c("box","dot","bar","violin","density","surface"), symbol="GeneSymbol", main=TRUE, na.rm=TRUE, group=NULL, subGroup=NULL, highlight=NULL, facet=NULL, shiny=FALSE, groupByGene=TRUE,...) {
  npOptions<-list(...)

  if(main==TRUE) {
    if(length(gene)>1) {
      main<-paste0(c("Gene Expression:",paste0(gene,collapse=", ")),collapse=" ")
    } else {
      main<-paste0(gene, " Expression")
    }
  }
  data<-getGeneData(x=x, gene=gene, plotType=plotType, symbol=symbol,group=group, subGroup=subGroup,highlight=highlight,facet=facet)
  if(is.null(subGroup)){
    subGroup<-FALSE
  } else {
    subGroup<-TRUE
  }
  if(is.null(highlight)){
    highlight<-FALSE
  } else {
    highlight<-TRUE
  }
  npOptions<-append(list(x=data$x,by=data$by,pointHighlights=highlight, subGroup=subGroup, facet=facet, na.rm=na.rm,main=main),npOptions)
  dataOut<-1
  if(plotType[1]=="box"){
    dataOut<-do.call("niceBox",npOptions)
  } else if (plotType[1]=="dot") {
    dataOut<-do.call("niceDots",npOptions)
  } else if (plotType[1]=="violin") {
    dataOut<-do.call("niceVio",npOptions)
  } else if (plotType[1]=="bar") {
    dataOut<-do.call("niceBar",npOptions)
  } else if (plotType[1]=="density") {
    dataOut<-do.call("niceDensity",npOptions)
  } else if (plotType[1]=="surface") {
    npOptions<- append(list(plotType="surface"),npOptions)
    dataOut<-do.call("niceDensity",npOptions)
  } else {
    stop("invalid plot type")
  }
  invisible(dataOut)
}

# genePlot.npOutput <- function(x, gene, by=NULL, plotType=c("box","dot","bar","violin","density"), symbol="GeneSymbol", main=TRUE, na.rm=FALSE, shiny=FALSE, groupByGene=TRUE, ...) {
#   npOptions<-list(...)
#   #NOT DONE
#   data<-1
#   if(grepl("ENST",rownames(exprs(x))[1])){warning("This appears to be isoform data which may cause problems. Please us isoPlot instead.", call.=FALSE)}
#   if(sum(gene %in% fData(x)[,symbol])>0){
#     if(length(gene)>1) {
#       data<-map(gene, function(g) exprs(x)[fData(x)[,symbol]==g,]) %>% as.data.frame()
#       colnames(data)<-gene
#     } else {
#       data<-exprs(x)[which(fData(x)[,symbol]==gene),]
#     }
#   } else if(sum(gene %in% rownames(exprs(x)))>0) {
#     if(length(gene)>1) {
#       data<-t(exprs(x)[gene,])
#     } else {
#       data<-exprs(x)[gene,]
#     }
#   }else {
#     stop("unable to identify gene listed in annotation or rownames of the data provided")
#   }
#   if(!is.null(by) & is.character(by)) {
#     by<-pData(x)[,by]
#   } else{
#     stopifnot((is.data.frame(by) | is.null(by) | is.factor(by)))
#   }
#   if(main==TRUE) {
#     if(length(gene)>1) {
#       main<-paste0(c("Gene Expression:",paste0(gene,collapse=", ")),collapse=" ")
#     } else {
#       main<-paste0(gene, " Expression")
#     }
#   }
#   npOptions<-append(list(x=data,by=by,na.rm=na.rm,main=main,groupByGene=groupByGene,plotType=plotType,shiny=shiny),npOptions)
#   dataOut<-do.call("genePlot", npOptions)
#   invisible(dataOut)
# }

