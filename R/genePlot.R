#' @title Plot Gene Expression Data
#' @description Vissualize gene expression data.
#'
#' @details
#' This is a major to do.
#'
#' @param x R data object; Most typically this is an \code{ExpressionSet} there is support for other datatypes as well.
#' @param gene character; Gene or vector of gene names. These can either be rownames from the gene expression data or looked up in the feature data.
#' @param by character vecter, or dataframe of factors; Will look up grouping factors in the phenotype data (\code{pData}) of the R data object \code{x}. Factors and dataframes of factors can also be passed directly.
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
genePlot <- function(x, gene, by=NULL, plotType=c("box","dot","bar","violin","density"), symbol="GeneSymbol", main=TRUE, na.rm=FALSE, shiny=FALSE, groupByGene=TRUE, useNormCounts=TRUE, ...) {UseMethod("genePlot",x)}

#' @importFrom purrr map
#' @importFrom Biobase exprs pData fData
#' @export
genePlot.default <- function(x, gene=NULL, by=NULL, plotType=c("box","dot","bar","violin","density"), symbol="GeneSymbol", main=TRUE, na.rm=FALSE, shiny=FALSE, groupByGene=TRUE,...) {
  npOptions<-list(...)

  if(main==TRUE) {
    if(length(gene)>1) {
      main<-paste0(c("Gene Expression:",paste0(gene,collapse=", ")),collapse=" ")
    } else {
      main<-paste0(gene, " Expression")
    }
  }
  npOptions<-append(list(x=x,by=by,na.rm=na.rm,main=main),npOptions)
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
  } else {
    stop("invalid plot type")
  }
  invisible(dataOut)
}

#' @import dplyr
#' @importFrom purrr map
#' @export
genePlot.numeric <- function(x, gene=NULL, by=NULL, plotType=c("box","dot","bar","violin","density"), symbol="GeneSymbol", main=TRUE, na.rm=FALSE, shiny=FALSE, groupByGene=TRUE,...) {
  npOptions<-list(...)
  groupings<-0
  gnames<-NULL
  if(!is.null(npOptions[["subGroup"]])){
    if(npOptions[["subGroup"]]==TRUE) {
      groupings<-1
      gnames="subGroup"
    }
  }
  if(!is.null(npOptions[["pointHighlights"]])){
    if(npOptions[["pointHighlights"]]==TRUE) {
      gnames=c(gnames,"pointHighlights")
      groupings<-groupings+1
    }
  }
  if (groupings==2) {
    gnames<- "subGroup and pointHighlights"
  }
  byDim<-0
  if(!is.null(dim(by))){
    byDim<-dim(by)[2]
    for(i in 1:dim(by)[2]) {
      by[,i]<-factor(by[,i])
    }
  }
  if((is.null(by) | is.factor(by) | byDim<=1) & groupings>0) {
    warning(paste0("By is empty with ",gnames," active.\nSetting grouping factors to FALSE.\nAdd grouping factors using the 'by' option to utilize the options.\n"), call.=FALSE)
    npOptions[["subGroup"]]<-FALSE
    npOptions[["pointHighlights"]]<-FALSE
    npOptions[["legend"]]<-FALSE
  } else if(groupings>0 & byDim<groupings+1) {
      if (groupings==2){gnames<-"pointHighlights"}
      warning(paste0("Not enough factors in 'by' to support grouping options.\nSetting ",gnames, " to FALSE.\n"),call. = FALSE)
      npOptions[[gnames]]<-FALSE
      if(is.null(npOptions[["legend"]])){
        npOptions[["legend"]]<-TRUE
      }
  } else if(groupings>0) {
    if(is.null(npOptions[["legend"]])){
      npOptions[["legend"]]<-TRUE
    }
  }

  if(main==TRUE & !is.null(gene)) {
    if(length(gene)>1) {
      main<-paste0(c("Gene Expression:",paste0(gene,collapse=", ")),collapse=" ")
    } else {
      main<-paste0(gene, " Expression")
    }
  }
  npOptions<-append(list(x=x,by=by,na.rm=na.rm,main=main),npOptions)
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
  } else {
    stop("invalid plot type")
  }
  invisible(dataOut)
}

#' @importFrom purrr map
#' @export
genePlot.data.frame <- function(x, gene=NULL, by=NULL, plotType=c("box","dot","bar","violin","density"), symbol="GeneSymbol", main=TRUE, na.rm=FALSE, shiny=FALSE, groupByGene=TRUE,...) {
  npOptions<-list(...)
  groupings<-0
  gnames<-NULL
  if(!is.null(gene)) {
    if(sum(gene %in% rownames(x))>0) {
      x<-t(x[gene,])
    } else {
      stop(paste0("Unable to identify gene(s) ",paste0(gene,collapse = ", "),"in the rownames of the data provided"))
    }
  }
  if(!is.null(npOptions[["subGroup"]])){
    if(npOptions[["subGroup"]]==TRUE) {
      groupings<-1
      gnames="subGroup"
    }
  }
  if(!is.null(npOptions[["pointHighlights"]])){
    if(npOptions[["pointHighlights"]]==TRUE) {
      gnames=c(gnames,"pointHighlights")
      groupings<-groupings+1
    }
  }
  if (groupings==2) {
    gnames<- "subGroup and pointHighlights"
  }
  byDim<-0
  if(!is.null(dim(by))){
    byDim<-dim(by)[2]
    for(i in 1:dim(by)[2]) {
      by[,i]<-factor(by[,i])
    }
  } else if (is.factor(by) | is.character(by)) {
    by<-factor(by)
    byDim<-1
  }
  if((is.null(by) |  byDim==0) & groupings>0) {
    warning(paste0("By is empty with ",gnames," active.\nSetting grouping factors to FALSE.\nAdd grouping factors using the 'by' option to utilize the options.\n"), call.=FALSE)
    npOptions[["subGroup"]]<-FALSE
    npOptions[["pointHighlights"]]<-FALSE
    npOptions[["legend"]]<-FALSE
  } else if(groupings>0 & byDim<groupings & plotType!="density") {
    if (groupings==2){gnames<-"pointHighlights"}
    warning(paste0("Not enough factors in 'by' to support grouping options.\nSetting ",gnames, " to FALSE.\n"), call.=FALSE)
    npOptions[[gnames]]<-FALSE
    if(is.null(npOptions[["legend"]])){
      npOptions[["legend"]]<-TRUE
    }
  } else if(groupings>0) {
    if(is.null(npOptions[["legend"]])){
      npOptions[["legend"]]<-TRUE
    }
  }
  if(groupings>0) {
    if(is.null(npOptions[["flipFacts"]])){
      npOptions[["flipFacts"]]<-groupByGene
    }
  }
  if(main==TRUE & !is.null(gene)) {
    if(length(gene)>1) {
      main<-paste0(c("Gene Expression:",paste0(gene,collapse=", ")),collapse=" ")
    } else {
      main<-paste0(gene, " Expression")
    }
  }
  npOptions<-append(list(x=x,by=by,na.rm=na.rm,main=main),npOptions)
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
  } else {
    stop("invalid plot type")
  }
  invisible(dataOut)
}

#' @importFrom purrr map
#' @importFrom Biobase exprs pData fData
#' @export
genePlot.ExpressionSet <- function(x, gene, by=NULL, plotType=c("box","dot","bar","violin","density"), symbol="GeneSymbol", main=TRUE, na.rm=FALSE, shiny=FALSE, groupByGene=TRUE, ...) {
  npOptions<-list(...)
  data<-1
  if(grepl("ENST",rownames(exprs(x))[1])){warning("This appears to be isoform data which may cause problems. Please us isoPlot instead.", call.=FALSE)}
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
  npOptions<-append(list(x=data,by=by,na.rm=na.rm,main=main,groupByGene=groupByGene,plotType=plotType,shiny=shiny),npOptions)
  dataOut<-do.call("genePlot", npOptions)
  invisible(dataOut)
}

#' @import tibble
#' @export
genePlot.tibble <- function(x, gene=NULL, by=NULL, plotType=c("box","dot","bar","violin","density"), symbol="GeneSymbol", main=TRUE, na.rm=FALSE, shiny=FALSE, groupByGene=TRUE, ...) {
  npOptions<-list(...)
  x<-as.data.frame(x)
  npOptions<-append(list(x=x,gene=gene,by=by,na.rm=na.rm,main=main,plotType=plotType,groupByGene=groupByGene),npOptions)
  dataOut<-do.call("genePlot", npOptions)
  invisible(dataOut)
}


#' @export
genePlot.matrix <- function(x, gene=NULL, by=NULL, plotType=c("box","dot","bar","violin","density"), symbol="GeneSymbol", main=TRUE, na.rm=FALSE, shiny=FALSE, groupByGene=TRUE, ...) {
  npOptions<-list(...)
  x<-data.frame(x)
  npOptions<-append(list(x=x,gene=gene,by=by,na.rm=na.rm,main=main,plotType=plotType,groupByGene=groupByGene),npOptions)
  dataOut<-do.call("genePlot", npOptions)
  invisible(dataOut)
}


#' @importClassesFrom EDASeq SeqExpressionSet
#' @importFrom purrr map
#' @importFrom EDASeq counts normCounts
#' @importFrom Biobase pData fData
#' @export
genePlot.SeqExpressionSet <- function(x, gene, by=NULL, plotType=c("box","dot","bar","violin","density"), symbol="GeneSymbol", main=TRUE, na.rm=FALSE, shiny=FALSE, groupByGene=TRUE, useNormCounts=TRUE, ...) {
  npOptions<-list(...)
  data<-1
  if(useNormCounts==TRUE & sum(is.na(normCounts(x)[1,]))>0) {
    warning("normCounts values in data set appear to be uninitialized as the first row return return as NAs\nProceeding with raw count data. Set useNormCounts=FALSE to avoid this warning.\n",call.=FALSE)
    useNormCounts<-FALSE
  }
  if(useNormCounts==FALSE) {
    warning("Raw count data lacks normalization and may not be a reliable metric.\nUse count based normalization protocols such as EDASeq or SVA to normalize or use data transformations such as voom (limma), cpm/rpkm (edgeR), regularized log (DESeq2) or variance stablizing transformation (DESeq2) for better results.", call. = FALSE)
  }
  if(grepl("ENST",rownames(counts(x))[1])){warning("This appears to be isoform data which may cause problems. Please us isoPlot instead.", call.=FALSE)}
  if(sum(gene %in% fData(x)[,symbol])>0){
    if(length(gene)>1) {
      if(useNormCounts==TRUE){
        data<-map(gene, function(g) normCounts(x)[fData(x)[,symbol]==g,]) %>% as.data.frame()
      } else {
        data<-map(gene, function(g) counts(x)[fData(x)[,symbol]==g,]) %>% as.data.frame()
      }
      colnames(data)<-gene
    } else {
      if(useNormCounts==TRUE){
        data<-normCounts(x)[which(fData(x)[,symbol]==gene),]
      } else {
        data<-counts(x)[which(fData(x)[,symbol]==gene),]
      }
    }
  } else if(sum(gene %in% rownames(counts(x)))>0) {
    if(length(gene)>1) {
      if(useNormCounts==TRUE){
        data<-t(normCounts(x)[gene,])
      } else {
        data<-t(counts(x)[gene,])
      }
    } else {
      if(useNormCounts==TRUE){
        data<-normCounts(x)[gene,]
      } else {
        data<-counts(x)[gene,]
      }
    }
  }else {
    stop("unable to identify gene listed in annotation or rownames of the data provided")
  }
  if(!is.null(by) & is.character(by)) {
    by<-pData(x)[,by]
  } else{
    stopifnot((is.data.frame(by) | is.null(by) | is.factor(by)))
  }
  if(main==TRUE) {
    if(length(gene)>1) {
      if(useNormCounts==TRUE) {
        main<-paste0(c("Normalize Gene Counts:",paste0(gene,collapse=", ")),collapse=" ")
    } else {
        main<-paste0(c("Raw Gene Counts:",paste0(gene,collapse=", ")),collapse=" ")
      }
    } else {
      if(useNormCounts==TRUE) {
        main<-paste0(gene, " Normalized Counts")
      } else {
        main<-paste0(gene, " Raw Counts")
      }
    }
  }
  npOptions<-append(list(x=data,by=by,na.rm=na.rm,main=main,groupByGene=groupByGene,plotType=plotType,shiny=shiny),npOptions)
  dataOut<-do.call("genePlot", npOptions)
  invisible(dataOut)
}



#' @importClassesFrom DESeq2 DESeqTransform
#' @importFrom purrr map
#' @importFrom SummarizedExperiment assay colData
#' @export
genePlot.DESeqTransform <- function(x, gene, by=NULL, plotType=c("box","dot","bar","violin","density"), symbol="GeneSymbol", main=TRUE, na.rm=FALSE, shiny=FALSE, groupByGene=TRUE, ...) {
  npOptions<-list(...)
  data<-1
  if(grepl("ENST",rownames(assay(x))[1])){warning("This appears to be isoform data which may cause problems. Please us isoPlot instead.", call.=FALSE)}
  if(sum(gene %in% rownames(assay(x)))>0) {
    if(length(gene)>1) {
      data<-t(assay(x)[gene,])
    } else {
      data<-assay(x)[gene,]
    }
  }else {
    stop("unable to identify gene listed in annotation or rownames of the data provided")
  }
  if(!is.null(by) & is.character(by)) {
    by<-colData(x)[,by]
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
  npOptions<-append(list(x=data,by=by,na.rm=na.rm,main=main,groupByGene=groupByGene,plotType=plotType,shiny=shiny),npOptions)
  dataOut<-do.call("genePlot", npOptions)
  invisible(dataOut)
}

#' @importClassesFrom limma EList
#' @importFrom purrr map
#' @export
genePlot.EList <- function(x, gene, by=NULL, plotType=c("box","dot","bar","violin","density"), symbol="GeneSymbol", main=TRUE, na.rm=FALSE, shiny=FALSE, groupByGene=TRUE, ...) {
  npOptions<-list(...)
  data<-1
  if(grepl("ENST",rownames(x$E)[1])){warning("This appears to be isoform data which may cause problems. Please us isoPlot instead.", call.=FALSE)}
  if(sum(gene %in% x$gene[,symbol])>0) {
    if(length(gene)>1) {
      data<-map(gene, function(g) x$E[rownames(x$gene)[x$gene[,symbol]==g],]) %>% as.data.frame()
      colnames(data)<-gene
    } else {
      data<-x$E[rownames(x$gene)[which(x$gene[,symbol]==gene)],]
    }
  } else if(sum(gene %in% rownames(x$E))>0) {
    if(length(gene)>1) {
      data<-t(x$E[gene,])
    } else {
      data<-x$E[gene,]
    }
  }else {
    stop("unable to identify gene listed in annotation or rownames of the data provided")
  }
  if(!is.null(by) & is.character(by)) {
    by<-x$design[,by]
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
  npOptions<-append(list(x=data,by=by,na.rm=na.rm,main=main,groupByGene=groupByGene,plotType=plotType,shiny=shiny),npOptions)
  dataOut<-do.call("genePlot", npOptions)
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

