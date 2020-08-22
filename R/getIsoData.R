#There are to major S3 functions at play here: showIsoforms and getIsoData.
#The showIsoforms function is used to extract and filter the isoform IDs based on user input with
#S3 variants are used to customize the data extraction to the data type.
#
#The getIsoData is a bit more lightweight than it's getGeneData counterpart due to the
#since the isoforms have already been identified. Moreover, once the data is is excrated
#we can call on getGeneData to finish the pre-processing for calling NicePlots. As with
#showIsoforms, S3 variants are used to customize the data extraction to the data type.
#
#Currently supported classes:
# ExpessionSet, SeqExpressionSet, DESeqTransform, SummarizedExperiment, EList

#' @title Display Isoform Annotation
#' @description Display and filter isoform annotation and find those associated with given genes.
#'
#' @details
#' This is a convenience function that will pull out isoform annotation from popular Bioconductor classes.
#' The isoforms can be selected by name and/or by thier association with a particular gene or genes.
#' If nothing else is specified, only the isoform IDs are returned. Isoforms can be filtered by \code{appris} or \code{transcriptType}.
#' If the \code{appris} option is set to \code{\link{TRUE}}, only isoforms with some type of appris annotation will be returned.
#' If it is set to a charater string such as 'principal', it will only return isoforms where that value is a substring of the appris tag.
#' The appris ID collumn, if it exisits, is identified by looking for 'appris' (case insensitive) in the annotation collumn names.
#' The other filter, \code{transcriptType}, returns only those isoforms where the transcript type tag contains the value as a subtring.
#' Finally, if \code{annotation} is set to \code{\link{TRUE}}, all annotation collumn are included with the output.
#' If a character or numeric vector is supplies, they will be used to filter the columns.
#'
#' Note that this can also be run with gene level data to return gene IDs based on gene symbol.
#'
#' @param x R data object with stored isoform annotation data; Most typically this is an \code{ExpressionSet} there is support for other datatypes as well.
#' @param isoforms character; A vector of isoform IDs to include in the output. Can be used in combination with with \code{genes}.
#' @param genes character; A vector of gene symbols. Will include isoforms associated with the genes listed in addition to the isoforms listed in the \code{isoforms} option. The \code{symbol} option controls the column used in the gene symbol look up.
#' @param annotation boolean or vector; If set to \code{\link{TRUE}} all annotation will be listed. Numeric or character vectors can be supplied to subset the anntation columns as desired. Default is \code{\link{FALSE}}
#' @param appris boolean or character; If set to TRUE, will return only isoforms with appris annotation. If set to a character string, will restirct isoforms to those with the character value matching a substring of the appris tag. Appris collumn is determined by the first collumn name to containing 'Appris' (case insenstive).
#' @param transcriptType character; Returns only those isoforms where the transcript type collumn has a substring that matches the character value supplied such as 'protein' in 'protein_coding'. The transcript type collumn is determined by the \code{ttype} option.
#' @param symbol character; Column name of the optional gene symbols column in the annotation. The default value is 'GeneSymbol'.
#' @param ttype character; Column name of the optional transcript type column in the annotation. The default value is 'transcript_type'.
#' @param ... additional parameters for S3 variants.
#'
#' @return a vector of isoform IDs or a dataframe of isoform IDs with requested annotation.
#'
#' @examples
#' ToDo<-1
#'
#' @export
#' @seealso \code{\link{isoPlot}}
showIsoforms <- function(x, isoforms=NULL, genes=NULL,annotation=FALSE, appris=NULL,transcriptType=NULL,symbol="GeneSymbol",ttype="transcript_type", ...) {UseMethod("showIsoforms",x)}

#' @export
showIsoforms.default <- function(x, isoforms=NULL, genes=NULL,annotation=FALSE, appris=NULL,transcriptType=NULL,symbol="GeneSymbol",ttype="transcript_type", ...) {
  IsoDat<-x
  foundSymbol<-TRUE
  if(!is.null(genes) & !(symbol %in% colnames(IsoDat))) {
    warning(paste0("Unable to locate gene symbol column ",symbol," in the isoform feature data.\nIsoform selection by gene symbol is disabled."),call.=FALSE)
    genes<-NULL
    foundSymbol<-FALSE
  }
  if(foundSymbol==TRUE & !is.null(genes)) {
    IsoDat<-IsoDat[which(rownames(IsoDat) %in% isoforms | IsoDat[,symbol] %in% genes), ]
  } else if (is.null(isoforms) & (foundSymbol==FALSE | is.null(genes))) {
    stop("Unable to identify valid isoforms from input.",call.=FALSE)
  } else {
    IsoDat<-IsoDat[which(rownames(IsoDat) %in% isoforms), ]
  }
  if(appris!=FALSE & !is.null(appris) & !is.na(appris)) {
    apprisLoc<-grep("appris",colnames(IsoDat),ignore.case = TRUE)
    if(length(apprisLoc)==0) {
      warning("Unable to find a column with 'appris' in the title.\nContinuing without appris filtering. See help for details.", call.=FALSE)
    } else {
      if(appris==TRUE){
        IsoDat<-IsoDat[!is.na(IsoDat[,apprisLoc]),]
      } else {
        IsoDat<-IsoDat[grepl(appris,IsoDat[,apprisLoc]),]
      }
    }
  }
  if(transcriptType!=FALSE & !is.null(transcriptType) & !is.na(transcriptType)) {
    IsoDat<-IsoDat[grep(transcriptType,IsoDat[,ttype], ignore.case = TRUE),]
  }
  if(annotation[1]==FALSE | is.null(annotation[1]) | is.na(annotation[1])) {
    IsoDat<-rownames(IsoDat)
  } else if (annotation[1]!=TRUE) {
    if(length(annotation)>1) {
      IsoDat<-IsoDat[,annotation]
    } else {
      myIsos<-rownames(IsoDat)
      IsoDat<-IsoDat[,annotation]
      if(is.factor(x = IsoDat)) {
        IsoDat<-as.character(IsoDat)
      }
      names(IsoDat)<-myIsos
    }
  }
  if(is.vector(IsoDat)) {
    if(length(IsoDat)==0) {
      stop("No valid isoforms selected after filtering.\nCheck options or see documentaiton for more details.",call.=FALSE)
    }
  } else if (dim(IsoDat)[1]==0) {
    stop("No valid isoforms selected after filtering.\nCheck options or see documentaiton for more details.",call.=FALSE)
  }
  IsoDat
}

#' @importClassesFrom Biobase ExpressionSet
#' @export
showIsoforms.ExpressionSet <- function(x, isoforms=NULL, genes=NULL,annotation=FALSE, appris=FALSE,transcriptType=FALSE,symbol="GeneSymbol",ttype="transcript_type", ...) {
  IsoDat<-fData(x)
  showIsoforms(IsoDat, isoforms=isoforms, genes=genes, annotation=annotation, appris=appris, transcriptType=transcriptType, symbol=symbol, ttype=ttype)
}

#' @importClassesFrom EDASeq SeqExpressionSet
#' @importFrom EDASeq counts normCounts
#' @importFrom Biobase fData
#' @export
showIsoforms.SeqExpressionSet <- function(x, isoforms=NULL, genes=NULL,annotation=FALSE, appris=FALSE,transcriptType=FALSE,symbol="GeneSymbol",ttype="transcript_type", ...) {
  IsoDat<-fData(x)
  showIsoforms(IsoDat, isoforms=isoforms, genes=genes, annotation=annotation, appris=appris, transcriptType=transcriptType, symbol=symbol, ttype=ttype)
}

#' @importClassesFrom DESeq2 DESeqTransform
#' @importFrom SummarizedExperiment rowData
#' @export
showIsoforms.DESeqTransform <- function(x, isoforms=NULL, genes=NULL,annotation=FALSE, appris=FALSE,transcriptType=FALSE,symbol="GeneSymbol",ttype="transcript_type", ...) {
  IsoDat<-rowData(x)
  showIsoforms(IsoDat, isoforms=isoforms, genes=genes, annotation=annotation, appris=appris, transcriptType=transcriptType, symbol=symbol, ttype=ttype)
}

#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment rowData
#' @export
showIsoforms.SummarizedExperiment <- function(x, isoforms=NULL, genes=NULL,annotation=FALSE, appris=FALSE,transcriptType=FALSE,symbol="GeneSymbol",ttype="transcript_type", ...) {
  IsoDat<-rowData(x)
  showIsoforms(IsoDat, isoforms=isoforms, genes=genes, annotation=annotation, appris=appris, transcriptType=transcriptType, symbol=symbol, ttype=ttype)
}

#' @importClassesFrom limma EList
#' @export
showIsoforms.EList <- function(x, isoforms=NULL, genes=NULL,annotation=FALSE, appris=FALSE,transcriptType=FALSE,symbol="GeneSymbol",ttype="transcript_type", ...) {
  IsoDat<-x$genes
  showIsoforms(IsoDat, isoforms=isoforms, genes=genes, annotation=annotation, appris=appris, transcriptType=transcriptType, symbol=symbol, ttype=ttype)
}


#' @title Extract Isoform Expression Data
#' @description Preprocess Isoform expression data for downstream analysis and plotting.
#'
#' @details
#' ToDO
#'
#' @param d R data object with stored isoform annotation data; Most typically this is an \code{ExpressionSet} there is support for other datatypes as well.
#' @param isoforms character; A vector of isoform IDs to include in the output. Can be used in combination with with \code{genes}.
#' @param plotType character; Can be set to "box", "violin, "dot", "bar" or "denisity" for boxplots, violin plots, dot plots, bar plots, and kernal desity plots, respectively.
#' @param group factor or name of factor to be exracted from \code{x} (e.g. \code{\link[Biobase]{pData}}). Used as the primary grouping factor.
#' @param subGroup factor or name of factor to be exracted from \code{x} (e.g. \code{\link[Biobase]{pData}}). Used to subgroup data unless multiple genes are selected in which case \code{subGroup} is ignored.
#' @param highlight factor or name of factor to be exracted from \code{x} (e.g. \code{\link[Biobase]{pData}}). Used to color data points by factor levels. Only valid for graphs with point overlays.
#' @param facet factor or name of factor to be exracted from \code{x} (e.g. \code{\link[Biobase]{pData}}). Split the data into multiple smaller graphs.
#' @param stack factor or name of factor to be exracted from \code{x} (e.g. \code{\link[Biobase]{pData}}). Used for stacked bar plots where both the individual and aggregate values are important. Valid only for bar plots.
#' @param useNormCounts logical; By default \code{genePlot} will try to use normCounts instead of counts in \code{SeqExpressionSets}. Set to FALSE to use raw counts instead, though this will generate a warning about useing non-normalized data.
#' @param appris boolean or character; If set to TRUE, will return only isoforms with appris annotation. If set to a character string, will restirct isoforms to those with the character value matching a substring of the appris tag. Appris collumn is determined by the first collumn name to containing 'Appris' (case insenstive).
#' @param transcriptType character; Returns only those isoforms where the transcript type collumn has a substring that matches the character value supplied such as 'protein' in 'protein_coding'. The transcript type collumn is determined by the \code{ttype} option.
#' @param symbol character; Column name of the optional gene symbols column in the annotation. The default value is 'GeneSymbol'.
#' @param ttype character; Column name of the optional transcript type column in the annotation. The default value is 'transcript_type'.
#' @param ... additional parameters for S3 variants.

#' @return ToDo
#'
#' @examples
#' ToDo<-1
#'
#' @seealso \code{\link{isoPlot}} \code{\link{showIsoforms}} \code{\link{getGeneData}}
getIsoData <- function(d,isoforms=NULL, plotType=plotType, group=NULL, subGroup=NULL, highlight=NULL, facet=NULL, stack=NULL, useNormCounts=TRUE, appris=NULL,transcriptType=NULL,symbol="GeneSymbol",ttype="transcript_type",  ...) {UseMethod("getIsoData",d)}

getIsoData.default <- function(d,isoforms=NULL, plotType=plotType, group=NULL, subGroup=NULL, highlight=NULL, facet=NULL, stack=NULL, useNormCounts=TRUE, appris=NULL,transcriptType=NULL,symbol="GeneSymbol",ttype="transcript_type", ...){
  idOptions<-list(...)
  if(is.vector(d) | is.factor(d)) {
    d<-as.numeric(as.character(d))
  } else {
    d<-as.data.frame(d)
    dataLength<-dim(d)[1]
    if(!is.null(isoforms)){
      if(sum(isoforms %in% rownames(d)) == length(isoforms)) {
        dx<-d[isoforms,]
      } else {
        stop(paste0("Some of the isoforms: ", paste0(isoforms,collapse=", "), " not found in the rownames of the data."),call. =FALSE)
      }
    }
    if(is.data.frame(d)){
      d<-data.frame(t(d))
    }
  }
  gdOptions<-append(list(x=d,plotType=plotType,group=group,subGroup=subGroup,highlight=highlight,facet=facet,stack=stack, rawData=FALSE), idOptions)
  return(do.call("getGeneData", gdOptions))
}

#' @importClassesFrom Biobase ExpressionSet
#' @importFrom purrr map
#' @importFrom magrittr %>%
#' @importFrom Biobase exprs pData fData
getIsoData.ExpressionSet <- function(d,isoforms=NULL, plotType=plotType, group=NULL, subGroup=NULL, highlight=NULL, facet=NULL, stack=NULL, useNormCounts=TRUE, appris=NULL,transcriptType=NULL,symbol="GeneSymbol",ttype="transcript_type", ...){
  idOptions<-list(...)
  isoDat<-exprs(d)[isoforms,]
  if(length(isoforms)>1) {
    isoDat<-t(isoDat)
  }

  #Start processing the factor options. Assuming values are either independent factors or colnames of pData
  if(!is.null(group) & length(group)==1 & any(group %in% colnames(pData(d)))) {
    group<-pData(d)[,group[1]]
  }
  if(!is.null(subGroup) & length(subGroup)==1 & any(subGroup %in% colnames(pData(d)))) {
    subGroup<-pData(d)[,subGroup[1]]
  }
  if(!is.null(highlight) & length(highlight)==1 & any(highlight %in% colnames(pData(d)))) {
    highlight<-pData(d)[,highlight[1]]
  }
  if(!is.null(facet) & length(facet)==1 & any(facet %in% colnames(pData(d)))) {
    facet<-pData(d)[,facet[1]]
  }
  if(!is.null(stack) & length(stack)==1 & any(stack %in% colnames(pData(d)))) {
    stack<-pData(d)[,stack[1]]
  }
  gdOptions<-list(x=isoDat,plotType=plotType,group=group,subGroup=subGroup,highlight=highlight,facet=facet,stack=stack, rawData=FALSE)
  gdOptions<-append(gdOptions,idOptions)
  return(do.call("getGeneData", gdOptions))
}


#' @importClassesFrom EDASeq SeqExpressionSet
#' @importFrom EDASeq counts normCounts
#' @importFrom Biobase pData fData
getIsoData.SeqExpressionSet<- function(d,isoforms=NULL, plotType=plotType, group=NULL, subGroup=NULL, highlight=NULL, facet=NULL, stack=NULL, useNormCounts=TRUE, appris=NULL,transcriptType=NULL,symbol="GeneSymbol",ttype="transcript_type", ...){
  idOptions<-list(...)
  if(sum(!is.na(normCounts(d)))==0 & useNormCounts==TRUE) {
    warning("The normCounts slot in this data set is empty.\nUsing raw count data instead...\n", call. = FALSE)
    useNormCounts<-FALSE
  }
  isoDat<-NULL
  if(useNormCounts==TRUE) {
    isoDat<-normCounts(d)[isoforms,]
  } else {
    isoDat<-counts(d)[isoforms,]
  }

  #Start processing the factor options. Assuming values are either independent factors or colnames of pData
  if (!is.null(group) & length(group) == 1 & any(group %in% colnames(pData(d)))) {
    group <- pData(d)[, group[1]]
  }
  if (!is.null(subGroup) & length(subGroup) == 1 & any(subGroup %in% colnames(pData(d)))) {
    subGroup <- pData(d)[, subGroup[1]]
  }
  if (!is.null(highlight) & length(highlight) == 1 & any(highlight %in% colnames(pData(d)))) {
    highlight <- pData(d)[, highlight[1]]
  }
  if (!is.null(facet) & length(facet) == 1 & any(facet %in% colnames(pData(d)))) {
    facet <- pData(d)[, facet[1]]
  }
  if(!is.null(stack) & length(stack)==1 & any(stack %in% colnames(pData(d)))) {
    stack<-pData(d)[,stack[1]]
  }
  gdOptions<-list(x=isoDat,plotType=plotType,group=group,subGroup=subGroup,highlight=highlight,facet=facet,stack=stack, rawData=FALSE)
  gdOptions<-append(gdOptions,idOptions)
  return(do.call("getGeneData", gdOptions))
}

#' @importClassesFrom DESeq2 DESeqTransform
#' @importFrom SummarizedExperiment assay colData rowData
getIsoData.DESeqTransform<- function(d,isoforms=NULL, plotType=plotType, group=NULL, subGroup=NULL, highlight=NULL, facet=NULL, stack=NULL, useNormCounts=TRUE, appris=NULL,transcriptType=NULL,symbol="GeneSymbol",ttype="transcript_type", ...){
  idOptions<-list(...)
  isoDat<-assay(d)[isoforms,]

#Start processing the factor options. Assuming values are either independent factors or colnames of pData
  if (!is.null(group) & length(group) == 1 & any(group %in% colnames(colData(d)))) {
    group <- colData(d)[, group[1]]
  }
  if (!is.null(subGroup) & length(subGroup) == 1 & any(subGroup %in% colnames(colData(d)))) {
    subGroup <- colData(d)[, subGroup[1]]
  }
  if (!is.null(highlight) & length(highlight) == 1 & any(highlight %in% colnames(colData(d)))) {
    highlight <- colData(d)[, highlight[1]]
  }
  if (!is.null(facet) & length(facet) == 1 & any(facet %in% colnames(colData(d)))) {
    facet <- colData(d)[, facet[1]]
  }
  if(!is.null(stack) & length(stack)==1 & any(stack %in% colnames(colData(d)))) {
    stack<-colData(d)[,stack[1]]
  }
  gdOptions<-list(x=isoDat,plotType=plotType,group=group,subGroup=subGroup,highlight=highlight,facet=facet,stack=stack, rawData=FALSE)
  gdOptions<-append(gdOptions,idOptions)
  return(do.call("getGeneData", gdOptions))
}

#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment assay colData rowData
getIsoData.SummarizedExperiment<- function(d,isoforms=NULL, plotType=plotType, group=NULL, subGroup=NULL, highlight=NULL, facet=NULL, stack=NULL, useNormCounts=TRUE, appris=NULL,transcriptType=NULL,symbol="GeneSymbol",ttype="transcript_type", ...){
  idOptions<-list(...)
  isoDat<-assay(d)[isoforms,]

  #Start processing the factor options. Assuming values are either independent factors or colnames of pData
  if (!is.null(group) & length(group) == 1 & any(group %in% colnames(colData(d)))) {
    group <- colData(d)[, group[1]]
  }
  if (!is.null(subGroup) & length(subGroup) == 1 & any(subGroup %in% colnames(colData(d)))) {
    subGroup <- colData(d)[, subGroup[1]]
  }
  if (!is.null(highlight) & length(highlight) == 1 & any(highlight %in% colnames(colData(d)))) {
    highlight <- colData(d)[, highlight[1]]
  }
  if (!is.null(facet) & length(facet) == 1 & any(facet %in% colnames(colData(d)))) {
    facet <- colData(d)[, facet[1]]
  }
  if(!is.null(stack) & length(stack)==1 & any(stack %in% colnames(colData(d)))) {
    stack<-colData(d)[,stack[1]]
  }
  gdOptions<-list(x=isoDat,plotType=plotType,group=group,subGroup=subGroup,highlight=highlight,facet=facet,stack=stack, rawData=FALSE)
  gdOptions<-append(gdOptions,idOptions)
  return(do.call("getGeneData", gdOptions))
}

#' @importClassesFrom limma EList
getIsoData.EList<- function(d,isoforms=NULL, plotType=plotType, group=NULL, subGroup=NULL, highlight=NULL, facet=NULL, stack=NULL, useNormCounts=TRUE, appris=NULL,transcriptType=NULL,symbol="GeneSymbol",ttype="transcript_type", ...){
  idOptions<-list(...)
  isoDat<-d$E[isoforms,]

  #Start processing the factor options. Assuming values are either independent factors or colnames of pData
  if (!is.null(group) & length(group) == 1 & any(group %in% colnames(d$targets) | group %in% colnames(d$design))) {
    if(any(group %in% colnames(d$design))){
      group <- d$design[,group[1]]
    } else {
      group <- d$targets[,group[1]]
    }
  }
  if (!is.null(subGroup) & length(subGroup) == 1 & any(subGroup %in% colnames(d$targets) | subGroup %in% colnames(d$design))) {
    if(any(subGroup %in% colnames(d$design))){
      subGroup <- d$design[,subGroup[1]]
    } else {
      subGroup <- d$targets[,subGroup[1]]
    }
  }
  if (!is.null(highlight) & length(highlight) == 1 & any(highlight %in% colnames(d$targets) | highlight %in% colnames(d$design))) {
    if(any(highlight %in% colnames(d$design))){
      highlight <- d$design[,highlight[1]]
    } else {
      highlight <- d$targets[,highlight[1]]
    }
  }
  if (!is.null(facet) & length(facet) == 1 & any(facet %in% colnames(d$targets) | facet %in% colnames(d$design))) {
    if(any(facet %in% colnames(d$design))){
      facet <- d$design[,facet[1]]
    } else {
      facet <- d$targets[,facet[1]]
    }
  }
  if(!is.null(stack) & length(stack)==1 & any(stack %in% colnames(d$targets) | stack %in% colnames(d$design))) {
    if(any(stack %in% colnames(d$design))){
      stack <- d$design[,stack[1]]
    } else {
      stack <- d$targets[,stack[1]]
    }
  }
  gdOptions<-list(x=isoDat,plotType=plotType,group=group,subGroup=subGroup,highlight=highlight,facet=facet,stack=stack, rawData=FALSE)
  gdOptions<-append(gdOptions,idOptions)
  return(do.call("getGeneData", gdOptions))
}
