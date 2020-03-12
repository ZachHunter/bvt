
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

#' @export
showIsoforms.EList <- function(x, isoforms=NULL, genes=NULL,annotation=FALSE, appris=FALSE,transcriptType=FALSE,symbol="GeneSymbol",ttype="transcript_type", ...) {
  IsoDat<-x$genes
  showIsoforms(IsoDat, isoforms=isoforms, genes=genes, annotation=annotation, appris=appris, transcriptType=transcriptType, symbol=symbol, ttype=ttype)
}


#' @title Extract Isoform Exrpression Data
#' @description Preprocess Isoform expression data for downstrean analysis and plotting.
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
  gdOptions<-list(x=d,plotType=plotType,group=group,subGroup=subGroup,highlight=highlight,facet=facet,stack=stack, rawData=FALSE)
  return(do.call("getGeneData", gdOptions))
}

getIsoData.ExpressionSet <- function(d,isoforms=NULL, plotType=plotType, group=NULL, subGroup=NULL, highlight=NULL, facet=NULL, stack=NULL, useNormCounts=TRUE, appris=NULL,transcriptType=NULL,symbol="GeneSymbol",ttype="transcript_type", ...){
  isoDat<-exprs(d)[isoforms,]
  if(length(isoforms)>1) {
    isoDat<-t(isoDat)
  }

  #Start processing the factor options. Assuming values are either independant factors or colnames of pData
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
  return(do.call("getGeneData", gdOptions))
}

