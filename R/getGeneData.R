#' @title Retieve Gene Expression and Factor Data
#' @description Returns selected gene expression and factor data based on data type.
#'
#' @details
#' Todo
#'
#' @param x R data object; Most typically this is an \code{\link[Biobase]{ExpressionSet}}, but there is support for other datatypes as well.
#' @param plotType character; Can be set to "box", "violin, "dot", "bar" or "denisity" for boxplots, violin plots, dot plots, bar plots, and kernal desity plots, respectively.
#' @param gene Gene name, rowname of an expression table actualy vector/matrix of expression. In the case of gene names, the feature annotation element indicated by \code{symbol} is search for matches prior to checking the rownames of the expression table (e.g. \code{\link[Biobase]{assayData}}).
#' @param group factor or name of factor to be exracted from \code{x} (e.g. \code{\link[Biobase]{pData}}). Used as the primary grouping factor.
#' @param subGroup factor or name of factor to be exracted from \code{x} (e.g. \code{\link[Biobase]{pData}}). Used to subgroup data unless multiple genes are selected in which case \code{subGroup} is ignored.
#' @param highlight factor or name of factor to be exracted from \code{x} (e.g. \code{\link[Biobase]{pData}}). Used to color data points by factor levels. Only valid for graphs with point overlays.
#' @param facet factor or name of factor to be exracted from \code{x} (e.g. \code{\link[Biobase]{pData}}). Split the data into multiple smaller graphs.
#' @param stack factor or name of factor to be exracted from \code{x} (e.g. \code{\link[Biobase]{pData}}). Used for stacked bar plots where both the individual and aggregate values are important. Valid only for bar plots.
#' @param symbol character; Colname of of gene symbols in the feature data of \code{x} (e.g. \code{\link[Biobase]{fData}}).
#' @param useNormCounts logical; By default \code{genePlot} will try to use normCounts instead of counts in \code{SeqExpressionSets}. Set to FALSE to use raw counts instead, though this will generate a warning about useing non-normalized data.
#' @param ... Any paramenter recognized by \code{NicePlots} functions.
#'
#' @examples
#' ToDo<-1
#'
#' @importFrom purrr map
#' @importFrom Biobase exprs pData fData
#' @seealso \code{\link[bvt]{genePlot}}, \code{\link[Biobase]{ExpressionSet}}, \code{\link[EDASeq]{SeqExpressionSet-class}}, \code{\link[limma]{EList-class}}, \code{\link[DESeq2]{DESeqTransform}}
getGeneData <- function(x, gene, plotType="box", group=NULL, subGroup=NULL, highlight=NULL, facet=NULL, stack=NULL, symbol="GeneSymbol", useNormCounts=TRUE, ...) {UseMethod("getGeneData",x)}

#strategy here is to have each special bioconductor data type preprocess the data and then call getGeneData again to organize the data in the default method.
#' @importFrom tibble is_tibble
#' @importFrom purrr map
getGeneData.default <- function(x, gene=NULL, plotType="box", group=NULL, subGroup=NULL, highlight=NULL, facet=NULL, stack=NULL, symbol="GeneSymbol", useNormCounts=TRUE, rawData=TRUE,...) {
  #rawData helps with the situation where a raw matrix, tibble or dataframe has been provided and is sent straight to default.
  if(rawData==TRUE) {
    if(is.vector(x) | is.factor(x)) {
      x<-as.numeric(as.character(x))
    } else {
      x<-as.data.frame(x)
      dataLength<-dim(x)[1]
      if(!is.null(gene)){
        if(sum(gene %in% rownames(x)) == length(gene)) {
          x<-x[gene,]
        } else {
          stop(paste0("Some of the genes: ", paste0(gene,collapse=", "), " not found in the rownames of the data."),call. =FALSE)
        }
      }
      if(is.data.frame(x)){
        x<-data.frame(t(x))
      }
    }
  }

  factorData<-list(group=group,subGroup=subGroup,stack,highlight=highlight)
  if(plotType[1]=="bar"){
    highlight<-NULL
  } else {
    stack<-NULL
  }
  sampleN<-0
  #since this is the generic catch all, lets try to standardize the input types.
  if(is.list(x) | is_tibble(x)) {x<-as.data.frame(x)}
  if(is.data.frame(x)){x<-as.matrix(x)}
  #Checking to make sure inputs are sane, first in the case of a matrix, then in the case of a vector.
  if(is.matrix(x)) {
    sampleN<-dim(x)[1]
    #subGrouping is not possible if multiple genes have been selected.
    if(dim(x)[2]>1 & !is.null(subGroup)) {
      if(is.null(group)){
        group<-subGroup
      }
      subGroup<-NULL
    }
  } else {
    sampleN<-length(x)
  }
  if(!is.null(group) & sampleN != length(group)) {
    stop(paste0("Factor data for 'group' is not the same length (",length(group),") of the gene data (",sampleN,").\nPlease check your data inputs.\n",call. = FALSE))
  }
  if(!is.null(subGroup) & sampleN != length(subGroup)) {
    stop(paste0("Factor data for 'subGroup' is not the same length (",length(subGroup),") of the gene data (",sampleN,").\nPlease check your data inputs.\n",call. = FALSE))
  }
  if(!is.null(highlight) & sampleN != length(highlight)) {
    stop(paste0("Factor data for 'highlight' is not the same length (",length(highlight),") of the gene data (",sampleN,").\nPlease check your data inputs.\n",call. = FALSE))
  }
  if(!is.null(stack) & sampleN != length(stack)) {
    stop(paste0("Factor data for 'stack' is not the same length (",length(highlight),") of the gene data (",sampleN,").\nPlease check your data inputs.\n",call. = FALSE))
  }
  NullNames<-FALSE
  if(is.null(group) & is.null(subGroup) ) {
    factorData[[1]]<-factor(rep("Data",sampleN))
    group=TRUE
    NullNames<-TRUE
    if(!is.vector(x) & (!is.null(highlight) | !is.null(stack))){
      temp<-factorData[[1]]
      factorIndex<-4
      if(!is.null(stack) & is.null(highlight)){
        factorIndex<-3
      }
      factorData[[1]]<-factorData[[factorIndex]]
      factorData[[factorIndex]]<-temp
    }
  }
  factorData<-data.frame(factorData[c(!is.null(group),!is.null(subGroup),!is.null(stack),!is.null(highlight))])
  return(list(x=x,by=factorData,facet=facet,NullNames=NullNames))
}

#' @importFrom purrr map
#' @importFrom magrittr %>%
#' @importFrom Biobase exprs pData fData
getGeneData.ExpressionSet <- function(x, gene, plotType="box", group=NULL, subGroup=NULL, highlight=NULL, facet=NULL, stack=NULL, symbol="GeneSymbol", useNormCounts=TRUE, ...) {
  gdOptions<-list(...)
  #issue a quick warning if the data looks like isoform data
  if(grepl("^ENST",rownames(exprs(x))[1]) | grepl("^[NX]{1}[MR]{1}_",rownames(exprs(x))[1])){warning("This appears to be isoform data which is fine as long as each isoform is treated like a gene.\nFor better isoform support, please us isoPlot instead.", call.=FALSE)}
  data<-1
  if(symbol %in% colnames(fData(x))) {
    if(sum(gene %in% fData(x)[,symbol])>0) {
      if(length(gene)>1) {
        data<-map(gene, function(g) exprs(x)[fData(x)[,symbol]==g,]) %>% as.data.frame()
        colnames(data)<-gene
      } else {
        data<-exprs(x)[which(fData(x)[,symbol]==gene),]
      }
    }
  } else if(sum(gene %in% rownames(exprs(x)))>0) {
    if(length(gene)>1) {
      data<-t(exprs(x)[gene,])
    } else {
      data<-exprs(x)[gene,]
    }
  } else {
    stop("unable to identify gene listed in annotation or rownames of the data provided")
  }

  #Start processing the factor options. Assuming values are either independant factors or colnames of pData
  if(!is.null(group) & length(group)==1 & any(group %in% colnames(pData(x)))) {
      group<-pData(x)[,group[1]]
  }
  if(!is.null(subGroup) & length(subGroup)==1 & any(subGroup %in% colnames(pData(x)))) {
    subGroup<-pData(x)[,subGroup[1]]
  }
  if(!is.null(highlight) & length(highlight)==1 & any(highlight %in% colnames(pData(x)))) {
    highlight<-pData(x)[,highlight[1]]
  }
  if(!is.null(facet) & length(facet)==1 & any(facet %in% colnames(pData(x)))) {
    facet<-pData(x)[,facet[1]]
  }
  if(!is.null(stack) & length(stack)==1 & any(stack %in% colnames(pData(x)))) {
    stack<-pData(x)[,stack[1]]
  }
  gdOptions<-append(list(x=data,gene=gene,plotType=plotType,group=group,subGroup=subGroup,highlight=highlight,facet=facet,stack=stack, rawData=FALSE),gdOptions)
  return(do.call("getGeneData", gdOptions))
}


#' @importClassesFrom EDASeq SeqExpressionSet
#' @importFrom purrr map
#' @importFrom magrittr %>%
#' @importFrom EDASeq counts normCounts
#' @importFrom Biobase pData fData
getGeneData.SeqExpressionSet <- function(x, gene, plotType="box", group=NULL, subGroup=NULL, highlight=NULL, facet=NULL, stack=NULL, symbol="GeneSymbol", useNormCounts=TRUE, ...) {
  gdOptions<-list(...)
  #Handling the situation where normCounts has not been initialized.
  if(sum(!is.na(normCounts(x)))==0 & useNormCounts==TRUE) {
    warning("The normCounts slot in this data set is empty.\nUsing raw count data instead...\n", call. = FALSE)
    useNormCounts<-FALSE
  }
  #issue a quick warning if the data looks like isoform data
  if(grepl("^ENST",rownames(counts(x))[1]) | grepl("^[NX]{1}[MR]{1}_",rownames(counts(x))[1])){warning("This appears to be isoform data which is fine as long as each isoform is treated like a gene.\nFor better isoform support, please us isoPlot instead.", call.=FALSE)}
  data<-1
  if(symbol %in% colnames(fData(x))) {
    if(sum(gene %in% fData(x)[,symbol])>0) {
      if(length(gene)>1) {
        if(useNormCounts==TRUE) {
          data<-map(gene, function(g) normCounts(x)[fData(x)[,symbol]==g,]) %>% as.data.frame()
        } else {
          data<-map(gene, function(g) counts(x)[fData(x)[,symbol]==g,]) %>% as.data.frame()
        }
        colnames(data)<-gene
      } else {
        if(useNormCounts==TRUE) {
          data<-normCounts(x)[which(fData(x)[,symbol]==gene),]
        } else {
          data<-counts(x)[which(fData(x)[,symbol]==gene),]
        }
      }
    }
  } else if(sum(gene %in% rownames(counts(x)))>0) {
    if(length(gene)>1) {
      if(useNormCounts==TRUE) {
        data<-t(normCounts(x)[gene,])
      } else {
        data<-t(counts(x)[gene,])
      }
    } else {
      if(useNormCounts==TRUE) {
        data<-normCounts(x)[gene,]
      } else {
        data<-counts(x)[gene,]
      }
    }
  } else {
    stop("unable to identify gene listed in annotation or rownames of the data provided")
  }
  #Start processing the factor options. Assuming values are either independant factors or colnames of pData
  if(!is.null(group) & length(group)==1 & any(group %in% colnames(pData(x)))) {
    group<-pData(x)[,group[1]]
  }
  if(!is.null(subGroup) & length(subGroup)==1 & any(subGroup %in% colnames(pData(x)))) {
    subGroup<-pData(x)[,subGroup[1]]
  }
  if(!is.null(highlight) & length(highlight)==1 & any(highlight %in% colnames(pData(x)))) {
    highlight<-pData(x)[,highlight[1]]
  }
  if(!is.null(facet) & length(facet)==1 & any(facet %in% colnames(pData(x)))) {
    facet<-pData(x)[,facet[1]]
  }
  if(!is.null(stack) & length(stack)==1 & any(stack %in% colnames(pData(x)))) {
    stack<-pData(x)[,stack[1]]
  }
  gdOptions<-append(list(x=data,gene=gene,plotType=plotType,group=group,subGroup=subGroup,highlight=highlight,facet=facet,stack=stack, rawData=FALSE),gdOptions)
  do.call("getGeneData", gdOptions)
}


#' @importClassesFrom DESeq2 DESeqTransform
#' @importFrom SummarizedExperiment assay colData rowData
getGeneData.DESeqTransform <- function(x, gene, plotType="box", group=NULL, subGroup=NULL, highlight=NULL, facet=NULL, stack=NULL, symbol="GeneSymbol", useNormCounts=TRUE, ...) {
  gdOptions<-list(...)
  #issue a quick warning if the data looks like isoform data
  if(grepl("^ENST",rownames(assay(x))[1]) | grepl("^[NX]{1}[MR]{1}_",rownames(assay(x))[1])){warning("This appears to be isoform data which is fine as long as each isoform is treated like a gene.\nFor better isoform support, please us isoPlot instead.", call.=FALSE)}
  data<-1
  if(!is.null(rowData(x)) & symbol %in% colnames(rowData(x))) {
    if(sum(gene %in% rowData(x)[,symbol])>0) {
      if(length(gene)>1) {
        data<-map(gene, function(g) assay(x)[rownames(rowData(x))[rowData(x)[,symbol]==g],]) %>% as.data.frame()
        colnames(data)<-gene
      } else {
        data<-assay(x)[which(rowData(x)[,symbol]==gene),]
      }
    }
  } else if(sum(gene %in% rownames(assay(x)))>0) {
    if(length(gene)>1) {
      data<-t(assay(x)[gene,])
    } else {
      data<-assay(x)[gene,]
    }
  } else {
    stop("unable to identify gene listed in annotation or rownames of the data provided")
  }

  #Start processing the factor options. Assuming values are either independant factors or colnames of pData
  if(!is.null(group) & length(group)==1 & any(group %in% colnames(colData(x)))) {
    group<-colData(x)[,group[1]]
  }
  if(!is.null(subGroup) & length(subGroup)==1 & any(subGroup %in% colnames(colData(x)))) {
    subGroup<-colData(x)[,subGroup[1]]
  }
  if(!is.null(highlight) & length(highlight)==1 & any(highlight %in% colnames(colData(x)))) {
    highlight<-colData(x)[,highlight[1]]
  }
  if(!is.null(facet) & length(facet)==1 & any(facet %in% colnames(colData(x)))) {
    facet<-colData(x)[,facet[1]]
  }
  if(!is.null(stack) & length(stack)==1 & any(stack %in% colnames(colData(x)))) {
    stack<-colData(x)[,stack[1]]
  }
  gdOptions<-append(list(x=data,gene=gene,plotType=plotType,group=group,subGroup=subGroup,highlight=highlight,facet=facet,stack=stack, rawData=FALSE),gdOptions)
  do.call("getGeneData", gdOptions)
}


#' @importFrom purrr map
#' @importFrom magrittr %>%
#' @importFrom SummarizedExperiment assay colData rowData
getGeneData.SummarizedExperiment <- function(x, gene, plotType="box", group=NULL, subGroup=NULL, highlight=NULL, facet=NULL, stack=NULL, symbol="GeneSymbol", useNormCounts=TRUE, ...) {
  gdOptions<-list(...)
  #issue a quick warning if the data looks like isoform data
  if(grepl("^ENST",rownames(assay(x))[1]) | grepl("^[NX]{1}[MR]{1}_",rownames(assay(x))[1])){warning("This appears to be isoform data which is fine as long as each isoform is treated like a gene.\nFor better isoform support, please us isoPlot instead.", call.=FALSE)}
  data<-1
  if(!is.null(rowData(x)) & symbol %in% colnames(rowData(x))) {
    if (sum(gene %in% rowData(x)[,symbol])>0) {
      if(length(gene)>1) {
        data<-map(gene, function(g) assay(x)[rownames(rowData(x))[rowData(x)[,symbol]==g],]) %>% as.data.frame()
        colnames(data)<-gene
      } else {
        data<-assay(x)[which(rowData(x)[,symbol]==gene),]
      }
    }
  } else if(sum(gene %in% rownames(assay(x)))>0) {
    if(length(gene)>1) {
      data<-t(assay(x)[gene,])
    } else {
      data<-assay(x)[gene,]
    }
  } else {
    stop("unable to identify gene listed in annotation or rownames of the data provided")
  }

  #Start processing the factor options. Assuming values are either independant factors or colnames of pData
  if(!is.null(group) & length(group)==1 & any(group %in% colnames(colData(x)))) {
    group<-colData(x)[,group[1]]
  }
  if(!is.null(subGroup) & length(subGroup)==1 & any(subGroup %in% colnames(colData(x)))) {
    subGroup<-colData(x)[,subGroup[1]]
  }
  if(!is.null(highlight) & length(highlight)==1 & any(highlight %in% colnames(colData(x)))) {
    highlight<-colData(x)[,highlight[1]]
  }
  if(!is.null(facet) & length(facet)==1 & any(facet %in% colnames(colData(x)))) {
    facet<-colData(x)[,facet[1]]
  }
  if(!is.null(stack) & length(stack)==1 & any(stack %in% colnames(colData(x)))) {
    stack<-colData(x)[,stack[1]]
  }
  gdOptions<-append(list(x=data,gene=gene,plotType=plotType,group=group,subGroup=subGroup,highlight=highlight,facet=facet,stack=stack, rawData=FALSE),gdOptions)
  do.call("getGeneData", gdOptions)

}

#' @importClassesFrom limma EList
#' @importFrom purrr map
#' @importFrom Biobase exprs pData fData
getGeneData.EList <- function(x, gene, plotType="box", group=NULL, subGroup=NULL, highlight=NULL, facet=NULL, stack=NULL, symbol="GeneSymbol", useNormCounts=TRUE, ...) {
  gdOptions<-list(...)
  #issue a quick warning if the data looks like isoform data
  if(grepl("^ENST",rownames(x$E)[1]) | grepl("^[NX]{1}[MR]{1}_",rownames(x$E)[1])){warning("This appears to be isoform data which is fine as long as each isoform is treated like a gene.\nFor better isoform support, please us isoPlot instead.", call.=FALSE)}
  data<-1
  if(!is.null(x$genes) & symbol %in% colnames(x$genes)) {
    if (sum(gene %in% x$genes[,symbol])>0) {
      if(length(gene)>1) {
        data<-map(gene, function(g) x$E[rownames(x$genes)[x$genes[,symbol]==g],]) %>% as.data.frame()
        colnames(data)<-gene
      } else {
        data<-x$E[which(x$genes[,symbol]==gene),]
      }
    }
  } else if(sum(gene %in% rownames(x$E))>0) {
    if(length(gene)>1) {
      data<-t(x$E[gene,])
    } else {
      data<-x$E[gene,]
    }
  } else {
    stop("unable to identify gene listed in annotation or rownames of the data provided")
  }

  #Start processing the factor options. Assuming values are either independant factors or colnames of pData
  if(!is.null(group) & length(group)==1 & (any(group %in% colnames(x$design)) | any(group %in% colnames(x$targets)))) {
    if (any(group %in% colnames(x$design))) {
      group<-x$design[,group[1]]
    } else {
      group<-x$targets[,group[1]]
    }
  }
  if(!is.null(subGroup) & length(subGroup)==1 & (any(subGroup %in% colnames(x$design)) | any(subGroup %in% colnames(x$targets)))) {
    if (any(subGroup %in% colnames(x$design))) {
      subGroup<-x$design[,subGroup[1]]
    } else {
      subGroup<-x$targets[,subGroup[1]]
    }
  }
  if(!is.null(highlight) & length(highlight)==1 & (any(highlight %in% colnames(x$design)) | any(highlight %in% colnames(x$targets)))) {
    if (any(highlight %in% colnames(x$design))) {
      highlight<-x$design[,highlight[1]]
    } else {
      highlight<-x$targets[,highlight[1]]
    }
  }
  if(!is.null(facet) & length(facet)==1 & (any(facet %in% colnames(x$design)) | any(facet %in% colnames(x$targets)))) {
    if (any(facet %in% colnames(x$design))) {
      facet<-x$design[,facet[1]]
    } else {
      facet<-x$targets[,facet[1]]
    }
  }
  if(!is.null(stack) & length(stack)==1 & (any(stack %in% colnames(x$design)) | any(stack %in% colnames(x$targets)))) {
    if (any(stack %in% colnames(x$design))) {
      stack<-x$design[,stack[1]]
    } else {
      stack<-x$targets[,stack[1]]
    }
  }
  gdOptions<-append(list(x=data,gene=gene,plotType=plotType,group=group,subGroup=subGroup,highlight=highlight,facet=facet,stack=stack, rawData=FALSE),gdOptions)
  do.call("getGeneData", gdOptions)
}

