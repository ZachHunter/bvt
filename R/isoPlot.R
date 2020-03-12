#' @title Plot Isoform Expression Data
#' @description Vissualize isoform expression data.
#'
#' @details
#' The \code{isoPlot} is designed to make vissualizatin of isoform expression data simple and easy for R novices and bioinformaticians alike.
#' The function is an S3 generic that accept various R and Bioconductor datasets as input and exracts the expression, factor and annotation data from them according to type.
#' The factors allow for spliting expression data from one or more genes into groups and for plot types with data point overlays. Points can be colored by factors levels as well.
#' If the input data is a Bioconductor data set such as an \code{\link[Biobase]{ExpressionSet}} and the \code{gene} option is used, \code{isoPlot} will attempt to look up the isoforms in the associated with the gene in the annotation data (e.g. \code{\link[Biobase]{fData}}) according to the data input type and look for the gene symbol column indicated by the \code{symbol} option (defaults to 'GeneSymbol').
#' If no matches are found the row names of are checked of the expression data are check for matches as well.
#' If charactar values are given for factor input, \code{isoPlot} will attempt to look up assocated phenotype data (e.g. \code{\link[Biobase]{pData}}).
#' One can also pass raw data vectors/data frames and/or factors to \code{isoPlots} to bypass this feature, which is critical for data sets and data formats where integrated phenotype and feature data is not available.
#' The \code{isoPlot} uses the \code{NicePlots} graphics library and any \code{NicePlots} option and/or theme can be used in conjuction with options detailed below.
#' The \code{plotType} options supported correspond to \code{NicePlots} functions and include box plots (\code{\link[NicePlots]{niceBox}}), dot plots (\code{\link[NicePlots]{niceDots}}), violin plots (\code{\link[NicePlots]{niceVio}}), bar plots (\code{\link[NicePlots]{niceBar}}) as well as both one/two dimentional kernal density plots (\code{\link[NicePlots]{niceDensity}}).
#' Supported data input types include: \code{\link[Biobase]{ExpressionSet}}, \code{\link[EDASeq]{SeqExpressionSet-class}}, \code{\link[limma]{EList-class}}, \code{\link[DESeq2]{DESeqTransform}}, as well as standard R data types such as \code{\link[base]{vector}}, \code{\link[base]{matrix}}, \code{\link[base]{data.frame}}, and \code{\link[tibble]{tibble}}.
#' \code{isoPlot} silently returns a list of class \code{npData} that conatains a summarized findings, p-values (if indicated), extracted plotting data, and plotting options.
#' All npData objects can be replotted using  the \code{\link[graphics]{plot}} function, \code{isoPlot} or any of the \code{NicePlots} functions.
#' Options passed to any of these, including \code{plotType} will override the options for the \code{npData} object.
#'
#' @param x R data object; Most typically this is an \code{ExpressionSet} there is support for other datatypes as well.
#' @param isoforms character; Isoform IDs or a vector of isoform IDS to plot.
#' @param gene character; Gene or vector of gene names. This is an optional setting that will return all of the isoforms associated with the gene.
#' @param appris boolean or character; If set to TRUE, will return only isoforms with appris annotation. If set to a character string, will restirct isoforms to those with the character value matching a substring of the appris tag. Appris collumn is determined by the first collumn name to containing 'Appris' (case insenstive).
#' @param transcriptType character; Returns only those isoforms where the transcript type collumn has a substring that matches the character value supplied such as 'protein' in 'protein_coding'. The transcript type collumn is determined by the \code{ttype} option.
#' @param asPercentage boolean; If set to \code{\link{TRUE}}, the isoform expression is given as a percetange of total gene expression (defaults to \code{\link{FALSE}})
#' @param group factor or name of factor to be exracted from \code{x} (e.g. \code{\link[Biobase]{pData}}). Used as the primary grouping factor.
#' @param subGroup factor or name of factor to be exracted from \code{x} (e.g. \code{\link[Biobase]{pData}}). Used to subgroup data unless multiple genes are selected in which case \code{subGroup} is ignored.
#' @param highlight factor or name of factor to be exracted from \code{x} (e.g. \code{\link[Biobase]{pData}}). Used to color data points by factor levels. Only valid for graphs with point overlays.
#' @param facet factor or name of factor to be exracted from \code{x} (e.g. \code{\link[Biobase]{pData}}). Split the data into multiple smaller graphs.
#' @param stack factor or name of factor to be exracted from \code{x} (e.g. \code{\link[Biobase]{pData}}). Used for stacked bar plots where both the individual and aggregate values are important. Valid only for bar plots.
#' @param plotType character; Can be set to "box", "violin, "dot", "bar" or "denisity" for boxplots, violin plots, dot plots, bar plots, and kernal desity plots, respectively.
#' @param main character; The main plot title. Defaults to true for automated generation.
#' @param symbol character; Colname of of gene symbols in the feature data of \code{x} (\code{fData}).
#' @param legend boolean or character; Draws a figure legend. Use to set the legend title which defaults to "Legend" if equals \code{\link{TRUE}}. Set to \code{\link{FALSE}} to disable.
#' @param na.rm logical; Removes \code{\link{NA}} values prior to ploting.
#' @param shiny logical; Use \code{\link[shiny]{shiny}} interfaces if available.
#' @param groupByGene logical; If more then one gene is listed and \code{grouByGene} is \code{TRUE}
#' @param useNormCounts logical; By default \code{genePlot} will try to use normCounts instead of counts in \code{SeqExpressionSets}. Set to FALSE to use raw counts instead, though this will generate a warning about useing non-normalized data.
#' @param ttype character; Column name of the optional transcript type column in the annotation. The default value is 'transcript_type'.
#' @param ... Any paramenter recognized by \code{NicePlots} functions.
#'
#' @return an list of class \code{npData}. This contains data necessary to regenerate the plot as well as summary statistcs.
#'
#' @examples
#' ToDo<-1
#'
#' @importFrom purrr map
#' @importFrom tidyr gather
#' @importFrom Biobase exprs pData fData
#' @export
#' @seealso \code{\link{genePlot}}, \code{\link{showIsoforms}}, \code{\link[NicePlots]{niceBox}}, \code{\link[NicePlots]{niceVio}}, \code{\link[NicePlots]{niceBar}}, \code{\link[NicePlots]{niceDots}}, \code{\link[NicePlots]{niceDensity}}
isoPlot <- function(x, isoforms=NULL, gene=NULL, plotType=c("box","dot","bar","violin","density","suface"), asPercentage=FALSE, symbol="GeneSymbol",legend=NULL, main=TRUE, na.rm=TRUE, group=NULL, subGroup=NULL, highlight=NULL, facet=NULL, stack=NULL, shiny=FALSE, groupByGene=FALSE, useNormCounts=TRUE, appris=FALSE, transcriptType=FALSE, ttype="transcript_type",...) {UseMethod("isoPlot",x)}

#' @importFrom purrr map
#' @importFrom tidyr gather
#' @importFrom NicePlots niceBox niceVio niceBar niceDensity
#' @importFrom Biobase exprs pData fData
#' @export
isoPlot.default <- function(x, isoforms=NULL, gene=NULL, plotType=c("box","dot","bar","violin","density","suface"), asPercentage=FALSE, symbol="GeneSymbol", legend=NULL, main=TRUE, na.rm=TRUE, group=NULL, subGroup=NULL, highlight=NULL, facet=NULL, stack=NULL, shiny=FALSE, groupByGene=FALSE, useNormCounts=TRUE, appris=FALSE, transcriptType=FALSE, ttype="transcript_type", ...) {

  npOptions<-list(...)
  #First lets handle the case that someone set something to FALSE or NA instead of just leaving it as NULL

  if(sum(isoforms==FALSE)==1 | sum(is.na(isoforms))==1) {isoforms<-NULL}
  if(sum(gene==FALSE)==1 | sum(is.na(gene))==1) {gene<-NULL}
  if((length(group)==1 & sum(group==FALSE)==1) | sum(is.na(group))==length(group)) {group<-NULL}
  if((length(subGroup)==1 & sum(subGroup==FALSE)==1) | sum(is.na(subGroup))==length(subGroup)) {subGroup<-NULL}
  if((length(stack)==1 & sum(stack==FALSE)==1) | sum(is.na(stack))==length(stack)) {stack<-NULL}
  if((length(highlight)==1 & sum(highlight==FALSE)==1) | sum(is.na(highlight))==length(highlight)) {highlight<-NULL}

  #Now lets get rid of incompatible options
  if(!is.null(highlight) & plotType[1]=="bar") {
    hightlight<-NULL
  }
  if(!is.null(stack) & plotType[1]!="bar") {
    stack<-NULL
  }

  #Looking up which isoforms meet the selected criteria
  isos<-NULL
  if(is.null(isoforms) & is.null(gene)){
    isos<-rownames(x)
  } else {
    isos<-showIsoforms(x,isoforms=isoforms,genes = gene, annotation = FALSE, appris = appris,transcriptType = transcriptType, symbol = symbol,ttype = ttype)
  }

  #Setting default title of nothing was given
  if(main==TRUE) {
     if(length(gene)>1) {
       main<-paste0(c(paste0(gene,collapse=", "),"Isoform Expression"),collapse=" ")
     } else if (!is.null(isoforms)) {
       main<-paste0(c(paste0(isoforms,collapse=", "),"Expression"),collapse=" ")
     } else if (!is.null(gene)) {
       main<-paste0(gene, " Expression")
     } else {
       "Isoform Expression"
     }
  }

  #Setting the legend to turn on automatically
  if(is.null(legend)){
    legend<-FALSE
    if(!is.null(subGroup) | !is.null(stack)| !is.null(highlight)) {
      legend<-"Legend"
    }
  }

  #Stack can be set to true if the goal is stacking isoforms.
  #This is handled separately from the getIsoData
  #This section of code just hides this use case from getIsoData.
  isoStack<-FALSE
  if(!is.null(stack)) {
    if(stack[1]==TRUE) {
      isoStack<-TRUE
      stack<-NULL
    }
  }

  #Quick test to see if searching by gene symbol is available
  SymbolFound<-FALSE
  if(symbol %in% colnames(showIsoforms(x, isoforms=isos, annotation = T))) {
    SymbolFound<-TRUE
  }
  #Collecting the expresion and factor data
  data<-getIsoData(d=x, isoforms=isos, plotType=plotType, symbol=symbol,group=group, subGroup=subGroup,highlight=highlight,facet=facet, stack=stack, useNormCounts=useNormCounts)

  #Convert isoforms as a percentage of gene expression.
  if(asPercentage==TRUE & SymbolFound==TRUE) {
    myGenes<-showIsoforms(x,isoforms=isos,symbol=symbol,annotation=symbol)
    uniGenes<-unique(myGenes)
    gexprs<-vector(mode = "list", length = length(uniGenes))
    names(gexprs)<-uniGenes
    for(cgene in uniGenes) {
      cisos<-showIsoforms(x,genes=cgene,annotation = FALSE)
      cDat<-getIsoData(d=x, isoforms=cisos, plotType=plotType, symbol=symbol,group=group, subGroup=subGroup,highlight=highlight,facet=facet, stack=stack, useNormCounts=useNormCounts)
      if(length(cisos)==1) {
        gexprs[[cgene]]<-cDat$x
      } else {
        gexprs[[cgene]]<-rowSums(cDat$x)
      }
      gexprs[[cgene]][which(gexprs[[cgene]]==0)]<-1 #Avoiding divide by zero errors. Iso will be zero anyway.
    }
    if (length(isos)==1) {
      data$x<-data$x/gexprs[[1]]*100
    } else {
      for (i in 1:dim(data$x)[2]) {
        data$x[,i]<-data$x[,i]/gexprs[[myGenes[i]]]*100
      }
    }
    npOptions<-append(npOptions,list("axisText"=c("","%")))
  }
  if(isoStack==TRUE){
    stack<-TRUE
  }
  if(is.null(stack)){
    stack<-FALSE
  } else {
    #If stack is TRUE we are assuming the user means to stack isoforms
    if(stack[1]==TRUE) {
      stackData<-NA
      if(is.vector(data$x)) {
        #There there is only one isoform then there is nothing to stack
        stack<-FALSE
      } else {
        stackData<-data.frame(data$x,data$by) %>%
          gather(key="isoforms",value="exprs",colnames(data$x))
        data$x<-stackData[,"exprs"]
        data$by<-stackData[,seq_len(dim(stackData)[2]-1)]
      }
    }
    stack<-TRUE
  }

  #if group or subgroup are left empty and gene symbol annotation is available
  #and if there are isoforms from more than one gene present, we will add a gene symbol factor level automatically
  if(SymbolFound==TRUE){
    iso2gene<-showIsoforms(x, isoforms = isos, symbol=symbol, annotation = symbol)
    if(length(unique(iso2gene))>1 & sum(c(is.null(group),is.null(subGroup)))>=1 & grepl("bar",plotType,ignore.case = TRUE) & isoStack==TRUE) {
      if("isoforms" %in% colnames(data$by) ){
        geneFact<-iso2gene[data$by$isoforms]
        print(head(data$by))
        if(sum(data$by$group =="data" | is.na(data$by$group)) == length(data$by$group)){
          data$by$group <- factor(geneFact)
        } else {
          data$by <- data.frame(geneFact,data$by)
        }
        subGroup<-TRUE
      }
    } else if (length(unique(iso2gene))>1 & sum(c(is.null(group),is.null(subGroup)))==2) {
      geneData<-data.frame(data$x,data$by) %>%
        gather(key="isoforms",value="exprs",colnames(data$x)) %>%
        select("isoforms", colnames(data$by),"exprs")
      data$x<-geneData[,"exprs"]
      geneFact<-iso2gene[geneData$isoforms]
      if(sum(data$by$group =="data", na.rm=TRUE) == length(data$by$group)-sum(is.na(data$by$group))){
        geneData$group <- factor(geneFact)
        data$by<-geneData[,seq_len(dim(geneData)[2]-1)]
      } else {
        data$by <- data.frame(geneFact,geneData[,seq_len(dim(geneData)[2]-1)])
      }
      print(head(data$by))
      subGroup<-TRUE
    }
  }

  #Now we convert the options to boolean TRUE/FALSE for compatibility with NicePlots
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

  if(!is.vector(data$x) & (!is.null(group) | !is.null(subGroup))) {
    subGroup<-TRUE
  }
  if(is.null(group) & subGroup==TRUE) {
    subGroup<-FALSE
  }
  if(plotType[1]=="density" & !is.null(group)) {
    subGroup<-TRUE
  }

  #Formatting options and adding new data
  npOptions<-append(list(x=data$x,by=data$by,pointHighlights=highlight,flipFacts=groupByGene, subGroup=subGroup, facet=facet,stack=stack, na.rm=na.rm,main=main, legend=legend),npOptions)
  if(groupByGene==TRUE & data$NullNames==TRUE) {
    if(is.factor(data$by)) {
      npOptions<-append(npOptions,list(subGroupLabels=rep("",length(levels(data$by)))))
    } else {
      npOptions<-append(npOptions,list(subGroupLabels=rep("",length(levels(data$by[,1])))))
    }
  }
  #Calling NicePlots
  dataOut<-1
  if(grepl("box", plotType[1], ignore.case = TRUE)){
    dataOut<-do.call("niceBox",npOptions)
  } else if (grepl("dot", plotType[1], ignore.case = TRUE)) {
    dataOut<-do.call("niceDots",npOptions)
  } else if (grepl("vio", plotType[1], ignore.case = TRUE)) {
    dataOut<-do.call("niceVio",npOptions)
  } else if (grepl("bar", plotType[1], ignore.case = TRUE)) {
    dataOut<-do.call("niceBar",npOptions)
  } else if (grepl("den",plotType[1], ignore.case = TRUE)) {
    dataOut<-do.call("niceDensity",npOptions)
  } else if (grepl("sur", plotType[1], ignore.case = TRUE)) {
    npOptions<- append(list(plotType="surface"),npOptions)
    dataOut<-do.call("niceDensity",npOptions)
  } else {
    stop("invalid plot type")
  }
  invisible(dataOut)
}



#' @importFrom purrr map
#' @importFrom NicePlots niceBox niceVio niceBar niceDensity
#' @importFrom Biobase exprs pData fData
#' @export
isoPlot.npData<-function(x, isoforms=NULL, gene=NULL, plotType=NULL, ...) {
  clOptions<-list(...)
  for(opt in names(clOptions)) {
    if(is.null(x$options[opt])){
      append(x$options,list(opt=clOptions[[opt]]))
    }else{
      x$options[[opt]]<-clOptions[[opt]]
    }
  }
  if(!is.null(x$options[["groupByGene"]])){
    if(x$options[["groupByGene"]]==TRUE) {
      x$options[["flipFacts"]]<-FALSE
    } else {
      x$options[["flipFacts"]]<-TRUE
    }
  }
  dataOut<-1
  if(grepl("box", plotType[1], ignore.case = TRUE)){
    dataOut<-do.call("niceBox",npOptions)
  } else if (grepl("dot", plotType[1], ignore.case = TRUE)) {
    dataOut<-do.call("niceDots",npOptions)
  } else if (grepl("vio", plotType[1], ignore.case = TRUE)) {
    dataOut<-do.call("niceVio",npOptions)
  } else if (grepl("bar", plotType[1], ignore.case = TRUE)) {
    dataOut<-do.call("niceBar",npOptions)
  } else if (grepl("den",plotType[1], ignore.case = TRUE)) {
    dataOut<-do.call("niceDensity",npOptions)
  } else if (grepl("sur", plotType[1], ignore.case = TRUE)) {
    npOptions<- append(list(plotType="surface"),npOptions)
    dataOut<-do.call("niceDensity",npOptions)
  } else {
    stop("invalid plot type")
  }

  invisible(dataOut)
}
