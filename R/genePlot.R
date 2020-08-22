#' @title Plot Gene Expression Data
#' @description Vissualize gene expression data.
#'
#' @details
#' The \code{genePlot} is designed to make vissualizatin of gene expression data simple and easy for R novices and bioinformaticians alike.
#' The function is an S3 generic that accept various R and Bioconductor datasets as input and exracts the expression, factor and annotation data from them according to type.
#' The factors allow for spliting expression data from one or more genes into groups and for plot types with data point overlays, points can be colored by factors levels as well.
#' If the input data is a Bioconductor data set such as an \code{\link[Biobase]{ExpressionSet}} and the \code{gene} option is used, \code{genePlot} will attempt to look up the genes in the associated feature annotation data (e.g. \code{\link[Biobase]{fData}}) according to the data input type and look for the gene symbol column indicated by the \code{symbol} option (defaults to 'GeneSymbol').
#' If no matches are found the row names of are checked of the expression data are check for matches as well.
#' If charactar values are given for factor input, \code{genePlot} will attempt to look up assocated phenotype data (e.g. \code{\link[Biobase]{pData}}).
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
#' @param stack factor or name of factor to be exracted from \code{x} (e.g. \code{\link[Biobase]{pData}}). Used for stacked bar plots where both the individual and aggregate values are important. Valid only for bar plots.
#' @param plotType character; Can be set to "box", "violin, "dot", "bar" or "denisity" for boxplots, violin plots, dot plots, bar plots, and kernal desity plots, respectively.
#' @param main character; The main plot title. Defaults to true for automated generation.
#' @param symbol character; Colname of of gene symbols in the feature data of \code{x} (\code{fData}).
#' @param legend boolean or character; Draws a figure legend. Use to set the legend title which defaults to "Legend" if equals \code{\link{TRUE}}. Set to \code{\link{FALSE}} to disable.
#' @param na.rm logical; Removes \code{\link{NA}} values prior to ploting.
#' @param shiny logical; Use \code{\link[shiny]{shiny}} interfaces if available.
#' @param groupByGene logical; If more then one gene is listed and \code{grouByGene} is \code{TRUE}
#' @param useNormCounts logical; By default \code{genePlot} will try to use normCounts instead of counts in \code{SeqExpressionSets}. Set to FALSE to use raw counts instead, though this will generate a warning about useing non-normalized data.
#' @param ... Any paramenter recognized by \code{NicePlots} functions.
#'
#' @return an list of class \code{npData}. This contains data necessary to regenerate the plot as well as summary statistcs.
#'
#' @examples
#' #While designed for use with bioconductor datasets,
#' #genePlot can also be used for generic data like iris.
#'
#' data(iris)
#'
#' #Using a vector of data
#' genePlot(iris$Sepal.Length,gene=NA, highlight=iris$Species, plotType="dot", pointSize=.75,
#' width=.5, pointShape=1, main="Distribution of Sepal Length")
#'
#' #Kernal denisty plots
#' genePlot(iris$Sepal.Length,gene=NA, group=iris$Species, plotType="density",
#' main="Distribution of Sepal Lengths by Species")
#'
#' #Plotting multiple collumns:
#' genePlot(t(iris[,1:4]),gene=c("Sepal.Length","Petal.Length"), plotType="violin",
#' highlight=iris$Species, pointShape=c(16:18), pointSize=.9)
#'
#' #Multiple collumns with grouping factors:
#' genePlot(t(iris[,1:4]),gene=c("Sepal.Width","Petal.Width"), plotType="bar", group=iris$Species)
#'
#' #Same with grouping order reveresed
#' genePlot(t(iris[,1:4]),gene=c("Sepal.Width","Petal.Width"), plotType="bar", group=iris$Species,
#' groupByGene=FALSE, theme=npColorTheme, errFun="t95ci", legend=TRUE)
#'
#' #2D distribution plotting
#' genePlot(t(iris[,1:4]),gene=c("Sepal.Width","Petal.Width"), plotType="density",
#' group=iris$Species, theme=npGGTheme)
#'
#' #Surface plotting of the above. Use rgl for interactive models.
#' genePlot(t(iris[,1:4]),gene=c("Sepal.Width","Petal.Width"), plotType="surface", legend=TRUE,
#' useRgl=FALSE, theta=60, phi=30, theme=npGGTheme)
#'
#' @importFrom purrr map
#' @importFrom Biobase exprs pData fData
#' @export
#' @seealso \code{\link[NicePlots]{niceBox}}, \code{\link[NicePlots]{niceVio}}, \code{\link[NicePlots]{niceBar}}, \code{\link[NicePlots]{niceDots}}, \code{\link[NicePlots]{niceDensity}}
genePlot <- function(x, gene=NULL, plotType=c("box","dot","bar","violin","density","suface"), symbol="GeneSymbol",legend=NULL, main=TRUE, na.rm=TRUE, group=NULL, subGroup=NULL, highlight=NULL, facet=NULL, stack=NULL, shiny=FALSE, groupByGene=TRUE, useNormCounts=TRUE, ...) {UseMethod("genePlot",x)}

#' @importFrom purrr map
#' @importFrom NicePlots niceBox niceVio niceBar niceDensity
#' @importFrom Biobase exprs pData fData
#' @export

genePlot.default <- function(x, gene=NULL, plotType=c("box","dot","bar","violin","density","surface"), symbol="GeneSymbol", legend=NULL, main=TRUE, na.rm=TRUE, group=NULL, subGroup=NULL, highlight=NULL, facet=NULL, stack=NULL, shiny=FALSE, groupByGene=TRUE, useNormCounts=TRUE, ...) {

  npOptions<-list(...)
  #First lets handle the case that someone set something to FALSE or NA instead of just leaving it as NULL
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

  #Setting default title of nothing was given
  if(main==TRUE) {
    if(length(gene)>1) {
      main<-paste0(c("Gene Expression:",paste0(gene,collapse=", ")),collapse=" ")
    } else {
      main<-paste0(gene, " Expression")
    }
  }

  #Setting the legend to turn on automatically
  if(is.null(legend)){
    legend<-FALSE
    if(!is.null(subGroup) | !is.null(stack)| !is.null(highlight)) {
      legend<-"Legend"
    }
  }
  #Collecting the expression and factor data
  data<-getGeneData(x=x, gene=gene, plotType=plotType, symbol=symbol,group=group, subGroup=subGroup,highlight=highlight,facet=facet, stack=stack, useNormCounts=useNormCounts)

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
  if(is.null(stack)){
    stack<-FALSE
  } else {
    stack<-TRUE
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

#' @importFrom purrr map
#' @importFrom NicePlots niceBox niceVio niceBar niceDensity
#' @importFrom Biobase exprs pData fData
#' @export
genePlot.npData<-function(x, gene=NULL, plotType=NULL, ...) {
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
    dataOut<-do.call("niceBox",x$options)
  } else if (grepl("dot", plotType[1], ignore.case = TRUE)) {
    dataOut<-do.call("niceDots",x$options)
  } else if (grepl("vio", plotType[1], ignore.case = TRUE)) {
    dataOut<-do.call("niceVio",x$options)
  } else if (grepl("bar", plotType[1], ignore.case = TRUE)) {
    dataOut<-do.call("niceBar",x$options)
  } else if (grepl("den",plotType[1], ignore.case = TRUE)) {
    dataOut<-do.call("niceDensity",x$options)
  } else if (grepl("sur", plotType[1], ignore.case = TRUE)) {
    x$options<- append(list(plotType="surface"),x$options)
    dataOut<-do.call("niceDensity",x$options)
  } else {
    stop("invalid plot type")
  }

  invisible(dataOut)
}

