#' @title Make a Gene Level Scatter Plot
#' @description Visualize gene expression data in 1, 2, and 3 dimension scatter/waterfall plots for exploratory data analysis.
#'
#' @details
#' Details will be forthcoming. Makes a 2D or 3D scatter plot.
#'
#' @param x R data object; Most typically this is an \code{ExpressionSet} there is support for other data types as well including \code{\link{matrix}} and \code{\link{data.frame}}.
#' @param genes character vector; Names of genes or or phenotype data that can be coerced into a numeric vector.
#' @param color vector/character; Should either be factor of values used to color points or the name of a gene/phenotype that can be used to construct one.
#' @param shape vector/character; Should either be factor of values used to control the shape (i.e. \code{pty}) of points or the name of a gene/phenotype that can be used to construct one.
#' @param size vector/character; Should either be factor of values used to control the shape (i.e. \code{pty}) of points or the name of a gene/phenotype that can be used to construct one.
#' @param trendline character; Valid options include \code{color}, \code{shape} or \code{density}.See details for more information. Setting to \code{\link{TRUE}} will cause the first grouping factor to be used.
#' @param facet factor or name of factor to be exacted from \code{x} (e.g. \code{\link[Biobase]{pData}}). Split the data into multiple smaller graphs.
#' @param main character; The main plot title. Defaults to true for automated generation.
#' @param symbol character; Column name of of gene symbols in the feature data of \code{x} (\code{fData}).
#' @param legend logical or character; Draws a figure legend. Use to set the legend title which defaults to "Legend" if equals \code{\link{TRUE}}. Set to \code{\link{FALSE}} to disable.
#' @param theme npTheme object; A valid npTheme object the controls default settings.
#' @param na.rm logical; Removes \code{\link{NA}} values prior to plotting.
#' @param shiny logical; Use \code{\link[shiny]{shiny}} interfaces if available.
#' @param isTidy logical; Transposes input data if set to \code{\link{TRUE}}. Biological expression data is often formatted with genes as rows and samples as columns which is the transpose of the more standard tidy data format. You can read more about tidy data at \code{vignette("tidy-data",package = "tidyr")}. Defaults to FALSE unless input data is a \code{\link[base]{matrix}}, \code{\link[base]{data.frame}}, or \code{\link[tibble]{tibble}}.
#' @param useNormCounts logical; By default \code{geneScatter} will try to use normCounts instead of counts in \code{SeqExpressionSets}. Set to FALSE to use raw counts instead, though this will generate a warning about using non-normalized data.
#' @param ... Any valid bvt plotting parameter that can be found in \code{\link{bvt_graphic_options}}.
#'
#' @return A list of class \code{\link{npData}}. This contains data necessary to regenerate the plot as well as summary statistics.
#'
#' @examples
#' #While bioinformatic data sets are the intended use case, bvt functions can be used with as regular
#' #plotting functions with data such as iris. As geneScatter is expecting expression data with
#' #patients as columns and genes as rows, it is necessary to take the transpose of the data set first.
#'
#' data("iris")
#' #basic usage
#' geneScatter(iris[,1:2], color=iris$Species, shape=iris$Species, size=iris$Petal.Length,
#'    legend=c("Species", "Petal Length"), main="Basic Scatter")
#'
#' #using adding a trendline
#' a<-geneScatter(iris[,3:4], color=iris$Species,  trendline=TRUE, theme=npGGTheme, verbose=TRUE,
#' corMethod="spearman", pointSize=.8,logScale=10, minorTick=3, minorGuides=TRUE, main="Using Trend Lines")
#' #to access the linear model or cor.test statics later:
#' a$stats
#'
#' #multiple trend lines
#' geneScatter(iris[,1:2], color=iris$Species=="setosa", shape=iris$Species,
#' trendline="color", theme=npColorTheme, legend=c("Is Setosa?","Species"), main="Multiple Trend Lines")
#'
#' #single variable plotting
#' geneScatter(iris[,3], color=iris$Species, main="Single Variable Plotting")
#'
#' #waterfall version of the above. Note type is the same as in base plotting \(i.e. "p","b","h","l"\)
#' orderedIris<-order(iris[,3], decreasing = TRUE)
#' geneScatter(iris[orderedIris,3], color=iris$Species[orderedIris], type="h", main="Waterfall Example")
#'
#' #3D plotting. You can set useRgl=TRUE for rgl based interactive graphics
#' geneScatter(iris[,1:3], color=iris$Species, logScale=2, size=iris[,4] ,pointSize=1)
#'
#' @importFrom purrr map
#' @importFrom Biobase exprs pData fData
#' @export
#' @seealso \code{\link[NicePlots]{niceScatter}}, \code{\link[NicePlots]{niceDensity}}
geneScatter<- function(x, genes=NULL, color=NULL, shape=NULL, size=NULL, trendline=FALSE, symbol="GeneSymbol",legend=NULL, main=TRUE, na.rm=TRUE, facet=NULL, theme=basicTheme, isTidy=if(is.data.frame(x) | is.matrix(x)){TRUE}else{FALSE},  shiny=FALSE, useNormCounts=TRUE, ...) {UseMethod("geneScatter",x)}

#' @importFrom purrr map
#' @importFrom NicePlots niceBox niceVio niceBar niceDensity
#' @importFrom Biobase exprs pData fData
#' @export

geneScatter.default <- function(x, genes=NULL, color=NULL, shape=NULL, size=NULL, trendline=FALSE, symbol="GeneSymbol",legend=NULL, main=TRUE, na.rm=TRUE, facet=NULL, theme=basicTheme, isTidy=if(is.data.frame(x) | is.matrix(x)){TRUE}else{FALSE},  shiny=FALSE, useNormCounts=TRUE, ...) {

  npOptions<-list(...)
  if(any(grepl("npTheme", class(theme)))) {
    npOptions$theme<-theme
  } else {
    warning("Selected theme is not of class 'npTheme'. See help for more details. Proceeding with default settings...",call. = FALSE)
  }

  #First lets handle the case that someone set something to FALSE or NA instead of just leaving it as NULL
  if(sum(genes==FALSE)==1 | sum(is.na(genes))==1) {genes<-NULL}
  if((length(color)==1 & sum(color==FALSE)==1) | sum(is.na(color))==length(color)) {color<-NULL}
  if((length(shape)==1 & sum(shape==FALSE)==1) | sum(is.na(shape))==length(shape)) {shape<-NULL}
  if((length(size)==1 & sum(size==FALSE)==1) | sum(is.na(size))==length(size)) {size<-NULL}
  if((length(facet)==1 & sum(facet==FALSE)==1) | sum(is.na(facet))==length(facet)) {facet<-NULL}

  #Setting default title of nothing was given
  if(main==TRUE) {
    if(length(genes)>1) {
      main<-paste0(c("Gene Expression:",paste0(genes,collapse=", ")),collapse=" ")
    } else {
      main<-paste0(genes, " Expression")
    }
  }

  #Setting the legend to turn on automatically
  if(is.null(legend)){
    legend<-FALSE
    if(!is.null(color) | !is.null(shape)| !is.null(size)) {
      legend<-"Legend"
    }
  }
  #Collecting the expression and factor data
  data<-getGeneData(x=x, gene=genes, plotType="box", symbol=symbol,group=color, subgroup=shape,highlight=size,facet=facet, stack=NULL, useNormCounts=useNormCounts, isTidy=isTidy)

  #Now we convert the options to boolean TRUE/FALSE for compatibility with NicePlots
  if(!is.null(color)){
    color<-data$by$group
  }
  if(!is.null(shape)){
    #This is an issue because we are being lazy and repurposing getGeneData for use here.
    #The subgroup variable is discarded if more than one gene is used as this effectively eats of one of the two possible grouping factors.
    #If this is the case then shape may be left as NULL. If this happens we just run it again as the grouping variable.
    if(is.null(data$by$subgroup)) {
      data2<-getGeneData(x=x, gene=genes, plotType="box", symbol=symbol,group=shape,subgroup=size, highlight=color,facet=facet, stack=NULL, useNormCounts=useNormCounts, isTidy=isTidy)
      shape<-data2$by$group
    } else {
      shape<-data$by$subgroup
    }
  }
  if(!is.null(size)){
    size<-data$by$highlight
  }

  #Formatting options and adding new data
  npOptions<-append(list(x=data$x,by=NULL,color=color, shape=shape, facet=facet,size=size, na.rm=na.rm,main=main, legend=legend, trendline=trendline),npOptions)

  #Calling NicePlots
  dataOut<-1
  dataOut<-do.call("niceScatter",npOptions)
  invisible(dataOut)
}

#' @importFrom purrr map
#' @importFrom NicePlots niceBox niceVio niceBar niceDensity
#' @importFrom Biobase exprs pData fData
#' @export
geneScatter.npData<-function(x, genes=NULL, ...) {
  clOptions<-list(...)
  for(opt in names(clOptions)) {
    if(is.null(x$options[opt])){
      append(x$options,list(opt=clOptions[[opt]]))
    }else{
      x$options[[opt]]<-clOptions[[opt]]
    }
  }
  dataOut<-1
  do.call("niceScatter",x$options)
  invisible(dataOut)
}

