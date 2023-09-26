#' @title Plot Gene Expression Data
#' @description Visualize gene expression data for exploratory data analysis
#'
#' @details
#' The \code{genePlot} function is designed to make visualization of gene expression data simple and easy for R novices and bioinformaticians alike.
#' The function is an S3 generic that accept various R and Bioconductor data sets as input and exacts the expression, factor and annotation data from them according to input data type.
#' The factors allow for splitting expression data from one or more genes into groups and for various plot types with data point overlays where points can be colored by factors levels as well.
#' If the input data is a Bioconductor data set such as an \code{\link[Biobase]{ExpressionSet}} and the \code{gene} option is used, \code{genePlot} will attempt to look up the genes in the associated
#' feature annotation data (e.g. \code{\link[Biobase]{fData}}) according to the data input type and look for the gene symbol column indicated by the \code{symbol} option (defaults to 'GeneSymbol').
#' If no matches are found the row names of are checked of the expression data are check for matches as well.
#' If character values are given for factor input, \code{genePlot} will attempt to look up associated phenotype data (e.g. \code{\link[Biobase]{pData}}).
#' One can also pass raw data vectors/data frames and/or factors to \code{genePlot} to bypass this feature, which is critical for data sets and data formats where integrated phenotype and feature data is not available.
#' The \code{genePlot} uses the \code{NicePlots} graphics library and any \code{NicePlots} option and/or theme can be used in conjunction with options detailed below.
#' The \code{plotType} options supported correspond to \code{NicePlots} functions and include box plots (\code{\link[NicePlots]{niceBox}}), dot plots (\code{\link[NicePlots]{niceDots}}), violin plots (\code{\link[NicePlots]{niceVio}}), bar plots (\code{\link[NicePlots]{niceBar}}) as well as both one/two dimensional kernel density plots (\code{\link[NicePlots]{niceDensity}}).
#' Supported data input types include: \code{\link[Biobase]{ExpressionSet}}, \code{\link[EDASeq]{SeqExpressionSet-class}}, \code{\link[SummarizedExperiment]{SummarizedExperiment}}, \code{\link[limma]{EList-class}}, \code{\link[DESeq2]{DESeqTransform}}, \code{DiffBind} \code{\link[DiffBind]{dba}} objects, as well as standard R data types such as \code{\link[base]{vector}}, \code{\link[base]{matrix}}, \code{\link[base]{data.frame}}, and \code{\link[tibble]{tibble}}.
#' \code{genePlot} silently returns a list of class \code{\link{npData}} that contains a summarized findings, p-values (if indicated), extracted plotting data, and plotting options.
#' All \code{\link{npData}} objects can be replotted using  the \code{\link[graphics]{plot}} function, \code{genePlot} or any of the \code{NicePlots} functions.
#' Options passed to any of these, including \code{plotType} will override the options for the \code{\link{npData}} object. For \code{\link[SummarizedExperiment]{SummarizedExperiment}} and \code{DiffBind} \code{\link[DiffBind]{dba}} objects, the \code{assayType} argument can be added to specify which assay to use if
#' multiple options are available.
#'
#' @param x R data object; Most typically this is an \code{ExpressionSet} there is support for other data types as well.
#' @param gene character; Gene or vector of gene names. These can either be rownames from the gene expression data or looked up in the feature data.
#' @param group factor or name of factor to be extracted from \code{x} (e.g. \code{\link[Biobase]{pData}}). Used as the primary grouping factor.
#' @param subgroup factor or name of factor to be extracted from \code{x} (e.g. \code{\link[Biobase]{pData}}). Used to subgroup data unless multiple genes are selected in which case \code{subgroup} is ignored.
#' @param highlight factor or name of factor to be extracted from \code{x} (e.g. \code{\link[Biobase]{pData}}). Used to color data points by factor levels. Only valid for graphs with point overlays.
#' @param facet factor or name of factor to be extracted from \code{x} (e.g. \code{\link[Biobase]{pData}}). Split the data into multiple smaller graphs.
#' @param stack factor or name of factor to be extracted from \code{x} (e.g. \code{\link[Biobase]{pData}}). Used for stacked bar plots where both the individual and aggregate values are important. Valid only for bar plots.
#' @param plotType character; Can be set to "box", "violin, "dot", "bar", "density" or "surface" for boxplots, violin plots, dot plots, bar plots, and kernel density plots, respectively.
#' @param main character; The main plot title. Defaults to true for automated generation.
#' @param symbol character; Column name of of gene symbols in the feature data of \code{x} (\code{fData}).
#' @param legend logical or character; Draws a figure legend. Use to set the legend title which defaults to "Legend" if equals \code{\link{TRUE}}. Set to \code{\link{FALSE}} to disable.
#' @param na.rm logical; Removes \code{\link{NA}} values prior to ploting.
#' @param shiny logical; Use \code{\link[shiny]{shiny}} interfaces if available.
#' @param groupByGene logical; If more then one gene is listed and \code{grouByGene} is \code{TRUE}
#' @param theme npTheme object; A valid npTheme object the controls default settings. Defaults to \code{basicTheme}.
#' @param useNormCounts logical; By default \code{genePlot} will try to use \code{normCounts()} instead of \code{counts()} in \code{SeqExpressionSets}. Set to FALSE to use raw counts instead, though this will generate a warning about using non-normalized data.
#' @param ... Any parameter recognized by \code{NicePlots} functions. A complete list of these options can be found in \code{\link{bvt_graphic_options}}.
#'
#' @return an list of class \code{\link{npData}}. This contains data necessary to regenerate the plot as well as summary statistics.
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
#' #Kernel density plots
#' genePlot(iris$Sepal.Length,gene=NA, group=iris$Species, plotType="density",
#' main="Distribution of Sepal Lengths by Species")
#'
#' #Plotting multiple columns:
#' genePlot(t(iris[,1:4]),gene=c("Sepal.Length","Petal.Length"), plotType="violin",
#' highlight=iris$Species, pointShape=c(16:18), pointSize=.9)
#'
#' #Multiple columns with grouping factors:
#' genePlot(t(iris[,1:4]),gene=c("Sepal.Width","Petal.Width"), plotType="bar", group=iris$Species)
#'
#' #Same with grouping order reversed
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
genePlot <- function(x, gene=NULL, plotType=c("box","dot","bar","violin","density","surface"), theme=basicTheme, symbol="GeneSymbol",legend=NULL, main=TRUE, na.rm=TRUE, group=NULL, subgroup=NULL, highlight=NULL, facet=NULL, stack=NULL, shiny=FALSE, groupByGene=TRUE, useNormCounts=TRUE, ...) {UseMethod("genePlot",x)}

#' @importFrom purrr map
#' @importFrom NicePlots niceBox niceVio niceBar niceDensity
#' @importFrom Biobase exprs pData fData
#' @export

genePlot.default <- function(x, gene=NULL, plotType=c("box","dot","bar","violin","density","surface"), theme=basicTheme, symbol="GeneSymbol", legend=NULL, main=TRUE, na.rm=TRUE, group=NULL, subgroup=NULL, highlight=NULL, facet=NULL, stack=NULL, shiny=FALSE, groupByGene=TRUE, useNormCounts=TRUE, ...) {
  dataOut<-1
  npOptions<-list(...)
  testenv<-new.env()
  if(any(grepl("npTheme", class(theme)))) {
    npOptions$theme<-theme
  } else {
    warning("Selected theme is not of class 'npTheme'. See help for more details. Proceeding with defaul settings...",call. = FALSE)
  }
  if(!is.null(npOptions$subtitle)) {
    npOptions$sub<-npOptions$subtitle
  }
  #First lets handle the case that someone set something to FALSE or NA instead of just leaving it as NULL
  if(sum(gene==FALSE)==1 | sum(is.na(gene))==1) {gene<-NULL}
  if((length(group)==1 & sum(group==FALSE)==1) | sum(is.na(group))==length(group)) {group<-NULL}
  if((length(subgroup)==1 & sum(subgroup==FALSE)==1) | sum(is.na(subgroup))==length(subgroup)) {subgroup<-NULL}
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
    if(!is.null(subgroup) | !is.null(stack)| !is.null(highlight)) {
      legend<-"Legend"
    }
  }
  if(legend==TRUE) {
    legend<-"Legend"
  }
  Tester<-1
  #Collecting the expression and factor data based on data type
  data<-getGeneData(x=x, gene=gene, plotType=plotType, symbol=symbol,group=group, subgroup=subgroup,highlight=highlight,facet=facet, stack=stack, useNormCounts=useNormCounts, assatType=npOptions$assayType)
  #Note this use of an alternative environment is due to some weird issues seen in RStudio with plotType take on strange values.
  #Unclear if this helped but left in place for now.
  assign("PT",plotType, envir = testenv)
  #Begin shiny GUI pre processing
  if(shiny[1]==TRUE) {
    shinyOpts<-append(list(plotType=plotType,highlight=highlight,groupByGene=groupByGene,group=group, gene=gene, subgroup=subgroup, facet=facet,stack=stack, na.rm=na.rm,main=main, legend=legend,symbol=symbol,useNormCounts=useNormCounts),npOptions)
    if(!is.null(shinyOpts$group)){
      if(length(shinyOpts$group)>1) {
        shinyOpts$group<-deparse(substitute(group))
      }
    }
    if(!is.null(shinyOpts$subgroup)){
      if(length(shinyOpts$subgroup)>1) {
        shinyOpts$subgroup<-deparse(substitute(subgroup))
      }
    }
    if(!is.null(shinyOpts$highlight)){
      if(length(shinyOpts$highlight)>1) {
        shinyOpts$highlight<-deparse(substitute(highlight))
      }
    }
    if(!is.null(shinyOpts$stack)){
      if(length(shinyOpts$stack)>1) {
        shinyOpts$stack<-deparse(substitute(stack))
      }
    }
    if(is.null(shinyOpts$logScale)){
      shinyOpts$logScale<-FALSE
    }
    if(is.null(shinyOpts$theme)){
      shinyOpts$theme<-theme
    }
    #Run Shiny Widget
    dataOut<-shinyGenePlot(data=x, geneList=data$geneList, factorList=data$factorList, gpOptions=shinyOpts, dbName=deparse(substitute(x)),themeName=deparse(substitute(theme)))
    #If RStudio is being used just return the npData object as ploting seems to cause issues.
    dataOut$options<-dataOut$options[names(dataOut$options)!="RSOverride"]
    if(Sys.getenv("RSTUDIO") == "1") {
      return(dataOut$npData)
    }
    #If not using RStudio, the options (vcommand from the shiny widget) is used to reprocesses the data
    #and then continues with the previously scheduled program.
    newOptions<-lapply(dataOut$options, function(o) eval(parse(text=o)))
    newOptions$x<-x
    shinyPlotType<-"box"
    if(is.null(newOptions$plotType)) {
      shinyPlotType<-"box"
    } else {
      shinyPlotType<-newOptions$plotType
    }
    data<-getGeneData(x=x, gene=newOptions$gene, plotType=shinyPlotType, symbol=symbol,group=newOptions$group, subgroup=newOptions$subgroup,highlight=newOptions$highlight,facet=newOptions$facet, stack=newOptions$stack, useNormCounts=useNormCounts)

    if(is.null(newOptions$main)==TRUE) {
      if(length(newOptions$gene)>1) {
        main<-paste0(c("Gene Expression:",paste0(newOptions$gene,collapse=", ")),collapse=" ")
      } else {
        main<-paste0(newOptions$gene, " Expression")
      }
    }

    #Setting the legend to turn on automatically
    if(is.null(newOptions$legend)){
      legend<-FALSE
      if(!is.null(newOptions$subgroup) | !is.null(newOptions$stack)| !is.null(newOptions$highlight)) {
        legend<-"Legend"
      }
    }
    if(legend==TRUE) {
      legend<-"Legend"
    }
    if(is.null(newOptions$groupByGene)) {
      groupByGene<-TRUE
    } else {
      groupByGene<-newOptions$groupByGene
    }
    snpOptions<-newOptions[!(names(newOptions) %in% c("x","gene","group","highlight","subgroup","stack","facet","normCounts","plotType","symbol","main","groupByGene"))]

    group<-newOptions$group
    highlight<-newOptions$highlight
    subgroup<-newOptions$subgroup
    stack<-newOptions$stack

    assign("PT",shinyPlotType, envir = testenv)
    assign("npo",snpOptions,envir = testenv)
  }
  if(shiny==FALSE) {
    assign("npo",npOptions,envir = testenv)
  }
  #Now we convert the options to logical TRUE/FALSE for compatibility with NicePlots
  if(is.null(subgroup)){
    subgroup<-FALSE
  } else {
    subgroup<-TRUE
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
  if(!is.vector(data$x) & (!is.null(group) | !is.null(subgroup))) {
    subgroup<-TRUE
  }
  if(is.null(group) & subgroup==TRUE) {
    subgroup<-FALSE
  }
  if(plotType[1]=="density" & !is.null(group)) {
    subgroup<-TRUE
  }

  #Formatting options and adding new data
  testenv$npo$main<-main
  testenv$npo$legend<-legend
  testenv$npo$flipFacts<-groupByGene
  testenv$npo$x<-data$x
  testenv$npo$by<-data$by
  testenv$npo$subgroup<-subgroup
  testenv$npo$pointHighlights<-highlight
  testenv$npo$facet<-FALSE
  testenv$npo$stack<-stack
  testenv$npo$na.rm<-na.rm

  if(groupByGene==TRUE & data$NullNames==TRUE) {
    if(is.factor(data$by)) {
      testenv$npo$subgroupLabels<-rep("",length(levels(data$by)))
    } else {
      testenv$npo$subgroupLabels<-rep("",length(levels(data$by[,1])))
    }
  }

  #Calling NicePlots
  dataOut<-1
  if(testenv$PT[1]=="box"){
    dataOut<-do.call("niceBox",testenv$npo)
  } else if (testenv$PT[1]=="dot") {
    dataOut<-do.call("niceDots",testenv$npo)
  } else if (testenv$PT[1]=="violin") {
    dataOut<-do.call("niceVio",testenv$npo)
  } else if (testenv$PT[1]=="bar") {
    dataOut<-do.call("niceBar",testenv$npo)
  } else if (testenv$PT[1]=="density") {
    dataOut<-do.call("niceDensity",testenv$npo)
  } else if (testenv$PT[1]=="surface") {
    testenv$npo$plotType<-"surface"
    dataOut<-do.call("niceDensity",testenv$npo)
  } else {
    stop("invalid plot type")
  }
  invisible(dataOut)
}

#' @importFrom purrr map
#' @importFrom NicePlots niceBox niceVio niceBar niceDensity
#' @importFrom Biobase exprs pData fData
#' @export
genePlot.npData<-function(x, gene=NULL, plotType=NULL,theme=basicTheme, ...) {
  #the big difference with the npData version of genePlot over the NicePlots equivalents is we are giving users the opportunity to add new factors to the npData object
  #With the generic plot version of npData, x and by are set but different options can be turned on and off

  clOptions<-list(...)
  if(!is.null(clOptions$subtitle)) {
    clOptions$sub<-clOptions$subtitle
  }

  if(is.null(plotType)) {
    plotType<-x$plotType
  }
  if(sum(c("group","stack","highlight","subgroup") %in% names(clOptions),na.rm=T)>0){
    #check for new group data
    if("group" %in% names(clOptions)) {
      if(is.vector(x$options$by)) {
        if(length(clOptions[["group"]]) == length(x$options$by)) {
          if(sum(x$options$stack[1],x$options$pointHighlights[1],x$options$subgroup[1],na.rm=TRUE)==0) {
            x$options$by<-clOptions[["group"]]
          } else {
            x$options$by<-data.frame(group=clOptions[["group"]],temp=x$options$by)
            if(x$options$subgroup[1]==TRUE) {
              colnames(x$options$by)[2]<- "subgroup"
            } else if (x$options$pointHighlights==TRUE){
              colnames(x$options$by)[2]<- "highlight"
            } else if (x$options$stack==TRUE){
              colnames(x$options$by)[2]<- "stack"
            }
          }
        } else {
          warning(paste0("Length of group (",length(clOptions[["group"]]),") is not the same as existing factor data (",length(x$options$by),").\nIgnoring input...\n\n"),call.=FALSE)
        }
        #Now handle the case that by is a data.frame for group
      } else {
        if(length(clOptions[["group"]]) == dim(x$options$by)[1]) {
          if(sum(x$options$stack[1],x$options$pointHighlights[1],x$options$subgroup[1],na.rm=TRUE) < dim(x$options$by)[2]) {
            x$options$by[,1] <-clOptions[["group"]]
            colnames(x$options$by)[1]<-"group"
          } else {
            x$options$by<-data.frame(group=clOptions[["group"]], x$options$by)
          }
        } else {
          warning(paste0("Length of group (",length(clOptions[["group"]]),") is not the same as the existing factor data (",dim(x$options$by)[1],").\nIgnoring input...\n\n"),call.=FALSE)
        }
      }
    }

    #Check for new subgroup data
    if ("subgroup" %in% names(clOptions)) {
      if (is.vector(x$options$by)) {
        if (length(clOptions[["subgroup"]]) == length(x$options$by)) {
          if ((is.null(plotType[1]) & x$plotType[1]=="bar") | plotType[1] == "bar")  {
            if("stack" %in% names(x$options)) {
              if(x$options$stack[1]==TRUE) {
                x$options$by<-data.frame(group=clOptions[["subgroup"]],stack=x$options$by)
              } else {
                x$options$by<-data.frame(group=x$options$by,subgroup=clOptions[["subgroup"]])
                x$options[["subgroup"]]<-TRUE
              }
            } else {
              x$options$by<-data.frame(group=x$options$by,subgroup=clOptions[["subgroup"]])
              x$options[["subgroup"]]<-TRUE
            }
          } else {
            if("pointHighlights" %in% names(x$options)) {
              if(x$options$pointHighlights[1]==TRUE) {
                x$options$by<-data.frame(group=clOptions[["subgroup"]],highlight=x$options$by)
              } else {
                x$options$by<-data.frame(group=x$options$by,subgroup=clOptions[["subgroup"]])
                x$options[["subgroup"]]<-TRUE
              }
            } else {
              x$options$by<-data.frame(group=x$options$by,subgroup=clOptions[["subgroup"]])
              x$options[["subgroup"]]<-TRUE
            }
          }
        } else {
          warning(paste0("Length of subgroup (",length(clOptions[["subgroup"]]),") is not the same as the existing factor data (",length(x$options$by),").\nIgnoring input...\n\n"),call.=FALSE)
        }
      } else {
        #by is a data.frame for subgroup
        if (length(clOptions[["subgroup"]]) == dim(x$options$by)[1]) {
          if ((is.null(plotType[1]) & x$plotType[1]=="bar") | plotType[1] == "bar")  {
            if("stack" %in% names(x$options)) {
              if(x$options$stack[1]==TRUE) {
                if(x$options$subgroup[1]==TRUE) {
                  if(dim(x$options$by)[2]>1) {
                    x$options$by[,2]<-clOptions[["subgroup"]]
                  } else {
                    x$options$by <- data.frame(group=clOptions[["subgroup"]],stack=x$options$by[,1])
                  }
                } else {
                  if(dim(x$options$by)[2]>1) {
                    x$options$by<-data.frame(group=x$options$by[,1],subgroup=clOptions[["subgroup"]],x$options$by[,-1])
                    x$options[["subgroup"]]<-TRUE
                  } else {
                    x$options$by <- data.frame(group=x$options$by[,1], subgroup=clOptions[["subgroup"]])
                    x$options[["subgroup"]]<-TRUE
                  }
                }
              } else {
                if(dim(x$options$by)[2]>1) {
                  x$options$by[,2]<-clOptions[["subgroup"]]
                  x$options[["subgroup"]]<-TRUE
                } else {
                  x$options$by <- data.frame(group=x$options$by[,1], subgroup=clOptions[["subgroup"]])
                  x$options[["subgroup"]]<-TRUE
                }
              }
            } else {
              if(dim(x$options$by)[2]>1) {
                x$options$by[,2]<-clOptions[["subgroup"]]
                x$options[["subgroup"]]<-TRUE
              } else {
                x$options$by <- data.frame(group=x$options$by[,1], subgroup=clOptions[["subgroup"]])
                x$options[["subgroup"]]<-TRUE
              }
            }
          } else {
            if("pointHighlights" %in% names(x$options)) {
              if(x$options$pointHighlights[1]==TRUE) {
                if(x$options$subgroup[1]==TRUE) {
                  if(dim(x$options$by)[2]>1) {
                    x$options$by[,2]<-data.frame(group=clOptions[["subgroup"]])
                  } else {
                    x$options$by<-data.frame(group=clOptions[["subgroup"]],highlight=by[,1])
                  }
                } else {
                  if(dim(x$options$by)[2]>1) {
                    x$options$by<-data.frame(group=x$options$by[,1],subgroup=clOptions[["subgroup"]],by[,-1])
                    x$options[["subgroup"]]<-TRUE
                  } else {
                    x$options$by<-data.frame(group=clOptions[["subgroup"]],highlight=by[,1])
                  }
                }
              } else {
                if(dim(x$options$by)[2]>1) {
                  x$options$by[,2]<-data.frame(group=clOptions[["subgroup"]])
                  x$options[["subgroup"]]<-TRUE
                } else {
                  x$options$by<-data.frame(group=x$options$by[,1],subgroup=clOptions[["subgroup"]])
                  x$options[["subgroup"]]<-TRUE
                }
              }
            } else {
              if(dim(x$options$by)[2]>1) {
                x$options$by[,2]<-data.frame(group=clOptions[["subgroup"]])
                x$options[["subgroup"]]<-TRUE
              } else {
                x$options$by<-data.frame(group=by[,1],subgroup=clOptions[["subgroup"]])
                x$options[["subgroup"]]<-TRUE
              }
            }
          }
        } else {
          warning(paste0("Length of subgroup (",length(clOptions[["subgroup"]]),") is not the same as the existing factor data (",length(x$options$by),").\nIgnoring input...\n\n"),call.=FALSE)
        }
      }
    }

    #Check for new stack data
    if("stack" %in% names(clOptions)) {
      if (is.vector(x$options$by)) {
        if (length(clOptions[["stack"]]) == length(x$options$by)) {
          if ((is.null(plotType[1]) & x$plotType[1]=="bar") | plotType[1] == "bar")  {
            x$options$by<-data.frame(group=x$options$by,stack=clOptions[["stack"]])
            x$options[["stack"]]<-TRUE
          }
        } else {
          warning(paste0("Length of stack (",length(clOptions[["stack"]]),") is not the same as the existing factor data (",length(x$options$by),").\nIgnoring input...\n\n"),call.=FALSE)
        }
      } else {
        #by is a data.frame for stack
        if (length(clOptions[["stack"]]) == dim(x$options$by)[1]) {
          if ((is.null(plotType[1]) & x$plotType[1]=="bar") | plotType[1] == "bar")  {
            if(dim(x$options$by)[2]>2) {
              if(x$options$subgroup[1]==TRUE) {
                x$options$by[,3]<-clOptions[["stack"]]
              } else {
                x$options$by[,2]<-clOptions[["stack"]]
              }
              x$options[["stack"]]<-TRUE
            } else {
              if(x$options$subgroup[1]==TRUE) {
                x$options$by<-data.frame(x$options$by,stack=clOptions[["stack"]])
              } else {
                x$options$by[,2]<-clOptions[["stack"]]
              }
              x$options[["stack"]]<-TRUE
            }
          }
        } else {
          warning(paste0("Length of stack (",length(clOptions[["stack"]]),") is not the same as the existing factor data (",length(x$options$by),").\nIgnoring input...\n\n"),call.=FALSE)
        }
      }
    }

    #Check for new highlight data
    if("highlight" %in% names(clOptions)) {
      if (is.vector(x$options$by)) {
        if (length(clOptions[["highlight"]]) == length(x$options$by)) {
          if (!((is.null(plotType[1]) & x$plotType[1]=="bar") | plotType[1] == "bar"))  {
            x$options$by<-data.frame(group=x$options$by,highlight=clOptions[["highlight"]])
            x$options[["pointHighlights"]]<-TRUE
          }
        } else {
          warning(paste0("Length of highlight (",length(clOptions[["highlight"]]),") is not the same as the existing factor data (",length(x$options$by),").\nIgnoring input...\n\n"),call.=FALSE)
        }
      } else {
        #by is a data.frame for highlight
        if (length(clOptions[["highlight"]]) == dim(x$options$by)[1]) {
          if (!((is.null(plotType[1]) & x$plotType[1]=="bar") | plotType[1] == "bar"))  {
            if(dim(x$options$by)[2]>2) {
              if(x$options$subgroup[1]==TRUE) {
                x$options$by[,3]<-clOptions[["highlight"]]
                x$options[["highlight"]]<-TRUE
              } else {
                x$options$by[,2]<-clOptions[["highlight"]]
                x$options[["pointHighlights"]]<-TRUE
              }
            } else {
              if(x$options$subgroup[1]==TRUE) {
                x$options$by<-data.frame(x$options$by,highlight=clOptions[["highlight"]])
                x$options[["pointHighlights"]]<-TRUE
               } else {
                x$options$by[,2]<-clOptions[["highlight"]]
                x$options[["pointHighlights"]]<-TRUE
              }
            }
          }
        } else {
          warning(paste0("Length of highlight (",length(clOptions[["highlight"]]),") is not the same as the existing factor data (",length(x$options$by),").\nIgnoring input...\n\n"),call.=FALSE)
        }
      }
    }
  }
  for(opt in names(clOptions)) {
    if(is.null(x$options[opt])){
      append(x$options,list(opt=clOptions[[opt]]))
    } else {
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

#' @name bvt_graphic_options
#' @title Plotting Options In BVT
#' @description Documenting plotting options that are universally available to BVT functions but not otherwise documented
#' @details
#' We can break these options into basic categories as there are a lot of them. This will cover \emph{Plot Enviroment} with includes
#' legend, axis, and title formatting; \emph{Plot Formatting} which includes how points are rendered,
#' determines what error bar measure, etc; \emph{Style Formatting} which controls color, line types, background, etc., and finally \emph{Advanced Settings}.
#' Much of this can be set interactively using the shiny interface and you can see the resulting code in the advanced tab. Another general note
#' is that all theme options can be passed as function arguments to override theme defaults. Finally, while we will mention again in the
#' \emph{Advanced Settings} section, use \code{RSOverride=TRUE} to stop bvt functions from resetting the graphic environment every time, wiping out
#' the plot history. This is particularly important when using layouts for multiple plots or in markdown scripts. This is due to a known bug
#' in RStudio where graphics devices don't always update properly. It is turned on by default when using RStudio. Note that the graphics issuse
#' will reset and appear in the proper fashion after resizing the plotting environment manually.
#'
#' \strong{Plot Enviroment}
#'
#' This section covers options that format axes, labels, titles, etc.
#' \describe{
#'   \item{\code{main}}{Sets the title of the plot, overriding defaults.}
#'   \item{\code{subtitle}}{Sets the subtitle for the plot. This slot is empty by default and used by the \code{showCalc} option if left blank}
#'   \item{\code{ylab}}{Sets the y-axis title. Note \code{xlab} and \code{zlab} can be used to set titles for the x and z axes where applicable.}
#'   \item{\code{yLim}}{Sets the lower bounds and upper bounds, respectively, on the y-axis using a length two numeric vector. Note that \code{xLim}
#' and \code{zLim} can also be used where applicable.}
#'   \item{\code{logScale}}{Transforms the data into log scale based on the number assigned. In order to stabilize low values, a small number, typically 1, is added to all data points be default when active. Set to false to disable.
#' Note that data labels represent the values prior to transformation but scale is still in log. In \emph{scatter plots} and 2 dimensional \emph{density plots} (i.e. contour plots),
#' a numeric vector can be supplied giving log base values (or FALSE for non-log) for the x, y, and z axis, respectively}
#'   \item{\code{expLabels}}{If \code{logScale} is active, data labels are written as powers of the logarithm base supplied when set to TRUE.}
#'   \item{\code{logAdjustment}}{The number added to each data point prior to log transformation. Default is 1. Can be set to 0.}
#'   \item{\code{axisText}}{Length two character vector that adds set text before and/or after the axis labels. eg. \code{c("","\%")} for percentages.
#' For plots with multiple data axis, these can be set individually by providing a list of the form \code{list(x=c("d","CT"),y=c("","\%"))}.}
#'   \item{\code{groupLabels}}{Replaces the group factor level names with the strings supplied in the character vector on the categorical axis.}
#'   \item{\code{subgroupLabels}}{Replaces the subgroup factor level names with the strings supplied in the character vector on the categorical axis.
#' This will update the legend levels when appropriate.}
#'   \item{\code{highlightLabels}}{Replaces the highlight factor level names with the strings supplied in the character. This will updated the legend levels.}
#'   \item{\code{sidePlot}}{When set to TRUE, rotates the entire plot 90 degrees clockwise.}
#'   \item{\code{rotateY}}{Rotates the data axis labels 90 degrees when set to TRUE.}
#'   \item{\code{rotateGroupLabels}}{Rotates the group labels by 90 degrees.}
#'   \item{\code{titleSize}}{Cex like scale factor that controls the size of the plot title set with \code{main}.}
#'   \item{\code{groupLabSize}}{Cex like scale factor that controls the size of the categorical axis group labels that cane be set with \code{groupLabels}.}
#'   \item{\code{subgroupLabSize}}{Cex like scale factor that controls the size fo the subgroup labels on the categorical axis. These labels can be set with \code{subgroupLabels}.}
#'   \item{\code{subSize}}{Cex like scale factor controlling the size of the plot subtitle. Subtitles can be set with \code{subtitle}.}
#'   \item{\code{axisLabelSize}}{Cex like scale factor controling the size of the axis titles in the plot. These can be set with }
#'   \item{\code{yAxisLabSize}}{Cex like scale factor that controls the size of the numeric data labels.}
#'   \item{\code{legendSize}}{Cex like scale factor that controls the size of the legend when applicable.}
#'   \item{\code{groupLabelSpacing}}{Controls the distance of the categorical axis group labels from the axis.}
#'   \item{\code{subgroupLabelSpacing}}{Controls the distance of the categorical axis subgroup labels from the axis.}
#'   \item{\code{legendSpacing}}{Controls the distance between line items in the legend.}
#'   \item{\code{minorTick}}{Sets the number of ticks marks used between major axis ticks. Can be set to 0 or FALSE to disable.}
#'   \item{\code{extendTicks}}{Determines if the minor ticks should be confined between the first and last major tick marks (default for R).
#' Defaults to TRUE which extends ticks to end of the axis.}
#' }
#'
#' \strong{Plot Formatting}
#'
#' This section are options that impact the actual plotting area rather than the axes, labels, titles, etc.
#' \describe{
#'   \item{\code{guides}}{Determines if guide lines should be drawn on the plot area to indicate the major ticks on the data axis.}
#'   \item{\code{minorGuides}}{Determines if minor guide lines are drawn to indicate minor ticks on the data axis. The number of minor ticks can be
#'   set using \code{minorTick}. This can be useful when \code{logScale} is active.}
#'   \item{\code{drawPoints}}{Determines if a data point overlay is drawn. Contours, boxes, violins, error bars, etc. will still be displayed.
#'   For \emph{box plots}, outlier data points (based on IQR distance from the innerquartile boundaries and set using the
#'   \code{outliers} parameter) are still drawn as per convention. Setting \code{outliers=0} will cause the whiskers to extend to the full range, eliminating
#'   the need for outlier plotting if that is desired.}
#'   \item{\code{pointMethod}}{Determines how the data point overlay is drawn. The options are: "linear" where all points drawn
#'   at the same position on the categorical axis; "jitter", where points are randomly assigned a categorical axis value within an interval;
#'    "distribution" where points are placed in order from highest to lowest on the categorical axis interval; and "beeswarm" where points are positioned to
#'    not overlap but pack together as close as possible. The width of the categorical axis interval where these methods can plot is determined
#'    by \code{pointLaneWidth}. The precise position of these points on the graph can be obtained in a \code{\link{npData}} object at \code{npData$options$xypos}.}
#'   \item{\code{swarmOverflow}}{If \code{pointMethod} is set to "beeswarm" and the points are unable to fit within the the region on the categorical
#'   axis determined by \code{pointLaneWidth}, this overflow can be an handled in different ways: "none" this ignores the \code{pointLaneWidth}
#'   setting and plots out the full swarm without restriction. Note that this can cause swarms from different factor levels to overlap; "gutter"
#'   causes the overflow points to be plotted along the outer edge of the allowed interval; "wrap" causes the swarm to just start over plotting
#'   the overflow points in a new swarm on top of the existing one; "random" causes overflow points to be randomly scatter on the categorical
#'   axis, and "omit" which causes the overflow points to be not be drawn.}
#'   \item{\code{aggFun}}{This options is used by \emph{dot} and \emph{bar plots} to determine the aggregate function used, i.e. either "mean" or "median".
#'   This will be represented as the central band or bar height, respectively.}
#'   \item{\code{errFun}}{This determines the method used for drawing error bars or confidence intervals in \emph{dot} and \emph{bar plots}. The options
#'   are: "se" this is the standard error of the mean and is the default setting; "sd" standard deviation; "range" gives the minimum and maximum
#'   values. Note that unlike \emph{box plots}, this is not impacted by the \code{outliers} setting; "t95ci" returns the 95\% confidence interval
#'   based on the t-distribution; and "boot95ci" which calculates empirical 95\% confidence intervals of the aggregate function determined by the
#'   \code{aggFun} setting. This uses resampling with replacement and can take some time with larger data sets. The number of bootstraps used is
#'   set with \code{curvePoints} though this is likely to be changed in a future update. Note that "t95ci", "se", and "sd" only make sense with
#'   \code{aggFun="mean"}. Also, "sd" and "se" have a default multiple of 2x but can be altered using \code{errorMultiple}.}
#'   \item{\code{errorMultiple}}{The number of standard deviations (sd) or standard errors of the mean (se) to be drawn. Default is 2.}
#'   \item{\code{errorBars}}{Controls if error bars are drawn for \emph{dot} and \emph{bar plots}. Note that mean/median bands and bars, respectively,
#'   will still be drawn.}
#'   \item{\code{barType}}{Used only for \emph{dot plots}, this option controls how the central mean/median point is indicated per factor levels: "bar"
#'   this is the default. Bar width is controls by \code{width}, and thickness by \code{lWidth}; "dot", draws a ball where the size is determined
#'   by \code{lWidth}; "none" the mean/median point is not drawn on the graph.}
#'   \item{\code{drawBox}}{Determines if the boxes are drawn for \emph{box} and \emph{violin plots}. Defaults to TRUE. Note that the median bar and
#'   the violin will still be drawn for \emph{box plots} and \emph{violin plots}, respectively.}
#'   \item{\code{trimViolins}}{Controls if the violins in \emph{violin plots} are truncated at the last data point or allowed to taper off.}
#'   \item{\code{trimCurves}}{Determines if single variable \emph{density plot} curves are truncated at the extreme data points or allowed to taper off.
#'   Default value is TRUE, trimming the curves.}
#'   \item{\code{trimTrendLines}}{Determines if trend lines in \code{\link{geneScatter}} stop at the extreme data points or continue to the edge of the graph.}
#'   \item{\code{showTrendConfidence}}{Determines if trend line confidence intervals are drawn in \code{\link{geneScatter}} plots.}
#'   \item{sizeScale}{Cex like scale that determines the maximum difference in point size when \code{size} is active in \code{\link{geneScatter}}.
#'   Note that the minimum point \code{size} is controlled by the \code{pointSize} setting.}
#'   \item{\code{sizeLevels}}{When \code{size} is given a continuous variable in \code{\link{geneScatter}} this value determines the number of factor levels
#'   shown in the figure legend.}
#'   \item{\code{drawRug}}{Draws a rug plot indicating values underneath the density curve of \emph{density plots} if set to TRUE.}
#'   \item{\code{nlevels}}{Determines the number of contour levels drawn in 2D \emph{density plots}.}
#' }
#'
#' \strong{Style Formatting}
#'
#' This section contains the graphic formatting options that are primarily stylistic. Many of these options are also theme settings.
#' Any theme setting can be overridden by providing the value as an argument to a \emph{bvt} function. Themes are named lists of settings that allow for different values
#' for different plot types by appending a either "BP", "VP", "DP", "Bar" or "2D" to the setting to indicate \emph{box plots}, \emph{violin plots},
#' \emph{dot plots}, \emph{bar plots}, or \emph{scatter plots} (including 2D \emph{density/contour plots}), respectively. Using these codes are unnecessary
#' for setting options on the fly. For example providing \code{width=1.25} to a \emph{violin plot}, will automatically override the theme value
#' for \code{widthVP} in the current theme. You can find the valid npThemes available in your environment using \code{\link[NicePlots]{npThemes}}
#' or create a new theme based off of an existing template with \code{\link[NicePlots]{newNPTheme}}. Colors in \emph{bvt} are controlled by
#' the \code{\link{plotColors}} list object. You can read more about how to set colors by following the link: \code{\link{plotColors}}. For the
#' many of the settings below, if only one value is given, then all factor levels, lines, points, etc. will have that same value.
#' If more than one is given, the value will be assigned based on a plotting factor in order of priority: \code{highlight}, \code{subgroup},
#' and \code{group}. Factor level assignments are based on position in the vector, so a factor with 5 levels should be associated with either
#' a length 1 or >= length 5 vector. Note that values are allowed to repeat and while colors will often automatically loop back to the begining
#' of the vector if it runs out of values the settings below will typically not plot past that factor level. By passing a vector containing some
#' \code{NA} values will cause selective drop out at those positions of the line, point, fill, etc. just for that factor level.
#'
#' \describe{
#'   \item{\code{width}}{Controls the maximum width of boxes, bars, vioins, etc.}
#'   \item{\code{pointLaneWidth}}{Controls the width on the categorical axis were data points can be drawn for each factor level.}
#'   \item{\code{vioBoxWidth}}{Determines the width of the inner box plot for \emph{violin plots}.}
#'   \item{\code{fontFamily}}{While future updates may support the selection of custom fonts, right now, \emph{bvt} supports font family selection
#'which includes the options "serif", "sans", and "mono".}
#'   \item{\code{pointSize}}{Cex scale to control the size of data points. This also sets the base point size when size scaling is active in \emph{scatter plots}.}
#'   \item{\code{pointShape}}{A numeric vector used to select the shape of points corresponding to base R's \code{pch} setting.  For \code{\link{geneScatter}}
#'   shape scales are directly associated with a factor using the \code{shape} option.}
#'   \item{\code{lWidth}}{Cex like factor controls the line width for plot features similar to \code{lwd}.}
#'   \item{\code{errorCapType}}{Controls how the end of the whiskers of \emph{box plots} and the tops of the error bars of \emph{dot} and \emph{bar plots}
#'   are drawn. Valid options are "ball", "bar" or "none". Size or width of this cap is determined by \code{errorBarCapWidth}.}
#'   \item{\code{errorBarCapWidth}}{Controls the size of the cap for error bars and whiskers. Setting the value to 1 will cause a bar to be the
#'   same width as the corresponding box or bar.}
#'   \item{\code{errorBarLineType}}{Controls the R line type (lty) used to draw error bars and whiskers.}
#'   \item{\code{vioBoxWidth}}{Controls the relative width of the box plot found at the center of violin plots.}
#'   \item{\code{plotColors}}{A named list of color settings used in \emph{bvt}. Information on how to use this list to make custom
#'   graphs can be found in the following link: \code{\link{plotColors}}.}
#' }
#'
#' \strong{Advanced Settings}
#'
#' These options are a bit more technical or at least don't fit cleanly in to one of the above categories.
#' \describe{
#'   \item{\code{outliers}}{This option controls how many inner quartile ranges (IQR) from the edge of the inner
#'   quartile a data point can be before being considered an outlier. The default value the standard 1.5 times the IQR.
#'   This principally impacts the whisker length on \emph{box plots} and can be set to 0 to automatically have the whiskers extend
#'   to the full range. Note that outlier points are plotted normally regardless.}
#'   \item{\code{useRgl}}{Determines if 3D \emph{scatter plots} and 2D density \emph{surface plots} are rendered interactively using the
#'   RGL library. While only listed as suggested for the package, its can be a helpful way to visualize the data.}
#'   \item{\code{RSOverride}}{There is a known bug in RStudio that causes plot environment not to updated properly. This make plots appear to overlap
#'   until manually resized, forcing a refresh. To get around this, \emph{bvt} will reset the plotting environment with every new graph when
#'   running on RStudio though this results in the loss of plot history. To disable this functionality, set \code{RSOverride=TRUE}. This is
#'   particularly important when constructing multi-panel figures and in rmarkdown files.}
#'   \item{\code{bandwidth}}{This controls the bandwidth setting for kernel density estimation used in \emph{violin plots} and \emph{density plots}.}
#'   \item{\code{curvePoints}}{This controls the number of line segments used to draw curves in \emph{bvt}. It also controls the number of bootstraps
#'   iterations run with using the "boot95ci" option with \code{errFun} though this is likely to change in future updates.}
#'   \item{\code{calcType}}{Statistical testing based on the \code{group} factor levels are automatically calculated in \code{\link{genePlot}} and \code{\link{isoPlot}}.
#'   The options are: "ttest" performs a pairwise t-test with  Holm–Bonferroni correction; "wilcox" performs a pairwise Wilcoxon rank sum test with
#'   Holm–Bonferroni correction; "anova" performs an analysis of variance test; "Tukey-HSD" performs the Tukey range test or honestly significant
#'   difference test; "none" disables statistical testing. Note that \code{subgroup} and \code{highlight} factor levels are ignored. The summery
#'   results can be reported automatically by setting \code{verbose=TRUE} or by looking at the stats section of an \code{\link{npData}} object returned by
#'   \code{\link{genePlot}} or \code{\link{isoPlot}}, i.e. \code{npData$stats}.}
#'   \item{\code{chowCalc}}{If set to TRUE and \code{subtitle} is NULL, the p-value from ANOVA or t-test/wilcoxon comparisons of two groups will be
#'   displayed below the graph.}
#'   \item{\code{corMethod}}{This is used by \code{\link{geneScatter}} to calculate correlation statistics for trend lines if present. The options are passed
#'   to \code{\link[stats]{cor.test}} and can be set to "pearson", "kendall", or "spearman". Full linear model and correlation data can be found in
#'   the returned \code{\link{npData}} object at \code{npData$stats}.}
#'   \item{\code{verbose}}{Causes statical testing results and descriptive statistics to print as output when set to TRUE.}
#' }
#' @seealso \code{\link{plotColors}}, \code{\link{genePlot}}, \code{\link{isoPlot}}, \code{\link{geneScatter}}
NULL

#' @name plotColors
#' @title Setting Color Options in BVT
#' @description How to use the \code{plotColors} vector to create custom plots and new themes in BVT
#' @details In BVT, colors are controlled by a named list of options. Each theme has a \code{plotColor} list associated with it.
#' Like theme options detailed in \code{\link{bvt_graphic_options}}, a \code{plotColors} list can be supplied directly to all BVT functions.
#' The full list does not need to be supplied, but valid list items will selectively override the corresponding \code{plotColors} setting in
#' the active theme. For instance \code{plotColors=list(points=c("red","gold","green"))} will override just the point color settings. For color
#' settings including "points", "lines", and "fill" that may be associated with factor data, these settings will get associated with the following
#' factors in order of priority: \code{highlight} (\code{points} only), \code{stack} (\emph{bar plots} only), \code{subgroup}, and \code{group}.
#' If a single color is given, the color will be applied uniformly (e.g. \code{plotColors=list(fill="blue")} will cause all fill colors to be blue
#' in the plot, while \code{plotColors=list(fill=c("blue","green"))} will cause the fill to alternate between blue and green based on \code{subgroup}
#' or \code{group} factor levels.). If a unique color is desired for each factor level, it is important that the coresponding \code{plotColors} setting
#' has a length equal to or greater than the number of factor levels. For \code{highlight} and \code{stack}, where alternating colors may be
#' critical to the plot and where "points" or "fill" are of length one, respectively, BVT will automatically check to see if "fill", "points", or "lines"
#' are length greater than one and will use the those colors as an alternative, only when "highlight" or "stack" is active.
#' Finally, the \code{\link[NicePlots]{setAlpha}} function provides a convenient way to assign an alpha transparency value to any valid R color including
#' named colors such as "steelblue" or "purple". Below is a list of all the valid \code{plotColors} settings:
#' \describe{
#'   \item{\code{bg}}{Sets the color of the plotting area or canvas. Must be of length 1. Setting this to "open" will make a transparent background
#'   with out the full bounding box (i.e. only the axes are drawn.)}
#'   \item{\code{marginBg}}{Sets the background color for the margins around the the plotting area were the axis labels, titles, etc. are drawn.
#'   Must be length 1. Note that if "bg" is set to "open", this will set the background for the entire plot.}
#'   \item{\code{guides}}{Color of the guide lines corresponding to the major tick marks. Providing more than one color will cause the guides
#'   to cycle through them in order.}
#'   \item{\code{minorGuides}}{Color of the guide lines corresponding to the minor tick marks. Providing more than one color will cause the guides
#'   to cycle through them in order.}
#'   \item{\code{lines}}{Color of line features in the plot. Providing more than one color will cause line color to alternate based on its associated factor as described above.}
#'   \item{\code{points}}{Color of the data point features in the plot. Providing more than one color will cause point colors to alternate based on its associated factor as described above.
#'   If this is length one and \code{highlight} is active, which needs more than one color to work, \code{fill} and \code{lines} color schemes may be used instead.}
#'   \item{\code{fill}}{Fill color for features in the plot. Providing more than one color will cause fill colors to alternate based on its associated factor as described above.
#'   If this is length one and \code{stack} is active, which needs more than one color to work, \code{points} and \code{lines} color schemes may be used instead.}
#'   \item{\code{axis}}{Sets the color of the axes.}
#'   \item{\code{majorTick}}{Color of the major tick marks. Currently, only the first color is used.}
#'   \item{\code{minorTick}}{Color of the minor tick marks. Currently, only the first color is used.}
#'   \item{\code{title}}{Color of the plot title. Currently, this must be a length one vector}
#'   \item{\code{numbers}}{Color of the numeric labels for the plots. Typically, this is the y-axis, but can include the x and z axes for \emph{scatter plots}.
#'   Currently, this must be a length one vector.}
#'   \item{\code{subtext}}{Color of the subtitle, if present. Currently, this must be a length one vector.}
#'   \item{\code{labels}}{Sets the color of the group factor labels on the categorical axis. If more than one color is given, the colors with cycle through the
#'   factor levels, repeating as necessary. Note that the text of the group levels can be customized using the \code{groupLabels} option.}
#'   \item{\code{subgroupLabels}}{Sets the color of the group factor labels on the categorical axis. If more than one color is given, the colors with cycle through the
#'   factor levels, repeating as necessary. Note that the text of the group levels can be customized using the \code{groupLabels} option.}
#'   \item{\code{axisLabels}}{Sets the color of the overall axis label/title. This is is distinct from numeric and group labels associated with the
#'   major tick marks. This must be a length one vector.}
#'   \item{\code{legendLineCol}}{Sets the color of the optional boundry box around the square color lables in the legend. Only the first color is used.}
#'   \item{\code{legendBG}}{Sets the background color for the legend. Default is \code{NA} which levels the background transparent. Appears to be broken in the current build.}
#'   \item{\code{vioBoxFill}}{Sets the fill color for the \emph{box plot} embedded in a \emph{violin plot}. The \code{fill} option will set the fill color for the violin overall,
#'   but this can be used to color the internal area of the box and whisker plot separately. Supports multiple colors and will cycle like the \code{fill} option with factor levels.}
#'   \item{\code{vioBoxLineCol}}{Sets the line color for the \emph{box plot} embedded in a \emph{violin plot}. The \code{lines} option will set the line color for the violin,
#'   but this can be used to set the line color of the box and whisker plot separately. Supports multiple colors and will cycle like the \code{lines} option with factor levels. }
#'   \item{\code{scaleDefaultColor}}{Place holder for color spectrum settings used when setting color to a continuos variable. Not currently supported.}
#' }
#' @examples
#' #Using plotColors on the fly
#' data(iris)
#' library(purrr)
#' new_colors<-map_chr(c("red","blue","green"), setAlpha, alpha=.5)
#' pcList<-list(
#'   lines=new_colors,
#'   points=new_colors,
#'   subgroupLabels=new_colors,
#'   bg="lightgrey",
#'   fill="white")
#'
#' genePlot(t(iris[,1:2]), group=iris$Species, plotColors=pcList, main="Color Example")
#'
#' #Making a new theme based on npGGTheme
#' customTheme<-newNPTheme(theme=npGGTheme, plotColors = pcList, errorBarLineType=1, fontFamily="sans")
#' genePlot(t(iris[,1:2]), group=iris$Species, theme=customTheme, main="Theme Example")
#' npThemes()
#'
#' @seealso \code{\link{bvt_graphic_options}}, \code{\link{genePlot}}, \code{\link{isoPlot}}, \code{\link{geneScatter}}, \code{\link[NicePlots]{newNPTheme}}
NULL

#' @name npData
#' @title The npData Object
#' @description Understanding the npData objects returned by BVT functions
#' @details An object of class \code{npData} is returned from all BCT plotting functions. The primary use of these is to allow for the graph
#' to be saved so it can be plotted later using \code{\link[base]{plot}} or a BVT plotting function. Modifications can also be made such as changing
#' a \emph{bar plot} to a \emph{box plot} or changing any of the \code{\link{bvt_graphic_options}} or \code{\link{plotColors}} settings. These
#' saved objects can also be passed to functions such as \code{link{pick_points}}, to further analyze the data or make new factor levels.
#' These objects also have some useful information stored in them stored in the a named list with the sections: \emph{summary}, \emph{stats}, \emph{plotType},
#' and \emph{options}. A breakdown of each section is below.
#'
#' \strong{npData$summary}
#'
#' This section is a containing the summary data and descriptive statistics for the plot. The data is broken down by \code{group}, \code{subgroup}, and
#' \code{highlight} factor levels if active and reflects the essential statistics of the plot. A \emph{bar plot} might have mean and 2x standard error of
#' the mean listed by factor level while a \emph{box plot} would have the quartile values. These summary statistics will change with \code{aggFun} and \code{errFun}
#' settings. Not yet implemented for \emph{scatter plots} generated with \code{\link{geneScatter}}.
#'
#' \strong{npData$stats}
#'
#' This section is a named list containing an overall p-value where applicable or \code{NA} (\code{npData$stats$p.value}), the name of the statistical
#' test performed (e.g. \code{\link[stats]{pairwise.wilcox.test}}; \code{npData$stats$test}), and the object returned by that test (\code{npData$stats$results}).
#' In the case of \code{\link{geneScatter}}, this section will provide a linear model from \code{\link[stats]{lm}}, and the output of \code{\link[stats]{cor.test}}
#' for each trend line present.
#'
#' \strong{npData$plotType}
#'
#' This is a text string that corresponds to the plot type. Valid options include \code{box}, \code{dot}, \code{bar}, \code{density}, and \code{scatter}.
#'
#' \strong{npData$options}
#'
#' This section has all the options and data needed to regenerate the plot. The three parts worth noting are \code{npData$options$x} which
#' as all the numeric data as a \code{\link[base]{vector}} or \code{\link[base]{data.frame}}, \code{npData$options$by} which contains a
#' \code{\link[base]{factor}} or a \code{\link[base]{data.frame}} of factors used to subset or highlight the data, and \code{npData$options$xypos}
#' which contains the actual x and y coordinates for all points on the graph including those subject to \code{jitter}, \code{distribution}, or \code{beeswarm}
#' settings from the \code{pointMethod} option. Note that \code{x} and \code{by} are the original data with \code{NA} still present while \code{xypos}
#' is post filtering and may be smaller and has the corresponding row numbers of the original \code{x} and \code{by} data to allow for mapping between the two
#' systems. This can be useful when making custom modifications to a graph as shown in the example below.
#'
#' @examples
#'
#' data(iris)
#' #plotting a saved graph with different options
#'
#' a<-genePlot(t(iris[,1:4]), group=iris$Species, plotType="bar", main="npData Example")
#' #Note that only the first two columns of data are used  below
#' #as this is the maximum allowed for the \emph{density plots}.
#' plot(a, plotType="density", theme=npColorTheme)
#'
#' #Example showing the use of npData$options$xypos with a toy data set
#'
#' salesData<-data.frame(Q1=c(1,4,2,2), Q2=c(2,3,4,2), Q3=c(4,4,2,3))
#' stores<-factor(c("Store 1","Store 2","Store 3","Store 4"))
#'
#' #Transposition is necessary since expression data is typical the transpose of tidy data
#' b<-genePlot(t(salesData), highlight=stores, plotType="dot",
#'    theme=npGGTheme, pointLaneWidth=2,
#'    plotColors=list(lines=setAlpha("black",.5)),
#'    legendSize=1, main="Toy Sales Data",
#'    pointMethod="jitter", ylab="Sales Volume (Millions)")
#'
#' #Drawing connecting lines
#' for(i in 1:4){
#'  lines(
#'     x=b$options$xypos[b$options$xypos$ID==i,"x"],
#'     y=b$options$xypos[b$options$xypos$ID==i,"y"],
#'     col=npGGTheme$plotColors$points[i]
#'  )
#' }
#'
#' @seealso \code{\link{genePlot}}, \code{\link{isoPlot}}, \code{\link{geneScatter}}
NULL
