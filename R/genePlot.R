#' @title Plot Gene Expression Data
#' @description Visualize gene expression data for exploratory data analysis
#'
#' @details
#' The \code{genePlot} is designed to make visualization of gene expression data simple and easy for R novices and bioinformaticians alike.
#' The function is an S3 generic that accept various R and Bioconductor data sets as input and exacts the expression, factor and annotation data from them according to type.
#' The factors allow for splitting expression data from one or more genes into groups and for plot types with data point overlays, points can be colored by factors levels as well.
#' If the input data is a Bioconductor data set such as an \code{\link[Biobase]{ExpressionSet}} and the \code{gene} option is used, \code{genePlot} will attempt to look up the genes in the associated feature annotation data (e.g. \code{\link[Biobase]{fData}}) according to the data input type and look for the gene symbol column indicated by the \code{symbol} option (defaults to 'GeneSymbol').
#' If no matches are found the row names of are checked of the expression data are check for matches as well.
#' If character values are given for factor input, \code{genePlot} will attempt to look up associated phenotype data (e.g. \code{\link[Biobase]{pData}}).
#' One can also pass raw data vectors/data frames and/or factors to \code{genePlots} to bypass this feature, which is critical for data sets and data formats where integrated phenotype and feature data is not available.
#' The \code{genePlot} uses the \code{NicePlots} graphics library and any \code{NicePlots} option and/or theme can be used in conjunction with options detailed below.
#' The \code{plotType} options supported correspond to \code{NicePlots} functions and include box plots (\code{\link[NicePlots]{niceBox}}), dot plots (\code{\link[NicePlots]{niceDots}}), violin plots (\code{\link[NicePlots]{niceVio}}), bar plots (\code{\link[NicePlots]{niceBar}}) as well as both one/two dimensional kernel density plots (\code{\link[NicePlots]{niceDensity}}).
#' Supported data input types include: \code{\link[Biobase]{ExpressionSet}}, \code{\link[EDASeq]{SeqExpressionSet-class}}, \code{\link[SummarizedExperiment]{SummarizedExperiment}}, \code{\link[limma]{EList-class}}, \code{\link[DESeq2]{DESeqTransform}}, \code{DiffBind} \code{\link[DiffBind]{dba}} objects, as well as standard R data types such as \code{\link[base]{vector}}, \code{\link[base]{matrix}}, \code{\link[base]{data.frame}}, and \code{\link[tibble]{tibble}}.
#' \code{genePlot} silently returns a list of class \code{npData} that contains a summarized findings, p-values (if indicated), extracted plotting data, and plotting options.
#' All npData objects can be replotted using  the \code{\link[graphics]{plot}} function, \code{genePlot} or any of the \code{NicePlots} functions.
#' Options passed to any of these, including \code{plotType} will override the options for the \code{npData} object. For \code{\link[SummarizedExperiment]{SummarizedExperiment}} and \code{DiffBind} \code{\link[DiffBind]{dba}} objects, the \code{assayType} argument can be added to specify which assay to use if
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
#' @param useNormCounts logical; By default \code{genePlot} will try to use \code{normCounts()} instead of \code{counts()} in \code{SeqExpressionSets}. Set to FALSE to use raw counts instead, though this will generate a warning about useing non-normalized data.
#' @param ... Any parameter recognized by \code{NicePlots} functions.
#'
#' @return an list of class \code{npData}. This contains data necessary to regenerate the plot as well as summary statistics.
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

