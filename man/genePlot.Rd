% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genePlot.R
\name{genePlot}
\alias{genePlot}
\title{Plot Gene Expression Data}
\usage{
genePlot(
  x,
  gene = NULL,
  plotType = c("box", "dot", "bar", "violin", "density", "surface"),
  theme = basicTheme,
  symbol = "GeneSymbol",
  legend = NULL,
  main = TRUE,
  na.rm = TRUE,
  group = NULL,
  subgroup = NULL,
  highlight = NULL,
  facet = NULL,
  stack = NULL,
  shiny = FALSE,
  groupByGene = TRUE,
  useNormCounts = TRUE,
  ...
)
}
\arguments{
\item{x}{R data object; Most typically this is an \code{ExpressionSet} there is support for other data types as well.}

\item{gene}{character; Gene or vector of gene names. These can either be rownames from the gene expression data or looked up in the feature data.}

\item{plotType}{character; Can be set to "box", "violin, "dot", "bar", "density" or "surface" for boxplots, violin plots, dot plots, bar plots, and kernel density plots, respectively.}

\item{theme}{npTheme object; A valid npTheme object the controls default settings. Defaults to \code{basicTheme}.}

\item{symbol}{character; Column name of of gene symbols in the feature data of \code{x} (\code{fData}).}

\item{legend}{logical or character; Draws a figure legend. Use to set the legend title which defaults to "Legend" if equals \code{\link{TRUE}}. Set to \code{\link{FALSE}} to disable.}

\item{main}{character; The main plot title. Defaults to true for automated generation.}

\item{na.rm}{logical; Removes \code{\link{NA}} values prior to ploting.}

\item{group}{factor or name of factor to be extracted from \code{x} (e.g. \code{\link[Biobase]{pData}}). Used as the primary grouping factor.}

\item{subgroup}{factor or name of factor to be extracted from \code{x} (e.g. \code{\link[Biobase]{pData}}). Used to subgroup data unless multiple genes are selected in which case \code{subgroup} is ignored.}

\item{highlight}{factor or name of factor to be extracted from \code{x} (e.g. \code{\link[Biobase]{pData}}). Used to color data points by factor levels. Only valid for graphs with point overlays.}

\item{facet}{factor or name of factor to be extracted from \code{x} (e.g. \code{\link[Biobase]{pData}}). Split the data into multiple smaller graphs.}

\item{stack}{factor or name of factor to be extracted from \code{x} (e.g. \code{\link[Biobase]{pData}}). Used for stacked bar plots where both the individual and aggregate values are important. Valid only for bar plots.}

\item{shiny}{logical; Use \code{\link[shiny]{shiny}} interfaces if available.}

\item{groupByGene}{logical; If more then one gene is listed and \code{grouByGene} is \code{TRUE}}

\item{useNormCounts}{logical; By default \code{genePlot} will try to use \code{normCounts()} instead of \code{counts()} in \code{SeqExpressionSets}. Set to FALSE to use raw counts instead, though this will generate a warning about using non-normalized data.}

\item{...}{Any parameter recognized by \code{NicePlots} functions. A complete list of these options can be found in \code{\link{bvt_graphic_options}}.}
}
\value{
an list of class \code{\link{npData}}. This contains data necessary to regenerate the plot as well as summary statistics.
}
\description{
Visualize gene expression data for exploratory data analysis
}
\details{
The \code{genePlot} function is designed to make visualization of gene expression data simple and easy for R novices and bioinformaticians alike.
The function is an S3 generic that accept various R and Bioconductor data sets as input and exacts the expression, factor and annotation data from them according to input data type.
The factors allow for splitting expression data from one or more genes into groups and for various plot types with data point overlays where points can be colored by factors levels as well.
If the input data is a Bioconductor data set such as an \code{\link[Biobase]{ExpressionSet}} and the \code{gene} option is used, \code{genePlot} will attempt to look up the genes in the associated
feature annotation data (e.g. \code{\link[Biobase]{fData}}) according to the data input type and look for the gene symbol column indicated by the \code{symbol} option (defaults to 'GeneSymbol').
If no matches are found the row names of are checked of the expression data are check for matches as well.
If character values are given for factor input, \code{genePlot} will attempt to look up associated phenotype data (e.g. \code{\link[Biobase]{pData}}).
One can also pass raw data vectors/data frames and/or factors to \code{genePlot} to bypass this feature, which is critical for data sets and data formats where integrated phenotype and feature data is not available.
The \code{genePlot} uses the \code{NicePlots} graphics library and any \code{NicePlots} option and/or theme can be used in conjunction with options detailed below.
The \code{plotType} options supported correspond to \code{NicePlots} functions and include box plots (\code{\link[NicePlots]{niceBox}}), dot plots (\code{\link[NicePlots]{niceDots}}), violin plots (\code{\link[NicePlots]{niceVio}}), bar plots (\code{\link[NicePlots]{niceBar}}) as well as both one/two dimensional kernel density plots (\code{\link[NicePlots]{niceDensity}}).
Supported data input types include: \code{\link[Biobase]{ExpressionSet}}, \code{\link[EDASeq]{SeqExpressionSet-class}}, \code{\link[SummarizedExperiment]{SummarizedExperiment}}, \code{\link[limma]{EList-class}}, \code{\link[DESeq2]{DESeqTransform}}, \code{DiffBind} \code{\link[DiffBind]{dba}} objects, as well as standard R data types such as \code{\link[base]{vector}}, \code{\link[base]{matrix}}, \code{\link[base]{data.frame}}, and \code{\link[tibble]{tibble}}.
\code{genePlot} silently returns a list of class \code{\link{npData}} that contains a summarized findings, p-values (if indicated), extracted plotting data, and plotting options.
All \code{\link{npData}} objects can be replotted using  the \code{\link[graphics]{plot}} function, \code{genePlot} or any of the \code{NicePlots} functions.
Options passed to any of these, including \code{plotType} will override the options for the \code{\link{npData}} object. For \code{\link[SummarizedExperiment]{SummarizedExperiment}} and \code{DiffBind} \code{\link[DiffBind]{dba}} objects, the \code{assayType} argument can be added to specify which assay to use if
multiple options are available.
}
\examples{
#While designed for use with bioconductor datasets,
#genePlot can also be used for generic data like iris.

data(iris)

#Using a vector of data
genePlot(iris$Sepal.Length,gene=NA, highlight=iris$Species, plotType="dot", pointSize=.75,
width=.5, pointShape=1, main="Distribution of Sepal Length")

#Kernel density plots
genePlot(iris$Sepal.Length,gene=NA, group=iris$Species, plotType="density",
main="Distribution of Sepal Lengths by Species")

#Plotting multiple columns:
genePlot(t(iris[,1:4]),gene=c("Sepal.Length","Petal.Length"), plotType="violin",
highlight=iris$Species, pointShape=c(16:18), pointSize=.9, main="Violins with Point Highlights")

#Multiple columns with grouping factors:
genePlot(t(iris[,1:4]),gene=c("Sepal.Width","Petal.Width"), plotType="bar", group=iris$Species,
main="Multi-Column Data With Grouping Factor")

#Same with grouping order reversed
genePlot(t(iris[,1:4]),gene=c("Sepal.Width","Petal.Width"), plotType="bar", group=iris$Species,
groupByGene=FALSE, theme=npColorTheme, errFun="t95ci", legend=TRUE, main="Multi-Column Data With Grouping Factor")

#2D distribution plotting
genePlot(t(iris[,1:4]),gene=c("Sepal.Width","Petal.Width"), plotType="density",
group=iris$Species, theme=npGGTheme, main="2D Density Example")

#Surface plotting of the above. Use rgl for interactive models.
genePlot(t(iris[,1:4]),gene=c("Sepal.Width","Petal.Width"), plotType="surface", legend=TRUE,
useRgl=FALSE, theta=60, phi=30, theme=npGGTheme)

}
\seealso{
\code{\link[NicePlots]{niceBox}}, \code{\link[NicePlots]{niceVio}}, \code{\link[NicePlots]{niceBar}}, \code{\link[NicePlots]{niceDots}}, \code{\link[NicePlots]{niceDensity}}
}
