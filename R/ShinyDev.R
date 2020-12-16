# library(shiny)
# library(miniUI)
# library(colourpicker)
# library(RColorBrewer)
# library(purrr)
#
# pick_points <- function(data) {
#   ui <- miniPage(
#     gadgetTitleBar(paste("Select points")),
#     miniContentPanel(padding = 0,
#                      plotOutput("plot1", height = "100%", brush = "brush")
#     ),
#     miniButtonBlock(
#       actionButton("add", "", icon = icon("thumbs-up")),
#       actionButton("sub", "", icon = icon("thumbs-down")),
#       actionButton("none", "" , icon = icon("ban")),
#       actionButton("all", "", icon = icon("refresh"))
#     )
#   )
#   server <- function(input, output) {
#     # For storing selected points
#     vals <- reactiveValues(keep = rep(TRUE, nrow(data$options$x)))
#
#     output$plot1 <- renderPlot({
#       # Plot the kept and excluded points as two separate data sets
#       #keep    <- data$options$xypos[ vals$keep, , drop = FALSE]
#       #keepFac<-keep$ID
#       # newFact<-NULL
#       # if(!is.vector(data$options$x)){
#       #   rowMod<- nrow(data$options$x)
#       #   rowSelector<-as.numeric(data$options$xypos$ID[vals$keep])
#       #   rowSelector<- rowSelector %% rowMod
#       #   rowSelector[which(rowSelector==0)]<-rowMod
#       #   newFact<- as.character(seq(rowMod)) %in% as.character(rowSelector)
#       #   # newFact<-rep(NA,rowMod)
#       #   # rowSelector<-as.numeric(data$options$xypos$ID) %% rowMod
#       #   # rowSelector[which(rowSelector==0)]<-rowMod
#       #   # rowSelector<-sort(unique(rowSelector))
#       #   # newFact[rowSelector]<-TRUE
#       #   # rowSelector<-as.numeric(data$options$xypos$ID[vals$keep]) %% rowMod
#       #   # rowSelector[which(rowSelector==0)]<-rowMod
#       #   # rowSelector<-sort(unique(rowSelector))
#       #   # newFact[rowSelector]<-FALSE
#       # } else {
#       #   rowMod<- length(data$options$x)
#       #   newFact<-rep(NA,rowMod)
#       #   newFact[as.numeric(data$options$xypos$ID)] <-TRUE
#       #   newFact[as.numeric(data$options$xypos$ID[vals$keep])] <- FALSE
#       # }
#
#       #exclude <- data[!vals$keep, , drop = FALSE]
#       #tester<-sum(newFact[!is.na(newFact)])
#       #if(tester==0 | tester == length(newFact[!is.na(newFact)])) {
#       if(sum(vals$keep)==length(vals$keep) | is.null(vals$keep)) {
#         genePlot(data, main="testing...",RSOveride=TRUE)
#       } else {
#         genePlot(data, main=paste0("True: ", sum(vals$keep,na.rm = TRUE),"; FALSE: ", sum(!vals$keep,na.rm = TRUE)), highlight=vals$keep, legend="Selected", RSOveride=TRUE)
#       }
#     })
#
#     # Update selected points
#     selected <- reactive({
#       bpoints<-brushedPoints(data$options$xypos, xvar = "x", yvar="y", input$brush, allRows = TRUE)$selected_
#       if(!is.vector(data$options$x)){
#         rowMod<- nrow(data$options$x)
#         rowSelector<-as.numeric(data$options$xypos$ID[bpoints])
#         rowSelector<- rowSelector %% rowMod
#         rowSelector[which(rowSelector==0)]<-rowMod
#         newFact<- as.character(seq(rowMod)) %in% as.character(rowSelector)
#       } else {
#         rowMod<- length(data$options$x)
#         newFact<-rep(NA,rowMod)
#         newFact[as.numeric(data$options$xypos$ID)] <-TRUE
#         newFact[as.numeric(data$options$xypos$ID[bpoints])] <- FALSE
#       }
#       newFact
#     })
#     observeEvent(input$add,  vals$keep <- vals$keep | selected())
#     observeEvent(input$sub,  vals$keep <- vals$keep & !selected())
#     observeEvent(input$all,  vals$keep <- rep(TRUE, nrow(data$options$xypos)))
#     observeEvent(input$none, vals$keep <- rep(FALSE, nrow(data$options$xypos)))
#
#     observeEvent(input$done, {
#       stopApp(data$options$x[!vals$keep,])
#     })
#     observeEvent(input$cancel, {
#       stopApp(NULL)
#     })
#
#   }
#
#   runGadget(ui, server)
# }



#' @title Shiny genePlot Widget
#' @description A shiny widget for bvt's genePlot functions
#'
#' @details
#' This is a the optional gui interface for using genePlot. It is an internal function that is not exported
#' The interface returns the options selected by the user.
#'data, genes, geneList, factors, factorList, theme, gpOptions
#' @param data \code{\link[Biobase]{ExpressionSet}} or other R data object containing the data and possibly feature or phenotype information. While intended for use with Bioconductor objects, \code{\link[base]{data.frame}}, \code{\link[base]{matrix}}, and \code{\link[tibble]{tibble}} can also be used.
#' @param genes character vector; names of genes already in use
#' @param geneList character vector; List of all possible gene names including expression data row names and gene symbols if in use.
#' @param factors active factors. Placeholder
#' @param factorList all possible factors from phenotype data.
#' @param gpOptions active options from the command line. Placeholder.
#'
#' @examples
#' ToDo<-1
#'
#' @importFrom NicePlots setAlpha basicTheme npThemes
#' @importFrom purrr map_chr
#' @importFrom RColorBrewer brewer.pal
#' @importFrom graphics abline par rect
#' @seealso \code{\link[bvt]{genePlot}}
shinyGenePlot <- function(data, genes, geneList, factors, factorList, gpOptions) {
  if(requireNamespace("shiny",quietly = TRUE) & requireNamespace("colourpicker",quietly = TRUE) & requireNamespace("miniUI",quietly = TRUE) == FALSE){
    stop("Missing required libraries for interactive shiny UI. Please install shiny, miniUI and colourpicker.")
  }
  factorList<-append(list("None"=FALSE),as.list(factorList))
  themes<-npThemes()
  ui <- miniUI::miniPage(
    miniUI::miniContentPanel(padding = 0,
     shiny::plotOutput("plot1", height = "100%")
    ),
    miniUI::miniTabstripPanel(
      miniUI::miniTabPanel("Select Data", icon = shiny::icon("database"),
        miniUI::miniContentPanel(padding = 2,
          shiny::fillRow(height = "80px",flex=c(3,4,1,4,1),
            shiny::fillCol(
              shiny::textInput("gene", label = "Genes Symbols/ID", value = "Enter text...")
            ),
            shiny::fillCol(
              style = "margin-top: 25px;",
              shiny::selectInput("geneChoiceList", label=NULL,choices=NULL)
            ),
            shiny::fillCol(
              style = "margin-top: 25px;",
              shiny::actionButton("addGene",label="",icon=shiny::icon("plus"))
            ),
            shiny::fillCol(
              shiny::selectizeInput("activeGenes",label="Selected Data", multiple=TRUE, choices=genes)
            ),
            shiny::fillCol(
              style = "margin-top: 25px;",
              shiny::actionButton("removeGene",label="",icon=shiny::icon("trash"))
            )
          ),
          shiny::fillRow(flex=c(7,1,5),
              shiny::fillCol(
              shiny::selectInput("groupFactorSelector", label="Select a Grouping Factor",choices=factorList),
              shiny::selectInput("subGroupFactorSelector", label="Select a Subgroup Factor",choices=factorList),
              shiny::selectInput("highlightFactorSelector", label="Select a Highlighting Factor",choices=factorList),
              shiny::selectInput("stackFactorSelector", label="Select a Stacking Factor",choices=factorList),
              shiny::checkboxInput("groupByGene",label="Group by Gene", value=TRUE)
            ),
            shiny::fillCol(
              shiny::br(),
              shiny::br(),
              shiny::br(),
              shiny::br(),
              shiny::br()
            ),
            shiny::fillCol(
              shiny::textInput("group", label = "Group", value = FALSE),
              shiny::textInput("subGroup", label = "Subgroup", value = FALSE),
              shiny::textInput("highlight", label = "Highlight", value = FALSE),
              shiny::textInput("stack", label = "Stack (Bar Plot Only)", value = FALSE),
              shiny::checkboxInput("asPercentage",label="As Percent (Stacked Bars)", value=FALSE)
            )
          )
        )
      ),
      miniUI::miniTabPanel("Plot Options", icon = shiny::icon("chart-area"),
        miniUI::miniContentPanel(padding = 5,
          shiny::fillCol(
            shiny::fillRow(
              shiny::selectInput("plotType", label = ("Select PlotType"),
                choices = list("Box Plot" = "box", "Dot Plot" = "dot", "Violin Plot" = "violin", "Bar Plot"="bar", "Density Plot"="density", "Surface Plot"="surface"),
                selected = 1,width="95%"),
              shiny::textInput("main",label="Plot Title",width="95%")
            ),
            shiny::fillRow(
              shiny::checkboxInput("plotPoints", label = "Plot Points", value = TRUE,width="95%"),
              shiny::selectInput("pointMethod", label = "Point Ploting Style",
                          choices = list("Jitter" = "jitter", "Beeswarm" = "beeswarm", "Distribution" = "distribution", "Linear"="linear"),
                          selected = 1,width="95%"),
              shiny::sliderInput("pointSize", label = "Point Size",
                          min=.05,max=3,value=.8,width="95%",step=.05)
            ),
            shiny::fillRow(
              shiny::checkboxInput("legend", label = "Plot Legend", value = FALSE,width="95%"),
              shiny::textInput("legendTitle", label = "Legend Title", value = "Legend",width="95%"),
              shiny::sliderInput("legendSize", label="Scale Legend", min=.25,max=2, value=.9,width="95%")
            ),
            shiny::fillRow(
              shiny::selectInput("logScaleY", label = "Log Scale",
                          choices = list("None"= FALSE, "2" = 2, "10" =10, "e" = exp(1)),
                          selected = 1,width="95%"),
              shiny::selectInput("logScaleX", label = "Log Scale X (Contour)",
                          choices =  list("None"= FALSE, "2" = 2, "10" =10, "e" = exp(1)),
                          selected = 1,width="95%"),
              shiny::sliderInput("logAdjustment", label="Log Adjustment", min=0,max=1, step=0.05, value=1,width="95%",)
            )

          ),
          shiny::fillCol(
            shiny::fillRow(
              shiny::textInput("ylab",label="Y-Axis Label",width="95%"),
              shiny::checkboxInput("expLabels", label = "Exp Labels (Log)", value = FALSE,width="95%"),
              shiny::checkboxInput("rotateY", label = "Rotate Data Labels", value = FALSE,width="95%")
            ),
            shiny::fillRow(
              shiny::textInput("groupNames",label="Group Labels",width="95%",value = NULL),
              shiny::textInput("subGroupNames",label="Subgroup Labels",width="95%",value = NULL),
              shiny::checkboxInput("rotateLabels", label = "Rotate Group Labels", value = FALSE,width="95%")
            ),
            shiny::fillRow(
              shiny::numericInput("minorTick",label="Minor Ticks",min=0, value= 4, step=1, width="95%"),
              shiny::checkboxInput("showMinorGuide", label = "Minor Guides", value = FALSE,width="95%"),
              shiny::checkboxInput("guides", label = "Major Guides", value = FALSE,width="95%")
            ),
            shiny::fillRow(
              shiny::selectInput("aggFun", label = "Central Tendancy",
                          choices = list("Median"= "median", "Mean" = "mean"),
                          selected = 1,width="95%"),
              shiny::selectInput("errFun", label = "Error Measure",
                          choices = list("Standard Deviation"= "sd", "Standard Error" = "se","Range"="range","95-CI T-Distribution"="t95ci","95-CI Bootstrap"="boot95ci"),
                          selected = 2,width="95%"),
              shiny::numericInput("errorMultiple", label="Error Multiplier", value=2,min=0,width="95%",step = 1)
            )
          ),
          shiny::fillCol(height = "50%",
              shiny::fillRow(
              shiny::checkboxInput("errorBars", label = "Draw Error Bars", value = TRUE, width="95%"),
              shiny::sliderInput("phi", label = "Phi Rotation (persp)", min=-180,max=180, value=30,width="95%"),
              shiny::sliderInput("theta", label="Theta Rotation (persp)", min=-180,max=180, value=30,width="95%")
            ),
            shiny::fillRow(
              shiny::checkboxInput("rug", label = "Density Rug Plot", value = FALSE,width="95%"),
              shiny::checkboxInput("trimCurves", label="Trim Curves to Data", value=TRUE, width="95%"),
              shiny::sliderInput("vioBoxWidth", label="Violin Box Width", min=0,max=1.5,step=.1,value=.4,width="95%")
            )
          )
        )
      ),
      miniUI::miniTabPanel("Format", icon = shiny::icon("sliders"),
        miniUI::miniContentPanel(padding = 5,
          shiny::fillCol(
            shiny::fillRow(
              shiny::h4("Themes and Color Selection") ,
              shiny::selectInput("theme", label = "Select a Theme",
                        choices = themes,
                        selected = themes[1],
                        width="95%"
              )
            ),
            shiny::fillRow(
              shiny::selectInput("colorBrewer", label="RColorBrewer Palettes",
                choices = c("BrBG (11)" ,"PiYG (11)", "PRGn (11)","PuOr (11)","RdBu (11)", "RdGy (11)",
                            "RdYlBu (11)", "RdYlGn (11)", "Spectral (11)", "Accent (8)", "Dark2 (8)", "Paired (12)",
                            "Pastel1 (9)", "Pastel2 (8)", "Set1 (9)", "Set2 (8)", "Set3 (12)", "Blues (9)" ,"BuGn (9)" ,
                            "BuPu (9)", "GnBu (9)", "Greens (9)", "Greys (9)", "Oranges (9)", "OrRd (9)", "PuBu (9)",
                            "PuBuGn (9)", "PuRd (9)", "Purples (9)", "RdPu (9)", "Reds (9)","YlGn (9)","YlGnBu (9)",
                            "YlOrBr (9)", "YlOrRd (9)"),
                width="95%"
              ),
              shiny::selectInput("nColors","# Colors",choices = list("11"=11,"10"=10,"9"=9,"8"=8, "7"=7,"6"=6,"5"=5,"4"=4,"3"=3), width="95%", selected ="11"),
              shiny::sliderInput("alphaSlider", label="Set Brewer Palette Alpha",
                          min = 0,max=1,step=.01,value=1,width = "95%"
              )
            ),
            shiny::fillRow(
              shiny::textInput("selectedColors",label="Selected Colors"),
              colourpicker::colourInput("colPicker", label= "Color Picker",
                          value = "#00FF0080",
                          allowTransparent = TRUE,
                          closeOnClick = FALSE
              )
            ),
            shiny::fillRow(
              shiny::plotOutput("colorPlot", height="40px")
            )
          ),
          shiny::hr(),
          shiny::h4("Points, Lines and Fills"),
          shiny::fillCol(
            shiny::fillRow(
              shiny::textInput("pointColors",label="Point Colors", width="95%"),
              shiny::textInput("pointShapes", label="Point Shapes", width="95%")
            ),
            shiny::fillRow(
              shiny::textInput("lineColors",label="Line Colors",  width="95%"),
              shiny::textInput("fillColors",label="Fill Colors", width="95%")
            ),
            shiny::fillRow(
              shiny::sliderInput("lwdSlider","Line Width", min=0,max=5,step=.2,value = 1,width="95%"),
              shiny::sliderInput("width","Lane Width", min=0,max=2,step=.1,value = .65,width="95%"),
              shiny::sliderInput("pointLaneWidth","Point Lane Width", min=0,max=2,step=.1,value = .5,width="95%")
            ),
            shiny::fillRow(
              shiny::textInput("errorBarLineType",label="Error Bar Line Type",width="95%"),
              shiny::selectInput("errorBarCapType", label="Error Bar Cap", choices=list("None"="none","Ball"="ball","Bar"="bar"),selected = "Bar",width="95%"),
              shiny::sliderInput("errorCapWidth", label="Error Cap Width", min=0,max=1,value=.2,step=.05, width="95%")
            )
          ),
          shiny::fillCol(height = "25%",
            shiny::fillRow(
              shiny::textInput("vioBoxFill", label="Violin Box Plot Color", width="95%"),
              shiny::textInput("vioBoxLineCol",label="Violin Box Line Color", width="95%"),
              shiny::selectInput("swarmOverflow", label="Swarm Overflow Strategy", choices=list(Random="random",Gutter="gutter",Wrap="wrap",Omit="Omit",None="none"),selected="Wrap",width="95%")
            )
          ),
          shiny::hr(),
          shiny::h4("Axis and Text"),
          shiny::fillCol(
            shiny::fillRow(
              shiny::selectInput("fontFamily", label="Font Family", choices = list("Sans Serif"="sans","Serif"="serif","Monotype"="mono"),selected="sans",width="95%"),
              shiny::textInput("titleCol", label="Title Color",width="95%"),
              shiny::sliderInput("titleSize",label="Title Size", min=.5, max=4, step = .1,value=1.5, width="95%")
            ),
            shiny::fillRow(
              shiny::textInput("labelCol", label="Group Label Color",  width="95%"),
              shiny::sliderInput("labelSize", label="Group Label Size", min=.2,max=2,step=0.02, value=.96, width="95%"),
              shiny::sliderInput("labelSpacing", label="Group Label Spacing", min=.4,max=2, step=0.02, value=.96, width="95%")
            ),
            shiny::fillRow(
              shiny::textInput("subGroupLabelCol", label="Subgroup Label Color", width="95%"),
              shiny::sliderInput("subGroupLabelSize", label="Subgroup Label Size", min=.2,max=2,step=0.02, value=.68, width="95%"),
              shiny::sliderInput("subGroupLabelSpacing", label="Subgroup Label Spacing", min=0.01, max=1.5, step=0.01, value=.26, width="95%")
            ),
            shiny::fillRow(
              shiny::textInput("dataAxisLablesCol",label="Data (Y-Axis) Label Color", width="95%"),
              shiny::textInput("axisCol",label="Axis Color", width="95%"),
              shiny::textInput("majorTickColor",label="Axis Tick Color", width="95%")
            )
          ),
          shiny::fillCol(height = "50%",
            shiny::fillRow(
              shiny::textInput("minorTickColor",label="Minor Tick Color", width="95%"),
              shiny::textInput("majorGuideColor",label="Guideline Color", width="95%"),
              shiny::textInput("minorGuidesColor",label="Minor Guide Color", width="95%")
            ),
            shiny::fillRow(
              shiny::textInput("axisPrepend", label="Prepend Axis Text",width="95%"),
              shiny::textInput("axisAppend", label="Append Axis Text",width="95%"),
              shiny::checkboxInput("extendTick", label="Extend Minor Ticks", value=TRUE, width="95%")
            )
          ),
          shiny::hr(),
          shiny::h4("Background and Legend"),
          shiny::fillCol(height = "50%",
            shiny::fillRow(
              shiny::textInput("canvasColor",label="Plot Area Color", width="95%"),
              shiny::textInput("marginCol",label="Background Color", width="95%"),
              shiny::textInput("subColor",label="Subtext Color", width="95%")
            ),
            shiny::fillRow(
              shiny::textInput("legendBorder",label="Legend Border Color", width="95%"),
              shiny::textInput("LegendLineCol", label="Legend Box Outline Color", width="95%"),
              shiny::textInput("LegendBG", label="Legend Background Color", width="95%")
            )
          )
        )
      ),
      miniUI::miniTabPanel("Advanced", icon = shiny::icon("code"),
        miniUI::miniContentPanel(padding = 5,
          shiny::h4("Misc Options"),
          shiny::fillCol(
            shiny::fillRow(
              shiny::checkboxInput("useRgl", label="Use RGL for 3D",value =TRUE, width="95%"),
              shiny::checkboxInput("sidePlot", label="Flip X/Y Axses", value=FALSE, width="95%"),
              shiny::sliderInput("curvePoints",label="Curve/Bootstrap Sampling",min=50,max=5000,step=50,value=1000, width="95%")
            ),
            shiny::fillRow(
              shiny::checkboxInput("verbose", label="Verbose Output", value=FALSE, width="95%"),
              shiny::selectInput("calcType", label="Statical Testing",choices=list(None="none",Wilcoxon="wilcox",T.Test="t.test",ANOVA="anvo",Tukey_HSD="Tukey"),selected = "None", width="95%"),
              shiny::textInput("bandwith", label="Bandwidth", value=NULL, width="95%")
            ),
            shiny::fillRow(
              shiny::checkboxInput("strictLimits",label="Strict Data Axis Limits", value=FALSE, width="95%"),
              shiny::sliderInput("outliers",label="Outlier IQR Threshold", min=0,max=5,step=.1,value=1.5, width="95%"),
              shiny::br()
            )
          ),
          shiny::h4("Actual genePlot Code"),
          shiny::fillCol(height="30%",
            shiny::fillRow(
              shiny::verbatimTextOutput("codeExample", placeholder = TRUE)
            )
          ),
          shiny::h4("Data Summary"),
          shiny::fillCol(
            shiny::fillRow(
              shiny::tableOutput("descriptive")
            )
          ),
          shiny::h4("Statistics"),
          shiny::fillCol(
            shiny::fillRow(
              shiny::verbatimTextOutput("stats", placeholder = TRUE)
            )
          )
        )
      )
    ),
    miniUI::miniButtonBlock(
      shiny::actionButton(inputId = "done",label = "Plot it!",width = "20%"),
      shiny::actionButton(inputId = "reset",label = "Reset",width = "20%"),
      shiny::actionButton(inputId = "cancel",label = "Cancel",width = "20%")
    )
  )

  server <- function(input, output, session) {

    aGenes<-shiny::reactiveValues(g=genes)
    RSOveride<-shiny::reactiveValues(rso=TRUE, npData=NULL)
    plotOptions<-shiny::reactiveValues(gene=genes,
                                plotType="box",
                                drawPoints=TRUE,
                                legend=FALSE,
                                minorTick=FALSE,
                                guides=TRUE,
                                minorGuides=FALSE,
                                pointMethod="jitter",
                                theme="basicTheme",
                                pointSize=.75,
                                main=TRUE,
                                logScale=FALSE,
                                expLabels=FALSE,
                                logAdjustment=1,
                                rotateLabels=FALSE,
                                ylab="",
                                main=TRUE,
                                rotateY=FALSE,
                                drawRug=FALSE,
                                errorMultiple=2,
                                aggFun="mean",
                                errFun="se",
                                errorBars=TRUE,
                                trimCurves=TRUE,
                                vioBoxWidth=.4,
                                normalize=FALSE,
                                groupByGene=TRUE,
                                LegendSize=.8,
                                plotColors=list(),
                                titleSize=basicTheme$titleSize,
                                groupLabSize=basicTheme$groupLabSize,
                                groupLabelSpacing=basicTheme$groupLabelSpacing,
                                fontFamily=basicTheme$fontFamily,
                                errorBarCapWidth=basicTheme$errorBarCapWidth,
                                errorBarLineType=1,
                                errorCapType=basicTheme$errorCapType,
                                subgroupLabelSpacing=basicTheme$subgroupLabelSpacing,
                                subGroupLabSize=basicTheme$subGroupLabSize,
                                axisText=c("",""),
                                extendTicks=TRUE,
                                strictLimits=FALSE,
                                sidePlot=FALSE,
                                verbose=FALSE,
                                curvePoints=basicTheme$curvePoints)


    output$plot1 <- shiny::renderPlot({
      RSOveride$npData<-genePlot(data,
                                 gene=plotOptions$gene,
                                 group=plotOptions$group,
                                 subGroup=plotOptions$subGroup,
                                 highlight=plotOptions$highlight,
                                 stack=plotOptions$stack,
                                 plotType=plotOptions$plotType,
                                 drawPoints=plotOptions$drawPoints,
                                 legend=plotOptions$legend,
                                 pointMethod=plotOptions$pointMethod,
                                 RSOveride=RSOveride$rso,
                                 theme=get(plotOptions$theme),
                                 normalize=plotOptions$asPercentage,
                                 groupByGene=plotOptions$groupByGene,
                                 main=plotOptions$main,
                                 pointSize=plotOptions$pointSize,
                                 LegendSize=plotOptions$legendSize,
                                 logScale=plotOptions$logScale,
                                 logAdjustment=plotOptions$logAdjustment,
                                 ylab=plotOptions$ylab,
                                 expLabels=plotOptions$expLabels,
                                 rotateY=plotOptions$rotateY,
                                 rotateLabels=plotOptions$rotateLabels,
                                 groupNames=plotOptions$groupNames,
                                 subGroupLabels=plotOptions$subGroupLabels,
                                 minorTick=plotOptions$minorTick,
                                 guides=plotOptions$guides,
                                 minorGuides=plotOptions$minorGuides,
                                 drawRug=plotOptions$drawRug,
                                 aggFun=plotOptions$aggFun,
                                 errFun=plotOptions$errFun,
                                 errorMultiple=plotOptions$errorMultiple,
                                 theta=plotOptions$theta,
                                 phi=plotOptions$phi,
                                 errorBars=plotOptions$errorBars,
                                 trimCurves=plotOptions$trimCurves,
                                 trimViolins=plotOptions$trimCurves,
                                 vioBoxWidth=plotOptions$vioBoxWidth,
                                 plotColors=plotOptions$plotColors,
                                 pointShape=plotOptions$pointShape,
                                 lWidth=plotOptions$lWidth,
                                 width=plotOptions$width,
                                 pointLaneWidth=plotOptions$pointLaneWidth,
                                 swarmOverflow=plotOptions$swarmOverflow,
                                 errorCapType=plotOptions$errorCapType,
                                 errorBarLineType=plotOptions$errorBarLineType,
                                 errorBarCapWidth=plotOptions$errorBarCapWidth,
                                 fontFamily=plotOptions$fontFamily,
                                 groupLabelSpacing=plotOptions$groupLabelSpacing,
                                 groupLabSize=plotOptions$groupLabSize,
                                 titleSize=plotOptions$titleSize,
                                 subGroupLabSize=plotOptions$subGroupLabSize,
                                 subgroupLabelSpacing=plotOptions$subgroupLabelSpacing,
                                 axisText=plotOptions$axisText,
                                 extendTicks=plotOptions$extendTick,
                                 verbose=plotOptions$verbose,
                                 curvePoints=plotOptions$curvePoints,
                                 sidePlot=plotOptions$sidePlot)

      if(RSOveride$rso==TRUE) {
        RSOveride$rso<-FALSE
      }
      output$descriptive<-shiny::renderTable(RSOveride$npData$summary)
      #utput$Statistics<-RSOveride$npData$stats
    })

    output$colorPlot<-shiny::renderPlot( {
      opar<-par(mar=c(0,0,0,0))
      plot(1,1,type="n", xlim=c(0,100), ylim=c(0,10),axes = FALSE, ylab="", xlab="")
      abline(h =5,lwd=5)
      sCol<-shiny::renderText(shiny::req( {input$selectedColors} ))
      if(!is.null(sCol()) & sCol()!="") {
        sCol<-trimws(unlist(strsplit(sCol(),split = ",")))
        sCol<-map_chr(sCol,function(x) if(!grepl("\\#",x)){setAlpha(x,input$alphaSlider)}else{x})
        inc<-100/length(sCol)
        for(i in seq(length(sCol))) {
          rect(inc*(i-1),-10,inc*i,20,col=sCol[i],lwd=0)
        }
      }
      par(opar)
    })


    ######Data Minitab
    shiny::observeEvent(input$gene, {
      geneChoices<-NULL
      cSymbol<-shiny::renderText({ input$gene })
      if(nchar(cSymbol())>=3) {
        geneChoices<-geneList[grepl(paste0("^",cSymbol()),geneList)]
        if(length(geneChoices)<=40) {
          shiny::updateSelectInput(session,"geneChoiceList",choices=geneChoices)
        }
      }
      if(nchar(cSymbol())<=1 & length(geneList)<=30) {
        shiny::updateSelectInput(session,"geneChoiceList",choices=geneList)
      }
    })

    shiny::observeEvent(input$addGene, {
      if(!is.null(shiny::req(input$geneChoiceList))) {
        cGene<-shiny::renderText({ shiny::req(input$geneChoiceList) })
        if(!is.null(cGene()) & !is.na(cGene()) & length(cGene())>0 & ! cGene() %in% aGenes$g & cGene()!="") {
          aGenes$g<-c(aGenes$g,cGene())
          shiny::updateSelectizeInput(session,"activeGenes",choices=aGenes$g)
          plotOptions$gene<-aGenes$g
        }
      }
    })

    shiny::observeEvent(input$removeGene, {
      if(!is.null(shiny::req( {input$activeGenes }))) {
        aGenes$g<-aGenes$g[! aGenes$g %in% shiny::req( {input$activeGenes} ) ]
        shiny::updateSelectizeInput(session,"activeGenes",choices=aGenes$g)
        plotOptions$gene<-aGenes$g
      }
    })

    shiny::observeEvent(input$groupFactorSelector, {
      if(!is.null(shiny::req( {input$groupFactorSelector }))) {
        shiny::updateTextInput(session, "group", value = input$groupFactorSelector )
      }
    })

    shiny::observeEvent(input$group, {
      if(!is.null(shiny::req( {input$group }))) {
        groupVal<-shiny::renderText( {input$group} )
        if (groupVal() == shiny::req( {input$groupFactorSelector} )) {
          plotOptions[["group"]]<-groupVal()
        } else {
          plotOptions[["group"]] <- tryCatch(
            { suppressWarnings(eval(parse(text=groupVal()))) },
            warning = function(w) {plotOptions[["group"]] <-groupVal()},
            error = function(e) {plotOptions[["group"]] <-FALSE}
          )
        }
      }
    })

    shiny::observeEvent(input$subGroupFactorSelector, {
      if(!is.null(shiny::req( {input$subGroupFactorSelector }))) {
        shiny::updateTextInput(session, "subGroup", value = input$subGroupFactorSelector )
      }
    })

    shiny::observeEvent(input$subGroup, {
      groupVal<-shiny::renderText( {shiny::req(input$group)} )
      if(!is.null(shiny::req( {input$subGroup })) & groupVal()!="FALSE") {
        subGroupVal<-shiny::renderText( {input$subGroup} )
        if(subGroupVal()!="") {
          if(subGroupVal() == shiny::req( {input$subGroupFactorSelector} )){
            plotOptions[["subGroup"]]<-subGroupVal()
          } else {
            plotOptions[["subGroup"]] <- tryCatch(
              { suppressWarnings(eval(parse(text=subGroupVal()))) },
              warning = function(w) {plotOptions[["subGroup"]] <-subGroupVal()},
              error = function(e) {plotOptions[["subGroup"]] <-FALSE}
            )
          }
        } else {
          plotOptions[["subGroup"]]<-FALSE
        }
      }
    })

    shiny::observeEvent(input$highlightFactorSelector, {
      if(!is.null(shiny::req( {input$highlightFactorSelector }))) {
        shiny::updateTextInput(session, "highlight", value = input$highlightFactorSelector )
      }
    })

    shiny::observeEvent(input$highlight, {
      if(!is.null(shiny::req( {input$highlight }))) {
        highlightVal<-shiny::renderText( {input$highlight} )
        if(highlightVal() == shiny::req( {input$highlightFactorSelector} )){
          plotOptions[["highlight"]]<-highlightVal()
        } else {
          plotOptions[["highlight"]] <- tryCatch(
            { suppressWarnings(eval(parse(text=highlightVal()))) },
            warning = function(w) {plotOptions[["highlight"]] <-highlightVal()},
            error = function(e) {plotOptions[["highlight"]] <-FALSE}
          )
        }
      }
    })

    shiny::observeEvent(input$stackFactorSelector, {
      if(!is.null(shiny::req( {input$stackFactorSelector }))) {
        shiny::updateTextInput(session, "stack", value = input$stackFactorSelector )
      }
    })

    shiny::observeEvent(input$stack, {
      if(!is.null(shiny::req( {input$stack }))) {
        stackVal<-shiny::renderText( {input$stack} )
        if(stackVal() == shiny::req( {input$stackFactorSelector} )){
          plotOptions[["stack"]]<-stackVal()
        } else {
          plotOptions[["stack"]] <- tryCatch(
            { suppressWarnings(eval(parse(text=stackVal()))) },
            warning = function(w) {plotOptions[["stack"]] <-stackVal()},
            error = function(e) {plotOptions[["stack"]] <-FALSE}
          )
        }
      }
    })

    shiny::observeEvent(input$groupByGene, {
      plotOptions$groupByGene<-input$groupByGene
    })

    shiny::observeEvent(input$asPercentage, {
      plotOptions$normalize<-input$asPercentage
    })


    ######Plot Options Minitab
    shiny::observeEvent(input$plotType, {
      pt<-shiny::renderText(shiny::req({input$plotType}))
      RSOveride$rso<-TRUE
      plotOptions$plotType<-pt()
      if(pt() != "density" & pt() != "surface") {
        plotOptions$logScale<-plotOptions$logScale[1]
      }
    })

    shiny::observeEvent(input$main, {
      titleText<-shiny::renderText( {input$main} )
      if(is.null(titleText()) | titleText() == "") {
        plotOptions$main<-TRUE
      } else {
        plotOptions$main<-titleText()
      }
    })

    shiny::observeEvent(input$plotPoints, {
      plotOptions$drawPoints<-input$plotPoints
    })

    shiny::observeEvent(input$pointMethod, {
      pm<-shiny::renderText(shiny::req({input$pointMethod}))
      plotOptions$pointMethod<-pm()
    })

    shiny::observeEvent(input$pointSize, {
      plotOptions$pointSize<-input$pointSize
    })

    shiny::observeEvent(input$legend, {
      groupVal<- shiny::renderText( {shiny::req(input$group)} )
      subGroupVal<- shiny::renderText( {shiny::req(input$subGroup)} )
      highlightVal<- shiny::renderText( {shiny::req(input$highlight)} )
      if(length(aGenes$g)>1) {
        if(groupVal() != "FALSE" | highlightVal() != FALSE) {
          legendVal<- shiny::renderText( {input$legend} )
          if(legendVal() == "TRUE") {
            plotOptions$legend <- "Legend"
          } else {
            plotOptions$legend <- FALSE
          }
        } else {
          plotOptions$legend<-FALSE
          shiny::updateCheckboxInput(session,"legend",value=FALSE)
        }
      } else {
        if(highlightVal() != FALSE | (groupVal() !=FALSE & subGroupVal() != FALSE)) {
          legendVal<- shiny::renderText( {input$legend} )
          if(legendVal() == "TRUE") {
            plotOptions$legend <- "Legend"
          } else {
            plotOptions$legend <- FALSE
          }
        } else {
          plotOptions$legend<-FALSE
          shiny::updateCheckboxInput(session,"legend",value=FALSE)
        }
      }
    })

    shiny::observeEvent(input$legendTitle, {
      if(shiny::req( {input$legend} ) == TRUE ) {
        legendTitleVal<-shiny::renderText(shiny::req( {input$legendTitle} ))
        if(is.null(legendTitleVal()) | legendTitleVal() == "") {
          plotOptions$legend<-"Legend"
        } else {
          plotOptions$legend<-legendTitleVal()
        }
      } else {
        plotOptions$legend <- FALSE
      }
    })

    shiny::observeEvent(input$legendSize, {
      plotOptions$legendSize<-input$legendSize
    })

    shiny::observeEvent(input$logScaleY, {
      logScaleYVal<-shiny::renderText(shiny::req( {input$logScaleY} ))
      plotTypeVal<-shiny::renderText(shiny::req( {input$plotType} ))
      if(plotTypeVal() == "density" | plotTypeVal() == "surface" ) {
        if(shiny::req( {input$logScaleX} ) != FALSE ) {
          if (input$logScaleY == FALSE) {
            plotOptions$logScale[1]<-FALSE
          } else {
            plotOptions$logScale[1]<-as.numeric(logScaleYVal())
          }
        } else {
          if (input$logScaleY == FALSE) {
            plotOptions$logScale<-FALSE
          } else {
            plotOptions$logScale<-as.numeric(logScaleYVal())
          }
        }
      } else {
        if (input$logScaleY == FALSE) {
          plotOptions$logScale<-FALSE
        } else {
          plotOptions$logScale<-as.numeric(logScaleYVal())
        }
      }
    })

    shiny::observeEvent(input$logScaleX, {
      logScaleXVal<-shiny::renderText(shiny::req( {input$logScaleX} ))
      plotTypeVal<-shiny::renderText(shiny::req( {input$plotType} ))
      if(plotTypeVal() == "density" | plotTypeVal() == "surface" ) {
        if(shiny::req( {input$logScaleX} ) != FALSE ) {
          plotOptions$logScale[2]<-as.numeric(logScaleXVal())
        } else {
          plotOptions$logScale[2]<-FALSE
        }
      }
    })

    shiny::observeEvent(input$logAdjustment, {
      plotOptions$logAdjustment <- input$logAdjustment
    })

    shiny::observeEvent(input$ylab, {
      ylabVal<-shiny::renderText( {input$ylab} )
      if(ylabVal() == "" | is.null(ylabVal())) {
        plotOptions$ylab<-NULL
      } else {
        plotOptions$ylab<-ylabVal()
      }
    })

    shiny::observeEvent(input$expLabels, {
      if( {input$logScaleX} !=FALSE | {input$logScaleY} !=FALSE ) {
        plotOptions$expLabels<- input$expLabels
      } else {
        plotOptions$expLabels<-FALSE
        shiny::updateCheckboxInput(session ,inputId = "expLabels",value = FALSE)
      }
    })

    shiny::observeEvent(input$rotateY, {
      plotOptions$rotateY<-input$rotateY
    })

    shiny:: observeEvent(input$groupNames, {
      groupNamesVal <- shiny::renderText( {input$groupNames} )
      if(is.null(groupNamesVal()) | groupNamesVal() == "" ){
        plotOptions$groupNames<-NULL
      } else {
        plotOptions$groupNames<-trimws(unlist(strsplit(groupNamesVal(),",")))
      }
    })

    shiny::observeEvent(input$subGroupNames, {
      subGroupNamesVal <- shiny::renderText( {input$subGroupNames} )
      if(is.null(subGroupNamesVal()) | subGroupNamesVal() == "" ){
        plotOptions$subGroupLabels<-NULL
      } else {
        plotOptions$subGroupLabels<-trimws(unlist(strsplit(subGroupNamesVal(),",")))
      }
    })

    shiny::observeEvent(input$rotateLabels, {
      plotOptions$rotateLabels<-input$rotateLabels
    })

    shiny::observeEvent(input$minorTick, {
      plotOptions$minorTick<-input$minorTick
    })

    shiny::observeEvent(input$showMinorGuide, {
      plotOptions$minorGuides<-input$showMinorGuide
    })

    shiny::observeEvent(input$guides, {
      plotOptions$guides=input$guides
    })

    shiny::observeEvent(input$aggFun, {
      aggFunVal<-shiny::renderText(shiny::req( {input$aggFun} ))
      plotOptions$aggFun<-aggFunVal()
    })

    shiny::observeEvent(input$errFun, {
      errFunVal<-shiny::renderText(shiny::req( {input$errFun} ))
      plotOptions$errFun<-errFunVal()
    })

    shiny::observeEvent(input$errorMultiple, {
      plotOptions$errorMultiple <- input$errorMultiple
    })

    shiny::observeEvent(input$errorBars,{
      plotOptions$errorBars<-input$errorBars
    })

    shiny::observeEvent(input$rug, {
      plotOptions$drawRug<-input$rug
    })

    shiny::observeEvent(input$phi, {
      plotOptions$phi<-input$phi
    })

    shiny::observeEvent(input$theta, {
      plotOptions$theta<-input$theta
    })

    shiny::observeEvent(input$trimCurves, {
      plotOptions$trimCurves<-input$trimCurves
      plotOptions$trimViolins<-input$trimCurves
    })

    shiny::observeEvent(input$vioBoxWidth, {
      plotOptions$vioBoxWidth<-input$vioBoxWidth
    })

    #####Format Minitab
    shiny::observeEvent(input$theme, {
      cTheme<-shiny::renderText(shiny::req({input$theme}))
      plotOptions$theme<-cTheme()
    })

    shiny::observeEvent(input$colorBrewer, {
      cbVal<-shiny::renderText(shiny::req({input$colorBrewer}))
      maxC<-as.numeric(gsub("(.*) \\((\\d+)\\)","\\2",cbVal()))
      pallete<-gsub("(.*) \\((\\d+)\\)","\\1",cbVal())
      nC<-maxC
      if(!is.null(pallete) & pallete !="" & is.numeric(maxC)) {
        if(input$nColors <= maxC) {
          nC<-input$nColors
          cbVal<-brewer.pal(name=pallete,n=nC)
        } else {
          cbVal<-brewer.pal(name=pallete,n=maxC)
        }
        cbVal<-map_chr(cbVal,setAlpha, alpha=input$alphaSlider)
        shiny::updateTextInput(session,inputId="selectedColors",value=paste0(cbVal,collapse=","))
        nVect<-rev(seq(3,maxC))
        names(nVect)<-nVect

        shiny::updateSelectInput(session, inputId="nColors", choices=as.list(nVect), selected = max(nVect))
      }
    })

    shiny::observeEvent(input$nColors, {
      nColVal<-shiny::renderText(shiny::req( {input$nColors} ))
      nColVal<-as.numeric( nColVal() )
      cbVal<-shiny::renderText(shiny::req({input$colorBrewer}))
      pallete<-gsub("(.*) \\((\\d+)\\)","\\1",cbVal())
      cbVal<-brewer.pal(name=pallete,n=nColVal)
      cbVal<-map_chr(cbVal,setAlpha, alpha=input$alphaSlider)
      shiny::updateTextInput(session,inputId="selectedColors",value=paste0(cbVal,collapse = ","))
    })

    shiny::observeEvent(input$alphaSlider, {
      scVal<-shiny::renderText({input$selectedColors})
      if(!is.null(scVal()) & scVal() != ""){
        scVal<- map_chr(trimws(unlist(strsplit( scVal(), split=","))), setAlpha, alpha=input$alphaSlider)
        shiny::updateTextInput(session,inputId="selectedColors",value=paste0(scVal,collapse = ","))
      }
    })

    shiny::observeEvent(input$pointColors, {
      cCol<-shiny::renderText(input$pointColors)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$points<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$points<- cVal
      }
    })

    shiny::observeEvent(input$pointShapes, {
      cPch<-shiny::renderText(input$pointShapes)
      if(is.null(cPch()) | cPch()=="") {
        plotOptions$pointShape<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cPch(),",")))
        if(!is.null(cVal[1]) & cVal[1]!=""){
          if(grepl("\\d+",cVal[1])) {
            cVal<-as.numeric(cVal)
          }
        }
        plotOptions$pointShape<- cVal
      }
    })



    shiny::observeEvent(input$lineColors, {
      cCol<-shiny::renderText(input$lineColors)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$lines<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$lines<- cVal
      }
    })

    shiny::observeEvent(input$fillColors, {
      cCol<-shiny::renderText(input$fillColors)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$fill<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$fill<- cVal
      }
    })

    shiny::observeEvent(input$lwdSlider, {
      plotOptions$lWidth<-shiny::req(input$lwdSlider)
    })

    shiny::observeEvent(input$width, {
      plotOptions$width<-shiny::req(input$width)
    })

    shiny::observeEvent(input$pointLaneWidth, {
      plotOptions$pointLaneWidth<-input$pointLaneWidth
    })

    # observeEvent(input$errorBarLineType, {
    #   cVal<-renderText(input$errorBarLineType)
    #   cTheme<-renderText({input$theme})
    #   if(is.null(cVal()) | cVal()=="") {
    #     plotOptions$errorBarLineType<-eval(parse(text=cTheme()))$errorBarLineTypeBP
    #   } else {
    #     cVal<-as.numeric(trimws(unlist(strsplit(cVal(), split=","))))
    #     plotOptions$errorBarLineType<-cVal
    #   }
    # })

    shiny::observeEvent(input$errorBarCapType, {
      cVal<-shiny::renderText(input$errorBarCapType)
      plotOptions$errorCapType<-cVal()
    })

    shiny::observeEvent(input$errorCapWidth, {
      plotOptions$errorBarCapWidth<-input$errorCapWidth
    })

    shiny::observeEvent(input$vioBoxFill, {
      cCol<-shiny::renderText(input$vioBoxFill)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$vioBoxFill<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$vioBoxFill<- cVal
      }
    })

    shiny::observeEvent(input$vioBoxLineCol, {
      cCol<-shiny::renderText(input$vioBoxLineCol)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$vioBoxLineCol<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$vioBoxLineCol<- cVal
      }
    })

    shiny::observeEvent(input$swarmOverflow, {
      cVal<-shiny::renderText(input$swarmOverflow)
      plotOptions$swarmOverflow<-cVal()
    })

    shiny::observeEvent(input$fontFamily, {
      cVal<-shiny::renderText(input$fontFamily)
      plotOptions$fontFamily<-cVal()
    })

    shiny::observeEvent(input$titleCol, {
      cCol<-shiny::renderText(input$titleCol)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$title<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$title<- cVal
      }
    })

    shiny::observeEvent(input$titleSize, {
      plotOptions$titleSize<-input$titleSize
    })

    shiny::observeEvent(input$labelCol, {
      cCol<-shiny::renderText(input$labelCol)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$labels<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$labels<- cVal
      }
    })

    shiny::observeEvent(input$labelSize, {
      plotOptions$groupLabSize<-input$labelSize
    })

    shiny::observeEvent(input$labelSpacing, {
      plotOptions$groupLabelSpacing<-input$labelSpacing
    })

    shiny::observeEvent(input$subGroupLabelCol, {
      cCol<-shiny::renderText(input$subGroupLabelCol)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$subGroupLabels<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$subGroupLabels<- cVal
      }
    })

    shiny::observeEvent(input$subGroupLabelSize, {
      plotOptions$subGroupLabSize<-input$subGroupLabelSize
    })

    shiny::observeEvent(input$subGroupLabelSpacing, {
      plotOptions$subgroupLabelSpacing<-input$subGroupLabelSpacing
    })

    shiny::observeEvent(input$dataAxisLablesCol, {
      cCol<-shiny::renderText(input$dataAxisLablesCol)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$numbers<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$numbers<- cVal
      }
    })

    shiny::observeEvent(input$axisCol, {
      cCol<-shiny::renderText(input$axisCol)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$axis<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$axis<- cVal
      }
    })

    shiny::observeEvent(input$majorTickColor, {
      cCol<-shiny::renderText(input$majorTickColor)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$majorTick<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$majorTick<- cVal
      }
    })

    shiny::observeEvent(input$minorTickColor, {
      cCol<-shiny::renderText(input$minorTickColor)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$minorTick<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$minorTick<- cVal
      }
    })

    shiny::observeEvent(input$majorGuideColor, {
      cCol<-shiny::renderText(input$majorGuideColor)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$guides<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$guides<- cVal
      }
    })

    shiny::observeEvent(input$minorGuidesColor, {
      cCol<-shiny::renderText(input$minorGuidesColor)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$minorGuides<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$minorGuides<- cVal
      }
    })

    shiny::observeEvent(input$axisPrepend, {
      preVal<-shiny::renderText(input$axisPrepend)
      aVal<-shiny::renderText(input$axisAppend)
      axisText<-c("","")
      if(is.null(preVal())) {
        axisText[1]<-""
      } else {
        axisText[1]<-preVal()
      }
      if(is.null(aVal())) {
        axisText[2]<-""
      } else {
        axisText[2]<-aVal()
      }
      plotOptions$axisText<-axisText
    })

    shiny::observeEvent(input$axisAppend, {
      preVal<-shiny::renderText(input$axisPrepend)
      aVal<-shiny::renderText(input$axisAppend)
      axisText<-c("","")
      if(is.null(preVal())) {
        axisText[1]<-""
      } else {
        axisText[1]<-preVal()
      }
      if(is.null(aVal())) {
        axisText[2]<-""
      } else {
        axisText[2]<-aVal()
      }
      plotOptions$axisText<-axisText
    })

    shiny::observeEvent(input$extendTick, {
      plotOptions$extendTick<-input$extendTick
    })

    shiny::observeEvent(input$canvasColor, {
      cCol<-shiny::renderText(input$canvasColor)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$bg<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$bg<- cVal
      }
    })

    shiny::observeEvent(input$marginCol, {
      cCol<-shiny::renderText(input$marginCol)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$marginBg<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$marginBg<- cVal
      }
    })

    shiny::observeEvent(input$subColor, {
      cCol<-shiny::renderText(input$subColor)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$subtext<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$subtext<- cVal
      }
    })

    shiny::observeEvent(input$legendBorder, {
      cCol<-shiny::renderText(input$legendBorder)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$LegendBorder<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$LegendBorder<- cVal
      }
    })

    shiny::observeEvent(input$LegendLineCol, {
      cCol<-shiny::renderText(input$LegendLineCol)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$LegendLineCol<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$LegendLineCol<- cVal
      }
    })

    shiny::observeEvent(input$LegendBG, {
      cCol<-shiny::renderText(input$LegendBG)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$LegendBG<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$LegendBG<- cVal
      }
    })


    #Advanced Tab
    shiny::observeEvent(input$useRgl, {
      plotOptions$LegendBG<-input$useRgl
    })

    shiny::observeEvent(input$sidePlot, {
      plotOptions$sidePlot<-input$sidePlot
    })

    shiny::observeEvent(input$curvePoints, {
      plotOptions$curvePoints<-input$curvePoints
    })

    shiny::observeEvent(input$verbose, {
      plotOptions$verbose<-input$verbose
    })


    #shiny::selectInput("calcType", label="Statical Testing",choices=list(None="none",Wilcoxon="wilcox",T.Test="t.test",ANOVA="anvo",Tukey_HSD="Tukey"),selected = "None", width="95%"),
  #shiny::textInput("bandwith", label="Bandwidth", value=NULL, width="95%")

    shiny::observeEvent(input$strictLimits, {
      plotOptions$strictLimits<-input$strictLimits
    })

  #shiny::sliderInput("outliers",label="Outlier IQR Threshold", min=0,max=5,step=.1,value=1.5, width="95%"),


    #Lower Button Panel
    shiny::observeEvent(input$done, {
      shiny::stopApp(RSOveride$npData)
    })
    shiny::observeEvent(input$cancel, {
      shiny::stopApp(NULL)
    })
  }
  viewer<-shiny::dialogViewer(dialogName = "Interactive genePlot UI",height = 2000)
  shiny::runGadget(ui, server, viewer=viewer)
}
