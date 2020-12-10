library(shiny)
library(miniUI)
library(colourpicker)
library(RColorBrewer)
library(purrr)

pick_points <- function(data) {
  ui <- miniPage(
    gadgetTitleBar(paste("Select points")),
    miniContentPanel(padding = 0,
                     plotOutput("plot1", height = "100%", brush = "brush")
    ),
    miniButtonBlock(
      actionButton("add", "", icon = icon("thumbs-up")),
      actionButton("sub", "", icon = icon("thumbs-down")),
      actionButton("none", "" , icon = icon("ban")),
      actionButton("all", "", icon = icon("refresh"))
    )
  )
  server <- function(input, output) {
    # For storing selected points
    vals <- reactiveValues(keep = rep(TRUE, nrow(data$options$x)))

    output$plot1 <- renderPlot({
      # Plot the kept and excluded points as two separate data sets
      #keep    <- data$options$xypos[ vals$keep, , drop = FALSE]
      #keepFac<-keep$ID
      # newFact<-NULL
      # if(!is.vector(data$options$x)){
      #   rowMod<- nrow(data$options$x)
      #   rowSelector<-as.numeric(data$options$xypos$ID[vals$keep])
      #   rowSelector<- rowSelector %% rowMod
      #   rowSelector[which(rowSelector==0)]<-rowMod
      #   newFact<- as.character(seq(rowMod)) %in% as.character(rowSelector)
      #   # newFact<-rep(NA,rowMod)
      #   # rowSelector<-as.numeric(data$options$xypos$ID) %% rowMod
      #   # rowSelector[which(rowSelector==0)]<-rowMod
      #   # rowSelector<-sort(unique(rowSelector))
      #   # newFact[rowSelector]<-TRUE
      #   # rowSelector<-as.numeric(data$options$xypos$ID[vals$keep]) %% rowMod
      #   # rowSelector[which(rowSelector==0)]<-rowMod
      #   # rowSelector<-sort(unique(rowSelector))
      #   # newFact[rowSelector]<-FALSE
      # } else {
      #   rowMod<- length(data$options$x)
      #   newFact<-rep(NA,rowMod)
      #   newFact[as.numeric(data$options$xypos$ID)] <-TRUE
      #   newFact[as.numeric(data$options$xypos$ID[vals$keep])] <- FALSE
      # }

      #exclude <- data[!vals$keep, , drop = FALSE]
      #tester<-sum(newFact[!is.na(newFact)])
      #if(tester==0 | tester == length(newFact[!is.na(newFact)])) {
      if(sum(vals$keep)==length(vals$keep) | is.null(vals$keep)) {
        genePlot(data, main="testing...",RSOveride=TRUE)
      } else {
        genePlot(data, main=paste0("True: ", sum(vals$keep,na.rm = TRUE),"; FALSE: ", sum(!vals$keep,na.rm = TRUE)), highlight=vals$keep, legend="Selected", RSOveride=TRUE)
      }
    })

    # Update selected points
    selected <- reactive({
      bpoints<-brushedPoints(data$options$xypos, xvar = "x", yvar="y", input$brush, allRows = TRUE)$selected_
      if(!is.vector(data$options$x)){
        rowMod<- nrow(data$options$x)
        rowSelector<-as.numeric(data$options$xypos$ID[bpoints])
        rowSelector<- rowSelector %% rowMod
        rowSelector[which(rowSelector==0)]<-rowMod
        newFact<- as.character(seq(rowMod)) %in% as.character(rowSelector)
      } else {
        rowMod<- length(data$options$x)
        newFact<-rep(NA,rowMod)
        newFact[as.numeric(data$options$xypos$ID)] <-TRUE
        newFact[as.numeric(data$options$xypos$ID[bpoints])] <- FALSE
      }
      newFact
    })
    observeEvent(input$add,  vals$keep <- vals$keep | selected())
    observeEvent(input$sub,  vals$keep <- vals$keep & !selected())
    observeEvent(input$all,  vals$keep <- rep(TRUE, nrow(data$options$xypos)))
    observeEvent(input$none, vals$keep <- rep(FALSE, nrow(data$options$xypos)))

    observeEvent(input$done, {
      stopApp(data$options$x[!vals$keep,])
    })
    observeEvent(input$cancel, {
      stopApp(NULL)
    })

  }

  runGadget(ui, server)
}



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
#'
#' @importFrom purrr map_chr
#' @importFrom RColorBrewer brewer.pal
#' @seealso \code{\link[bvt]{genePlot}}
shinyGenePlot <- function(data, genes, geneList, factors, factorList, gpOptions) {
  if(require("shiny",quietly = TRUE) & require("colourpicker",quietly = TRUE) & require("miniUI",quietly = TRUE) == FALSE){
    stop("Missing required libraries for interactive shiny UI. Please install shiny, miniUI and colourpicker.")
  }
  factorList<-append(list("None"=FALSE),as.list(factorList))
  themes<-npThemes()
  ui <- miniPage(
    miniContentPanel(padding = 0,
     plotOutput("plot1", height = "100%")
    ),
    miniTabstripPanel(
      miniTabPanel("Select Data", icon = icon("database"),
        miniContentPanel(padding = 2,
          fillRow(height = "80px",flex=c(3,4,1,4,1),
            fillCol(
              textInput("gene", label = "Genes Symbols/ID", value = "Enter text...")
            ),
            fillCol(
              style = "margin-top: 25px;",
              selectInput("geneChoiceList", label=NULL,choices=NULL)
            ),
            fillCol(
              style = "margin-top: 25px;",
              actionButton("addGene",label="",icon=icon("plus"))
            ),
            fillCol(
              selectizeInput("activeGenes",label="Selected Data", multiple=TRUE, choices=genes)
            ),
            fillCol(
              style = "margin-top: 25px;",
              actionButton("removeGene",label="",icon=icon("trash"))
            )
          ),
          fillRow(flex=c(7,1,5),
            fillCol(
              selectInput("groupFactorSelector", label="Select a Grouping Factor",choices=factorList),
              selectInput("subGroupFactorSelector", label="Select a Subgroup Factor",choices=factorList),
              selectInput("highlightFactorSelector", label="Select a Highlighting Factor",choices=factorList),
              selectInput("stackFactorSelector", label="Select a Stacking Factor",choices=factorList),
              checkboxInput("groupByGene",label="Group by Gene", value=TRUE)
            ),
            fillCol(
              br(),
              br(),
              br(),
              br(),
              br()
            ),
            fillCol(
              textInput("group", label = "Group", value = FALSE),
              textInput("subGroup", label = "Subgroup", value = FALSE),
              textInput("highlight", label = "Highlight", value = FALSE),
              textInput("stack", label = "Stack (Bar Plot Only)", value = FALSE),
              checkboxInput("asPercentage",label="As Percent (Stacked Bars)", value=FALSE)
            )
          )
        )
      ),
      miniTabPanel("Plot Options", icon = icon("chart-area"),
        miniContentPanel(padding = 5,
          fillCol(
            fillRow(
              selectInput("plotType", label = ("Select PlotType"),
                choices = list("Box Plot" = "box", "Dot Plot" = "dot", "Violin Plot" = "violin", "Bar Plot"="bar", "Density Plot"="density", "Surface Plot"="surface"),
                selected = 1,width="95%"),
              textInput("main",label="Plot Title",width="95%")
            ),
            fillRow(
              checkboxInput("plotPoints", label = "Plot Points", value = TRUE,width="95%"),
              selectInput("pointMethod", label = "Point Ploting Style",
                          choices = list("Jitter" = "jitter", "Beeswarm" = "beeswarm", "Distribution" = "distribution", "Linear"="linear"),
                          selected = 1,width="95%"),
              sliderInput("pointSize", label = "Point Size",
                          min=.05,max=3,value=.8,width="95%",step=.05)
            ),
            fillRow(
              checkboxInput("legend", label = "Plot Legend", value = FALSE,width="95%"),
              textInput("legendTitle", label = "Legend Title", value = "Legend",width="95%"),
              sliderInput("legendSize", label="Scale Legend", min=.25,max=2, value=.9,width="95%")
            ),
            fillRow(
              selectInput("logScaleY", label = "Log Scale",
                          choices = list("None"= FALSE, "2" = 2, "10" =10, "e" = exp(1)),
                          selected = 1,width="95%"),
              selectInput("logScaleX", label = "Log Scale X (Contour)",
                          choices =  list("None"= FALSE, "2" = 2, "10" =10, "e" = exp(1)),
                          selected = 1,width="95%"),
              sliderInput("logAdjustment", label="Log Adjustment", min=0,max=1, step=0.05, value=1,width="95%",)
            )

          ),
          fillCol(
            fillRow(
              textInput("ylab",label="Y-Axis Label",width="95%"),
              checkboxInput("expLabels", label = "Exp Labels (Log)", value = FALSE,width="95%"),
              checkboxInput("rotateY", label = "Rotate Data Labels", value = FALSE,width="95%")
            ),
            fillRow(
              textInput("groupNames",label="Group Labels",width="95%",value = NULL),
              textInput("subGroupNames",label="Subgroup Labels",width="95%",value = NULL),
              checkboxInput("rotateLabels", label = "Rotate Group Labels", value = FALSE,width="95%")
            ),
            fillRow(
              numericInput("minorTick",label="Minor Ticks",min=0, value= 4, step=1, width="95%"),
              checkboxInput("showMinorGuide", label = "Minor Guides", value = FALSE,width="95%"),
              checkboxInput("guides", label = "Major Guides", value = FALSE,width="95%")
            ),
            fillRow(
              selectInput("aggFun", label = "Central Tendancy",
                          choices = list("Median"= "median", "Mean" = "mean"),
                          selected = 1,width="95%"),
              selectInput("errFun", label = "Error Measure",
                          choices = list("Standard Deviation"= "sd", "Standard Error" = "se","Range"="range","95-CI T-Distribution"="t95ci","95-CI Bootstrap"="boot95ci"),
                          selected = 2,width="95%"),
              numericInput("errorMultiple", label="Error Multiplier", value=2,min=0,width="95%",step = 1)
            )
          ),
          fillCol(height = "50%",
            fillRow(
              checkboxInput("errorBars", label = "Draw Error Bars", value = TRUE, width="95%"),
              sliderInput("phi", label = "Phi Rotation (persp)", min=-180,max=180, value=30,width="95%"),
              sliderInput("theta", label="Theta Rotation (persp)", min=-180,max=180, value=30,width="95%")
            ),
            fillRow(
              checkboxInput("rug", label = "Density Rug Plot", value = FALSE,width="95%"),
              checkboxInput("trimCurves", label="Trim Curves to Data", value=TRUE, width="95%"),
              sliderInput("vioBoxWidth", label="Violin Box Width", min=0,max=1.5,step=.1,value=.4,width="95%")
            )
          )
        )
      ),
      miniTabPanel("Format", icon = icon("sliders"),

        miniContentPanel(padding = 5,

          fillCol(
            fillRow(
              h4("Themes and Color Selection") ,
              selectInput("theme", label = "Select a Theme",
                        choices = themes,
                        selected = themes[1],
                        width="95%"
              )
            ),
            fillRow(
              selectInput("colorBrewer", label="RColorBrewer Palettes",
                choices = c("BrBG (11)" ,"PiYG (11)", "PRGn (11)","PuOr (11)","RdBu (11)", "RdGy (11)",
                            "RdYlBu (11)", "RdYlGn (11)", "Spectral (11)", "Accent (8)", "Dark2 (8)", "Paired (12)",
                            "Pastel1 (9)", "Pastel2 (8)", "Set1 (9)", "Set2 (8)", "Set3 (12)", "Blues (9)" ,"BuGn (9)" ,
                            "BuPu (9)", "GnBu (9)", "Greens (9)", "Greys (9)", "Oranges (9)", "OrRd (9)", "PuBu (9)",
                            "PuBuGn (9)", "PuRd (9)", "Purples (9)", "RdPu (9)", "Reds (9)","YlGn (9)","YlGnBu (9)",
                            "YlOrBr (9)", "YlOrRd (9)"),
                width="95%"
              ),
              selectInput("nColors","# Colors",choices = list("11"=11,"10"=10,"9"=9,"8"=8, "7"=7,"6"=6,"5"=5,"4"=4,"3"=3), width="95%", selected ="11"),
              sliderInput("alphaSlider", label="Set Brewer Palette Alpha",
                          min = 0,max=1,step=.01,value=1,width = "95%"
              )
            ),
            fillRow(
              textInput("selectedColors",label="Selected Colors"),
              colourInput("colPicker", label= "Color Picker",
                          value = "#00FF0080",
                          allowTransparent = TRUE,
                          closeOnClick = FALSE
              )
            ),
            fillRow(
              plotOutput("colorPlot", height="40px")
            )
          ),
          hr(),
          h4("Points, Lines and Fills"),
          fillCol(
            fillRow(
              textInput("pointColors",label="Point Colors", width="95%"),
              textInput("pointShapes", label="Point Shapes", width="95%")
            ),
            fillRow(
              textInput("lineColors",label="Line Colors",  width="95%"),
              textInput("fillColors",label="Fill Colors", width="95%")
            ),
            fillRow(
              sliderInput("lwdSlider","Line Width", min=0,max=5,step=.2,value = 1,width="95%"),
              sliderInput("width","Lane Width", min=0,max=2,step=.1,value = .65,width="95%"),
              sliderInput("pointLaneWidth","Point Lane Width", min=0,max=2,step=.1,value = .5,width="95%")
            ),
            fillRow(
              textInput("errorBarLineType",label="Error Bar Line Type",width="95%"),
              selectInput("errorBarCapType", label="Error Bar Cap", choices=list("None"="none","Ball"="ball","Bar"="bar"),selected = "Bar",width="95%"),
              sliderInput("errorCapWidth", label="Error Cap Width", min=0,max=1,value=.2,step=.05, width="95%")
            )
          ),
          fillCol(height = "25%",
            fillRow(
              textInput("vioBoxFill", label="Violin Box Plot Color", width="95%"),
              textInput("vioBoxLineCol",label="Violin Box Line Color", width="95%"),
              selectInput("swarmOverflow", label="Swarm Overflow Strategy", choices=list(Random="random",Gutter="gutter",Wrap="wrap",Omit="Omit",None="none"),selected="Wrap",width="95%")
            )
          ),
          hr(),
          h4("Axis and Text"),
          fillCol(
            fillRow(
              selectInput("fontFamily", label="Font Family", choices = list("Sans Serif"="sans","Serif"="serif","Monotype"="mono"),selected="sans",width="95%"),
              textInput("titleCol", label="Title Color",width="95%"),
              sliderInput("titleSize",label="Title Size", min=.5, max=4, step = .1,value=1.5, width="95%")
            ),
            fillRow(
              textInput("labelCol", label="Group Label Color",  width="95%"),
              sliderInput("labelSize", label="Group Label Size", min=.2,max=2,step=0.02, value=.96, width="95%"),
              sliderInput("labelSpacing", label="Group Label Spacing", min=.4,max=2, step=0.02, value=.96, width="95%")
            ),
            fillRow(
              textInput("subGroupLabelCol", label="Subgroup Label Color", width="95%"),
              sliderInput("subGroupLabelSize", label="Subgroup Label Size", min=.2,max=2,step=0.02, value=.68, width="95%"),
              sliderInput("subGroupLabelSpacing", label="Subgroup Label Spacing", min=0.01, max=1.5, step=0.01, value=.26, width="95%")
            ),
            fillRow(
              textInput("dataAxisLablesCol",label="Data (Y-Axis) Label Color", width="95%"),
              textInput("axisCol",label="Axis Color", width="95%"),
              textInput("majorTickColor",label="Axis Tick Color", width="95%")
            )
          ),
          fillCol(height = "50%",
            fillRow(
              textInput("minorTickColor",label="Minor Tick Color", width="95%"),
              textInput("majorGuideColor",label="Guideline Color", width="95%"),
              textInput("minorGuidesColor",label="Minor Guide Color", width="95%")
            ),
            fillRow(
              textInput("axisPrepend", label="Prepend Axis Text",width="95%"),
              textInput("axisAppend", label="Append Axis Text",width="95%"),
              checkboxInput("extendTick", label="Extend Minor Ticks", value=TRUE, width="95%")
            )
          ),
          hr(),
          h4("Background and Legend"),
          fillCol(height = "50%",
            fillRow(
              textInput("canvasColor",label="Plot Area Color", width="95%"),
              textInput("marginCol",label="Background Color", width="95%"),
              textInput("subColor",label="Subtext Color", width="95%")
            ),
            fillRow(
              textInput("legendBorder",label="Legend Border Color", width="95%"),
              textInput("LegendLineCol", label="Legend Box Outline Color", width="95%"),
              textInput("LegendBG", label="Legend Background Color", width="95%")
            )
          )
        )
      ),
      miniTabPanel("Advanced", icon = icon("code"),
        miniContentPanel(padding = 5,
          h4("Misc Options"),
          fillCol(
            fillRow(
              checkboxInput("useRgl", label="Use RGL for 3D",value =TRUE, width="95%"),
              checkboxInput("sidePlot", label="Flip X/Y Axses", value=FALSE, width="95%"),
              sliderInput("curvePoints",label="Curve/Bootstrap Sampling",min=50,max=5000,step=50,value=1000, width="95%")
            ),
            fillRow(
              checkboxInput("verbose", label="Verbose Output", value=FALSE, width="95%"),
              selectInput("calcType", label="Statical Testing",choices=list(None="none",Wilcoxon="wilcox",T.Test="t.test",ANOVA="anvo",Tukey_HSD="Tukey"),selected = "None", width="95%"),
              textInput("bandwith", label="Bandwidth", value=NULL, width="95%")
            ),
            fillRow(
              checkboxInput("strictLimits",label="Strict Data Axis Limits", value=FALSE, width="95%"),
              sliderInput("outliers",label="Outlier IQR Threshold", min=0,max=5,step=.1,value=1.5, width="95%"),
              br()
            )
          ),
          h4("Actual genePlot Code"),
          fillCol(height="30%",
            fillRow(
              verbatimTextOutput("codeExample", placeholder = TRUE)
            )
          ),
          h4("Data Summary"),
          fillCol(
            fillRow(
              tableOutput("descriptive")
            )
          ),
          h4("Statistics"),
          fillCol(
            fillRow(
              verbatimTextOutput("stats", placeholder = TRUE)
            )
          )
        )
      )
    ),
    miniButtonBlock(
      actionButton(inputId = "done",label = "Plot it!",width = "20%"),
      actionButton(inputId = "reset",label = "Reset",width = "20%"),
      actionButton(inputId = "cancel",label = "Cancel",width = "20%")
    )
  )

  server <- function(input, output, session) {

    aGenes<-reactiveValues(g=genes)
    RSOveride<-reactiveValues(rso=TRUE, npData=NULL)
    plotOptions<-reactiveValues(gene=genes,
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
                                extendTicks=TRUE)


    output$plot1 <- renderPlot({
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
                                 extendTicks=plotOptions$extendTicks
                                 )

      if(RSOveride$rso==TRUE) {
        RSOveride$rso<-FALSE
      }
      output$descriptive<-renderTable(RSOveride$npData$summary)
      #utput$Statistics<-RSOveride$npData$stats
    })

    output$colorPlot<-renderPlot( {
      opar<-par(mar=c(0,0,0,0))
      plot(1,1,type="n", xlim=c(0,100), ylim=c(0,10),axes = FALSE, ylab="", xlab="")
      abline(h =5,lwd=5)
      sCol<-renderText(req( {input$selectedColors} ))
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
    observeEvent(input$gene, {
      geneChoices<-NULL
      cSymbol<-renderText({ input$gene })
      if(nchar(cSymbol())>=3) {
        geneChoices<-geneList[grepl(paste0("^",cSymbol()),geneList)]
        if(length(geneChoices)<=40) {
          updateSelectInput(session,"geneChoiceList",choices=geneChoices)
        }
      }
      if(nchar(cSymbol())<=1 & length(geneList)<=30) {
        updateSelectInput(session,"geneChoiceList",choices=geneList)
      }
    })

    observeEvent(input$addGene, {
      if(!is.null(req(input$geneChoiceList))) {
        cGene<-renderText({ req(input$geneChoiceList) })
        if(!is.null(cGene()) & !is.na(cGene()) & length(cGene())>0 & ! cGene() %in% aGenes$g & cGene()!="") {
          aGenes$g<-c(aGenes$g,cGene())
          updateSelectizeInput(session,"activeGenes",choices=aGenes$g)
          plotOptions$gene<-aGenes$g
        }
      }
    })

    observeEvent(input$removeGene, {
      if(!is.null(req( {input$activeGenes }))) {
        aGenes$g<-aGenes$g[! aGenes$g %in% req( {input$activeGenes} ) ]
        updateSelectizeInput(session,"activeGenes",choices=aGenes$g)
        plotOptions$gene<-aGenes$g
      }
    })

    observeEvent(input$groupFactorSelector, {
      if(!is.null(req( {input$groupFactorSelector }))) {
        updateTextInput(session, "group", value = input$groupFactorSelector )
      }
    })

    observeEvent(input$group, {
      if(!is.null(req( {input$group }))) {
        groupVal<-renderText( {input$group} )
        if (groupVal() == req( {input$groupFactorSelector} )) {
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

    observeEvent(input$subGroupFactorSelector, {
      if(!is.null(req( {input$subGroupFactorSelector }))) {
        updateTextInput(session, "subGroup", value = input$subGroupFactorSelector )
      }
    })

    observeEvent(input$subGroup, {
      groupVal<-renderText( {req(input$group)} )
      if(!is.null(req( {input$subGroup })) & groupVal()!="FALSE") {
        subGroupVal<-renderText( {input$subGroup} )
        if(subGroupVal()!="") {
          if(subGroupVal() == req( {input$subGroupFactorSelector} )){
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

    observeEvent(input$highlightFactorSelector, {
      if(!is.null(req( {input$highlightFactorSelector }))) {
        updateTextInput(session, "highlight", value = input$highlightFactorSelector )
      }
    })

    observeEvent(input$highlight, {
      if(!is.null(req( {input$highlight }))) {
        highlightVal<-renderText( {input$highlight} )
        if(highlightVal() == req( {input$highlightFactorSelector} )){
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

    observeEvent(input$stackFactorSelector, {
      if(!is.null(req( {input$stackFactorSelector }))) {
        updateTextInput(session, "stack", value = input$stackFactorSelector )
      }
    })

    observeEvent(input$stack, {
      if(!is.null(req( {input$stack }))) {
        stackVal<-renderText( {input$stack} )
        if(stackVal() == req( {input$stackFactorSelector} )){
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

    observeEvent(input$groupByGene, {
      plotOptions$groupByGene<-input$groupByGene
    })

    observeEvent(input$asPercentage, {
      plotOptions$normalize<-input$asPercentage
    })


    ######Plot Options Minitab
    observeEvent(input$plotType, {
      pt<-renderText(req({input$plotType}))
      RSOveride$rso<-TRUE
      plotOptions$plotType<-pt()
      if(pt() != "density" & pt() != "surface") {
        plotOptions$logScale<-plotOptions$logScale[1]
      }
    })

    observeEvent(input$main, {
      titleText<-renderText( {input$main} )
      if(is.null(titleText()) | titleText() == "") {
        plotOptions$main<-TRUE
      } else {
        plotOptions$main<-titleText()
      }
    })

    observeEvent(input$plotPoints, {
      plotOptions$drawPoints<-input$plotPoints
    })

    observeEvent(input$pointMethod, {
      pm<-renderText(req({input$pointMethod}))
      plotOptions$pointMethod<-pm()
    })

    observeEvent(input$pointSize, {
      plotOptions$pointSize<-input$pointSize
    })

    observeEvent(input$legend, {
      groupVal<- renderText( {req(input$group)} )
      subGroupVal<- renderText( {req(input$subGroup)} )
      highlightVal<- renderText( {req(input$highlight)} )
      if(length(aGenes$g)>1) {
        if(groupVal() != "FALSE" | highlightVal() != FALSE) {
          legendVal<- renderText( {input$legend} )
          if(legendVal() == "TRUE") {
            plotOptions$legend <- "Legend"
          } else {
            plotOptions$legend <- FALSE
          }
        } else {
          plotOptions$legend<-FALSE
          updateCheckboxInput(session,"legend",value=FALSE)
        }
      } else {
        if(highlightVal() != FALSE | (groupVal() !=FALSE & subGroupVal() != FALSE)) {
          legendVal<- renderText( {input$legend} )
          if(legendVal() == "TRUE") {
            plotOptions$legend <- "Legend"
          } else {
            plotOptions$legend <- FALSE
          }
        } else {
          plotOptions$legend<-FALSE
          updateCheckboxInput(session,"legend",value=FALSE)
        }
      }
    })

    observeEvent(input$legendTitle, {
      if(req( {input$legend} ) == TRUE ) {
        legendTitleVal<-renderText(req( {input$legendTitle} ))
        if(is.null(legendTitleVal()) | legendTitleVal() == "") {
          plotOptions$legend<-"Legend"
        } else {
          plotOptions$legend<-legendTitleVal()
        }
      } else {
        plotOptions$legend <- FALSE
      }
    })

    observeEvent(input$legendSize, {
      plotOptions$legendSize<-input$legendSize
    })

    observeEvent(input$logScaleY, {
      logScaleYVal<-renderText(req( {input$logScaleY} ))
      plotTypeVal<-renderText(req( {input$plotType} ))
      if(plotTypeVal() == "density" | plotTypeVal() == "surface" ) {
        if(req( {input$logScaleX} ) != FALSE ) {
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

    observeEvent(input$logScaleX, {
      logScaleXVal<-renderText(req( {input$logScaleX} ))
      plotTypeVal<-renderText(req( {input$plotType} ))
      if(plotTypeVal() == "density" | plotTypeVal() == "surface" ) {
        if(req( {input$logScaleX} ) != FALSE ) {
          plotOptions$logScale[2]<-as.numeric(logScaleXVal())
        } else {
          plotOptions$logScale[2]<-FALSE
        }
      }
    })

    observeEvent(input$logAdjustment, {
      plotOptions$logAdjustment <- input$logAdjustment
    })

    observeEvent(input$ylab, {
      ylabVal<-renderText( {input$ylab} )
      if(ylabVal() == "" | is.null(ylabVal())) {
        plotOptions$ylab<-NULL
      } else {
        plotOptions$ylab<-ylabVal()
      }
    })

    observeEvent(input$expLabels, {
      if( {input$logScaleX} !=FALSE | {input$logScaleY} !=FALSE ) {
        plotOptions$expLabels<- input$expLabels
      } else {
        plotOptions$expLabels<-FALSE
        updateCheckboxInput(session ,inputId = "expLabels",value = FALSE)
      }
    })

    observeEvent(input$rotateY, {
      plotOptions$rotateY<-input$rotateY
    })

    observeEvent(input$groupNames, {
      groupNamesVal <- renderText( {input$groupNames} )
      if(is.null(groupNamesVal()) | groupNamesVal() == "" ){
        plotOptions$groupNames<-NULL
      } else {
        plotOptions$groupNames<-trimws(unlist(strsplit(groupNamesVal(),",")))
      }
    })

    observeEvent(input$subGroupNames, {
      subGroupNamesVal <- renderText( {input$subGroupNames} )
      if(is.null(subGroupNamesVal()) | subGroupNamesVal() == "" ){
        plotOptions$subGroupLabels<-NULL
      } else {
        plotOptions$subGroupLabels<-trimws(unlist(strsplit(subGroupNamesVal(),",")))
      }
    })

    observeEvent(input$rotateLabels, {
      plotOptions$rotateLabels<-input$rotateLabels
    })

    observeEvent(input$minorTick, {
      plotOptions$minorTick<-input$minorTick
    })

    observeEvent(input$showMinorGuide, {
      plotOptions$minorGuides<-input$showMinorGuide
    })

    observeEvent(input$guides, {
      plotOptions$guides=input$guides
    })

    observeEvent(input$aggFun, {
      aggFunVal<-renderText(req( {input$aggFun} ))
      plotOptions$aggFun<-aggFunVal()
    })

    observeEvent(input$errFun, {
      errFunVal<-renderText(req( {input$errFun} ))
      plotOptions$errFun<-errFunVal()
    })

    observeEvent(input$errorMultiple, {
      plotOptions$errorMultiple <- input$errorMultiple
    })

    observeEvent(input$errorBars,{
      plotOptions$errorBars<-input$errorBars
    })

    observeEvent(input$rug, {
      plotOptions$drawRug<-input$rug
    })

    observeEvent(input$phi, {
      plotOptions$phi<-input$phi
    })

    observeEvent(input$theta, {
      plotOptions$theta<-input$theta
    })

    observeEvent(input$trimCurves, {
      plotOptions$trimCurves<-input$trimCurves
      plotOptions$trimViolins<-input$trimCurves
    })

    observeEvent(input$vioBoxWidth, {
      plotOptions$vioBoxWidth<-input$vioBoxWidth
    })

    #####Format Minitab
    observeEvent(input$theme, {
      cTheme<-renderText(req({input$theme}))
      plotOptions$theme<-cTheme()
    })

    observeEvent(input$colorBrewer, {
      cbVal<-renderText(req({input$colorBrewer}))
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
        updateTextInput(session,inputId="selectedColors",value=paste0(cbVal,collapse=","))
        nVect<-rev(seq(3,maxC))
        names(nVect)<-nVect

        updateSelectInput(session, inputId="nColors", choices=as.list(nVect), selected = max(nVect))
      }
    })

    observeEvent(input$nColors, {
      nColVal<-renderText(req( {input$nColors} ))
      nColVal<-as.numeric( nColVal() )
      cbVal<-renderText(req({input$colorBrewer}))
      pallete<-gsub("(.*) \\((\\d+)\\)","\\1",cbVal())
      cbVal<-brewer.pal(name=pallete,n=nColVal)
      cbVal<-map_chr(cbVal,setAlpha, alpha=input$alphaSlider)
      updateTextInput(session,inputId="selectedColors",value=paste0(cbVal,collapse = ","))
    })

    observeEvent(input$alphaSlider, {
      scVal<-renderText({input$selectedColors})
      if(!is.null(scVal()) & scVal() != ""){
        scVal<- map_chr(trimws(unlist(strsplit( scVal(), split=","))), setAlpha, alpha=input$alphaSlider)
        updateTextInput(session,inputId="selectedColors",value=paste0(scVal,collapse = ","))
      }
    })

    observeEvent(input$pointColors, {
      cTheme<-renderText({input$theme})
      cVal<-paste0(eval(parse(text=cTheme()))$plotColors$points,collapse = ",")
      cCol<-renderText(input$pointColors)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$points<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$points<- cVal
      }
    })

    observeEvent(input$pointShapes, {
      cTheme<-renderText({input$theme})
      cVal<-paste0(eval(parse(text=cTheme()))$pointShape,collapse = ",")
      cPch<-renderText(input$pointShapes)
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



    observeEvent(input$lineColors, {
      cTheme<-renderText({input$theme})
      cVal<-paste0(eval(parse(text=cTheme()))$plotColors$lines,collapse = ",")
      cCol<-renderText(input$lineColors)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$lines<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$lines<- cVal
      }
    })

    observeEvent(input$fillColors, {
      cTheme<-renderText({input$theme})
      cVal<-paste0(eval(parse(text=cTheme()))$plotColors$fill,collapse = ",")
      cCol<-renderText(input$fillColors)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$fill<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$fill<- cVal
      }
    })

    observeEvent(input$lwdSlider, {
      plotOptions$lWidth<-req(input$lwdSlider)
    })

    observeEvent(input$width, {
      plotOptions$width<-req(input$width)
    })

    observeEvent(input$pointLaneWidth, {
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

    observeEvent(input$errorBarCapType, {
      cVal<-renderText(input$errorBarCapType)
      plotOptions$errorCapType<-cVal()
    })

    observeEvent(input$errorCapWidth, {
      plotOptions$errorBarCapWidth<-input$errorCapWidth
    })

    observeEvent(input$vioBoxFill, {
      cTheme<-renderText({input$theme})
      cVal<-paste0(eval(parse(text=cTheme()))$plotColors$vioBoxFill,collapse = ",")
      cCol<-renderText(input$vioBoxFill)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$vioBoxFill<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$vioBoxFill<- cVal
      }
    })

    observeEvent(input$vioBoxLineCol, {
      cTheme<-renderText({input$theme})
      cVal<-paste0(eval(parse(text=cTheme()))$plotColors$vioBoxLineCol,collapse = ",")
      cCol<-renderText(input$vioBoxLineCol)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$vioBoxLineCol<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$vioBoxLineCol<- cVal
      }
    })

    observeEvent(input$swarmOverflow, {
      cVal<-renderText(input$swarmOverflow)
      plotOptions$swarmOverflow<-cVal()
    })

    observeEvent(input$fontFamily, {
      cVal<-renderText(input$fontFamily)
      plotOptions$fontFamily<-cVal()
    })

    observeEvent(input$titleCol, {
      cTheme<-renderText({input$theme})
      cVal<-paste0(eval(parse(text=cTheme()))$plotColors$title,collapse = ",")
      cCol<-renderText(input$titleCol)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$title<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$title<- cVal
      }
    })

    observeEvent(input$titleSize, {
      plotOptions$titleSize<-input$titleSize
    })

    observeEvent(input$labelCol, {
      cTheme<-renderText({input$theme})
      cVal<-paste0(eval(parse(text=cTheme()))$plotColors$labels,collapse = ",")
      cCol<-renderText(input$labelCol)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$labels<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$labels<- cVal
      }
    })

    observeEvent(input$labelSize, {
      plotOptions$groupLabSize<-input$labelSize
    })

    observeEvent(input$labelSpacing, {
      plotOptions$groupLabelSpacing<-input$labelSpacing
    })

    observeEvent(input$subGroupLabelCol, {
      cTheme<-renderText({input$theme})
      cVal<-paste0(eval(parse(text=cTheme()))$plotColors$subGroupLabels,collapse = ",")
      cCol<-renderText(input$subGroupLabelCol)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$subGroupLabels<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$subGroupLabels<- cVal
      }
    })

    observeEvent(input$subGroupLabelSize, {
      plotOptions$subGroupLabSize<-input$subGroupLabelSize
    })

    observeEvent(input$subGroupLabelSpacing, {
      plotOptions$subgroupLabelSpacing<-input$subGroupLabelSpacing
    })

    observeEvent(input$dataAxisLablesCol, {
      cTheme<-renderText({input$theme})
      cVal<-paste0(eval(parse(text=cTheme()))$plotColors$numbers,collapse = ",")
      cCol<-renderText(input$dataAxisLablesCol)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$numbers<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$numbers<- cVal
      }
    })

    observeEvent(input$axisCol, {
      cTheme<-renderText({input$theme})
      cVal<-paste0(eval(parse(text=cTheme()))$plotColors$axis,collapse = ",")
      cCol<-renderText(input$axisCol)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$axis<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$axis<- cVal
      }
    })

    observeEvent(input$majorTickColor, {
      cTheme<-renderText({input$theme})
      cVal<-paste0(eval(parse(text=cTheme()))$plotColors$majorTick,collapse = ",")
      cCol<-renderText(input$majorTickColor)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$majorTick<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$majorTick<- cVal
      }
    })

    observeEvent(input$minorTickColor, {
      cTheme<-renderText({input$theme})
      cVal<-paste0(eval(parse(text=cTheme()))$plotColors$minorTick,collapse = ",")
      cCol<-renderText(input$minorTickColor)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$minorTick<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$minorTick<- cVal
      }
    })

    observeEvent(input$majorGuideColor, {
      cTheme<-renderText({input$theme})
      cVal<-paste0(eval(parse(text=cTheme()))$plotColors$guides,collapse = ",")
      cCol<-renderText(input$majorGuideColor)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$guides<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$guides<- cVal
      }
    })

    observeEvent(input$minorGuidesColor, {
      cTheme<-renderText({input$theme})
      cVal<-paste0(eval(parse(text=cTheme()))$plotColors$minorGuides,collapse = ",")
      cCol<-renderText(input$minorGuidesColor)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$minorGuides<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$minorGuides<- cVal
      }
    })

    observeEvent(input$axisPrepend, {
      preVal<-renderText(input$axisPrepend)
      aVal<-renderText(input$axisAppend)
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

    observeEvent(input$axisAppend, {
      preVal<-renderText(input$axisPrepend)
      aVal<-renderText(input$axisAppend)
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

    observeEvent(input$extendTick, {
      plotOptions$extendTick<-input$extendTick
    })

    observeEvent(input$canvasColor, {
      cTheme<-renderText({input$theme})
      cVal<-paste0(eval(parse(text=cTheme()))$plotColors$bg,collapse = ",")
      cCol<-renderText(input$canvasColor)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$bg<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$bg<- cVal
      }
    })

    observeEvent(input$marginCol, {
      cTheme<-renderText({input$theme})
      cVal<-paste0(eval(parse(text=cTheme()))$plotColors$marginBg,collapse = ",")
      cCol<-renderText(input$marginCol)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$marginBg<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$marginBg<- cVal
      }
    })

    observeEvent(input$subColor, {
      cTheme<-renderText({input$theme})
      cVal<-paste0(eval(parse(text=cTheme()))$plotColors$subtext,collapse = ",")
      cCol<-renderText(input$subColor)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$subtext<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$subtext<- cVal
      }
    })

    observeEvent(input$legendBorder, {
      cTheme<-renderText({input$theme})
      cVal<-paste0(eval(parse(text=cTheme()))$LegendBorder,collapse = ",")
      cCol<-renderText(input$legendBorder)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$LegendBorder<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$LegendBorder<- cVal
      }
    })

    observeEvent(input$LegendLineCol, {
      cTheme<-renderText({input$theme})
      cVal<-paste0(eval(parse(text=cTheme()))$LegendLineCol,collapse = ",")
      cCol<-renderText(input$LegendLineCol)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$LegendLineCol<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$LegendLineCol<- cVal
      }
    })

    observeEvent(input$LegendBG, {
      cTheme<-renderText({input$theme})
      cVal<-paste0(eval(parse(text=cTheme()))$LegendBG,collapse = ",")
      cCol<-renderText(input$LegendBG)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$LegendBG<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$LegendBG<- cVal
      }
    })

    #Lower Button Panel
    observeEvent(input$done, {
      stopApp(RSOveride$npData)
    })
    observeEvent(input$cancel, {
      stopApp(NULL)
    })
  }
  viewer<-dialogViewer(dialogName = "Interactive genePlot UI",height = 2000)
  runGadget(ui, server, viewer=viewer)
}
