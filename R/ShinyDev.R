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
#' @description A shiny widget for bvt's genePlot
#'
#' @details
#' This is a the optional gui interface for using genePlot. It is an internal function that is not exported
#' The interface returns the options selected by the user. The GUI allows for extensive customization and returns
#' the effective plotting command and npData object for plotting. This is due to the way that RStudio handles graphics a
#' and seems to have issues with immediately ploting the data. For this reason the data is plotted in base R
#' but only the npData object is returned in RStudio. This can to plotted using plot or genePlot or any nicePlots function.
#'
#' @param data \code{\link[Biobase]{ExpressionSet}} or other R data object containing the data and possibly feature or phenotype information. While intended for use with Bioconductor objects, \code{\link[base]{data.frame}}, \code{\link[base]{matrix}}, and \code{\link[tibble]{tibble}} can also be used.
#' @param geneList character vector; List of all possible gene names including expression data row names and gene symbols if in use.
#' @param factorList all possible factors from phenotype data.
#' @param gpOptions active options from the command line.
#' @param dbName character; Name of the database passed to \code{\link[bvt]{genePlot}}. If using the \code{\link[datasets]{iris}} data set this should be set to "iris".
#' @param themeName character; Name of of the selected npTheme as a character string. Used to set the default value in shiny. As \code{\link[base]{eval}} and \code{\link[base]{parse}} are used to access the theme settings any npTheme object in the enviroment will do, including custom themes.
#'
#' @examples
#' ToDo<-1
#'
#' @importFrom NicePlots setAlpha basicTheme npThemes
#' @importFrom purrr map_chr
#' @importFrom RColorBrewer brewer.pal
#' @importFrom graphics abline par rect
#' @seealso \code{\link[bvt]{genePlot}}
shinyGenePlot <- function(data, geneList, factorList, gpOptions,dbName="data",themeName="basicTheme") {
  if(requireNamespace("shiny",quietly = TRUE) & requireNamespace("colourpicker",quietly = TRUE) & requireNamespace("miniUI",quietly = TRUE) == FALSE){
    stop("Missing required libraries for interactive shiny UI. Please install shiny, miniUI and colourpicker.")
  }
  iTheme<-basicTheme
  if(!is.null(gpOptions$theme)) {
    iTheme<-gpOptions$theme
  }
  IPT<-"BP"
  if(!is.null(gpOptions$plotType[1])){
    if(grepl("box",ignore.case = T, x = gpOptions$plotType[1])) {
      IPT<-"BP"
    }
    else if(grepl("vio",ignore.case = T, x = gpOptions$plotType[1])) {
      IPT<-"VP"
    }
    else if(grepl("dot",ignore.case = T, x = gpOptions$plotType[1])) {
      IPT<-"DP"
    }
    else if(grepl("bar",ignore.case = T, x = gpOptions$plotType[1])) {
      IPT<-"Bar"
    }
    else {
      IPT<-"2D"
    }
  }
  #The above code is reading the plot type and theme information so we can access the right theme settings
  #for determining defaults. This is a very long function with a lot of UI element so a quick guide to the
  #code is probably in order. This section defines the UI elements using miniUI::miniContentPanel,
  #The elements below the plot preview window are broken into tabs (ie miniTabstripPanel) and are listed
  #more or less in order of the tabs from left to right and from the top of the tab page to the bottom.
  factorList<-append(list("None"=FALSE),as.list(factorList))
  themes<-npThemes()
  ui <- miniUI::miniPage(
    miniUI::miniContentPanel(padding = 0,
     shiny::plotOutput("plot1", height = "100%")
    ),
    miniUI::miniTabstripPanel(
      ############## Begin Data Tab Section ####################3
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
              shiny::selectizeInput("activeGenes",label="Selected Data", multiple=TRUE, choices=gpOptions$gene)
            ),
            shiny::fillCol(
              style = "margin-top: 25px;",
              shiny::actionButton("removeGene",label="",icon=shiny::icon("trash"))
            )
          ),
          shiny::fillRow(flex=c(7,1,5),
              shiny::fillCol(
              shiny::selectInput("groupFactorSelector", label="Select a Grouping Factor",choices=factorList,selected = if(is.null(gpOptions$group)){"None"} else if (gpOptions$group %in% factorList){gpOptions$group}else{"None"}),
              shiny::selectInput("subgroupFactorSelector", label="Select a Subgroup Factor",choices=factorList,selected = if(is.null(gpOptions$subgroup)){"None"} else if (gpOptions$subgroup %in% factorList){gpOptions$subgroup}else{"None"}),
              shiny::selectInput("highlightFactorSelector", label="Select a Highlighting Factor",choices=factorList,selected = if(is.null(gpOptions$highlight)){"None"} else if (gpOptions$highlight %in% factorList){gpOptions$highlight}else{"None"}),
              shiny::selectInput("stackFactorSelector", label="Select a Stacking Factor",choices=factorList,selected = if(is.null(gpOptions$stack)){"None"} else if (gpOptions$stack %in% factorList){gpOptions$stack}else{"None"}),
              shiny::checkboxInput("groupByGene",label="Group by Gene", value=if(is.null(gpOptions$groupByGene)){TRUE}else{gpOptions$groupByGene})
            ),
            shiny::fillCol(
              shiny::br(),
              shiny::br(),
              shiny::br(),
              shiny::br(),
              shiny::br()
            ),
            shiny::fillCol(
              shiny::textInput("group", label = "Group", value = gpOptions$group),
              shiny::textInput("subgroup", label = "Subgroup", value = gpOptions$subgroup),
              shiny::textInput("highlight", label = "Highlight", value = gpOptions$highlight),
              shiny::textInput("stack", label = "Stack (Bar Plot Only)", value = gpOptions$stack),
              shiny::checkboxInput("asPercentage",label="As Percent (Stacked Bars)", value= if(is.null(gpOptions$asPercentage)){FALSE}else{gpOptions$asPercentage})
            )
          )
        )
      ),
      ############### Begin Plot Options Tab Section ###########################
      miniUI::miniTabPanel("Plot Options", icon = shiny::icon("chart-area"),
        miniUI::miniContentPanel(padding = 5,
          shiny::fillCol(height = "85%",
            shiny::fillRow(
              shiny::selectInput("plotType", label = ("Select PlotType"),
                choices = list("Box Plot" = "box", "Dot Plot" = "dot", "Violin Plot" = "violin", "Bar Plot"="bar", "Density Plot"="density", "Surface Plot"="surface"),
                selected = if(is.null(gpOptions$plotType)){"box"}else if(any(grepl(gpOptions$plotType[1], c("box","dot","violin","bar","density","surface"), ignore.case = TRUE))) {c("box","dot","violin","bar","density","surface")[grep(gpOptions$plotType[1],ignore.case = TRUE, x = c("box","dot","violin","bar","density","surface"))[1]]}else {"box"},width="95%"),
              shiny::textInput("main",label="Plot Title",width="95%",value = if(is.null(gpOptions$main)){NULL}else if (gpOptions$main==paste0("Gene Expression: ",paste0(gpOptions$gene,collapse=", "))){NULL}else{gpOptions$main})
            ),
            shiny::fillRow(
              shiny::checkboxInput("plotPoints", label = "Plot Points", value = if(is.null(gpOptions$plotPoints)){TRUE}else{gpOptions$plotPoints},width="95%"),
              shiny::selectInput("pointMethod", label = "Point Ploting Style",
                          choices = list("Jitter" = "jitter", "Beeswarm" = "beeswarm", "Distribution" = "distribution", "Linear"="linear"),
                          selected = if(is.null(gpOptions$pointMethod)){if(IPT %in% c("2D","Bar")){iTheme$pointMethodBP} else {iTheme[[paste0("pointMethod",IPT)]]}} else if(any(grepl(gpOptions$pointMethod ,c("jitter","beeswarm","distribution","linear"), ignore.case = TRUE))){c("jitter","beeswarm","distribution","linear")[grep(gpOptions$pointMethod ,c("jitter","beeswarm","distribution","linear"), ignore.case = TRUE)[1]]}else{"jitter"},width="95%"),
              shiny::sliderInput("pointSize", label = "Point Size",
                          min=.05,max=3,value=if(is.null(gpOptions$pointSize)){if(IPT=="Bar"){iTheme$pointSizeBP}else{iTheme[[paste0("pointSize",IPT)]]}}else{gpOptions$pointSize},width="95%",step=.05)
            ),
            shiny::fillRow(
              shiny::checkboxInput("legend", label = "Plot Legend", value = if(is.null(gpOptions$legend)){FALSE}else{if(gpOptions$legend==FALSE){FALSE}else{TRUE}},width="95%"),
              shiny::textInput("legendTitle", label = "Legend Title", value = if(is.null(gpOptions$legend)){"Legend"}else{if(gpOptions$legend==FALSE){"Legend"}else{gpOptions$legend}},width="95%"),
              shiny::sliderInput("legendSize", label="Scale Legend", min=.25,max=2, value=if(is.null(gpOptions$legendSize)){iTheme$legendSize}else{gpOptions$legendSize},width="95%")
            )
          ),
          shiny::fillCol(height = "85%",
            shiny::fillRow(
              shiny::selectInput("logScaleY", label = "Log Scale",
                                 choices = list("None"= FALSE, "2" = 2, "10" =10, "e" = exp(1)),
                                 selected = if(is.null(gpOptions$logScale)){"None"} else {if(gpOptions$logScale==2){"2"} else if(gpOptions$logScale==10){"10"} else if(gpOptions$logScale=="e"){"e"}else{"None"}},width="95%"),
              shiny::selectInput("logScaleX", label = "Log Scale X (Contour)",
                                 choices =  list("None"= FALSE, "2" = 2, "10" =10, "e" = exp(1)),
                                 selected = FALSE,width="95%"),
              shiny::sliderInput("logAdjustment", label="Log Adjustment", min=0,max=1, step=0.05, value=if(is.null(gpOptions$logAdjustment)){1}else{gpOptions$logAdjustment},width="95%",)
            ),
            shiny::fillRow(
              shiny::textInput("ylab",label="Y-Axis Label",
                               value = if(is.null(gpOptions$ylab)){NULL}else if (gpOptions$ylab==""){NULL} else {gpOptions$ylab},
                               width="95%"),
              shiny::textInput("subtitle", value = if(is.null(gpOptions$subtitle)){NULL}else if (gpOptions$subtitle==""){NULL} else {gpOptions$subtitle}, label="Subtitle",width="95%"),
              shiny::checkboxInput("rotateY", label = "Rotate Data Labels", value = if(is.null(gpOptions$rotateY)){FALSE} else {if(gpOptions$rotateY==TRUE){TRUE}else{FALSE}},width="95%")
            ),
            shiny::fillRow(
              shiny::textInput("groupNames",value=if(is.null(gpOptions$groupNames)){NULL}else if (gpOptions$groupNames==""){NULL} else {paste0(gpOptions$groupNames, collapse = ",")}, label="Group Labels",width="95%"),
              shiny::textInput("subgroupNames",label="Subgroup Labels",value=if(is.null(gpOptions$subgroupNames)){NULL}else if (gpOptions$subgroupNames==""){NULL} else {paste0(gpOptions$subgroupNames, collapse = ", ")},width="95%"),
              shiny::checkboxInput("rotateLabels", label = "Rotate Group Labels",value=if(is.null(gpOptions$rotateLabels)){FALSE} else {if(gpOptions$rotateLabels==TRUE){TRUE}else{FALSE}}, width="95%")
            )
          ),
          shiny::fillCol(height = "85%",
            shiny::fillRow(
              shiny::numericInput("minorTick",label="Minor Ticks",min=0,value=if(!is.null(gpOptions$minorTick)){as.numeric(gpOptions$minorTick)} else if(gpOptions$logScale == FALSE) {as.numeric(iTheme$minorTick)}else{as.numeric(iTheme$minorTickLS)}, step=1, width="95%"),
              shiny::checkboxInput("showMinorGuide", label = "Minor Guides", value = if(is.null(gpOptions$minorGuides)){if(is.null(gpOptions$guides)){iTheme$guides}else{gpOptions$guides}}else{gpOptions$minorGuides},width="95%"),
              shiny::checkboxInput("guides", label = "Major Guides", value = if(is.null(gpOptions$guides)){iTheme$guides}else{gpOptions$guides},width="95%")
            ),
            shiny::fillRow(
              shiny::selectInput("aggFun", label = "Central Tendancy",
                      choices = list("Median"= "median", "Mean" = "mean"),
                      selected = if(is.null(gpOptions$aggFun)){"mean"} else if (grepl(gpOptions$aggFun, ignore.case = TRUE,"median")){"median"}else{"mean"},
                      width="95%"),
              shiny::selectInput("errFun", label = "Error Measure",
                      choices = list("Standard Deviation"= "sd", "Standard Error" = "se","Range"="range","95-CI T-Distribution"="t95ci","95-CI Bootstrap"="boot95ci"),
                      selected = if(is.null(gpOptions$errFun)){"se"}else if (any(grepl(gpOptions$errFun,ignore.case = TRUE, c("sd","se","range","t95ci","boot95ci")))){grep(gpOptions$errFun,ignore.case = TRUE, c("sd","se","range","t95ci","boot95ci"),value = TRUE)[1]}else{"se"},
                      width="95%"),
              shiny::numericInput("errorMultiple", label="Error Multiplier", value=if(is.null(gpOptions$errorMultiple)){2}else{gpOptions$errorMultiple},min=0,width="95%",step = 1)
            ),
            shiny::fillRow(
              shiny::checkboxInput("errorBars", label = "Draw Error Bars", value = if(is.null(gpOptions$errorBars)){TRUE}else{gpOptions$errorBars}, width="95%"),
              shiny::sliderInput("phi", label = "Phi Rotation (persp)", min=-180,max=180, value=if(is.null(gpOptions$phi)){30}else{gpOptions$phi},width="95%"),
              shiny::sliderInput("theta", label="Theta Rotation (persp)", min=-180,max=180, value=if(is.null(gpOptions$theta)){30}else{gpOptions$theta},width="95%")
            )
          ),
          shiny::fillCol(height = "85%",
            shiny::fillRow(
              shiny::checkboxInput("rug", label = "Density Rug Plot", value = if(is.null(gpOptions$drawRug)){FALSE}else{gpOptions$drawRug},width="95%"),
              shiny::checkboxInput("trimCurves", label="Trim Curves to Data", value=if(is.null(gpOptions$trimCurves) & is.null(gpOptions$trimViolins)){TRUE}else if(is.null(gpOptions$trimCurves)){gpOptions$trimViolins}else{gpOptions$trimCurves}, width="95%"),
              shiny::sliderInput("vioBoxWidth", label="Violin Box Width", min=0,max=1.5,step=.1,value=if(is.null(gpOptions$vioBoxWidth)){iTheme$vioBoxWidth}else{gpOptions$vioBoxWidth},width="95%")
            ),
            shiny::fillRow(
              shiny::textInput("yLim",label="Data Limits (min,max)",value = if(is.null(gpOptions$yLim)){NULL}else{paste0(gpOptions$yLim,collapse = ", ")}, width="95%"),
              shiny::textInput("xLim",label="Y Axis Limits",value = if(is.null(gpOptions$xLim)){NULL}else{paste0(gpOptions$xLim,collapse = ", ")}, width="95%"),
              shiny::checkboxInput("sidePlot", label="Flip X/Y Axses", value=if(is.null(gpOptions$sidePlot)){FALSE}else{gpOptions$sidePlot}, width="95%")
            ),
            shiny::fillRow(
              shiny::checkboxInput("expLabels", label = "Exp Labels (Log)", value = if(is.null(gpOptions$expLabels)){FALSE}else{gpOptions$expLabels},width="95%"),
              shiny::selectInput("barType", label="Bar Type (Dot Plot)",
                                   choices = list("Bar"= "bar", "Dot" = "dot","None"="none"),
                                   selected = if(is.null(gpOptions$barType)){"bar"}else if(grepl(gpOptions$barType, ignore.case = TRUE,"none")){3} else if(grepl(gpOptions$barType, ignore.case = TRUE,"dot")){2}else{1},
                                   width="95%"),
              shiny::checkboxInput("drawBox", label="Draw Box (Box & Vio Plots)", value= if(is.null(gpOptions$drawBox)){TRUE}else{gpOptions$drawBox}, width="95%")
            )
          )
        )
      ),
      ################## Begin Format Tab Section ########################################
      miniUI::miniTabPanel("Format", icon = shiny::icon("sliders"),
        miniUI::miniContentPanel(padding = 5,
          shiny::fillCol(
            shiny::fillRow(
              shiny::h4("Themes and Color Selection") ,
              shiny::selectInput("theme", label = "Select a Theme",
                        choices = themes,
                        selected = if(themeName %in% themes){themeName} else {"basicTheme"},
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
          shiny::fillCol(height="85%",
            shiny::fillRow(
              shiny::textInput("pointColors",label="Point Colors", value=if(is.null(gpOptions$plotColors)){NULL}else if(all(gpOptions$plotColors$points==iTheme$plotColors$points)){NULL}else{paste0(gpOptions$plotColors$points, collapse = ", ")}, width="95%"),
              shiny::textInput("pointShapes", label="Point Shapes", value=if(is.null(gpOptions$pointShape)){NULL} else if(IPT=="Bar"){NULL} else if(all(gpOptions$pointShape==iTheme[[paste0("pointShape",IPT)]])){NULL}else{paste0(gpOptions$pointShape, collapse = ", ")}, width="95%")
            ),
            shiny::fillRow(
              shiny::textInput("lineColors",label="Line Colors", value=if(is.null(gpOptions$plotColors)){NULL}else if(all(gpOptions$plotColors$lines==iTheme$plotColors$lines)){NULL}else{paste0(gpOptions$plotColors$lines, collapse = ", ")},  width="95%"),
              shiny::textInput("fillColors",label="Fill Colors", value=if(is.null(gpOptions$plotColors)){NULL}else if(all(gpOptions$plotColors$f==iTheme$plotColors$fill)){NULL}else{paste0(gpOptions$plotColors$fill, collapse = ", ")}, width="95%")
            ),
            shiny::fillRow(
              shiny::sliderInput("lWidth","Line Width", min=0,max=5,step=.2,value = if(is.null(gpOptions$lWidth)){iTheme[[paste0("lWidth",IPT)]]}else{gpOptions$lWidth},width="95%"),
              shiny::sliderInput("width","Lane Width", min=0,max=2,step=.1,value = if(is.null(gpOptions$width)){if(IPT == "2D"){iTheme$widthBP}else{iTheme[[paste0("width",IPT)]]}}else{gpOptions$width}, width="95%"),
              shiny::sliderInput("pointLaneWidth","Point Lane Width", min=0,max=5,step=.1,value = if(is.null(gpOptions$pointLaneWidth)){if(IPT %in% c("Bar", "2D")){NULL}else{iTheme[[paste0("pointLaneWidth",IPT)]]}}else{gpOptions$pointLaneWidth},width="95%")
            )

          ),
          shiny::fillCol(height = "55%",
            shiny::fillRow(
              shiny::textInput("errorBarLineType",label="Error Bar Line Type",value = if(is.null(gpOptions$errorBarLineType)){iTheme[[paste0("errorBarLineType",IPT)]]}else{gpOptions$errorBarLineType}, width="95%"),
              shiny::selectInput("errorCapType", label="Error Bar Cap",
                                 choices=list("None"="none","Ball"="ball","Bar"="bar"),
                                 selected = if(is.null(gpOptions$errorCapType)){iTheme$errorCapType}else if (any(grepl(gpOptions$errorCapType, c("none","bar","ball"), ignore.case = TRUE))) {grep(gpOptions$errorCapType, c("none", "ball","bar"), ignore.case = TRUE, value=TRUE)[1]}else {iTheme$errorCapType},
                                 width="95%"),
              shiny::sliderInput("errorCapWidth", label="Error Cap Width", min=0,max=1,value=if(is.null(gpOptions$errorBarCapWidth)){if(IPT=="2D"){iTheme$errorBarCapWidthBP}else{iTheme[[paste0("errorBarCapWidth",IPT)]]}}else{gpOptions$errorBarCapWidth},step=.05, width="95%")
            ),
            shiny::fillRow(
              shiny::textInput("vioBoxFill", label="Violin Box Plot Color", value=if(is.null(gpOptions$plotColors)){NULL}else if(all(gpOptions$plotColors$vioBoxFill==iTheme$plotColors$vioBoxFill)){NULL}else{paste0(gpOptions$plotColors$vioBoxFill, collapse = ", ")}, width="95%"),
              shiny::textInput("vioBoxLineCol",label="Violin Box Line Color", value=if(is.null(gpOptions$plotColors)){NULL}else if(all(gpOptions$plotColors$vioBoxLineCol==iTheme$plotColors$vioBoxLineCol)){NULL}else{paste0(gpOptions$plotColors$vioBoxLineCol, collapse = ", ")}, width="95%"),
              shiny::selectInput("swarmOverflow", label="Swarm Overflow Strategy",
                                 choices=list(Random="random",Gutter="gutter",Wrap="wrap",Omit="omit",None="none"),
                                 selected=if(is.null(gpOptions$swarmOverflow)){iTheme$swarmOverflow} else if (grepl(gpOptions$swarmOverflow, c("random","gutter","wrap","omit","none"), ignore.case = TRUE)){grep(gpOptions$swarmOverflow, c("random","gutter","wrap","omit","none"),ignore.case = TRUE)[1]}else{iTheme$swarmOverflow},
                                 width="95%")
            )
          ),
          shiny::hr(),
          shiny::h4("Axis and Text"),
          shiny::fillCol(height = "85%",
            shiny::fillRow(
              shiny::selectInput("fontFamily", label="Font Family",
                                 choices = list("Sans Serif"="sans","Serif"="serif","Monotype"="mono"),
                                 selected=if(is.null(gpOptions$fontFamily)){iTheme$fontFamily} else if (any(grepl(gpOptions$fontFamily, c("sans","serif","mono"), ignore.case = TRUE))){grep(gpOptions$fontFamily, c("sans","serif","mono"), ignore.case = TRUE,value=TRUE)}else{iTheme$fontFamily},
                                 width="95%"),
              shiny::textInput("titleCol", label="Title Color",value=if(is.null(gpOptions$plotColors)){NULL}else if(is.null(gpOptions$plotColors$title)){NULL}else if(gpOptions$plotColors$title[1]==iTheme$plotColors$title[1]){NULL}else{gpOptions$plotColors$title[1]},width="95%"),
              shiny::sliderInput("titleSize",label="Title Size", min=.5, max=4, step = .1,value=if(is.null(gpOptions$titleSize)){iTheme$titleSize}else{gpOptions$titleSize}, width="95%")
            ),
            shiny::fillRow(
              shiny::textInput("labelCol", label="Group Label Color", value=if(is.null(gpOptions$plotColors)){NULL} else if(is.null(gpOptions$plotColors$labels)){NULL} else {paste0(gpOptions$plotColors$labels, collapse = ", ")}, width="95%"),
              shiny::sliderInput("groupLabSize", label="Group Label Size", min=.2,max=2,step=0.02, value=if(is.null(gpOptions$groupLabSize)){iTheme$groupLabSize}else{gpOptions$groupLabSize}, width="95%"),
              shiny::sliderInput("groupLabelSpacing", label="Group Label Spacing", min=.4,max=2, step=0.02, value=if(is.null(gpOptions$groupLabelSpacing)){iTheme$groupLabelSpacing}else{gpOptions$groupLabelSpacing}, width="95%")
            ),
            shiny::fillRow(
              shiny::textInput("subgroupLabelCol", label="Subgroup Label Color",value=if(is.null(gpOptions$plotColors)){NULL} else if(is.null(gpOptions$plotColors$subgroupLabels)){NULL} else {paste0(gpOptions$plotColors$subgroupLabels, collapse = ", ")}, width="95%"),
              shiny::sliderInput("subgroupLabSize", label="Subgroup Label Size", min=.2,max=2,step=0.02, value=if(is.null(gpOptions$subgroupLabSize)){iTheme$subgroupLabSize}else{gpOptions$subgroupLabSize}, width="95%"),
              shiny::sliderInput("subgroupLabelSpacing", label="Subgroup Label Spacing", min=0.01, max=1.5, step=0.01, value=if(is.null(gpOptions$subgroupLabelSpacing)){iTheme$subgroupLabelSpacing}else{gpOptions$subgroupLabelSpacing}, width="95%")
            )
          ),
          shiny::fillCol(height = "85%",
            shiny::fillRow(
              shiny::textInput("yAxisLabelCol",label="Y-Axis Label Color",value=if(is.null(gpOptions$plotColors)){NULL} else if(is.null(gpOptions$plotColors$axisLabels)){NULL} else {gpOptions$plotColors$axisLabels[1]}, width="95%"),
              shiny::sliderInput("axisLabelSize",label = "Y-Axis Label Size", max = 3, min = 0 ,value=if(is.null(gpOptions$axisLabelSize)){iTheme$axisLabelSize}else{gpOptions$axisLabelSize},step = 0.1, width="95%"),
              shiny::textInput("axisPrepend", label="Prepend Axis Text", value=if(is.null(gpOptions$axisText)){NULL}else{gpOptions$axisText[1]},width="95%")
            ),
            shiny::fillRow(
              shiny::textInput("subtitleCol",label="Subtitle Color",value=if(is.null(gpOptions$plotColors)){NULL} else if(is.null(gpOptions$plotColors$subtext)){NULL} else {gpOptions$plotColors$subtext[1]}, width="95%"),
              shiny::sliderInput("subSize",label = "Subtitle Size", max = 3, min = 0 ,value=if(is.null(gpOptions$subSize)){iTheme$subSize}else{gpOptions$subSize},step = 0.1, width="95%"),
              shiny::textInput("axisAppend", label="Append Axis Text",value=if(is.null(gpOptions$axisText)){NULL}else{gpOptions$axisText[2]},width="95%")
            ),
            shiny::fillRow(
              shiny::textInput("dataAxisLablesCol",label="Data Label Color", value=if(is.null(gpOptions$plotColors)){NULL} else if(is.null(gpOptions$plotColors$numbers)){NULL} else {paste0(gpOptions$plotColors$numbers,collapse = ", ")}, width="95%"),
              shiny::sliderInput("yAxisLabSize", label="Data Label Size", width="95%", min=0, max=3,step=0.1, value=if(is.null(gpOptions$yAxisLabSize)){iTheme$yAxisLabSize}else{gpOptions$yAxisLabSize}),
              shiny::textInput("axisCol",label="Axis Color",value=if(is.null(gpOptions$plotColors)){NULL} else if(is.null(gpOptions$plotColors$axis)){NULL} else {gpOptions$plotColors$axis}, width="95%")
            )
          ),
          shiny::fillCol(height = "55%",
            shiny::fillRow(
              shiny::textInput("minorTickColor",label="Minor Tick Color",value=if(is.null(gpOptions$plotColors)){NULL} else if(is.null(gpOptions$plotColors$minorTick)){NULL} else {paste0(gpOptions$plotColors$minorTick, collapse = ", ")}, width="95%"),
              shiny::textInput("majorGuideColor",label="Guideline Color",value=if(is.null(gpOptions$plotColors)){NULL} else if(is.null(gpOptions$plotColors$guides)){NULL} else {paste0(gpOptions$plotColors$guides, collapse = ", ")}, width="95%"),
              shiny::textInput("minorGuidesColor",label="Minor Guide Color",value=if(is.null(gpOptions$plotColors)){NULL} else if(is.null(gpOptions$plotColors$minorGuides)){NULL} else {paste0(gpOptions$plotColors$minorGuides, collapse = ", ")}, width="95%")
            ),
            shiny::fillRow(
              shiny::textInput("majorTickColor",label="Axis Tick Color", value=if(is.null(gpOptions$plotColors)){NULL} else if(is.null(gpOptions$plotColors$majorTick)){NULL} else {paste0(gpOptions$plotColors$majorTick, collapse = ", ")}, width="95%"),
              shiny::checkboxInput("extendTick", label="Extend Minor Ticks", value=if(is.null(gpOptions$extendTicks)){TRUE}else{gpOptions$extendTicks},width="95%"),
              shiny::br()
            )
          ),
          shiny::hr(),
          shiny::h4("Background and Legend"),
          shiny::fillCol(height = "55%",
            shiny::fillRow(
              shiny::textInput("canvasColor",label="Plot Area Color",value=if(is.null(gpOptions$plotColors)){NULL} else if(is.null(gpOptions$plotColors$bg)){NULL} else {paste0(gpOptions$plotColors$bg, collapse = ", ")}, width="95%"),
              shiny::textInput("marginCol",label="Background Color",value=if(is.null(gpOptions$plotColors)){NULL} else if(is.null(gpOptions$plotColors$marginBg)){NULL} else {paste0(gpOptions$plotColors$marginBg, collapse = ", ")}, width="95%"),
              shiny::textInput("LegendBG", label="Legend Background Color", value=if(is.null(gpOptions$plotColors)){NULL} else if(is.null(gpOptions$plotColors$legendBG)){NULL} else {paste0(gpOptions$plotColors$legendBG, collapse = ", ")}, width="95%")
            ),
            shiny::fillRow(
              shiny::textInput("legendBorder",label="Legend Border Color", value=if(is.null(gpOptions$plotColors)){NULL} else if(is.null(gpOptions$plotColors$legendBorder)){NULL} else {paste0(gpOptions$plotColors$legendBorder, collapse = ", ")}, width="95%"),
              shiny::textInput("LegendLineCol", label="Legend Box Outline Color",value=if(is.null(gpOptions$plotColors)){NULL} else if(is.null(gpOptions$plotColors$legendLineCol)){NULL} else {paste0(gpOptions$plotColors$legendLineCol, collapse = ", ")}, width="95%"),
              shiny::br()
            )
          )
        )
      ),
      ################ Begin Advanced Tab Section ##################################
      miniUI::miniTabPanel("Advanced", icon = shiny::icon("code"),
        miniUI::miniContentPanel(padding = 5,
          shiny::h4("Misc Options"),
          shiny::fillCol(height="85%",
            shiny::fillRow(
              shiny::checkboxInput("useRgl", label="Use RGL for 3D",value = if(is.null(gpOptions$useRgl)){TRUE}else{gpOptions$useRgl}, width="95%"),
              shiny::numericInput("nlevels",label="# of Contour Levels", value = if(is.null(gpOptions$nLevels)){10}else{gpOptions$nlevels},min = 3,max=100,step = 1, width="95%"),
              shiny::sliderInput("curvePoints",label="Curve/Bootstrap Sampling",min=20,max=1000,step=20,value=if(is.null(gpOptions$curvePoints)){iTheme$curvePoints}else{gpOptions$curvePoints}, width="95%")
            ),
            shiny::fillRow(
              shiny::checkboxInput("verbose", label="Verbose Output", value=if(is.null(gpOptions$verbose)){FALSE}else{gpOptions$verbose}, width="95%"),
              shiny::selectInput("calcType", label="Statistical Testing",
                                 choices=list(None="none",Wilcoxon="wilcox",T.Test="t.test",ANOVA="anova",Tukey_HSD="Tukey"),
                                 selected = if(is.null(gpOptions$calcType)){"none"} else if(any(grepl(gpOptions$calcType,c("none","wilcox","t.test","anova","Tukey"), ignore.case = TRUE))) {grep(gpOptions$calcType,c("none","wilcox","t.test","anova","Tukey"), ignore.case = TRUE, value=TRUE)}else{"none"},
                                 width="95%"),
              shiny::textInput("bandwidth", label="Bandwidth", value=if(is.null(gpOptions$bandwidth)){NULL}else{gpOptions$bandwidth}, width="95%")
            ),
            shiny::fillRow(
              shiny::checkboxInput("strictLimits",label="Strict Data Axis Limits", value=if(is.null(gpOptions$strictLimits)){FALSE}else{gpOptions$strictLimits}, width="95%"),
              shiny::sliderInput("outliers",label="Outlier IQR Threshold", min=0,max=5,step=.1,value=if(is.null(gpOptions$outliers)){1.5}else{gpOptions$outliers}, width="95%"),
              shiny::br()
            )
          ),
          shiny::h4("Actual genePlot Code"),
          shiny::fillCol(height="55%",
            shiny::fillRow(
              shiny::textOutput("codeExample",)
            )
          ),
          shiny::h4("Statistics"),
          shiny::fillCol(height="65%",
            shiny::fillRow(
              shiny::verbatimTextOutput("stats", placeholder = TRUE)
            )
          ),
          shiny::h4("Data Summary"),
          shiny::fillCol(
            shiny::fillRow(
              shiny::tableOutput("descriptive")
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


  #This is the beginning or the server sections where the reactive values and element response code resides
  #RSOveride is a element that stores the npData output from the plot preview and the command line code
  #needed to generate the plot. The plotOptions reactive value is a list of all current settings and is the variable
  #UI elements update based on user input. The plotOptions settings are then passed to genePlot to create the
  #plot preview and the output is stored in RSOveride as previously mentioned.
  #
  #The following sections examines each UI element to update the vcommond list variable. Unlike plot options which
  #stores all current settings, this tries to contain only minimal set of options to generate the current code.
  #This is used to generate the code equivalent example in the advanced tab and is stored in RSOveride to be returned
  #to genePlot for use in non-RStudio environments (graphic output weirdness RStudio prevents direct plotting from the shiny app for some reason.)
  #The final section is all of the UI observation code to respond to the user and update plotOptions. The gpOptions
  #variable is used to remember the initial settings for initial defaults and also options reset.
  #The observation code does try to update all options set to default in the plot type or theme changes.
  #Elements not set to default values are left alone. Note some settings don't exsist for some plot types
  #such as point options for bar plots. In these cases box plot options are used as an alternative default.
  server <- function(input, output, session) {
    aGenes<-shiny::reactiveValues(g=gpOptions$gene)
    cTheme<-"basicTheme"
    PS<-"BP"
    #vcommand<-shiny::reactiveValues(x=dbName)
    RSOveride<-shiny::reactiveValues(rso=TRUE, npData=NULL, options=NULL)
    shiny::isolate(
      plotOptions<-shiny::reactiveValues(gene=gpOptions$gene,
                                group=if(is.null(gpOptions$group)){NULL}else{gpOptions$group},
                                subgroup=if(is.null(gpOptions$subgroup)){NULL}else{gpOptions$subgroup},
                                highlight=if(is.null(gpOptions$highlight)){NULL}else{gpOptions$highlight},
                                stack=if(is.null(gpOptions$stack)){NULL}else{gpOptions$stack},
                                plotType=if(is.null(gpOptions$plotType[1])){"violin"}else{gpOptions$plotType[1]},
                                drawPoints=if(is.null(gpOptions$drawPoints)){TRUE}else{gpOptions$drawPoints},
                                legend=if(is.null(gpOptions$legend)){FALSE}else{gpOptions$legend},
                                logScale=if(is.null(gpOptions$logScale)){FALSE} else {gpOptions$logScale},
                                minorTick=if(is.null(gpOptions$minorTick)){if(gpOptions$logScale==FALSE){iTheme$minorTick}else{iTheme$minorTickLS}}else{gpOptions$minorTick},
                                guides=if(is.null(gpOptions$guides)){iTheme$guides}else{gpOptions$guides},
                                minorGuides=if(is.null(gpOptions$minorGuides)){iTheme$minorGuides}else{gpOptions$minorGuides},
                                pointMethod=if(is.null(gpOptions$pointMethod)){if(IPT %in% c("Bar","2D")){"jitter"}else{iTheme[[paste0("pointMethod",IPT)]]}}else{gpOptions$pointMethod},
                                theme=themeName,
                                pointSize=if(is.null(gpOptions$pointSize)){if(IPT=="Bar"){iTheme$pointSizeBP}else{iTheme[[paste0("pointSize",IPT)]]}}else{gpOptions$pointSize},
                                main=if(is.null(gpOptions$main)) {TRUE} else {gpOptions$main},
                                expLabels=if(is.null(gpOptions$expLabels)){FALSE} else {gpOptions$expLabels},
                                logAdjustment=if(is.null(gpOptions$logAdjustment)){1} else {gpOptions$logAdjustment},
                                rotateLabels=if(is.null(gpOptions$rotateLabels)){FALSE} else {gpOptions$rotateLabels},
                                ylab=if(is.null(gpOptions$ylab)){""} else {gpOptions$ylab},
                                rotateY=if(is.null(gpOptions$rotateY)){FALSE} else {gpOptions$rotateY},
                                drawRug=if(is.null(gpOptions$drawRug)){FALSE} else {gpOptions$drawRug},
                                errorMultiple=if(is.null(gpOptions$errorMultiple)){2} else {gpOptions$errorMultiple},
                                aggFun=if(is.null(gpOptions$aggFun)){"mean"} else {gpOptions$aggFun},
                                errFun=if(is.null(gpOptions$errFun)){"se"} else {gpOptions$errFun},
                                errorBars=if(is.null(gpOptions$errorBars)){TRUE} else {gpOptions$errorBars},
                                trimCurves=if(is.null(gpOptions$trimCurves)){TRUE} else {gpOptions$trimCurves},
                                vioBoxWidth=if(is.null(gpOptions$vioBoxWidth)){iTheme[["vioBoxWidth"]]}else{gpOptions$vioBoxWidth},
                                normalize=if(is.null(gpOptions$normalize)){FALSE} else {gpOptions$normalize},
                                groupByGene=if(is.null(gpOptions$groupByGene)){TRUE} else {gpOptions$groupByGene},
                                legendSize=if(is.null(gpOptions$legendSize)){iTheme[["legendSize"]]}else{gpOptions$legendSize},
                                plotColors=if(is.null(gpOptions$plotColors)){list()}else{gpOptions$plotColors},
                                titleSize=if(is.null(gpOptions$titleSize)){iTheme[["titleSize"]]}else{gpOptions$titleSize},
                                groupLabSize=if(is.null(gpOptions$groupLabSize)){iTheme[["groupLabSize"]]}else{gpOptions$groupLabSize},
                                groupLabelSpacing=if(is.null(gpOptions$groupLabelSpacing)){iTheme[["groupLabelSpacing"]]}else{gpOptions$groupLabelSpacing},
                                fontFamily=if(is.null(gpOptions$fontFamily)){iTheme[["fontFamily"]]}else{gpOptions$fontFamily},
                                errorBarCapWidth=if(is.null(gpOptions$errorBarCapWidth)){iTheme[[paste0("errorBarCapWidth",IPT)]]}else{gpOptions$errorBarCapWidth},
                                errorBarLineType=if(is.null(gpOptions$errorBarLineType)){iTheme[[paste0("errorBarLineType",IPT)]]}else{gpOptions$errorBarLineType},
                                errorCapType=if(is.null(gpOptions$errorCapType)){iTheme[["errorCapType"]]}else{gpOptions$errorCapType},
                                subgroupLabelSpacing=if(is.null(gpOptions$subgroupLabelSpacing)){iTheme[["subgroupLabelSpacing"]]}else{gpOptions$subgroupLabelSpacing},
                                subgroupLabSize=if(is.null(gpOptions$subgroupLabSize)){iTheme[["subgroupLabSize"]]}else{gpOptions$subgroupLabSize},
                                axisText=if(is.null(gpOptions$axisText)){c("","")}else{gpOptions$axisText},
                                extendTicks=if(is.null(gpOptions$extendTicks)){TRUE}else{gpOptions$extendTicks},
                                strictLimits=if(is.null(gpOptions$strictLimits)){FALSE}else{gpOptions$strictLimits},
                                sidePlot=if(is.null(gpOptions$sidePlot)){FALSE}else{gpOptions$sidePlot},
                                verbose=if(is.null(gpOptions$verbose)){FALSE}else{gpOptions$verbose},
                                curvePoints=if(is.null(gpOptions$curvePoints)){iTheme[["curvePoints"]]}else{gpOptions$curvePoints},
                                nlevels=if(is.null(gpOptions$nlevels)){10}else{gpOptions$nlevels},
                                bandwidth=if(is.null(gpOptions$bandwidth)){NULL}else{gpOptions$bandwidth},
                                outliers=if(is.null(gpOptions$bandwidth)){1.5}else{gpOptions$bandwidth},
                                subtitle=if(is.null(gpOptions$subtitle)){NULL}else{gpOptions$subtitle},
                                drawBox=if(is.null(gpOptions$drawBox)){TRUE}else{gpOptions$drawBox},
                                barType=if(is.null(gpOptions$barType)){"bar"}else{gpOptions$barType},
                                calcType=if(is.null(gpOptions$calcType)){"none"}else{gpOptions$calcType},
                                axisLabelSize=if(is.null(gpOptions$axisLabelSize)){iTheme[["axisLabelSize"]]}else{gpOptions$axisLabelSize},
                                yAxisLabSize=if(is.null(gpOptions$yAxisLabSize)){iTheme[["yAxisLabSize"]]}else{gpOptions$yAxisLabSize},
                                subSize=if(is.null(gpOptions$subSize)){iTheme[["subSize"]]}else{gpOptions$subSize},
                                symbol=gpOptions$symbol,
                                useNormCounts=gpOptions$useNormCounts)
      )

      output$plot1 <- shiny::renderPlot({
      RSOveride$npData<-genePlot(data,
                                 gene=plotOptions$gene,
                                 group=plotOptions$group,
                                 subgroup=plotOptions$subgroup,
                                 highlight=plotOptions$highlight,
                                 stack=plotOptions$stack,
                                 plotType=plotOptions$plotType[1],
                                 drawPoints=plotOptions$drawPoints,
                                 legend=plotOptions$legend,
                                 pointMethod=plotOptions$pointMethod,
                                 RSOveride=RSOveride$rso,
                                 theme=get(plotOptions$theme),
                                 normalize=plotOptions$asPercentage,
                                 groupByGene=plotOptions$groupByGene,
                                 main=plotOptions$main,
                                 pointSize=plotOptions$pointSize,
                                 legendSize=plotOptions$legendSize,
                                 logScale=plotOptions$logScale,
                                 logAdjustment=plotOptions$logAdjustment,
                                 ylab=plotOptions$ylab,
                                 expLabels=plotOptions$expLabels,
                                 rotateY=plotOptions$rotateY,
                                 rotateLabels=plotOptions$rotateLabels,
                                 groupNames=plotOptions$groupNames,
                                 subgroupLabels=plotOptions$subgroupLabels,
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
                                 subgroupLabSize=plotOptions$subgroupLabSize,
                                 subgroupLabelSpacing=plotOptions$subgroupLabelSpacing,
                                 axisText=plotOptions$axisText,
                                 extendTicks=plotOptions$extendTick,
                                 verbose=plotOptions$verbose,
                                 curvePoints=plotOptions$curvePoints,
                                 sidePlot=plotOptions$sidePlot,
                                 xLim=plotOptions$xLim,
                                 yLim=plotOptions$yLim,
                                 nlevels=plotOptions$nlevels,
                                 bandwidth=plotOptions$bandwidth,
                                 outliers=plotOptions$outliers,
                                 subtitle=plotOptions$subtitle,
                                 drawBox=plotOptions$drawBox,
                                 barType=plotOptions$barType,
                                 calcType=plotOptions$calcType,
                                 yAxisLabSize=plotOptions$yAxisLabSize,
                                 subSize=plotOptions$subSize,
                                 axisLabelSize=plotOptions$axisLabelSize,
                                 symbol=plotOptions$symbol,
                                 useNormCounts=plotOptions$useNormCounts)

      #######Begin section constructing minimal options needed using vcommand ############
      #This is used to generate the virtual code example in the advanced tab
      #Outside of RStudio this command is actually used to generate the the plot after shiny is finished
      #If using RStudio, due to some unexplained issues, the npData object generated from the preview plot
      #and plotOptions is returned and can be plotted later.
      vcommand<-list(x=dbName)
      if(!is.null(plotOptions$gene[1])) {
        if(plotOptions$gene[1]!="") {
          if(length(plotOptions$gene)==1) {
            vcommand$gene<-paste0("\"",plotOptions$gene,"\"")
          } else {
            vcommand$gene<-paste0("c(",paste0("\"",plotOptions$gene,"\"",collapse = ","),")")
          }
        }
      }
      if(!is.null(plotOptions$group[1])) {
        groupVal<-shiny::renderText( {input$group} )
        if(!is.null(groupVal())) {
          if(groupVal() != "" & groupVal() != FALSE){
            if(groupVal() == input$groupFactorSelector) {
              vcommand$group<-paste0("\"",groupVal(),"\"")
            } else {
              vcommand$group<-groupVal()
            }
          }
        }
      }
      if(!is.null(plotOptions$subgroup[1])) {
        subgroupVal<-shiny::renderText( {input$subgroup} )
        if(!is.null(subgroupVal())) {
          if(subgroupVal() != "" & subgroupVal() != FALSE){
            if(subgroupVal() == input$subgroupFactorSelector) {
              vcommand$subgroup<-paste0("\"",subgroupVal(),"\"")
            } else {
              vcommand$subgroup<-subgroupVal()
            }
          }
        }
      }
      if(!is.null(plotOptions$highlight[1])) {
        highlightVal<-shiny::renderText( {input$highlight} )
        if(!is.null(highlightVal())) {
          if(highlightVal() != "" & highlightVal() != FALSE){
            if(highlightVal() == input$highlightFactorSelector) {
              vcommand$highlight<-paste0("\"",highlightVal(),"\"")
            } else {
              vcommand$highlight<-highlightVal()
            }
          }
        }
      }
      if(!is.null(plotOptions$stack[1])) {
        stackVal<-shiny::renderText( {input$stack} )
        if(!is.null(stackVal())) {
          if(stackVal() != "" & stackVal() != FALSE){
            if(stackVal() == input$stackFactorSelector) {
              vcommand$stack<-paste0("\"",stackVal(),"\"")
            } else {
              vcommand$stack<-stackVal()
            }
          }
        }
      }
      if(!is.null(plotOptions$stack[1])) {
        if(plotOptions$stack[1]!=FALSE) {
          if(is.character(plotOptions$stack[1]) & length(plotOptions$stack)==1) {
            vcommand$stack<-paste0("\"",plotOptions$stack,"\"")
          } else {
            stackVal<-shiny::renderText( {input$stack} )
            vcommand$stack<-stackVal()
          }
        }
      }
      PS<-"BP"
      if(plotOptions$plotType[1] != "box") {
        vcommand$plotType<-paste0("\"",plotOptions$plotType[1],"\"")
        if(plotOptions$plotType[1]=="violin") {
          PS<-"VP"
        } else if(plotOptions$plotType[1]=="dot") {
          PS<-"DP"
        } else if(plotOptions$plotType[1]=="bar") {
          PS<-"Bar"
        } else {
          PS<-"2D"
        }
      }
      if(plotOptions$drawPoints==FALSE){
        vcommand$drawPoints<-FALSE
      }
      if(!is.null(plotOptions$legend) & plotOptions$legend!=FALSE){
        lt<-'"Legend"'
        legendTitleVal<-shiny::renderText(shiny::req( {input$legendTitle} ))
        if(!is.null(legendTitleVal) & legendTitleVal()!=""){
          lt<-paste0("\"",legendTitleVal(),"\"")
        }
        vcommand$legend<-lt
      }
      cTheme<-eval(parse(text=plotOptions$theme))
      if(plotOptions$theme!="basicTheme") {
        vcommand$theme<-plotOptions$theme
      }
      if(plotOptions$plotType[1] %in% c("box","dot","violin") ) {
        if(plotOptions$pointMethod!=cTheme[[paste0("pointMethod",PS)]]) {
          vcommand$pointMethod<-paste0("\"",plotOptions$pointMethod,"\"")
        }

        if(!is.null(plotOptions$pointLaneWidth)) {
          if(cTheme[[paste0("pointLaneWidth",PS)]]!=plotOptions$pointLaneWidth & plotOptions$plotType %in% c("box","dot","violin")) {
            vcommand$pointLaneWidth<-plotOptions$pointLaneWidth
          }
        }
      }
      if(plotOptions$plotType[1] != "bar") {
        if(plotOptions$pointSize!=cTheme[[paste0("pointSize",PS)]]) {
          vcommand$pointSize<-plotOptions$pointSize
        }
      }
      if(!is.null(plotOptions$asPercentage) & !is.null(plotOptions$stack)) {
        if(PS=="Bar" & plotOptions$asPercentage==TRUE & plotOptions$stack!=FALSE) {
          vcommand$normalize<-TRUE
        }
      }
      if(plotOptions$groupByGene==FALSE) {
        vcommand$groupByGene<-FALSE
      }

      if(!is.null(plotOptions$main)) {
        if(plotOptions$main!=TRUE & plotOptions$main != ""){
          if(length(plotOptions$gene)>1) {
            if(plotOptions$main[1] != paste0(c("Gene Expression:",paste0(plotOptions$gene,collapse=", ")),collapse=" ")) {
              vcommand$main<-paste0("\"",plotOptions$main,"\"")
            }
          } else if (length(plotOptions$gene)==1 & plotOptions$main[1] != paste0(plotOptions$gene, " Expression")) {
            vcommand$main<-paste0("\"",plotOptions$main,"\"")
          }
        }
      }

      if(plotOptions$legendSize!=cTheme$legendSize) {
        vcommand$legendSize<-plotOptions$legendSize
      }
      if(sum(plotOptions$logScale)>0) {
        vcommand$logScale<-plotOptions$logScale
      }
      if(plotOptions$logAdjustment!=1 & sum(plotOptions$logScale)>0) {
        vcommand$logAdjustment<-plotOptions$logAdjustment
      }
      if(!is.null(plotOptions$ylab)) {
        if(plotOptions$ylab!=""){
          vcommand$ylab<-paste0("\"",plotOptions$ylab,"\"")
        }
      }
      if(plotOptions$expLabels==TRUE) {
        vcommand$expLabels<-plotOptions$expLabels
      }
      if(plotOptions$rotateY==TRUE) {
        vcommand$rotateY<-plotOptions$rotateY
      }
      if(plotOptions$rotateLabels==TRUE) {
        vcommand$rotateLabels<-plotOptions$rotateLabels
      }
      if(!is.null(plotOptions$groupNames[1])) {
        if(plotOptions$groupNames[1]!="") {
          vcommand$groupNames<-plotOptions$groupNames
        }
      }
      if(!is.null(plotOptions$subgroupLabels[1])){
        if(plotOptions$subgroupLabels[1]!="") {
          vcommand$subgroupLabels<-plotOptions$subgroupLabels
        }
      }
      if((sum(plotOptions$logScale)>0 & cTheme$minorTick!=plotOptions$minorTick) | (sum(plotOptions$logScale)>0 & cTheme$minorTickLS!=plotOptions$minorTick)){
        vcommand$minorTick<-plotOptions$minorTick
      }
      if(cTheme$guides!=plotOptions$guides) {
        vcommand$guides<-plotOptions$guides
      }
      if(!is.null(plotOptions$minorGuides) & !is.null(plotOptions$guides)){
        if(plotOptions$minorGuides!=plotOptions$guides) {
          vcommand$minorGuides<-plotOptions$minorGuides
        }
      }
      if(plotOptions$drawRug==TRUE) {
        vcommand$drawRug<-plotOptions$drawRug
      }
      if(PS %in% c("DP","Bar")){
        if(plotOptions$aggFun!="mean") {
          vcommand$aggFun<-paste0("\"",plotOptions$aggFun,"\"")
        }
        if(plotOptions$errFun!="se") {
          vcommand$errFun<-paste0("\"",plotOptions$errFun,"\"")
        }
        if(plotOptions$errorMultiple != 2) {
          vcommand$errorMultiple<-plotOptions$errorMultiple
        }
        if(!is.null(plotOptions$errorBars)) {
          if(plotOptions$errorBars != TRUE) {
            vcommand$errorBars<-plotOptions$errorBars
          }
        }
      }
      # RSOveride$options<-vcommand
      if(plotOptions$plotType[1]=="surface") {
        if(!is.null(plotOptions$theta)) {
          vcommand$theta<-plotOptions$theta
        }
        if(!is.null(plotOptions$phi)) {
          vcommand$phi<-plotOptions$phi
        }
      }
      if(plotOptions$plotType[1]=="density"){
        if(!is.null(plotOptions$trimCurves)){
          vcommand$trimCurves<-plotOptions$trimCurves
        }
      }
      if(PS=="VP"){
        if(!is.null(plotOptions$trimCurves)){
          if(plotOptions$trimCurves!=TRUE) {
            vcommand$trimViolins<-plotOptions$trimCurves
          }
        }
        if(!is.null(plotOptions$vioBoxWidth)) {
          if(cTheme$vioBoxWidth != plotOptions$vioBoxWidth) {
            vcommand$vioBoxWidth<-plotOptions$vioBoxWidth
          }
        }
      }
      if(!is.null(plotOptions$pointShape)) {
        if(sum(cTheme[[paste0("pointShape",PS)]]==plotOptions$pointShape)!=length(plotOptions$pointShape)) {
          if(length(plotOptions$pointShape)>1) {
            vcommand$pointShape<-paste0("c(",paste0(plotOptions$pointShape,collapse=","),")")
          } else {
            vcommand$pointShape<-plotOptions$pointShape
          }
        }
      }
      if(!is.null(plotOptions$lWidth)) {
        if(cTheme[[paste0("lWidth",PS)]]!=plotOptions$lWidth) {
          vcommand$lWidth<-plotOptions$lWidth
        }
      }
      if(plotOptions$plotType[1] != "surface" & plotOptions$plotType[1] != "density"){
        if(!is.null(plotOptions$width)) {
          if(cTheme[[paste0("width",PS)]]!=plotOptions$width) {
            vcommand$width<-plotOptions$width
          }
        }
      }

      if(!is.null(plotOptions$swarmOverflow)) {
        cpm<-shiny::renderText(plotOptions$pointMethod)
        if(cTheme$swarmOverflow!=plotOptions$swarmOverflow & cpm()=="beeswarm") {
          vcommand$swarmOverflow<-paste0("\"",plotOptions$swarmOverflow,"\"")
        }
      }
      # plotColors=plotOptions$plotColors,
      if(!is.null(plotOptions$errorCapType)) {
        if(cTheme$errorCapType!=plotOptions$errorCapType & plotOptions$plotType %in% c("bar","box","violin")) {
          vcommand$errorCapType<-paste0("\"",plotOptions$errorCapType,"\"")
        }
      }
      if(!is.null(plotOptions$errorBarLineType)) {
        if(!identical(cTheme[[paste0("errorBarLineType",PS)]],plotOptions$errorBarLineType)) {
          vcommand$errorBarLineType<-plotOptions$errorBarLineType
        }
      }
      if(plotOptions$plotType[1] != "density" & plotOptions$plotType[1] !="surface")
      if(!is.null(plotOptions$errorBarCapWidth)) {
        if(cTheme[[paste0("errorBarCapWidth",PS)]]!=plotOptions$errorBarCapWidth) {
          vcommand$errorBarCapWidth<-plotOptions$errorBarCapWidth
        }
      }
      if(!is.null(plotOptions$fontFamily)) {
        if(cTheme$fontFamily!=plotOptions$fontFamily) {
          vcommand$fontFamily<-paste0("\"",plotOptions$fontFamily,"\"")
        }
      }
      if(!is.null(plotOptions$groupLabelSpacing)) {
        if(cTheme$groupLabelSpacing!=plotOptions$groupLabelSpacing) {
          vcommand$groupLabelSpacing<-plotOptions$groupLabelSpacing
        }
      }
      if(!is.null(plotOptions$groupLabSize)) {
        if(cTheme$groupLabSize!=plotOptions$groupLabSize) {
          vcommand$groupLabSize<-plotOptions$groupLabSize
        }
      }
      if(!is.null(plotOptions$titleSize)) {
        if(cTheme$titleSize!=plotOptions$titleSize) {
          vcommand$titleSize<-plotOptions$titleSize
        }
      }
      if(!is.null(plotOptions$axisLabelSize)) {
        if(cTheme$axisLabelSize!=plotOptions$axisLabelSize) {
          vcommand$axisLabelSize<-plotOptions$axisLabelSize
        }
      }
      if(!is.null(plotOptions$yAxisLabSize)) {
        if(cTheme$yAxisLabSize!=plotOptions$yAxisLabSize) {
          vcommand$yAxisLabSize<-plotOptions$yAxisLabSize
        }
      }
      if(!is.null(plotOptions$subgroupLabSize)) {
        if(cTheme$subgroupLabSize!=plotOptions$subgroupLabSize) {
          vcommand$subgroupLabSize<-plotOptions$subgroupLabSize
        }
      }
      if(!is.null(plotOptions$subgroupLabelSpacing)) {
        if(cTheme$subgroupLabelSpacing!=plotOptions$subgroupLabelSpacing) {
          vcommand$subgroupLabelSpacing<-plotOptions$subgroupLabelSpacing
        }
      }
      if(!is.null(plotOptions$subgroupLabelSpacing)) {
        if(cTheme$subgroupLabelSpacing!=plotOptions$subgroupLabelSpacing) {
          vcommand$subgroupLabelSpacing<-plotOptions$subgroupLabelSpacing
        }
      }
      if(!is.null(plotOptions$axisText)) {
        if(is.list(plotOptions$axisText)) {
          if (sum(plotOptions$axisText$x == c("",""))!=2 & sum(plotOptions$axisText$y == c("",""))!=2) {
            vcommand$axisText<-plotOptions$axisText
          }
        } else {
          if(sum(plotOptions$axisText==c("",""))!=2) {
            vcommand$axisText<-plotOptions$axisText
          }
        }
      }
      if(!is.null(plotOptions$extendTick)){
        if(plotOptions$extendTick==FALSE) {
          vcommand$extendTick<-plotOptions$extendTick
        }
      }
      if(!is.null(plotOptions$verbose)){
        if(plotOptions$verbose==TRUE) {
          vcommand$verbose<-plotOptions$verbose
        }
      }
      if(!is.null(plotOptions$curvePoints)) {
        if(cTheme$curvePoints!=plotOptions$curvePoints) {
          vcommand$curvePoints<-plotOptions$curvePoints
        }
      }
      if(!is.null(plotOptions$sidePlot)){
        if(plotOptions$sidePlot==TRUE) {
          vcommand$sidePlot<-plotOptions$sidePlot
        }
      }

      if(!is.null(plotOptions$nlevels)){
        if(plotOptions$nlevels!=10 & plotOptions$plotType[1] =="density" & length(plotOptions$gene)>1) {
          vcommand$nlevels<-plotOptions$nlevels
        }
      }

      if(!is.null(plotOptions$yLim)) {
        vcommand$yLim<-plotOptions$yLim
      }
      if(!is.null(plotOptions$xLim) & plotOptions$plotType[1] =="density" & length(plotOptions$gene)>1) {
        vcommand$xLim<-plotOptions$xLim
      }

      #RSOveride=RSOveride$rso,

      if(!is.null(plotOptions$bandwidth)){
        if(plotOptions$plotType[1] =="violin" | plotOptions$plotType[1] =="density" | plotOptions$plotType[1] =="surface") {
          vcommand$bandwidth<-plotOptions$bandwidth
        }
      }

      if(!is.null(plotOptions$outliers)){
        if(plotOptions$outliers!=1.5 & plotOptions$plotType[1] !="density" & plotOptions$plotType[1] !="surface") {
          vcommand$outliers<-plotOptions$outliers
        }
      }

      if(!is.null(plotOptions$subtitle)) {
        if(plotOptions$subtitle != ""){
          vcommand$subtitle<-paste0("\"",plotOptions$subtitle,"\"")
        }
      }

      if(!is.null(plotOptions$subSize)){
        if(cTheme$subSize != plotOptions$subSize) {
          vcommand$subSize<-plotOptions$subSize
        }
      }

      if(!is.null(plotOptions$drawBox) & (plotOptions$plotType[1] =="box" | plotOptions$plotType[1] =="violin")) {
        if(plotOptions$drawBox == FALSE) {
          vcommand$drawBox<-plotOptions$drawBox
        }
      }

      if(!is.null(plotOptions$barType) & plotOptions$plotType[1] =="dot") {
        if(plotOptions$barType != "bar"){
          vcommand$barType<-paste0("\"",plotOptions$barType,"\"")
        }
      }

      if(!is.null(plotOptions$calcType)) {
        if(plotOptions$calcType!="none") {
          vcommand$calcType<-paste0("\"",plotOptions$calcType,"\"")
        }
      }

      if(!is.null(input$useRgl)) {
        if(input$useRgl==TRUE & plotOptions$plotType[1] =="surface"){
          vcommand$useRgL<-TRUE
        }
      }

      if(length(plotOptions$plotColors)>0) {
        nPCL<-length(plotOptions$plotColors)
        nPlotColors<-vector(mode="character",length=nPCL)
        nPCKeep<-rep(TRUE,nPCL)
        names(nPlotColors)<-names(plotOptions$plotColors)
        for(i in seq(nPCL)) {
          cCol<-plotOptions$plotColors[[i]]
          if(length(cCol)>1) {
            nPlotColors[i]<-paste0("c(",paste0('"',cCol,'"',collapse = ","),")")
          } else {
            if(is.na(cCol[1]) | cCol[1]=="NA") {
              nPCKeep[i]<-FALSE
            } else {
              nPlotColors[i]<-paste0('"',cCol,'"')
            }
          }
        }
        nPlotColors<-nPlotColors[nPCKeep]
        if(sum(nPCKeep)>0){
          nPlotColors<-paste0(names(nPlotColors),"=",nPlotColors,collapse=", ")
          nPlotColors<-paste0("list(",nPlotColors,")")
          vcommand$plotColors<-nPlotColors
        }
      }
      #colorTranslator<-c("points","lines","fill","vioBoxFill","vioBoxLineCol",)
      #names(colorTranslator)<-c("pointColors", "lineColors", "fillColors", "vioBoxFill", "vioBoxLineCol", "titleCol", "labelCol","subgroupLabelCol","yAxisLabelCol","subtitleCol","dataAxisLablesCol","axisCol","minorTickColor","majorGuideColor","minorGuidesColor","majorTickColor","canvasColor","LegendBG","legendBorder","LegendLineCol")
      if(gpOptions$symbol!="GeneSymbol"){
        vcommand$symbol<-gpOptions$symbol
      }
      if(gpOptions$useNormCounts==FALSE) {
        vcommand$useNormCounts<-FALSE
      }
      RSOveride$options<-vcommand

      #This needs to change after the first pass. RSOveride set to TRUE nukes the graphics environment.
      #if RStudio is being used to prevent ghost like over ploting issues caused by it refusing update
      #Within the shiny widget it needs to be run on the first pass only.
      if(RSOveride$rso==TRUE) {
        RSOveride$rso<-FALSE
      }
      output$descriptive<-shiny::renderTable(RSOveride$npData$summary)
      output$Statistics<-shiny::renderText(RSOveride$npData$stats)
      output$codeExample<-shiny::renderText(paste0("genePlot(",paste(names(vcommand),vcommand,sep="=",collapse = ", "),")"),quoted = FALSE)
    })

      #Code for the color plot preview in the format tab
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

    ################ Begin the UI Observation Section ##############################

    #################### Data Minitab ########################
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
      if(!is.null(shiny::req( {input$groupFactorSelector })) & RSOveride$rso==FALSE) {
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

    shiny::observeEvent(input$subgroupFactorSelector, {
      if(!is.null(shiny::req( {input$subgroupFactorSelector })) & RSOveride$rso==FALSE) {
        shiny::updateTextInput(session, "subgroup", value = input$subgroupFactorSelector )
      }
    })

    shiny::observeEvent(input$subgroup, {
      groupVal<-shiny::renderText( {shiny::req(input$group)} )
      if(!is.null(shiny::req( {input$subgroup })) & groupVal()!="FALSE") {
        subgroupVal<-shiny::renderText( {input$subgroup} )
        if(subgroupVal()!="") {
          if(subgroupVal() == shiny::req( {input$subgroupFactorSelector} )){
            plotOptions[["subgroup"]]<-subgroupVal()
          } else {
            plotOptions[["subgroup"]] <- tryCatch(
              { suppressWarnings(eval(parse(text=subgroupVal()))) },
              warning = function(w) {plotOptions[["subgroup"]] <-subgroupVal()},
              error = function(e) {plotOptions[["subgroup"]] <-FALSE}
            )
          }
        } else {
          plotOptions[["subgroup"]]<-FALSE
        }
      }
    })

    shiny::observeEvent(input$highlightFactorSelector, {
      if(!is.null(shiny::req( {input$highlightFactorSelector })) & RSOveride$rso==FALSE) {
        shiny::updateTextInput(session, "highlight", value = input$highlightFactorSelector )
      }
    })

    shiny::observeEvent(input$highlight, {
      if(!is.null(shiny::req( {input$highlight }))) {
        highlightVal<-shiny::renderText( {input$highlight} )
        if(highlightVal()!="") {
          if(highlightVal() == shiny::req( {input$highlightFactorSelector} )){
            plotOptions[["highlight"]]<-highlightVal()
          } else {
            plotOptions[["highlight"]] <- tryCatch(
              { suppressWarnings(eval(parse(text=highlightVal()))) },
              warning = function(w) {plotOptions[["highlight"]] <-highlightVal()},
              error = function(e) {plotOptions[["highlight"]] <-FALSE}
            )
          }
        } else {
          plotOptions[["highlight"]]<-FALSE
        }
      }
    })

    shiny::observeEvent(input$stackFactorSelector, {
      if(!is.null(shiny::req( {input$stackFactorSelector })) & RSOveride$rso==FALSE) {
        shiny::updateTextInput(session, "stack", value = input$stackFactorSelector )
      }
    })

    shiny::observeEvent(input$stack, {
      if(!is.null(shiny::req( {input$stack }))) {
        stackVal<-shiny::renderText( {input$stack} )
        if(stackVal()!="") {
          if(stackVal() == shiny::req( {input$stackFactorSelector} )){
            plotOptions[["stack"]]<-stackVal()
          } else {
            plotOptions[["stack"]] <- tryCatch(
              { suppressWarnings(eval(parse(text=stackVal()))) },
              warning = function(w) {plotOptions[["stack"]] <-stackVal()},
              error = function(e) {plotOptions[["stack"]] <-FALSE}
            )
          }
        } else {
          plotOptions[["stack"]]<-FALSE
        }
      }
    })

    shiny::observeEvent(input$groupByGene, {
      plotOptions$groupByGene<-input$groupByGene
    })

    shiny::observeEvent(input$asPercentage, {
      plotOptions$normalize<-input$asPercentage
    })


    ######Plot Options Minitab#######################################

    shiny::observeEvent(input$plotType, {
      pt<-shiny::renderText({input$plotType})
      RSOveride$rso<-TRUE

      preVal<-shiny::renderText(input$axisPrepend)
      aVal<-shiny::renderText(input$axisAppend)
      axisText<-c("","")
      oPT<-"2D"
      if(plotOptions$plotType=="dot"){
        oPT<-"DP"
      } else if (plotOptions$plotType=="violin") {
        oPT<-"VP"
      } else if (plotOptions$plotType=="bar") {
        oPT<-"Bar"
      } else if (plotOptions$plotType=="box") {
        oPT<-"BP"
      }
      cPT<-"2D"
      if (pt() == "dot") {
        cPT<-"DP"
      } else if (pt() == "violin"){
        cPT<-"VP"
      } else if (pt() == "bar") {
        cPT<-"Bar"
      } else if (pt() == "box") {
        cPT<-"BP"
      }

      myTheme<-get(plotOptions$theme)
      if(pt() != "density" & pt() != "surface") {
        plotOptions$logScale<-plotOptions$logScale[1]
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
        if(oPT=="2D"){oPT<-"BP"} #BP is the most likely setting with these options prior to a 2D plot as it is the defaul
        if(input$lWidth[1] == myTheme[[paste0("lWidth",oPT)]]){
          shiny::updateSliderInput(session,"lWidth",value = myTheme[[paste0("lWidth",cPT)]])
        }
        if(input$errorBarLineType[1] == myTheme[[paste0("errorBarLineType",oPT)]]){
          shiny::updateTextInput(session,"errorBarLineType",value = myTheme[[paste0("errorBarLineType",cPT)]])
        }
        if(input$errorCapWidth[1] == myTheme[[paste0("errorBarCapWidth",oPT)]]){
          shiny::updateSliderInput(session,"errorCapWidth",value = myTheme[[paste0("errorBarCapWidth",cPT)]])
        }
        if(input$width[1] == myTheme[[paste0("width",oPT)]]){
          shiny::updateSliderInput(session,"width",value = myTheme[[paste0("width",cPT)]])
        }
        if(cPT != "Bar" &  oPT != "Bar") {
          if(input$pointSize[1] == myTheme[[paste0("pointSize",oPT)]]){
            shiny::updateSliderInput(session,"pointSize",value = myTheme[[paste0("pointSize",cPT)]])
          }
          if(identical(input$pointShapes,myTheme[[paste0("pointShape",oPT)]])){
            shiny::updateTextInput(session,"pointShapes",value = myTheme[[paste0("pointShape",cPT)]])
          }
          if(input$pointLaneWidth[1] == myTheme[[paste0("pointLaneWidth",oPT)]]){
            shiny::updateSliderInput(session,"pointLaneWidth",value = myTheme[[paste0("pointLaneWidth",cPT)]])
          }
          if(input$pointMethod[1] == myTheme[[paste0("pointMethod",oPT)]]){
            shiny::updateSelectInput(session,"pointMethod",selected  = myTheme[[paste0("pointMethod",cPT)]])
          }
        } else if (oPT != "Bar") {
          cPT<-"BP" #Again, the most likely previously used default.
          if(input$pointSize[1] == myTheme[[paste0("pointSize",cPT)]]){
            shiny::updateSliderInput(session,"pointSize",value = myTheme[[paste0("pointSize",cPT)]])
          }
          if(identical(input$pointShapes,myTheme[[paste0("pointShape",oPT)]])){
            shiny::updateTextInput(session,"pointShapes",value = myTheme[[paste0("pointShape",cPT)]])
          }
          if(input$pointLaneWidth[1] == myTheme[[paste0("pointLaneWidth",oPT)]]){
            shiny::updateSliderInput(session,"pointLaneWidth",value = myTheme[[paste0("pointLaneWidth",cPT)]])
          }
          if(input$pointMethod[1] == myTheme[[paste0("pointMethod",oPT)]]){
            shiny::updateSelectInput(session,"pointMethod",selected  = myTheme[[paste0("pointMethod",cPT)]])
          }
        }

      } else {
        if(pt() != "density" | length(plotOptions$gene)>1){
          plotOptions$logScale<-plotOptions$logScale
        }
        axisText<-list(x=axisText,y=axisText)
        if(is.null(preVal())) {
          axisText$x[1]<-""
        } else {
          axisText$x[1]<-preVal()
        }
        if(is.null(aVal())) {
          axisText$x[2]<-""
        } else {
          axisText$x[2]<-aVal()
        }
        axisText$y<-axisText$x
        if(oPT != "Bar"){
          if(input$pointSize[1] == myTheme[[paste0("pointSize",oPT)]]){
            shiny::updateSliderInput(session,"pointSize",value = myTheme[[paste0("pointSize2D")]])
          }
          if(identical(input$pointShapes, myTheme[[paste0("pointShape",oPT)]])){
            shiny::updateTextInput(session,"pointShapes",value = myTheme[[paste0("pointShape2D")]])
          }
        }
        if(input$lWidth == myTheme[[paste0("lWidth",oPT)]]){
          shiny::updateSliderInput(session,"lWidth",value = myTheme[[paste0("lWidth2D")]])
        }
      }
      plotOptions$axisText<-axisText
      plotOptions$plotType<-pt()
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
      subgroupVal<- shiny::renderText( {shiny::req(input$subgroup)} )
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
        if(highlightVal() != FALSE | (groupVal() !=FALSE & subgroupVal() != FALSE)) {
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
      myTheme<-get(plotOptions$theme)
      logScaleYVal<-shiny::renderText(shiny::req( {input$logScaleY} ))
      plotTypeVal<-shiny::renderText(shiny::req( {input$plotType} ))
      if(plotTypeVal() == "density" | plotTypeVal() == "surface" ) {
        if(shiny::req( {input$logScaleX} ) != FALSE ) {
          if (input$logScaleY == FALSE) {
            if(plotOptions$logScale[1]!=FALSE & myTheme$minorTickLS == input$minorTick) {
              shiny::updateNumericInput(session,"minorTick",value=as.numeric(myTheme$minorTick))
            }
            plotOptions$logScale[1]<-FALSE
          } else {
            if(plotOptions$logScale[1]==FALSE & myTheme$minorTick == input$minorTick) {
              shiny::updateNumericInput(session,"minorTick",value=as.numeric(value=myTheme$minorTickLS))
            }
            plotOptions$logScale[1]<-as.numeric(logScaleYVal())
          }
        } else {
          if (input$logScaleY == FALSE) {
            if(plotOptions$logScale[1]!=FALSE & myTheme$minorTickLS == input$minorTick) {
              shiny::updateNumericInput(session,"minorTick",value=as.numeric(myTheme$minorTick))
            }
            plotOptions$logScale<-FALSE
          } else {
            if(plotOptions$logScale[1]==FALSE & myTheme$minorTick == input$minorTick) {
              shiny::updateNumericInput(session,"minorTick",value=as.numeric(myTheme$minorTickLS))
            }
            plotOptions$logScale<-as.numeric(logScaleYVal())
          }
        }
      } else {
        if (input$logScaleY == FALSE) {
          if(plotOptions$logScale[1]!=FALSE & myTheme$minorTickLS == input$minorTick) {
            shiny::updateNumericInput(session,"minorTick",value=as.numeric(myTheme$minorTick))
          }
          plotOptions$logScale<-FALSE
        } else {
          if(plotOptions$logScale[1]==FALSE & myTheme$minorTick == as.numeric(input$minorTick)) {
            shiny::updateNumericInput(session,"minorTick",value=as.numeric(myTheme$minorTickLS))
          }
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

    shiny::observeEvent(input$subtitle, {
      subTitleVal<-shiny::renderText( {input$subtitle} )
      if(subTitleVal() == "" | is.null(subTitleVal())) {
        plotOptions$subtitle<-NULL
      } else {
        plotOptions$subtitle<-subTitleVal()
      }
    })

    shiny::observeEvent(input$subSize, {
      if(!is.null(input$subSize)) {
        plotOptions$subSize<-input$subSize
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

    shiny::observeEvent(input$subgroupNames, {
      subgroupNamesVal <- shiny::renderText( {input$subgroupNames} )
      if(is.null(subgroupNamesVal()) | subgroupNamesVal() == "" ){
        plotOptions$subgroupLabels<-NULL
      } else {
        plotOptions$subgroupLabels<-trimws(unlist(strsplit(subgroupNamesVal(),",")))
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

    shiny::observeEvent(input$xLim, {
      limVal<-shiny::renderText(shiny::req( input$xLim ))
      xlimVal<-as.numeric(trimws(unlist(strsplit(gsub("[\\(\\)]","",limVal()),split = ","))))
      if(length(xlimVal) < 2 | is.null(limVal()) | limVal()=="" | sum(is.na(xlimVal))>0) {
        plotOptions$xLim<- NULL
      } else {
        plotOptions$xLim<-xlimVal[1:2]
      }
    })

    shiny::observeEvent(input$yLim, {
      limVal<-shiny::renderText(shiny::req( input$yLim ))
      ylimVal<-as.numeric(trimws(unlist(strsplit(gsub("[\\(\\)]","",limVal()),split = ","))))
      if(length(ylimVal) < 2 | is.null(limVal()) | limVal()=="" | sum(is.na(ylimVal))>0) {
        plotOptions$yLim<- NULL
      } else {
        plotOptions$yLim<-ylimVal[1:2]
      }
    })

    shiny::observeEvent(input$sidePlot, {
      plotOptions$sidePlot<-input$sidePlot
    })

    shiny::observeEvent(input$barType, {
      cbt<-shiny::renderText(input$barType)
      plotOptions$barType<-cbt()
    })

    shiny::observeEvent(input$drawBox, {
      plotOptions$drawBox<-input$drawBox
    })

    ################# Format Minitab #################################
    shiny::observeEvent(input$theme, {
      cTheme<-shiny::renderText(shiny::req({input$theme}))
      oTheme<-get(plotOptions$theme)
      plotOptions$theme<-cTheme()
      myTheme<-get(cTheme())
      cPT<-"2D"
      if(plotOptions$plotType=="dot"){
        cPT<-"DP"
      } else if (plotOptions$plotType=="violin") {
        cPT<-"VP"
      } else if (plotOptions$plotType=="bar") {
        cPT<-"Bar"
      } else if (plotOptions$plotType=="box") {
        cPT<-"BP"
      }
      #If current settings match the old theme values, they are updated to the settings of the new theme
      #updating sliderInputs
      for(i in c("legendSize","width","lWidth","pointLaneWidth","errorCapWidth","pointSize","vioBoxWidth","titleSize","axisLabelSize","curvePoints","groupLabelSpacing","groupLabSize","subgroupLabSize","subgroupLabelSpacing","yAxisLabSize","subSize")) {
        #Some options have plotType specific settings that need to be handled usng th cPT variable.
        if(i %in% c("width","lWidth","pointLaneWidth","errorCapWidth","pointSize")) {
          if(i == "pointSize" | i == "pointLaneWidth") {
            if(cPT == "Bar") { #As bar plots don't have point overlay options, updating against the default boxplot instead
              if (identical(oTheme[[paste0(i,"BP")]],plotOptions[[i]])) {
                shiny::updateSliderInput(session, i, value=myTheme[[paste0(i,"BP")]])
              }
            } else {
              if (identical(oTheme[[paste0(i,cPT)]], plotOptions[[i]])) {
                shiny::updateSliderInput(session, i, value=myTheme[[paste0(i,cPT)]])
              }
            }
          } else {
            if (identical(oTheme[[paste0(i,cPT)]], plotOptions[[i]])) {
              shiny::updateSliderInput(session, i, value=myTheme[[paste0(i,cPT)]])
            }
          }
        } else if(identical(oTheme[[i]],plotOptions[[i]])) {
          shiny::updateSliderInput(session,i,value=myTheme[[i]])
        }
      }
      #updating selectInputs
      for(i in c("errorCapType","fontFamily","swarmOverflow")) {
        if(identical(oTheme[[i]],plotOptions[[i]])) {
         shiny:: updateSelectInput(session,i,selected = myTheme[[i]])
        }
      }
      if(cPT == "Bar") {
        if(oTheme[["pointMethodBP"]] == plotOptions[["pointMethod"]]) {
          shiny::updateSelectInput(session, "pointMethod", selected = myTheme[["pointMethodBP"]])
        }
        if(identical(oTheme[["pointShapesBP"]], plotOptions$pointShapes)) {
          shiny::updateTextInput(session, "pointShapes", value = myTheme[["pointShapesBP"]] )
        }
      } else {
        if(cPT=="2D") {
          if(oTheme[["pointMethodBP"]] == plotOptions$pointMethod) {
            shiny::updateSelectInput(session, "pointMethod", selected = myTheme[["pointMethodBP"]])
          }
        } else {
          if(oTheme[[paste0("pointMethod",cPT)]] == plotOptions$pointMethod) {
            shiny::updateSelectInput(session, "pointMethod", selected = myTheme[[paste0("pointMethod",cPT)]])
          }
        }
        if(identical(oTheme[[paste0("pointShapes",cPT)]], plotOptions$pointShapes)) {
          shiny::updateTextInput(session, "pointShapes", value = myTheme[[paste0("pointShapes",cPT)]])
        }
      }
      if (identical(oTheme[[paste0("errorBarLineType",cPT)]], plotOptions$errorBarLineType)) {
        shiny::updateTextInput(session, "errorBarLineType", value = myTheme[[paste0("errorBarLineType",cPT)]])
      }
      cLs<-shiny::renderText(input$logScaleY)
      if(cLs() == "None") {
        if(oTheme$minorTick == plotOptions$minorTick) {
          if(myTheme$minorTick==FALSE) {
            shiny::updateNumericInput(session,"minorTick",value=0)
          } else {
            shiny::updateNumericInput(session,"minorTick",value=myTheme$minorTick)
          }
        }
      } else {
        if(oTheme$minorTickLS == plotOptions$minorTick) {
          if(myTheme$minorTickLS ==  FALSE) {
            shiny::updateNumericInput(session,"minorTick",value=0)
          } else {
            shiny::updateNumericInput(session,"minorTick",value=myTheme$minorTickLS)
          }
        }
      }
      shiny::updateCheckboxInput(session,"guides",value=myTheme$guides)
      shiny::updateCheckboxInput(session,"showMinorGuide",value=myTheme$guides)
   })

    shiny::observeEvent(input$colorBrewer, {
      cbVal<-shiny::renderText(shiny::req({input$colorBrewer}))
      maxC<-as.numeric(gsub("(.*) \\((\\d+)\\)","\\2",cbVal()))
      pallete<-gsub("(.*) \\((\\d+)\\)","\\1",cbVal())
      nCT<-shiny::renderText(input$nColors)
      nC<-maxC
      if(!is.null(pallete) & pallete !="" & is.numeric(maxC)) {
        if(as.numeric(nCT()) <= maxC) {
          nC<-as.numeric(nCT())
        }
        cbVal<-brewer.pal(name=pallete,n=nC)
        cbVal<-map_chr(cbVal,setAlpha, alpha=input$alphaSlider)
        shiny::updateTextInput(session,inputId="selectedColors",value=paste0(cbVal,collapse=","))
        nVect<-rev(seq(maxC))
        names(nVect)<-nVect
        shiny::updateSelectInput(session, inputId="nColors", choices=as.list(nVect), selected = nC)
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

    shiny::observeEvent(input$lWidth, {
      plotOptions$lWidth<-shiny::req(input$lWidth)
    })

    shiny::observeEvent(input$width, {
      plotOptions$width<-shiny::req(input$width)
    })

    shiny::observeEvent(input$pointLaneWidth, {
      plotOptions$pointLaneWidth<-input$pointLaneWidth
    })

    shiny::observeEvent(input$errorBarLineType, {
      cVal<-shiny::renderText(input$errorBarLineType)
      if(is.null(cVal()) | cVal()=="") {
        plotOptions$errorBarLineType<-NULL
      } else {
        cVal<-as.numeric(trimws(unlist(strsplit(cVal(), split=","))))
        plotOptions$errorBarLineType<-cVal
      }
    })

    shiny::observeEvent(input$errorCapType, {
      cVal<-shiny::renderText(input$errorCapType)
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

    shiny::observeEvent(input$groupLabSize, {
      plotOptions$groupLabSize<-input$groupLabSize
    })

    shiny::observeEvent(input$groupLabelSpacing, {
      plotOptions$groupLabelSpacing<-input$groupLabelSpacing
    })

    shiny::observeEvent(input$subgroupLabelCol, {
      cCol<-shiny::renderText(input$subgroupLabelCol)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$subgroupLabels<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$subgroupLabels<- cVal
      }
    })

    shiny::observeEvent(input$axisLabelSize, {
      plotOptions$axisLabelSize<-input$axisLabelSize
    })

    shiny::observeEvent(input$subgroupLabSize, {
      plotOptions$subgroupLabSize<-input$subgroupLabSize
    })

    shiny::observeEvent(input$subgroupLabelSpacing, {
      plotOptions$subgroupLabelSpacing<-input$subgroupLabelSpacing
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

    shiny::observeEvent(input$yAxisLabSize, {
      if(!is.null(input$yAxisLabSize)) {
        plotOptions$yAxisLabSize<- input$yAxisLabSize
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
      ptVal<-shiny::renderText(input$plotType)
      if(ptVal() == "density" | ptVal()=="surface") {
        axisText<-list(x=axisText,y=axisText)
        if(is.null(preVal())) {
          axisText$x[1]<-""
        } else {
          axisText$x[1]<-preVal()
        }
        if(is.null(aVal())) {
          axisText$x[2]<-""
        } else {
          axisText$x[2]<-aVal()
        }
        axisText$y<-axisText$x
      } else {
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
      }
      plotOptions$axisText<-axisText
    })

    shiny::observeEvent(input$axisAppend, {
      preVal<-shiny::renderText(input$axisPrepend)
      aVal<-shiny::renderText(input$axisAppend)
      axisText<-c("","")
      ptVal<-shiny::renderText(input$plotType)
      if(ptVal() == "density" | ptVal()=="surface") {
        axisText<-list(x=axisText,y=axisText)
        if(is.null(preVal())) {
          axisText$x[1]<-""
        } else {
          axisText$x[1]<-preVal()
        }
        if(is.null(aVal())) {
          axisText$x[2]<-""
        } else {
          axisText$x[2]<-aVal()
        }
        axisText$y<-axisText$x
      } else {
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
        plotOptions$plotColors$bg<- cVal[1]
      }
    })

    shiny::observeEvent(input$marginCol, {
      cCol<-shiny::renderText(input$marginCol)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$marginBg<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$marginBg<- cVal[1]
      }
    })

    shiny::observeEvent(input$subtitleCol, {
      cCol<-shiny::renderText(input$subtitleCol)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$subtext<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$subtext<- cVal[1]
      }
    })

    shiny::observeEvent(input$yAxisLabelCol, {
      cCol<-shiny::renderText(input$yAxisLabelCol)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$axisLabels<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$axisLabels<- cVal[1]
      }
    })

    shiny::observeEvent(input$legendBorder, {
      cCol<-shiny::renderText(input$legendBorder)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$LegendBorder<-NULL
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$legendBorder<- cVal
      }
    })

    shiny::observeEvent(input$LegendLineCol, {
      cCol<-shiny::renderText(input$LegendLineCol)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$legendLineCol<-NA
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$legendLineCol<- cVal
      }
    })

    shiny::observeEvent(input$LegendBG, {
      cCol<-shiny::renderText(input$LegendBG)
      if(is.null(cCol()) | cCol()=="") {
        plotOptions$plotColors$legendBG<-NA
      } else {
        cVal<-trimws(unlist(strsplit(cCol(),",")))
        plotOptions$plotColors$legendBG<- cVal
      }
    })


    ########### Advanced Tab ######################3
    shiny::observeEvent(input$useRgl, {
      plotOptions$useRgl<-input$useRgl
    })

    shiny::observeEvent(input$nlevels, {
      plotOptions$nlevels<-input$nlevels
    })

    shiny::observeEvent(input$curvePoints, {
      plotOptions$curvePoints<-input$curvePoints
    })

    shiny::observeEvent(input$verbose, {
      plotOptions$verbose<-input$verbose
    })

    shiny::observeEvent(input$bandwidth, {
      cbw<-shiny::renderText(input$bandwidth)
      if(!is.null(cbw()) &  cbw()!="") {
        plotOptions$bandwidth<-as.numeric(cbw())
      } else {
        plotOptions$bandwidth<-NULL
      }
    })

    shiny::observeEvent(input$strictLimits, {
      plotOptions$strictLimits<-input$strictLimits
    })

    shiny::observeEvent(input$outliers, {
      plotOptions$outliers<-input$outliers
    })

    shiny::observeEvent(input$calcType, {
      cStat<- shiny::renderText(input$calcType)
      if(is.null(cStat())) {
        plotOptions$calcType<-"none"
      } else {
        plotOptions$calcType<-cStat()
      }
    })

    #Lower Button Panel
    shiny::observeEvent(input$done, {
      vcomList<-shiny::reactiveValuesToList(RSOveride)
      vcomList$options<-vcomList$options[names(vcomList$options)!="RSOveride"]
      #vcomList<-vcomList$options
      shiny::stopApp(returnValue = vcomList)
    })
    shiny::observeEvent(input$reset, {
      shiny::updateSelectizeInput(session,"activeGenes",choices=gpOptions$gene)
      shiny::updateCheckboxInput(session,"groupByGene", value=if(is.null(gpOptions$groupByGene)){TRUE}else{gpOptions$groupByGene})
      shiny::updateTextInput(session,"group",value = gpOptions$group)
      shiny::updateTextInput(session,"subgroup",value = gpOptions$subgroup)
      shiny::updateTextInput(session,"highlight",value = gpOptions$highlight)
      shiny::updateTextInput(session,"stack",value = gpOptions$stack)
      shiny::updateTextInput(session,"group",value = gpOptions$group)
      shiny::updateCheckboxInput(session, "asPercentage", value= if(is.null(gpOptions$asPercentage)){FALSE}else{gpOptions$asPercentage})
      shiny::updateSelectInput(session,"plotType",selected = if(is.null(gpOptions$plotType)){"box"}else if(any(grepl(gpOptions$plotType[1], c("box","dot","violin","bar","density","surface"), ignore.case = TRUE))) {c("box","dot","violin","bar","density","surface")[grep(gpOptions$plotType[1],ignore.case = TRUE, x = c("box","dot","violin","bar","density","surface"))[1]]}else {"box"})
      shiny::updateTextInput(session,"main",value = if(is.null(gpOptions$main)){NULL}else if (gpOptions$main==paste0("Gene Expression: ",paste0(gpOptions$gene,collapse=", "))){NULL}else{gpOptions$main})
      shiny::updateCheckboxInput(session, "plotPoints",value = if(is.null(gpOptions$plotPoints)){TRUE}else{gpOptions$plotPoints})
      shiny::updateSelectInput(session, "pointMethod",selected = if(is.null(gpOptions$pointMethod)){if(IPT %in% c("2D","Bar")){iTheme$pointMethodBP} else {iTheme[[paste0("pointMethod",IPT)]]}} else if(any(grepl(gpOptions$pointMethod ,c("jitter","beeswarm","distribution","linear"), ignore.case = TRUE))){c("jitter","beeswarm","distribution","linear")[grep(gpOptions$pointMethod ,c("jitter","beeswarm","distribution","linear"), ignore.case = TRUE)[1]]}else{"jitter"})
      shiny::updateSliderInput(session, "pointSize", value=if(is.null(gpOptions$pointSize)){if(IPT=="Bar"){iTheme$pointSizeBP}else{iTheme[[paste0("pointSize",IPT)]]}}else{gpOptions$pointSize})
      shiny::updateCheckboxInput(session, "legend",value = if(is.null(gpOptions$legend)){FALSE}else{if(gpOptions$legend==FALSE){FALSE}else{TRUE}})
      shiny::updateTextInput(session, "legendTitle", value = if(is.null(gpOptions$legend)){"Legend"}else{if(gpOptions$legend==FALSE){"Legend"}else{gpOptions$legend}})
      shiny::updateSliderInput(session, "legendSize", value = if(is.null(gpOptions$legendSize)){iTheme$legendSize}else{gpOptions$legendSize})
      shiny::updateSelectInput(session, "logScaleY", selected = if(is.null(gpOptions$logScale)){FALSE} else {if(gpOptions$logScale==2){"2"} else if(gpOptions$logScale==10){"10"} else if(gpOptions$logScale=="e"){"e"}else{FALSE}})
      #shiny::updateSelectInput(session, "logScaleX", )
      shiny::updateSliderInput(session, "logAdjustment",value=if(is.null(gpOptions$logAdjustment)){1}else{gpOptions$logAdjustment})
      shiny::updateTextInput(session, "ylab", value = if(is.null(gpOptions$ylab)){""}else if (gpOptions$ylab==""){""} else {gpOptions$ylab})
      shiny::updateTextInput(session, "subtitle",  value = if(is.null(gpOptions$subtitle)){""}else if (gpOptions$subtitle==""){""} else {gpOptions$subtitle})
      shiny::updateCheckboxInput(session, "rotateY", value = if(is.null(gpOptions$rotateY)){FALSE} else {if(gpOptions$rotateY==TRUE){TRUE}else{FALSE}})
      shiny::updateTextInput(session, "groupNames", value=if(is.null(gpOptions$groupNames)){""}else if (gpOptions$groupNames==""){} else {paste0(gpOptions$groupNames, collapse = ",")})
      shiny::updateTextInput(session, "subgroupNames", value=if(is.null(gpOptions$subgroupNames)){""}else if (gpOptions$subgroupNames==""){""} else {paste0(gpOptions$subgroupNames, collapse = ", ")})
      shiny::updateCheckboxInput(session, "rotateLabels", value=if(is.null(gpOptions$rotateLabels)){FALSE} else {if(gpOptions$rotateLabels==TRUE){TRUE}else{FALSE}})
      shiny::updateNumericInput(session, "minorTick", value=if(!is.null(gpOptions$minorTick)){as.numeric(gpOptions$minorTick)} else if(gpOptions$logScale==FALSE) {as.numeric(iTheme$minorTick)}else{as.numeric(iTheme$minorTickLS)})
      shiny::updateCheckboxInput(session, "showMinorGuide",value = if(is.null(gpOptions$minorGuides)){if(is.null(gpOptions$guides)){iTheme$guides}else{gpOptions$guides}}else{gpOptions$minorGuides})
      shiny::updateCheckboxInput(session, "guides", value = if(is.null(gpOptions$guides)){iTheme$guides}else{gpOptions$guides})
      shiny::updateSelectInput(session, "aggFun", selected = if(is.null(gpOptions$aggFun)){"mean"} else if (grepl(gpOptions$aggFun, ignore.case = TRUE,"median")){"median"}else{"mean"})
      shiny::updateSelectInput(session, "errFun", selected = if(is.null(gpOptions$errFun)){"se"}else if (any(grepl(gpOptions$errFun,ignore.case = TRUE, c("sd","se","range","t95ci","boot95ci")))){grep(gpOptions$errFun,ignore.case = TRUE, c("sd","se","range","t95ci","boot95ci"),value = TRUE)[1]}else{"se"})
      shiny::updateNumericInput(session, "errorMultiple", value=if(is.null(gpOptions$errorMultiple)){2}else{gpOptions$errorMultiple})
      shiny::updateCheckboxInput(session, "errorBars", value = if(is.null(gpOptions$errorBars)){TRUE}else{gpOptions$errorBars})
      shiny::updateSliderInput(session, "phi", value=if(is.null(gpOptions$phi)){30}else{gpOptions$phi})
      shiny::updateSliderInput(session, "theta", value=if(is.null(gpOptions$theta)){30}else{gpOptions$theta})
      shiny::updateCheckboxInput(session, "rug", value = if(is.null(gpOptions$drawRug)){FALSE}else{gpOptions$drawRug})
      shiny::updateCheckboxInput(session, "trimCurves",  value=if(is.null(gpOptions$trimCurves) & is.null(gpOptions$trimViolins)){TRUE}else if(is.null(gpOptions$trimCurves)){gpOptions$trimViolins}else{gpOptions$trimCurves})
      shiny::updateSliderInput(session, "vioBoxWidth", value=if(is.null(gpOptions$vioBoxWidth)){iTheme$vioBoxWidth}else{gpOptions$vioBoxWidth})
      shiny::updateTextInput(session, "yLim", value = if(is.null(gpOptions$yLim)){""}else{paste0(gpOptions$yLim,collapse = ", ")})
      shiny::updateTextInput(session, "xLim", value = if(is.null(gpOptions$xLim)){""}else{paste0(gpOptions$xLim,collapse = ", ")})
      shiny::updateCheckboxInput(session, "sidePlot", value=if(is.null(gpOptions$sidePlot)){FALSE}else{gpOptions$sidePlot})
      shiny::updateCheckboxInput(session, "expLabels", value = if(is.null(gpOptions$expLabels)){FALSE}else{gpOptions$expLabels})
      shiny::updateSelectInput(session, "barType", selected = if(is.null(gpOptions$barType)){"bar"}else if(grepl(gpOptions$barType, ignore.case = TRUE,"none")){3} else if(grepl(gpOptions$barType, ignore.case = TRUE,"dot")){2}else{1})
      shiny::updateCheckboxInput(session, "drawBox", value= if(is.null(gpOptions$drawBox)){TRUE}else{gpOptions$drawBox})
      shiny::updateSelectInput(session, "theme", selected = if(themeName %in% themes){themeName} else {"basicTheme"})
      shiny::updateTextInput(session, "pointColors", value=if(is.null(gpOptions$plotColors)){""}else if(all(gpOptions$plotColors$points==iTheme$plotColors$points)){""}else{paste0(gpOptions$plotColors$points, collapse = ", ")})
      shiny::updateTextInput(session, "pointShapes", value=if(is.null(gpOptions$pointShape)){""} else if(IPT=="Bar"){""} else if(all(gpOptions$pointShape==iTheme[[paste0("pointShape",IPT)]])){}else{paste0(gpOptions$pointShape, collapse = ", ")})
      shiny::updateTextInput(session, "lineColors", value=if(is.null(gpOptions$plotColors)){""}else if(all(gpOptions$plotColors$lines==iTheme$plotColors$lines)){""}else{paste0(gpOptions$plotColors$lines, collapse = ", ")})
      shiny::updateTextInput(session, "fillColors", value=if(is.null(gpOptions$plotColors)){""}else if(all(gpOptions$plotColors$fill==iTheme$plotColors$fill)){""}else{paste0(gpOptions$plotColors$fill, collapse = ", ")})
      shiny::updateSliderInput(session, "lWidth", value=if(is.null(gpOptions$lWidth)){iTheme[[paste0("lWidth",IPT)]]}else{gpOptions$lWidth})
      shiny::updateSliderInput(session, "width", value = if(is.null(gpOptions$width)){if(IPT == "2D"){iTheme$widthBP}else{iTheme[[paste0("width",IPT)]]}}else{gpOptions$width})
      shiny::updateSliderInput(session, "pointLaneWidth", value = if(is.null(gpOptions$pointLaneWidth)){if(IPT %in% c("Bar", "2D")){NULL}else{iTheme[[paste0("pointLaneWidth",IPT)]]}}else{gpOptions$pointLaneWidth})
      shiny::updateTextInput(session, "errorBarLineType", value = if(is.null(gpOptions$errorBarLineType)){iTheme[[paste0("errorBarLineType",IPT)]]}else{gpOptions$errorBarLineType})
      shiny::updateSelectInput(session, "errorCapType", selected = if(is.null(gpOptions$errorCapType)){iTheme$errorCapType}else if (any(grepl(gpOptions$errorCapType, c("none","bar","ball"), ignore.case = TRUE))) {grep(gpOptions$errorCapType, c("none", "ball","bar"), ignore.case = TRUE, value=TRUE)[1]}else {iTheme$errorCapType})
      shiny::updateSliderInput(session, "errorCapWidth", value=if(is.null(gpOptions$errorBarCapWidth)){if(IPT=="2D"){iTheme$errorBarCapWidthBP}else{iTheme[[paste0("errorBarCapWidth",IPT)]]}}else{gpOptions$errorBarCapWidth})
      shiny::updateTextInput(session, "vioBoxFill", value=if(is.null(gpOptions$plotColors)){""}else if(all(gpOptions$plotColors$vioBoxFill==iTheme$plotColors$vioBoxFill)){""}else{paste0(gpOptions$plotColors$vioBoxFill, collapse = ", ")})
      shiny::updateTextInput(session, "vioBoxLineCol", value=if(is.null(gpOptions$plotColors)){""}else if(all(gpOptions$plotColors$vioBoxLineCol==iTheme$plotColors$vioBoxLineCol)){""}else{paste0(gpOptions$plotColors$vioBoxLineCol, collapse = ", ")})
      shiny::updateSelectInput(session, "swarmOverflow", selected=if(is.null(gpOptions$swarmOverflow)){iTheme$swarmOverflow} else if (grepl(gpOptions$swarmOverflow, c("random","gutter","wrap","omit","none"), ignore.case = TRUE)){grep(gpOptions$swarmOverflow, c("random","gutter","wrap","omit","none"),ignore.case = TRUE)[1]}else{iTheme$swarmOverflow})
      shiny::updateSelectInput(session, "fontFamily", selected=if(is.null(gpOptions$fontFamily)){iTheme$fontFamily} else if (any(grepl(gpOptions$fontFamily, c("sans","serif","mono"), ignore.case = TRUE))){grep(gpOptions$fontFamily, c("sans","serif","mono"), ignore.case = TRUE,value=TRUE)}else{iTheme$fontFamily})
      shiny::updateTextInput(session, "titleCol", value=if(is.null(gpOptions$plotColors)){""}else if(is.null(gpOptions$plotColors$title)){""}else if(gpOptions$plotColors$title[1]==iTheme$plotColors$title[1]){""}else{gpOptions$plotColors$title[1]})
      shiny::updateSliderInput(session, "titleSize", value=if(is.null(gpOptions$titleSize)){iTheme$titleSize}else{gpOptions$titleSize})
      shiny::updateTextInput(session, "labelCol", value=if(is.null(gpOptions$plotColors)){""}else if(is.null(gpOptions$plotColors$title)){""}else if(gpOptions$plotColors$title[1]==iTheme$plotColors$title[1]){""}else{gpOptions$plotColors$title[1]})
      shiny::updateSliderInput(session, "groupLabSize", value=if(is.null(gpOptions$groupLabSize)){iTheme$groupLabSize}else{gpOptions$groupLabSize})
      shiny::updateSliderInput(session, "groupLabelSpacing", value=if(is.null(gpOptions$groupLabelSpacing)){iTheme$groupLabelSpacing}else{gpOptions$groupLabelSpacing})
      shiny::updateTextInput(session, "subgroupLabelCol", value=if(is.null(gpOptions$plotColors)){""} else if(is.null(gpOptions$plotColors$subgroupLabels)){""} else {paste0(gpOptions$plotColors$subgroupLabels, collapse = ", ")})
      shiny::updateSliderInput(session, "subgroupLabSize", value=if(is.null(gpOptions$subgroupLabSize)){iTheme$subgroupLabSize}else{gpOptions$subgroupLabSize})
      shiny::updateSliderInput(session, "subgroupLabelSpacing", value=if(is.null(gpOptions$subgroupLabelSpacing)){iTheme$subgroupLabelSpacing}else{gpOptions$subgroupLabelSpacing})
      shiny::updateTextInput(session, "yAxisLabelCol", value=if(is.null(gpOptions$plotColors)){""} else if(is.null(gpOptions$plotColors$axisLabels)){""} else {gpOptions$plotColors$axisLabels[1]})
      shiny::updateSliderInput(session, "axisLabelSize", value=if(is.null(gpOptions$axisLabelSize)){iTheme$axisLabelSize}else{gpOptions$axisLabelSize})
      shiny::updateTextInput(session, "axisPrepend", value=if(is.null(gpOptions$axisText)){""}else{gpOptions$axisText[1]})
      shiny::updateTextInput(session, "subtitleCol", value=if(is.null(gpOptions$plotColors)){""} else if(is.null(gpOptions$plotColors$subtext)){""} else {gpOptions$plotColors$subtext[1]})
      shiny::updateSliderInput(session, "subSize", value=if(is.null(gpOptions$subSize)){iTheme$subSize}else{gpOptions$subSize})
      shiny::updateTextInput(session, "axisAppend", value=if(is.null(gpOptions$axisText)){""}else{gpOptions$axisText[2]})
      shiny::updateTextInput(session, "dataAxisLablesCol", value=if(is.null(gpOptions$plotColors)){""} else if(is.null(gpOptions$plotColors$numbers)){""} else {gpOptions$plotColors$numbers})
      shiny::updateSliderInput(session, "yAxisLabSize", value=if(is.null(gpOptions$yAxisLabSize)){iTheme$yAxisLabSize}else{gpOptions$yAxisLabSize})
      shiny::updateTextInput(session, "axisCol", value=if(is.null(gpOptions$plotColors)){""} else if(is.null(gpOptions$plotColors$axis)){""} else {gpOptions$plotColors$axis})
      shiny::updateTextInput(session, "minorTickColor", value=if(is.null(gpOptions$plotColors)){""} else if(is.null(gpOptions$plotColors$minorTick)){""} else {paste0(gpOptions$plotColors$minorTick, collapse = ", ")})
      shiny::updateTextInput(session, "majorGuideColor", value=if(is.null(gpOptions$plotColors)){""} else if(is.null(gpOptions$plotColors$guides)){""} else {paste0(gpOptions$plotColors$guides, collapse = ", ")})
      shiny::updateTextInput(session, "minorGuidesColor", value=if(is.null(gpOptions$plotColors)){""} else if(is.null(gpOptions$plotColors$minorGuides)){""} else {paste0(gpOptions$plotColors$minorGuides, collapse = ", ")})
      shiny::updateTextInput(session, "majorTickColor", value=if(is.null(gpOptions$plotColors)){""} else if(is.null(gpOptions$plotColors$majorTick)){""} else {paste0(gpOptions$plotColors$majorTick, collapse = ", ")})
      shiny::updateCheckboxInput(session, "extendTick", value=if(is.null(gpOptions$extendTicks)){TRUE}else{gpOptions$extendTicks})
      shiny::updateTextInput(session, "canvasColor", value=if(is.null(gpOptions$plotColors)){""} else if(is.null(gpOptions$plotColors$bg)){""} else {paste0(gpOptions$plotColors$bg, collapse = ", ")})
      shiny::updateTextInput(session, "marginCol", value=if(is.null(gpOptions$plotColors)){""} else if(is.null(gpOptions$plotColors$marginBg)){""} else {paste0(gpOptions$plotColors$marginBg, collapse = ", ")})
      shiny::updateTextInput(session, "LegendBG", value=if(is.null(gpOptions$plotColors)){""} else if(is.null(gpOptions$plotColors$legendBG)){""} else {paste0(gpOptions$plotColors$legendBG, collapse = ", ")})
      shiny::updateTextInput(session, "legendBorder", value=if(is.null(gpOptions$plotColors)){""} else if(is.null(gpOptions$plotColors$legendBorder)){""} else {paste0(gpOptions$plotColors$legendBorder, collapse = ", ")})
      shiny::updateTextInput(session, "LegendLineCol", value=if(is.null(gpOptions$plotColors)){""} else if(is.null(gpOptions$plotColors$legendLineCol)){""} else {paste0(gpOptions$plotColors$legendLineCol, collapse = ", ")})
      shiny::updateCheckboxInput(session, "useRgl", value = if(is.null(gpOptions$useRgl)){TRUE}else{gpOptions$useRgl})
      shiny::updateNumericInput(session, "nlevels",value = if(is.null(gpOptions$nLevels)){10}else{gpOptions$nlevels})
      shiny::updateSliderInput(session, "curvePoints", value=if(is.null(gpOptions$curvePoints)){iTheme$curvePoints}else{gpOptions$curvePoints})
      shiny::updateCheckboxInput(session, "verbose", value=if(is.null(gpOptions$verbose)){FALSE}else{gpOptions$verbose})
      shiny::updateSelectInput(session, "calcType", selected = if(is.null(gpOptions$calcType)){"none"} else if(any(grepl(gpOptions$calcType,c("none","wilcox","t.test","anova","Tukey"), ignore.case = TRUE))) {grep(gpOptions$calcType,c("none","wilcox","t.test","anova","Tukey"), ignore.case = TRUE, value=TRUE)}else{"none"})
      shiny::updateTextInput(session, "bandwidth", value=if(is.null(gpOptions$bandwidth)){""}else{gpOptions$bandwidth})
      shiny::updateCheckboxInput(session, "strictLimits",  value=if(is.null(gpOptions$strictLimits)){FALSE}else{gpOptions$strictLimits})
      shiny::updateSliderInput(session, "outliers", value=if(is.null(gpOptions$outliers)){1.5}else{gpOptions$outliers})
    })
    shiny::observeEvent(input$cancel, {
      shiny::stopApp(NULL)
    })
  }
  viewer<-shiny::dialogViewer(dialogName = "Interactive genePlot UI",height = 2000)
  shiny::runGadget(ui, server, viewer=viewer)
}
