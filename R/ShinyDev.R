library(shiny)
library(miniUI)

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
