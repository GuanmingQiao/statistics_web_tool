#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(plyr)
require('plyr')
# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(
  
  sidebarLayout(
    sidebarPanel(
      fileInput("csv_file", label = "Upload a file", accept = c('.csv')),
      
      uiOutput("selectdatatype"),
      
      uiOutput("selectdeletetype"),
      
      uiOutput("selectcriter"),
      
      uiOutput("selectpredic"),
      
      uiOutput("bootstrap"),
      
      actionButton("submit", "Submit")
    ),
    
    mainPanel(
      tags$div("Relative Weight",id="Weightmessage"),
      tags$div(""),
      dataTableOutput("results"),
      tags$div(""),
      tags$div("Confidence Interval",id="Confmessage"),
      tags$div(""),
      dataTableOutput("CI"),
      tags$div(""),
      tags$div("Significance",id="Sigmessage"),
      tags$div(""),
      dataTableOutput("Sig"),
      tags$div("Predictor Comparison--(Only Calculated if User Choose to Calculate Predictor's Weight)",id="Predmassage"),
      tags$div(""),
      dataTableOutput("Pred"),
      
      conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                       tags$div("Loading...",id="loadmessage"))
    )
  )
  
))

multRegress <- function(mydata) {
  numVar<<-NCOL(mydata)
  Variables<<- names(mydata)[2:numVar]
  #if input is correlation matrix, this line is unnecessary
  mydata<-cor(mydata, use="pairwise.complete.obs")
  
  RXX<-mydata[2:numVar,2:numVar]
  RXY<-mydata[2:numVar,1]
  
  RXX.eigen<-eigen(RXX)
  D<-diag(RXX.eigen$val)
  delta<-sqrt(D)
  
  lambda<-RXX.eigen$vec%*%delta%*%t(RXX.eigen$vec)
  lambdasq<-lambda^2
  beta<-solve(lambda)%*%RXY
  rsquare<<-sum(beta^2)
  
  RawWgt<-lambdasq%*%beta^2
  import<-(RawWgt/rsquare)*100
  
  result<<-data.frame(Variables, Raw.RelWeight=RawWgt, Rescaled.RelWeight=import)
}

# Define server logic required to draw a histogram
server <- shinyServer(function(input, output) {
  
  #Create data variable here so we can access it later outside of observeEvent()
  data <- reactiveValues(result = NULL, raw_data = NULL, thedata = NULL, CI.Results = NULL, CI.Significance = NULL, CI.Predictor.Comparison=NULL)
  
  output$results <- renderDataTable({
    
    inFile <- input$csv_file
    
    if (is.null(inFile))
      return(NULL)
    
    data$raw_data <- read.csv(inFile$datapath)
    
    #Make the dropdown list input for predictors and criterion
    list_variables <- names(data$raw_data)
    
    output$selectdatatype <- renderUI({
      radioButtons("datatype", "Is the data...", choices = 
                     list("A correlation matrix" = 1, "Raw data" = 2), selected = 1, inline = T)
      
    })
    
    output$selectdeletetype <- renderUI({
      conditionalPanel(
        condition = "input.datatype == 2",
        radioButtons("deletetype", "Missing data options", choices =
                       list("Pairwise deletion" = 1, "Listwise deletion" = 2), selected = 1, inline = T)
      )
    })
    
    real_list_variables <- as.list(list_variables)
    output$selectcriter <- renderUI({
      selectInput("criterion_select", "Choose criterion", real_list_variables)
    })
    
    predic_list_variables <- real_list_variables
    
    output$selectpredic <- renderUI({
      selectInput("predictor_select", "Choose predictors", predic_list_variables, multiple = TRUE)
    })
    
    observeEvent(input$submit, {
      
      criter <- input$criterion_select
      predic <- input$predictor_select
      
      if (is.null(criter) || is.null(predic))
        return(NULL)
      data$thedata<-data$raw_data[,c(criter,predic)]
      
      Labels<-names(data$thedata)[2:length(data$thedata)]
      
      data$result <- multRegress(data$thedata)
    })
    
    data$result
  })
  
  output$bootstrap <- renderUI({
    
    predic <- input$predictor_select
    variables <- names(data$raw_data)
    
    conditionalPanel(
      condition = "input.datatype == 2",
      hr(),
      h5("Bootstrapping"),
      numericInput("replications", "Number of bootstrap replications", 10000, min = NA, max = NA),
      numericInput("alpha", "Alpha value", .05, min = NA, max = NA),
      checkboxInput("testweight", "Test if a predictor's weight is significantly different 
                    from the weight of all other predictors", value = FALSE),
      conditionalPanel(
        condition = "input.testweight == true",
        selectInput("comparepredic", "Predictor to be compared", predic, multiple = FALSE)
      ),
      checkboxInput("checksign", "Test for statistical significance between 2 groups",
                    value = FALSE),
      conditionalPanel(
        condition = "input.checksign == true",
        selectInput("selectgroupvar", "Grouping variable", variables, multiple = FALSE)
      )
    )
  })
  
  output$Sig <-renderDataTable({
    if (is.null(data$thedata))
      return(NULL)
    CI.Significance <- NULL
    
    #only execute if raw data is selected
    if (input$datatype == 2) {
      Labels <-names(data$thedata)[2:length(data$thedata)]
      progress <- shiny::Progress$new(min=0,max=3)
      progress$set(message = "Initializing Output", value = 0)
      on.exit(progress$close())
      
      myRbootci<-function(x){
        boot.ci(multRBoot,conf=0.95,type="bca",index=x)
      }
      
      runRBoot<-function(num){
        INDEX<-1:num
        test<-lapply(INDEX,FUN=myRbootci)
        test2<-t(sapply(test,'[[',i=4))
        CIresultSig<<-data.frame(Labels,Significance.Lower.Bound=test2[,4],Significance.Upper.Bound=test2[,5])
      }
      
      multBootrand<-function(mydata, indices){
        mydata<-mydata[indices,]
        multRWeights<-multRegress(mydata)
        multReps<-multRWeights$Raw.RelWeight
        randWeight<-multReps[length(multReps)]
        randStat<-multReps[-(length(multReps))]-randWeight
        return(randStat)
      }
      
      if (!require("boot", character.only = TRUE)) {
        install.packages("boot")
        library(boot)
      }
      
      #Bootstrapped Confidence interval tests of Significance
      #Please be patient -- This can take a few minutes to run
      randVar<-rnorm(length(data$thedata[,1]),0,1)
      randData<-cbind(data$thedata,randVar)
      progress$set(message = "Bootstrapping Significance", value = 1)
      multRBoot<-boot(randData,multBootrand, 10000)
      progress$set(message = "Bootstrapping Significance CI", value = 2)
      multRci<-boot.ci(multRBoot,conf=0.95, type="bca")
      runRBoot(length(randData[,2:(numVar-1)]))
      progress$set(message = "Writing to data", value = 3)
      data$CI.Significance<-CIresultSig
    }
    if (is.null(data$CI.Significance))
      return(NULL)
    else
      data$CI.Significance
  })
  
  
  output$CI <- renderDataTable({
    
    if (is.null(data$thedata))
      return(NULL)
    
    CI.Results <- NULL
    #only execute if raw data is selected
    if (input$datatype == 2) {
      Labels <-names(data$thedata)[2:length(data$thedata)]
      progress <- shiny::Progress$new(min=0,max=5)
      progress$set(message = "Initializing Output", value = 0)
      on.exit(progress$close())
      
      multBootstrap<-function(mydata, indices){
        mydata<-mydata[indices,]
        multWeights<-multRegress(mydata)
        return(multWeights$Raw.RelWeight)
      }
      
      mybootci<-function(x){
        boot.ci(multBoot,conf=0.95, type="bca", index=x)
      }
      
      
      runBoot<-function(num){
        INDEX<-1:num
        test<-lapply(INDEX, FUN=mybootci)
        test2<-t(sapply(test,'[[',i=4)) #extracts confidence interval
        CIresult<<-data.frame(Variables, CI.Lower.Bound=test2[,4],CI.Upper.Bound=test2[,5])
      }
      
      progress$set(message = "Calculating Confidence Interval", value = 1)
      
      if (!require("boot", character.only = TRUE)) {
        install.packages("boot")
        library(boot)
      }
      #Normal Multiple Bootstrap
      progress$set(message = "Running Bootstrap", value = 2)
      multBoot<-boot(data$thedata, multBootstrap, 10000)
      progress$set(message = "Bootstrapping Confidence Interval", value = 3)
      multci<-boot.ci(multBoot,conf=0.95, type="bca")
      progress$set(message = "Finalizing Output Length", value = 4)
      runBoot(length(data$thedata[,2:numVar]))
      progress$set(message = "Outputting CI", value = 5)
      data$CI.Results<-CIresult
    }
    if (is.null(data$CI.Results))
      return(NULL)
    else
      data$CI.Results
  })
  
  output$Pred <-renderDataTable({
    if (is.null(data$thedata))
      return(NULL)
    CI.Predictor.Comparison <- NULL
    
    #only execute if raw data is selected and testweight is selected
    if (input$datatype == 2 & input$testweight == TRUE) {
      Labels <-names(data$thedata)[2:length(data$thedata)]
      progress <- shiny::Progress$new(min=0,max=6)
      progress$set(message = "Initializing Output", value = 0)
      on.exit(progress$close())
      
      progress$set(message = "Performing Multiple Bootstraps", value = 1)
      multBootcomp<-function(mydata, indices){
        mydata<-mydata[indices,]
        multCWeights<-multRegress(mydata)
        multCeps<-multCWeights$Raw.RelWeight
        comp2Stat<-multCeps-multCeps[1]
        comp2Stat<-comp2Stat[-1]
        Labels2<<-Labels[-1]
        return(comp2Stat)
      }
      
      myCbootci<-function(x){
        boot.ci(multC2Boot,conf=0.95,type="bca",index=x)
      }
      
      progress$set(message = "Evaluating Comparative Bootstrap", value = 2)
      runCBoot<-function(num){
        INDEX<-1:num
        test<-lapply(INDEX,FUN=myCbootci)
        test2<-t(sapply(test,'[[',i=4))
        CIresultComp<<-data.frame(Labels2,Comparative.Importance.Lower.Bound=test2[,4],Comparative.Importance.Upper.Bound=test2[,5])
      }
    
      if (!require("boot", character.only = TRUE)) {
        install.packages("boot")
        library(boot)
      }
      
      progress$set(message = "Finalizing Comparative Output", value = 3)
      multC2Boot<-boot(data$thedata, multBootcomp, 10000)
      progress$set(message = "4", value = 4)
      multC2ci<-boot.ci(multC2Boot,conf=0.95, type="bca")
      progress$set(message = "5", value = 5)
      runCBoot(numVar-2)
      progress$set(message = "6", value = 6)
      data$CI.Predictor.Comparison<-CIresultComp
      
    }
    if (is.null(data$CI.Predictor.Comparison))
      return(NULL)
    else
      data$CI.Predictor.Comparison
  })
  
})

# Run the application 
shinyApp(ui = ui, server = server)