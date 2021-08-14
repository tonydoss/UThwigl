
if (!require("shiny")) install.packages("shiny")
if (!require("shinycssloaders")) install.packages("shinycssloaders")
if (!require("DT")) install.packages("DT")
if (!require("devtools")) install.packages("devtools")
if (!require("UThwigl")) devtools::install_github("benmarwick/UThwigl")
if (!require("waiter")) install.packages("waiter")
if (!require("rhandsontable")) install.packages("rhandsontable")

require(devtools)
require(shiny)
require(shinycssloaders)
require(UThwigl)
require(DT)
require(cowplot)
require(waiter)
require(rhandsontable)

waiting_screen <- tagList(
  spin_solar(),
  HTML("<br><br><br><br><br>Generating model output...")
) 

# Define UI 
ui <- bootstrapPage(
  mainPanel(
    use_waiter(),
    # Application title
    titlePanel("UThwigl::csUTh : compute closed-system uranium-thorium ages"),
    
    tabsetPanel(   id = "inTabset",
      
      tabPanel("Load the data", 
               p("Before uploading, check that your CSV file contains columns with these names:"),
               HTML("
               <li> <b>Sample_ID</b>: sample identification code
               <li> <b>U234_U238_CORR</b>: [<sup>234</sup>U/<sup>238</sup>U] activity ratios 
               <li> <b>U234_U238_CORR_Int2SE</b>: the 2 sigma errors of the activity ratios
               <li> <b>Th230_U238_CORR</b>: [<sup>230</sup>Th/<sup>238</sup>U] activity ratios 
               <li> <b>Th230_U238_CORR_Int2SE</b>: the 2 sigma errors of the activity ratios
               <li> <b>Th232_U238_CORR</b>: [<sup>232</sup>Th/<sup>238</sup>U] activity ratios 
               <li> <b>Th232_U238_CORR_Int2SE</b>:   the 2 sigma errors of the activity ratios 
              "),
               tags$hr(),
               fileInput("file1", 
                         "Choose CSV file", 
                         accept = c("text/csv", 
                                    "text/comma-separated-values,text/plain", 
                                    ".csv")),
               actionButton('gotoinspect', 'Go to inspect the data')          
               ), # end of tab
      
      tabPanel("Inspect the data",
               value = "inspectthedata",
               p("Here is the raw data from the CSV file"),
               DT::dataTableOutput('contents'),
               actionButton('gotosetmodel', 'Go to set the model parameters')   
      ), # end of tab
      
      tabPanel("Set model parameters",
               value = "setmodelparameters",
               fluidRow(
                 column(4, 
                        # defaults for Pan2018
                        textInput(inputId = "sample_name", label = "Sample name", value = "YP003"),
                        numericInput("nbitchoice", "Number of iterations:", 100, min = 1, max = 1e6),
                        numericInput("detcorrectionchoice", "Do a detrital correction? (1 = yes, 0 = no):", 1, min = 0, max = 1),
                        numericInput("keepfiltereddata", "Save filtered data on which an outlier test was performed? (1 = yes, 0 = no):", 0, min = 0, max = 1),                        

                        
                 ),
                 column(4,
                        numericInput("R28det", "(232Th/238U) activity ratio of the detritus:", 0.8, min = 0.01, max = 10),
                        numericInput("R28det_err", "Error on the (232Th/238U) activity ratio of the detritus:", 0.08, min = 0.01, max = 10),
                        numericInput("R08det", "(230Th/238U) activity ratio of the detritus:", 1, min = 0.01, max = 10),
                        numericInput("R08det_err", "Error on the (230Th/238U) activity ratio of the detritus:", 0.05, min = 0.01, max = 10),
                        numericInput("R48det", "(234U/238U) activity ratio of the detritus:", 1, min = 0.5, max = 50),
                        numericInput("R48det_err", "Error on the (234U/238U) activity ratio of the detritus:", 0.02, min = 0.01, max = 50)

                 )
               ),
               actionButton("run", label = "Run simulation and visualise the output")
               
      ),  # end of tab
      
      tabPanel("Visualise the model",
               value = "visualise",
               HTML("<p><b>Plot legend</b><p>
                <b>A.</b> Closed-system ages<p>
                <b>B.</b> Initial [<sup>234</sup>U/<sup>238</sup>U] activity ratios for each sample analysis<p>
                <p>"),
                         tags$hr(),
                         # defaults for Pan2018
                         # show a spinner while we wait for the plots to draw
                         #  withSpinner()
                       plotOutput("plots", 
                                  width = "100%", 
                                  height = "600px"),
                                  color="blue", 
                                  size = 5,
      actionButton('gotooutput', 'Go to inspect the output')
      ), # end of tab
                
                tabPanel("Inspect the model", 
                         value = "modeloutput",
                         # defaults for Pan2018
                         tableOutput("print_summary"),
                         tags$hr(),
                         tableOutput("model_results_table"),
                         tags$hr()
                         
                ) # end of tab
                
                
                # end  tabset
    )))


server <- function(input, output, session) {
  
  # activate the buttons to move between tabs
  observeEvent(input$gotoinspect, {
    updateTabsetPanel(session, "inTabset",
                      selected = "inspectthedata")
  })
  
  observeEvent(input$gotosetmodel, {
    updateTabsetPanel(session, "inTabset",
                      selected = "setmodelparameters"  )
  })
  
  observeEvent(input$run, {
    updateTabsetPanel(session, "inTabset",
                      selected = "visualise"  )
  })
  
  observeEvent(input$gotooutput, {
    updateTabsetPanel(session, "inTabset",
                      selected = "modeloutput"  )
  })
  
  fname = tempfile(fileext = ".csv")
  
  observe({
    # remove button and isolate to update file automatically
    # after each table change
    input$saveBtn
    hot = isolate(input$hot)
    if (!is.null(hot)) {
      write.csv(hot_to_r(input$hot), fname)
      print(fname)
    }
  })
  
  output$hot = renderRHandsontable({
    if (!is.null(input$hot)) {
      DF <<-  hot_to_r(input$hot)
    } else {
      DF <<-  read.csv(system.file("extdata/input", "Pan2018.csv", package = "UThwigl"), stringsAsFactors = FALSE) 
    }
    
    rhandsontable(DF) %>%
      hot_table(highlightCol = TRUE, highlightRow = TRUE)
  })


  
  output$contents <- DT::renderDataTable({
    
    
    inFile <- input$file1
    
    if (is.null(inFile))
      return(DF)
    
    read.csv(inFile$datapath)
  })
  
  observeEvent(input$run, {
    updateTabsetPanel(session, inputId = "csUTh", selected = "vis")
  })
  
  
  observeEvent(input$run, {
    model_output <-  reactive({
      
      waiter_show(
        html = waiting_screen,
        color = "#333e48",
        id = 'plots',
        hide_on_render = TRUE
      )
      
      showNotification("Starting model run...")

      inFile <- input$file1
      if (is.null(inFile)) {
        input_data <- DF
        } else { 
        input_data <- read.csv(inFile$datapath)
          }

      output <- 
        csUTh(input_data,
              sample_name = input$sample_name,
              nbitchoice = input$nbitchoice,
              detcorrectionchoice = input$detcorrectionchoice,
              keepfiltereddata = input$keepfiltereddata,
              R28det = input$R28det,
              R28det_err = input$R28det_err,
              R08det = input$R08det,
              R08det_err = input$R08det_err,
              R48det = input$R48det,
              print_summary = FALSE,
              with_plots = FALSE)
      
      showNotification("Model run complete.")
      output
    })
    
    # get some of the output from the function to display
    
    output$model_results_table <- renderTable({ model_output()$results })
    output$print_summary <- renderText({ model_output()$print_summary })
    
    # draw the plots
    
    output$plots <- renderPlot({ 
      
      the_data <- model_output()
      
      library(ggplot2)
      p2 <- initial_234U_238U_plot(the_data) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
      
      # plot ages
      p1 <- ages_plot(the_data) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
      
      # draw plots in a panel
      p3 <- cowplot::plot_grid(p1, p2, ncol = 2)
      
      p3
      
      
      
    })
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
