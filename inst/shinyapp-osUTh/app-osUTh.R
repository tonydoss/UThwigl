
if (!require("shiny")) install.packages("shiny")
if (!require("shinycssloaders")) install.packages("shinycssloaders")
if (!require("DT")) install.packages("DT")
if (!require("devtools")) install.packages("devtools")
if (!require("UThwigl")) devtools::install_github("tonydoss/UThwigl")
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
    titlePanel("UThwigl::osUTh : compute open-system uranium-thorium ages using the diffusion-adsorption-decay (DAD) model"),
    
    tabsetPanel(   id = "abset",
      
      tabPanel("Load the data", 
               p("Before uploading, check that your CSV file contains columns with these names:"),
               HTML("
               <li> <b>U234_U238</b>: [<sup>234</sup>U/<sup>238</sup>U] activity ratios 
               <li> <b>U234_U238_2SE</b>: the 2 sigma errors of the activity ratios
               <li> <b>Th230_U238</b>: [<sup>230</sup>Th/<sup>238</sup>U] activity ratios 
               <li> <b>Th230_U238_2SE</b>: the 2 sigma errors of the activity ratios
               <li> <b>U_ppm</b>: uranium concentrations (in ppm)
               <li> <b>U_ppm_2SE</b>: the 2 sigma errors of the uranium concentrations
               <li> <b>x</b>: x coordinates (in mm) of the analyses, and the outer and inner surfaces of the sample 
               <li> <b>y</b>: y coordinates (in mm) of the analyses, and the outer and inner surfaces of the sample 
               <li> <b>Comments</b>: two rows must show 'outer surface' and 'inner surface' (with the corresponding x and y coordinates) 

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
                        # defaults for Hobbit_MH2T
                        numericInput("nbit", "Number of iterations:", 100, min = 1, max = 1e6),
                        numericInput("fsumtarget", "Value of squared sum", 0.01, min = 0.001, max = 100),
                        # numericInput("l", "Thickness of sample (cm):", 5.35, min = 0.01, max = 10),
                        numericInput("U_0", "Uranium concentration at the sample surface (ppm):", 25, min = 0.01, max = 500),
                        numericInput("U48_0_min", "Min [234U/238U] at the surface:", 1.265, min = 0.5, max = 50),
                        numericInput("U48_0_max", "Max [234U/238U] at the surface:", 1.275, min = 0.5, max = 50)
                 ),
                 column(4,
                        numericInput("T_min", "Age min (yr):", 1e3, min = 1e3, max = 500e3),
                        numericInput("T_max", "Age max (yr):", 20e3, min = 1e3, max = 500e3),
                        numericInput("K_min", "Min U diffusion coefficient:", 1e-13, min = 1e-16, max = 1e-9),
                        numericInput("K_max", "Max U diffusion coefficient:", 1e-11, min = 1e-16, max = 1e-9)
                 )
               ),
               actionButton("run", label = "Run simulation and visualise the output")
               
      ),  # end of tab
      
      tabPanel("Visualise the model",
               value = "visualise",
               HTML("<p><b>Plot legend</b><p>
                <b>A.</b> A histogram of the solution ages. <p>
                <b>B.</b> The U concentrations in the sample as a function of the relative distance from the center. <p>
                <b>C.</b> The measured (in blue) and modelled (in red) [<sup>234</sup>U/<sup>238</sup>U] activity ratios as a function of the relative distance from the center, and <p>
                <b>D.</b> The measured (in blue) and modelled (in red) [<sup>230</sup>Th/<sup>238</sup>U] activity ratios as a function of the relative distance from the center
                <p>"),
                         tags$hr(),
                         # defaults for Hobbit_MH2T
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
                         # defaults for Hobbit_MH2T
                         tableOutput("model_results_table"),
                         tags$hr(),
                         tableOutput("model_results_items"),
                         tags$hr(),
                         tableOutput("U48_0_final"),
                         tags$hr(),
                         tableOutput("output_data")
                         
                ) # end of tab
                
                
                # end  tabset
    )))


server <- function(input, output, session) {
  
  # activate the buttons to move between tabs
  observeEvent(input$gotoinspect, {
    updateTabsetPanel(session, "abset",
                      selected = "inspectthedata")
  })
  
  observeEvent(input$gotosetmodel, {
    updateTabsetPanel(session, "abset",
                      selected = "setmodelparameters"  )
  })
  
  observeEvent(input$run, {
    updateTabsetPanel(session, "abset",
                      selected = "visualise"  )
  })
  
  observeEvent(input$gotooutput, {
    updateTabsetPanel(session, "abset",
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
      pr(fname)
    }
  })
  
  output$hot = renderRHandsontable({
    if (!is.null(input$hot)) {
      DF <<-  hot_to_r(input$hot)
    } else {
      DF <<-  read.csv(system.file("extdata/input", "Hobbit_MH2T_for_iDAD.csv", package = "UThwigl"), stringsAsFactors = FALSE) 
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
    updateTabsetPanel(session, inputId = "osUTh", selected = "vis")
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
        osUTh(input_data,
                 nbit = input$nbit,
                 fsum_target = input$fsumtarget,
                 U48_0_min = input$U48_0_min,
                 U48_0_max = input$U48_0_max,
                 U_0 = input$U_0,
                 K_min = input$K_min,
                 K_max = input$K_max,
                 T_min = input$T_min,
                 T_max = input$T_max,
                 print_age = FALSE,
                 with_plots = FALSE)
      
      showNotification("Model run complete.")
      output
    })
    
    # get some of the output from the function to display
    
    output$model_results_table <- renderTable({ model_output()$results })
    
    output$model_results_items <- renderTable({ data.frame(model_output()[2:5]) })
    
    output$U48_0_final <- renderTable({ model_output()$U48_0_final })
    
    output$output_data <- renderTable({ model_output()$output_data })
    
    # draw the plots
    
    output$plots <- renderPlot({ 
      
      the_data <- model_output()
      
      big_size = 18
      less_big_size = 14
      point_size = 4
      
      T_sol_plot <- T_sol_plot(the_data, 
                               big_size = big_size, 
                               less_big_size = less_big_size)
      
      u_conc_t2_plot_out <- u_conc_profile_plot(the_data, 
                                                big_size = big_size, 
                                                less_big_size = less_big_size,
                                                point_size = point_size)
      
      u234_u238_ratio_plot_out <- u234_u238_ratio_plot(the_data, 
                                                       big_size = big_size, 
                                                       less_big_size = less_big_size, 
                                                       point_size = point_size)
      
      th230_u238_ratio_plot_out <- th230_u238_ratio_plot(the_data, 
                                                         big_size = big_size, 
                                                         less_big_size = less_big_size,
                                                         point_size = point_size)
      
      # combine the plots o one panel
      p1 <-
        cowplot::plot_grid(T_sol_plot,
                           u_conc_t2_plot_out,
                           u234_u238_ratio_plot_out,
                           th230_u238_ratio_plot_out,
                           labels = "AUTO",
                           ncol = 2)
      
      # pr the plots
      p1
      
      
    })
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
