library(shiny)
library(iDADwigl)

# Define UI for application that draws a histogram
ui <- fluidPage(
  # Application title
  titlePanel("iDADwigl"),
  fileInput("fileinput", "Choose CSV file", accept = c(
    "text/csv",
    "text/comma-separated-values,text/plain",
    ".csv")
  ),
  numericInput("nbit", "Number of iterations:", 1, min = 1, max = 1e6),
  numericInput("fsumtarget", "Value of squared sum", 0.05, min = 0.001, max = 100),
  numericInput("U48_0_min", "Min (234U/238U) at the surface:", 1, min = 0.5, max = 50),
  numericInput("U48_0_max", "Max (234U/238U) at the surface:", 1, min = 0.5, max = 50),
  numericInput("l", "Thickness of sample (cm):", 1, min = 0.01, max = 10),
  numericInput("U_0", "Uranium concentration at the sample surface (ppm):", 1, min = 0.01, max = 500),
  numericInput("K_min", "Min U diffusion coefficient:", 1e-12, min = 1e-16, max = 1e-9),
  numericInput("K_max", "Max U diffusion coefficient:", 1e-12, min = 1e-16, max = 1e-9),
  numericInput("T_min", "Age min (yr):", 1e4, min = 1e3, max = 500e3),
  numericInput("T_max", "Age max (yr):", 1e4, min = 1e3, max = 500e3),
  verbatimTextOutput("value")
)


server <- function(input, output) {
  # output <- iDADwigl(input$fileinput,
  #                    input$nbit,
  #                    input$fsum_target,
  #                    input$U48_0_min,
  #                    input$U48_0_max,
  #                    input$l,
  #                    input$U_0,
  #                    input$K_min,
  #                    input$K_max,
  #                    input$T_min,
  #                    input$T_max)
}

# Run the application 
shinyApp(ui = ui, server = server)

