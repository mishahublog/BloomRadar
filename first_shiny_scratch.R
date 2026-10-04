library(shiny)

ui <- fluidPage(
  titlePanel("Reactive Expression Demo"),
  sliderInput("threshold", "Select Minimum Value:", min = 10, max = 50, value = 20),
  sidebarLayout(
    sidebarPanel(
      # User selects a threshold value
      
      # A secondary UI input that shouldn't trigger data recalculation
      textInput("plot_title", "Custom Plot Title:", value = "My Data Plot")
    ),
    mainPanel(
      plotOutput("data_plot"),
      tableOutput("data_summary")
    )
  )
)

server <- function(input, output) {
  
  # 1. THE REACTIVE EXPRESSION
  # This runs ONCE when input$threshold changes. 
  # It does NOT re-run if input$plot_title changes!
  filtered_data <- reactive({
    message("Running heavy data processing step...") # Look at your R console to see when this triggers
    
    # Simulating a dataset (e.g., environmental observations)
    base_data <- data.frame(
      x = 1:100,
      y = rnorm(100, mean = 30, sd = 10)
    )
    
    # Filter data based on the slider input
    subset(base_data, y > input$threshold)
  })
  
  # 2. Plot Output
  output$data_plot <- renderPlot({
    # We call our reactive expression like a function: filtered_data()
    df <- filtered_data() 
    
    plot(df$x, df$y, main = input$plot_title, xlab = "Index", ylab = "Value", 
         col = "royalblue", pch = 16)
  })
  
  # 3. Summary Table Output
  output$data_summary <- renderTable({
    # Pulls the cached data instantly! Does NOT re-run the subsetting code.
    df <- filtered_data() 
    
    # Return a quick summary statistics table
    data.frame(
      Total_Points = nrow(df),
      Mean_Value = mean(df$y),
      Max_Value = max(df$y)
    )
  })
}

shinyApp(ui, server)