library(shiny)
library(bslib)
library(leaflet)
library(leaflet.extras)
library(shinyWidgets)

# 1. Flatten/prepare your dataset
# (Assumes simulation_results is already loaded in your environment)
combined_data <- simulation_results

combined_data$simulation_step <- as.integer(combined_data$simulation_step)
combined_data$x <- as.numeric(combined_data$x)
combined_data$y <- as.numeric(combined_data$y)
combined_data$age_hours <- as.numeric(combined_data$age_hours)

unique_steps <- sort(unique(combined_data$simulation_step))

# Pre-calculate max age globally for stable color scaling across frames
max_age_val <- max(combined_data$age_hours, na.rm = TRUE)
if(is.infinite(max_age_val) || max_age_val == 0) max_age_val <- 1

# 2. Shiny UI Setup 
ui <- page_fluid(
  titlePanel("Simulation Step Heatmap Animation"),
  
  card(
    # Custom CSS positioning overrides
    tags$head(
      tags$style(HTML("
        /* Create a relative container so the absolute layout binds tightly to the card bounds */
        .map-wrapper-box {
          position: relative;
          width: 100%;
          height: 650px;
        }
        
        /* Force the Leaflet map to take up the full wrapper dimensions */
        #heatmap_plot {
          width: 100% !important;
          height: 100% !important;
        }
        
        /* absolute overlay sitting flush on the bottom edge */
        .absolute-full-timeline {
          position: absolute;
          bottom: 0;
          left: 0;
          width: 100% !important;
          z-index: 1000; /* Forces the bar to float cleanly above map tiles */
          background-color: rgba(255, 255, 255, 0.25) !important; /* Low opacity */
          backdrop-filter: blur(5px); /* Frosted glass aesthetic */
          padding: 15px 40px 10px 40px; /* Generates internal slider breathing room */
          border-top: 1px solid rgba(0,0,0,0.1);
        }
        
        /* Style text variables so they pop cleanly above map backgrounds */
        .absolute-full-timeline label, .absolute-full-timeline .control-label, .absolute-full-timeline #value5 {
          font-weight: bold;
          color: #111111;
          text-shadow: 1px 1px 2px rgba(255,255,255,0.8);
        }
        
        /* Adjust Leaflet built-in controls upward so they don't get blocked */
        .leaflet-bottom {
          bottom: 110px !important;
        }
      "))
    ),
    card_header("Simulation Map View - Zoom-Responsive Tighter Clusters"),
    
    # Combined Layout Box
    div(class = "map-wrapper-box",
        leafletOutput("heatmap_plot"),
        
        # Absolute timeline overlay
        div(class = "absolute-full-timeline",
            sliderTextInput(
              inputId = "slider5",
              label = "Simulation Timeline Control",
              choices = unique_steps,       
              selected = min(unique_steps),   
              grid = TRUE,                 
              animate = animationOptions(interval = 100, loop = TRUE),
              width = "100%"
            ),
            textOutput("value5")
        )
    )
  )
)

# 3. Shiny Server Setup
server <- function(input, output, session) {
  
  # Render the base map
  output$heatmap_plot <- renderLeaflet({
    leaflet() %>%
      addTiles() %>%
      setView(lng = mean(combined_data$x, na.rm = TRUE), 
              lat = mean(combined_data$y, na.rm = TRUE), 
              zoom = 5)
  })
  
  # Text status engine with safeguard indexing
  output$value5 <- renderText({
    req(input$slider5)
    hours_val <- unique(combined_data$elapsed_hours[combined_data$simulation_step == input$slider5])
    paste0("Step: ", input$slider5, " | Time: ", hours_val[1], " hours")
  })
  
  # Map updater pipeline reacting to both timeline slider changes AND map zooming
  observe({
    req(input$slider5)
    
    # Capture current zoom level (defaults to 5 if map hasn't registered yet)
    current_zoom <- if (!is.null(input$heatmap_plot_zoom)) input$heatmap_plot_zoom else 5
    
    # Dynamically scale down radius and blur as zoom increases 
    # This keeps points tight and crisp instead of blowing up into giant blobs when zooming in
    dynamic_radius <- max(6, 30 - (current_zoom - 5) * 2.5)
    dynamic_blur   <- max(3, dynamic_radius * 0.4)
    
    step_slice <- combined_data[combined_data$simulation_step == input$slider5, ]
    step_slice <- step_slice[!is.na(step_slice$x) & !is.na(step_slice$y) & !is.na(step_slice$age_hours), ]
    
    # Normalize age for heatmap intensity scaling (0 to 1)
    step_slice$norm_age <- pmin(1, step_slice$age_hours / max_age_val)
    
    leafletProxy("heatmap_plot", data = step_slice) %>%
      clearHeatmap() %>%  
      addHeatmap(
        lng = ~x, 
        lat = ~y, 
        radius = dynamic_radius, # Scales tighter on zoom-in
        blur = dynamic_blur,     # Scales sharper on zoom-in
        max = 1,
        minOpacity = 0.4,
        intensity = ~norm_age, 
        gradient = c(
          "0.0" = "blue",   # Young cells
          "0.5" = "green",  # Maturing cells
          "1.0" = "red"     # Older cells
        )
      )
  })
}

shinyApp(ui, server)