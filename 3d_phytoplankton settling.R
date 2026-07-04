library(plotly)
library(colorRamps)

# 1. Ensure data is ordered sequentially by cell ID and time step
plot_data <- master_export_df[order(master_export_df$id, master_export_df$simulation_step), ]

# 2. Invert the Z axis value for standard oceanographic mapping 
# This ensures 0m is at the surface and higher numbers drop downwards
plot_data$z_plot <- -abs(plot_data$z)

# 3. Initialize the interactive 3D canvas
fig <- plot_ly(
  data = plot_data, 
  x = ~x, 
  y = ~y, 
  z = ~z_plot, 
  type = 'scatter3d', 
  mode = 'lines+markers',
  # Color lines by cell ID to distinguish individual particle tracks
  color = ~as.factor(id), 
  colors = colorRamps::matlab.like(length(unique(plot_data$id))),
  marker = list(
    size = 3.5, 
    opacity = 0.7,
    # Color markers by age to show maturity progression along the track
    color = ~age_hours, 
    colorscale = 'Viridis',
    colorbar = list(title = "Cell Age (hrs)", len = 0.5)
  ),
  line = list(width = 3)
)

# 4. Configure clean 3D axis titles and bounding layout box
fig <- fig %>% layout(
  title = "3D Phytoplankton Particle Trajectories & Mixed Layer Drift",
  scene = list(
    xaxis = list(title = "Longitude (°E)"),
    yaxis = list(title = "Latitude (°N)"),
    zaxis = list(title = "Depth Profile (m)")
  ),
  margin = list(l = 0, r = 0, b = 0, t = 50)
)

# Render the interactive plot window
fig
