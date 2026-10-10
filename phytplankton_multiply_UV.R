############################################################
# PHYTOPLANKTON PATCH VIEWER (Shiny + Leaflet + leaflet-velocity)
#
#  - Everything comes from phytoplankton_simulation_results.csv
#  - Flow : the u / v columns of the table (current at each cell)
#           interpolated to a grid per step (IDW), because the table
#           only holds u/v where cells are. The flow fades to zero
#           away from the cells.
#  - Patch: centroid track (+ drift speed), fading cell trails, points
############################################################

# install.packages(c("shiny", "leaflet", "leaflet.extras2", "jsonlite"))

library(shiny)
library(leaflet)
library(leaflet.extras2)
library(jsonlite)

############################################################
# 1. LOAD DATA
############################################################

SIM_FILE <- "phytoplankton_simulation_results.csv"

if (!file.exists(SIM_FILE)) {
  stop("Could not find ", SIM_FILE, " in the app folder.")
}

##########################################################
# 1a. SIMULATION TABLE (the patch)
##########################################################

sim <- read.csv(SIM_FILE)

needed <- c("simulation_step", "elapsed_hours", "id", "x", "y", "z", "u", "v")
missing_cols <- setdiff(needed, names(sim))
if (length(missing_cols) > 0) {
  stop("Table is missing columns: ", paste(missing_cols, collapse = ", "),
       "\nRe-run the model with the u/v edit to create u and v.")
}

sim <- sim[is.finite(sim$x) & is.finite(sim$y) &
             is.finite(sim$u) & is.finite(sim$v), ]
if (nrow(sim) == 0) stop("The simulation table has no cells with positions and u/v.")

steps   <- sort(unique(sim$simulation_step))
n_steps <- length(steps)

# Centroid (mean position) of the patch for every step
cent <- do.call(rbind, lapply(split(sim, sim$simulation_step), function(d) {
  data.frame(
    step  = d$simulation_step[1],
    hours = d$elapsed_hours[1],
    x     = mean(d$x),
    y     = mean(d$y),
    n     = nrow(d)
  )
}))
cent <- cent[order(cent$step), ]

##########################################################
# 1b. FLOW FIELD FROM THE TABLE (u, v of the cells)
##########################################################

# The table only has u/v at the cell positions, so a continuous field is
# built per step by inverse-distance weighting (IDW) onto a regular grid.
# The field fades to zero with distance from the nearest cell
# ("influence radius"), so it is only meaningful near the patch.

FLOW_PAD <- 2          # degrees of padding around all positions in the table
nx <- 70
ny <- 55
x_min <- min(sim$x) - FLOW_PAD
x_max <- max(sim$x) + FLOW_PAD
y_min <- min(sim$y) - FLOW_PAD
y_max <- max(sim$y) + FLOW_PAD
dx <- (x_max - x_min) / (nx - 1)
dy <- (y_max - y_min) / (ny - 1)

gx_vec <- seq(x_min, x_max, length.out = nx)
gy_vec <- seq(y_max, y_min, length.out = ny)   # north -> south
grid_x <- rep(gx_vec, times = ny)              # row-major: W -> E in each row
grid_y <- rep(gy_vec, each  = nx)

MAX_POINTS <- 1500   # cap for speed if the population grows large

make_field <- function(d, power, radius) {
  
  if (nrow(d) > MAX_POINTS) d <- d[sample(nrow(d), MAX_POINTS), ]
  
  d2 <- outer(grid_x, d$x, "-")^2 + outer(grid_y, d$y, "-")^2
  
  w <- 1 / (d2^(power / 2) + 1e-9)
  w_sum <- rowSums(w)
  
  u <- as.vector(w %*% d$u) / w_sum
  v <- as.vector(w %*% d$v) / w_sum
  
  # fade towards zero away from the cells
  dmin <- sqrt(apply(d2, 1, min))
  fade <- exp(-(dmin / radius)^2)
  
  list(u = u * fade, v = v * fade)
}

# leaflet-velocity JSON (row-major, N -> S, W -> E). addVelocity() takes a
# JSON object / file path / data.frame (not a plain list) via `content`.
make_velocity <- function(u, v) {
  header <- function(param_number) {
    list(
      parameterCategory = 2,
      parameterNumber   = param_number,
      lo1 = x_min,
      lo2 = x_max,
      la1 = y_max,
      la2 = y_min,
      dx  = dx,
      dy  = dy,
      nx  = nx,
      ny  = ny,
      refTime = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    )
  }
  
  u[!is.finite(u)] <- 0
  v[!is.finite(v)] <- 0
  
  jsonlite::toJSON(
    list(
      list(header = header(2), data = u),
      list(header = header(3), data = v)
    ),
    auto_unbox = TRUE,
    digits = 5
  )
}

############################################################
# 2. PATCH HELPERS
############################################################

hav_km <- function(lon1, lat1, lon2, lat2) {
  r <- 6371
  p <- pi / 180
  a <- sin((lat2 - lat1) * p / 2)^2 +
    cos(lat1 * p) * cos(lat2 * p) * sin((lon2 - lon1) * p / 2)^2
  2 * r * asin(sqrt(a))
}

# Covariance of the cells in local kilometres (x scaled by cos(lat))
local_km_cov <- function(d) {
  if (nrow(d) < 3) return(NULL)
  kx <- 111.32 * cos(mean(d$y) * pi / 180)
  ky <- 110.57
  X  <- cbind((d$x - mean(d$x)) * kx, (d$y - mean(d$y)) * ky)
  list(S = cov(X), kx = kx, ky = ky)
}

# RMS patch radius in km
patch_rms_km <- function(d) {
  lc <- local_km_cov(d)
  if (is.null(lc)) return(NA_real_)
  sqrt(sum(diag(lc$S)))
}

############################################################
# 3. MAP BOUNDS (all positions in the table, padded)
############################################################

D_PAD <- 2   # degrees of padding around all positions in the table

d_xlo <- min(sim$x) - D_PAD
d_xhi <- max(sim$x) + D_PAD
d_ylo <- min(sim$y) - D_PAD
d_yhi <- max(sim$y) + D_PAD

# Tight frame around the data (all steps) for the "Reset view" button
r_padx <- max(0.1, 0.1 * (max(sim$x) - min(sim$x)))
r_pady <- max(0.1, 0.1 * (max(sim$y) - min(sim$y)))
r_xlo <- min(sim$x) - r_padx
r_xhi <- max(sim$x) + r_padx
r_ylo <- min(sim$y) - r_pady
r_yhi <- max(sim$y) + r_pady

############################################################
# 3b. VARIABLES FOR THE DISTRIBUTION PLOTS
############################################################

var_info <- list(
  age_hours            = list(label = "Cell age (h)",             scale = 1),
  settling_velocity_ms = list(label = "Settling velocity (m/day)", scale = 86400),
  diameter             = list(label = "Cell diameter (\u00b5m)",   scale = 1),
  z                    = list(label = "Depth (m)",                scale = 1)
)
var_info <- var_info[names(var_info) %in% names(sim)]

dist_choices <- setNames(
  names(var_info),
  vapply(var_info, function(i) i$label, character(1))
)

############################################################
# 4. UI
############################################################

ui <- fluidPage(
  
  tags$head(tags$style(HTML("
    html, body { height: 100%; }
    .container-fluid { padding: 0; }
    .panel-box {
      background: rgba(255,255,255,0.94);
      padding: 14px 16px;
      border-radius: 8px;
      box-shadow: 0 1px 6px rgba(0,0,0,0.3);
    }
    .panel-box h5 { margin: 6px 0 4px 0; font-weight: 600; }
    .info-line { font-size: 12.5px; line-height: 1.5; }
  "))),
  
  leafletOutput("map", width = "100%", height = "100vh"),
  
  absolutePanel(
    top = 82, left = 10, width = "auto", fixed = TRUE,
    actionButton("reset_view", "Reset view", icon = icon("crosshairs"))
  ),
  
  absolutePanel(
    class = "panel-box",
    top = 12, right = 12, width = 320, draggable = TRUE,
    style = "max-height: 95vh; overflow-y: auto;",
    
    h4("Phytoplankton patch"),
    
    sliderInput(
      "step_idx", "Simulation step",
      min = 1, max = n_steps, value = 1, step = 1,
      animate = animationOptions(interval = 2000, loop = TRUE)
    ),
    uiOutput("info"),
    
    hr(),
    h5("Patch"),
    checkboxInput("show_track", "Centroid track", TRUE),
    sliderInput("trail_len", "Cell trails (steps)",
                min = 0, max = 10, value = 3, step = 1),
    checkboxInput("show_cells", "Show points (individual cells)", TRUE),
    
    hr(),
    h5("Distributions"),
    selectInput("dist_var", "Variable",
                choices = dist_choices, selected = dist_choices[1]),
    radioButtons("dist_scope", NULL,
                 choices = c("Current step" = "step",
                             "All steps"    = "all"),
                 selected = "step", inline = TRUE),
    sliderInput("dist_bins", "Bins",
                min = 5, max = 50, value = 20, step = 1),
    plotOutput("dist_plot", height = "200px"),
    verbatimTextOutput("dist_stats"),
    
    hr(),
    h5("Currents (u, v from the table)"),
    checkboxInput("show_flow", "Show flow animation", TRUE),
    sliderInput("radius", "Influence radius (deg)",
                min = 0.5, max = 10, value = 2, step = 0.5),
    sliderInput("power", "IDW power",
                min = 1, max = 4, value = 2, step = 0.5),
    sliderInput("velocity_scale", "Velocity scale",
                min = 0.1, max = 5, value = 1, step = 0.1),
    sliderInput("particle_mult", "Particle density",
                min = 0.5, max = 10, value = 3, step = 0.5),
    sliderInput("line_width", "Line width",
                min = 0.5, max = 4, value = 2, step = 0.5)
  )
)

############################################################
# 5. SERVER
############################################################

server <- function(input, output, session) {
  
  cur_step  <- reactive(steps[input$step_idx])
  step_data <- reactive(sim[sim$simulation_step == cur_step(), ])
  
  ##########################################################
  # INFO PANEL
  ##########################################################
  
  output$info <- renderUI({
    d <- step_data()
    s <- cur_step()
    ci <- which(cent$step == s)
    
    lines <- sprintf("Step %s &middot; %s h elapsed", s, cent$hours[ci])
    lines <- c(lines, sprintf("Cells: <b>%d</b>", nrow(d)))
    lines <- c(lines, sprintf("Centroid: %.3f&deg;E, %.3f&deg;N",
                              cent$x[ci], cent$y[ci]))
    
    rms <- patch_rms_km(d)
    if (is.finite(rms)) {
      lines <- c(lines, sprintf("RMS patch radius: %.1f km", rms))
    }
    
    if (ci > 1) {
      dist_km <- hav_km(cent$x[ci - 1], cent$y[ci - 1], cent$x[ci], cent$y[ci])
      dt_s <- (cent$hours[ci] - cent$hours[ci - 1]) * 3600
      if (dt_s > 0) {
        lines <- c(lines, sprintf("Centroid drift: %.1f km (%.2f m/s)",
                                  dist_km, dist_km * 1000 / dt_s))
      }
    }
    
    if (nrow(d) > 0) {
      lines <- c(lines, sprintf("Mean current at cells: u=%.3f, v=%.3f m/s",
                                mean(d$u), mean(d$v)))
    }
    
    div(class = "info-line", HTML(paste(lines, collapse = "<br>")))
  })
  
  ##########################################################
  # BASE MAP (OpenStreetMap-based)
  ##########################################################
  
  output$map <- renderLeaflet({
    leaflet(options = leafletOptions(preferCanvas = TRUE)) |>
      addProviderTiles(providers$OpenTopoMap,   group = "Terrain (OpenTopoMap)") |>
      addProviderTiles(providers$OpenStreetMap, group = "OpenStreetMap") |>
      fitBounds(d_xlo, d_ylo, d_xhi, d_yhi) |>
      addLayersControl(
        baseGroups = c("Terrain (OpenTopoMap)", "OpenStreetMap"),
        options = layersControlOptions(collapsed = TRUE),
        position = "topleft"
      )
  })
  
  ##########################################################
  # FLOW LAYER (u, v FROM THE TABLE, IDW-INTERPOLATED PER STEP)
  ##########################################################
  
  observe({
    proxy <- leafletProxy("map") |> removeVelocity(group = "velocity")
    if (!isTRUE(input$show_flow)) return()
    
    d <- step_data()
    if (nrow(d) == 0) return()
    
    f <- make_field(d, input$power, input$radius)
    
    max_speed <- max(sqrt(f$u^2 + f$v^2), na.rm = TRUE)
    if (!is.finite(max_speed) || max_speed <= 0) max_speed <- 1
    
    proxy |>
      addVelocity(
        content = make_velocity(f$u, f$v),
        group   = "velocity",
        layerId = "velocity",
        options = velocityOptions(
          speedUnit          = "m/s",
          minVelocity        = 0,
          maxVelocity        = max_speed,
          # tuned for ocean currents (~0.1-1 m/s); the default is for wind
          velocityScale      = (0.05 / max_speed) * input$velocity_scale,
          particleMultiplier = input$particle_mult / 1000,
          lineWidth          = input$line_width,
          opacity            = 0.97,
          colorScale         = c("#1d3557", "#2a6f97", "#6a4c93",
                                 "#c1121f", "#7b0828")
        )
      )
  })
  
  ##########################################################
  # CENTROID TRACK
  ##########################################################
  
  observe({
    proxy <- leafletProxy("map") |> clearGroup("Track")
    if (!isTRUE(input$show_track)) return()
    
    tr <- cent[cent$step <= cur_step(), ]
    if (nrow(tr) == 0) return()
    
    if (nrow(tr) > 1) {
      proxy <- proxy |>
        addPolylines(
          lng = tr$x, lat = tr$y,
          color = "#000000", weight = 3, opacity = 0.9,
          group = "Track"
        )
    }
    
    last <- tr[nrow(tr), ]
    
    proxy |>
      addCircleMarkers(
        lng = tr$x, lat = tr$y,
        radius = 3, weight = 1, color = "#ffffff",
        fillColor = "#000000", fillOpacity = 1,
        group = "Track"
      ) |>
      addCircleMarkers(
        lng = last$x, lat = last$y,
        radius = 8, weight = 2, color = "#ffffff",
        fillColor = "#e63946", fillOpacity = 1,
        label = sprintf("Patch centroid, step %s (%d cells)",
                        last$step, last$n),
        group = "Track"
      )
  })
  
  ##########################################################
  # CELL TRAILS
  ##########################################################
  
  observe({
    proxy <- leafletProxy("map") |> clearGroup("Trails")
    
    k <- input$trail_len
    if (k < 1) return()
    
    s <- cur_step()
    d <- sim[sim$simulation_step <= s & sim$simulation_step >= s - k, ]
    
    # only cells that exist in the current step
    now_ids <- step_data()$id
    d <- d[d$id %in% now_ids, ]
    if (nrow(d) == 0) return()
    
    # cap for speed
    ids <- unique(d$id)
    if (length(ids) > 300) ids <- sample(ids, 300)
    d <- d[d$id %in% ids, ]
    d <- d[order(d$id, d$simulation_step), ]
    
    # one NA-separated vector -> one polyline per cell
    pieces <- split(d, d$id)
    pieces <- pieces[vapply(pieces, nrow, integer(1)) > 1]
    if (length(pieces) == 0) return()
    
    lng <- unlist(lapply(pieces, function(p) c(p$x, NA)))
    lat <- unlist(lapply(pieces, function(p) c(p$y, NA)))
    
    proxy |>
      addPolylines(
        lng = lng, lat = lat,
        color = "#2a6f97", weight = 1.5, opacity = 0.7,
        group = "Trails"
      )
  })
  
  ##########################################################
  # CELLS (OVERLAPPING CELLS MERGED INTO ONE MARKER WITH A COUNT)
  #
  # Daughter cells are created ~100 m from their parent and the patch
  # cohesion step pulls them back together, so they overlap on the map.
  # Cells closer than ~MERGE_PX screen pixels are merged into one marker
  # whose size grows with the number of cells. Zoom in and they separate.
  ##########################################################
  
  MERGE_PX <- 8
  
  observe({
    proxy <- leafletProxy("map") |> clearGroup("Cells")
    
    d <- step_data()
    if (!isTRUE(input$show_cells) || nrow(d) == 0) return()
    
    zoom <- input$map_zoom
    if (is.null(zoom)) zoom <- 6
    
    # size of one screen pixel in degrees at this zoom level
    deg_per_px <- 360 / (256 * 2^zoom)
    bin <- MERGE_PX * deg_per_px
    
    key <- paste(floor(d$x / bin), floor(d$y / bin))
    grp <- match(key, unique(key))
    n   <- as.integer(tabulate(grp))
    
    m <- data.frame(
      x = as.vector(rowsum(d$x, grp)) / n,
      y = as.vector(rowsum(d$y, grp)) / n,
      z = as.vector(rowsum(d$z, grp)) / n,
      n = n
    )
    
    ids_txt <- vapply(split(d$id, grp), function(i) {
      s <- paste(head(i, 8), collapse = ", ")
      if (length(i) > 8) s <- paste0(s, ", ...")
      s
    }, character(1))
    
    popup_txt <- sprintf(
      "<b>%d cell%s here</b><br>ids: %s<br>Mean depth: %.1f m",
      m$n, ifelse(m$n == 1, "", "s"), ids_txt, m$z
    )
    
    proxy <- proxy |>
      addCircleMarkers(
        lng = m$x, lat = m$y,
        radius = pmin(22, 5 + 3 * log2(m$n)),
        weight = 1.5,
        color = "#ffffff",
        fillColor = ifelse(m$n > 1, "#1b9e5a", "#2ecc71"),
        fillOpacity = 0.9,
        popup = popup_txt,
        group = "Cells"
      )
    
    # count labels only on merged markers
    mm <- m[m$n > 1, , drop = FALSE]
    if (nrow(mm) > 0) {
      proxy |>
        addLabelOnlyMarkers(
          lng = mm$x, lat = mm$y,
          label = as.character(mm$n),
          labelOptions = labelOptions(
            noHide = TRUE, direction = "center", textOnly = TRUE,
            style = list("font-weight" = "bold", "font-size" = "11px",
                         "color" = "#ffffff",
                         "text-shadow" = "0 0 3px #0b3d20")
          ),
          group = "Cells"
        )
    }
  })
  
  # re-draw the merged cells when the zoom changes (input$map_zoom is
  # read inside the observer above, so this happens automatically)
  
  ##########################################################
  # RESET VIEW: zoom to the extent of the data only
  ##########################################################
  
  observeEvent(input$reset_view, {
    leafletProxy("map") |>
      fitBounds(r_xlo, r_ylo, r_xhi, r_yhi)
  })
  
  ##########################################################
  # DISTRIBUTIONS (age, settling velocity, diameter, depth)
  ##########################################################
  
  dist_values <- reactive({
    req(input$dist_var)
    info <- var_info[[input$dist_var]]
    d <- if (identical(input$dist_scope, "all")) sim else step_data()
    v <- d[[input$dist_var]] * info$scale
    v[is.finite(v)]
  })
  
  output$dist_plot <- renderPlot({
    req(input$dist_var)
    info <- var_info[[input$dist_var]]
    v <- dist_values()
    
    # fixed x range over ALL steps so histograms are comparable in time
    rng <- range(sim[[input$dist_var]] * info$scale, na.rm = TRUE)
    if (diff(rng) == 0) rng <- rng + c(-0.5, 0.5)
    breaks <- seq(rng[1], rng[2], length.out = input$dist_bins + 1)
    
    scope_txt <- if (identical(input$dist_scope, "all")) {
      "all steps"
    } else {
      paste("step", cur_step())
    }
    
    op <- par(mar = c(4, 4, 3, 1))
    on.exit(par(op))
    
    if (length(v) == 0) {
      plot.new()
      title(main = paste0(info$label, " - ", scope_txt))
      text(0.5, 0.5, "No data")
      return(invisible())
    }
    
    hist(
      v, breaks = breaks,
      col = "#2a6f97", border = "white",
      main = sprintf("%s - %s (n = %d)", info$label, scope_txt, length(v)),
      xlab = info$label,
      ylab = if (identical(input$dist_scope, "all")) "Cell records" else "Cells",
      cex.main = 0.9
    )
    abline(v = median(v), col = "#e63946", lty = 2, lwd = 2)
    legend("topright", legend = "median", col = "#e63946",
           lty = 2, lwd = 2, bty = "n", cex = 0.8)
  })
  
  output$dist_stats <- renderText({
    v <- dist_values()
    if (length(v) == 0) return("No data")
    sprintf(
      "n=%d  mean=%.4g  median=%.4g\nsd=%.4g  min=%.4g  max=%.4g",
      length(v), mean(v), median(v),
      if (length(v) > 1) sd(v) else NA_real_, min(v), max(v)
    )
  })
}

shinyApp(ui, server)