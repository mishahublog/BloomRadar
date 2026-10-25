############################################################
# PHYTOPLANKTON PATCH VIEWER (Shiny + Leaflet + leaflet-velocity)
#
#  - Everything comes from phytoplankton_simulation_results.csv
#  - Flow : the u / v columns of the table (current at each cell)
#           interpolated to a grid per step (IDW), because the table
#           only holds u/v where cells are. The flow fades to zero
#           away from the cells.
#  - Patch: HEAT MAP of cell density (+ cell counts as numbers),
#           centroid track (+ drift speed), fading cell trails
############################################################

# install.packages(c("shiny", "leaflet", "leaflet.extras", "leaflet.extras2", "jsonlite"))

library(shiny)
library(leaflet)
library(leaflet.extras2)    # addVelocity()
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

local_km_cov <- function(d) {
  if (nrow(d) < 3) return(NULL)
  kx <- 111.32 * cos(mean(d$y) * pi / 180)
  ky <- 110.57
  X  <- cbind((d$x - mean(d$x)) * kx, (d$y - mean(d$y)) * ky)
  list(S = cov(X), kx = kx, ky = ky)
}

patch_rms_km <- function(d) {
  lc <- local_km_cov(d)
  if (is.null(lc)) return(NA_real_)
  sqrt(sum(diag(lc$S)))
}

# Push points apart until none are closer than `min_dist` (same units as
# x / y). Used so the count numbers and depth labels do not pile up.
repel_xy <- function(x, y, min_dist, iters = 40) {
  n <- length(x)
  if (n < 2 || min_dist <= 0) return(list(x = x, y = y))
  x <- x + sin(seq_len(n) * 12.9898) * 1e-3 * min_dist   # separates exact overlaps
  y <- y + cos(seq_len(n) * 78.233)  * 1e-3 * min_dist
  for (k in seq_len(iters)) {
    ddx <- outer(x, x, "-")
    ddy <- outer(y, y, "-")
    dm  <- sqrt(ddx^2 + ddy^2)
    diag(dm) <- Inf
    close <- dm < min_dist
    if (!any(close)) break
    push <- ifelse(close, (min_dist - dm) / (2 * dm), 0)
    x <- x + rowSums(push * ddx)
    y <- y + rowSums(push * ddy)
  }
  list(x = x, y = y)
}

############################################################
# 3. MAP BOUNDS
############################################################

D_PAD <- 2

d_xlo <- min(sim$x) - D_PAD
d_xhi <- max(sim$x) + D_PAD
d_ylo <- min(sim$y) - D_PAD
d_yhi <- max(sim$y) + D_PAD

r_padx <- max(0.1, 0.1 * (max(sim$x) - min(sim$x)))
r_pady <- max(0.1, 0.1 * (max(sim$y) - min(sim$y)))
r_xlo <- min(sim$x) - r_padx
r_xhi <- max(sim$x) + r_padx
r_ylo <- min(sim$y) - r_pady
r_yhi <- max(sim$y) + r_pady

############################################################
# 3a. HEAT MAP GRADIENT (low density -> high density)
############################################################

# Weather-radar style colour scale (weak -> strong), like a reflectivity map
radar_cols <- c("#4a5fd0", "#00b4f0", "#00e000", "#009000",
                "#ffff00", "#ffa500", "#ff0000", "#a00000", "#ff00ff")

# Fixed grid over the whole data extent. Cell counts per grid box are
# pre-computed for every step; smoothing (a Gaussian) is applied in the app.
RAD_PADX <- max(0.2, 0.3 * (max(sim$x) - min(sim$x)))
RAD_PADY <- max(0.2, 0.3 * (max(sim$y) - min(sim$y)))
rad_xlo <- min(sim$x) - RAD_PADX
rad_xhi <- max(sim$x) + RAD_PADX
rad_ylo <- min(sim$y) - RAD_PADY
rad_yhi <- max(sim$y) + RAD_PADY
GNX <- 140
GNY <- 110
rdx <- (rad_xhi - rad_xlo) / GNX
rdy <- (rad_yhi - rad_ylo) / GNY

count_grid <- function(d) {
  ix <- pmin(GNX, pmax(1, floor((d$x - rad_xlo) / rdx) + 1))
  iy <- pmin(GNY, pmax(1, floor((d$y - rad_ylo) / rdy) + 1))
  matrix(tabulate((ix - 1) * GNY + iy, nbins = GNX * GNY), nrow = GNY, ncol = GNX)
}
count_list <- lapply(steps, function(s) count_grid(sim[sim$simulation_step == s, ]))

gauss_mat <- function(n, sigma) {
  idx <- seq_len(n)
  exp(-outer(idx, idx, "-")^2 / (2 * sigma^2))
}
smooth_grid <- function(C, sigma) {
  gauss_mat(GNY, sigma) %*% C %*% gauss_mat(GNX, sigma)
}

z_dom <- range(sim$z[is.finite(sim$z)])
if (diff(z_dom) == 0) z_dom <- z_dom + c(-0.5, 0.5)

# depth colours: shallow = yellow, deep = dark purple
depth_pal <- colorBin(
  palette = "viridis", domain = z_dom, bins = 6, pretty = TRUE,
  reverse = TRUE
)

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

    .bottom-bar {
      position: fixed; left: 50%; bottom: 14px;
      transform: translateX(-50%);
      width: min(720px, 86vw);
      padding: 0 6px;
      background: transparent;
      z-index: 1000;
    }
    .bottom-bar .form-group { margin-bottom: 0; width: 100%; }
    .bottom-bar .irs { width: 100%; }
    .bottom-bar label.control-label,
    .bottom-bar .irs-min, .bottom-bar .irs-max, .bottom-bar .irs-grid-text {
      color: #111 !important; font-weight: 600;
      text-shadow: 0 0 3px #fff, 0 0 3px #fff, 0 0 6px #fff !important;
    }
    .bottom-bar .irs-min, .bottom-bar .irs-max {
      background: rgba(255,255,255,0.7) !important;
    }
    .bottom-bar .irs-line { box-shadow: 0 0 4px rgba(255,255,255,0.9); }
    .bottom-bar .slider-animate-container { text-align: left; margin-top: 2px; }
    .bottom-bar .slider-animate-button {
      background: rgba(255,255,255,0.85); color: #111;
      border-radius: 4px; padding: 0 8px; text-shadow: none;
    }

    .leaflet-bottom { bottom: 0; }
    @media (max-width: 1150px) { .leaflet-bottom { bottom: 84px; } }

    .settings-modal {
      position: fixed; top: 0; left: 0; right: 0; bottom: 0;
      z-index: 3000;
      display: flex; align-items: center; justify-content: center;
      background: rgba(0,0,0,0.45);
      visibility: hidden; opacity: 0; pointer-events: none;
      transition: opacity .15s ease, visibility .15s;
    }
    .settings-modal.open {
      visibility: visible; opacity: 1; pointer-events: auto;
    }
    .settings-dialog {
      background: #fff; border-radius: 10px;
      width: min(760px, 94vw); max-height: 90vh;
      display: flex; flex-direction: column;
      box-shadow: 0 8px 30px rgba(0,0,0,0.35);
    }
    .settings-head {
      display: flex; justify-content: space-between; align-items: center;
      padding: 12px 18px; border-bottom: 1px solid #ddd;
    }
    .settings-head h4 { margin: 0; }
    .settings-close {
      border: 0; background: transparent; font-size: 26px; line-height: 1;
      color: #666; cursor: pointer; padding: 0 4px;
    }
    .settings-close:hover { color: #000; }
    .settings-body {
      overflow-y: auto; padding: 14px 18px 6px 18px;
      display: grid;
      grid-template-columns: repeat(2, minmax(240px, 1fr));
      gap: 0 32px;
    }
    @media (max-width: 760px) {
      .settings-body { grid-template-columns: 1fr; }
    }
    .settings-body h5 {
      margin: 4px 0 8px 0; font-weight: 600;
      border-bottom: 1px solid #eee; padding-bottom: 4px;
    }
    .settings-foot {
      padding: 10px 18px; border-top: 1px solid #ddd; text-align: right;
    }
    .panel-head {
      display: flex; justify-content: space-between; align-items: center;
      gap: 8px;
    }
    .panel-head h4 { margin: 0; }
  "))),
  
  tags$script(HTML("
    function openSettings() {
      document.getElementById('settings_modal').classList.add('open');
      setTimeout(function() {
        window.dispatchEvent(new Event('resize'));
      }, 60);
    }
    function closeSettings() {
      document.getElementById('settings_modal').classList.remove('open');
    }
    document.addEventListener('keydown', function(e) {
      if (e.key === 'Escape') closeSettings();
    });
  ")),
  
  leafletOutput("map", width = "100%", height = "100vh"),
  
  absolutePanel(
    top = 82, left = 10, width = "auto", fixed = TRUE,
    actionButton("reset_view", "Reset view", icon = icon("crosshairs"))
  ),
  
  div(
    class = "bottom-bar",
    sliderInput(
      "step_idx", "Simulation step",
      min = 1, max = n_steps, value = 1, step = 1,
      width = "100%",
      animate = animationOptions(interval = 800, loop = TRUE)
    )
  ),
  
  absolutePanel(
    class = "panel-box",
    top = 12, right = 12, width = 320, draggable = TRUE,
    style = "max-height: calc(100vh - 130px); overflow-y: auto;",
    
    div(
      class = "panel-head",
      h4("Phytoplankton patch"),
      tags$button(
        type = "button", class = "btn btn-default btn-sm",
        onclick = "openSettings()",
        icon("cog"), " Settings"
      )
    ),
    hr(style = "margin: 8px 0;"),
    uiOutput("info"),
    
    hr(style = "margin: 8px 0;"),
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
    verbatimTextOutput("dist_stats")
  ),
  
  div(
    id = "settings_modal", class = "settings-modal",
    onclick = "if (event.target === this) closeSettings();",
    
    div(
      class = "settings-dialog",
      
      div(
        class = "settings-head",
        h4("Settings"),
        tags$button(
          type = "button", class = "settings-close",
          onclick = "closeSettings()", `aria-label` = "Close",
          HTML("&times;")
        )
      ),
      
      div(
        class = "settings-body",
        
        # left column: patch
        div(
          h5("Patch"),
          checkboxInput("show_heat", "Radar-style density field", TRUE),
          sliderInput("heat_sigma", "Smoothing (grid boxes)",
                      min = 1, max = 15, value = 6, step = 0.5),
          sliderInput("heat_min", "Hide below (fraction of max)",
                      min = 0, max = 0.5, value = 0.03, step = 0.01),
          sliderInput("heat_opacity", "Field opacity",
                      min = 0.2, max = 1, value = 0.85, step = 0.05),
          sliderInput("repel_px", "Spot repel distance (px, 0 = off)",
                      min = 0, max = 80, value = 30, step = 5),
          checkboxInput("show_cells", "Show numbers / hover info", TRUE),
          checkboxInput("show_depth_labels", "Write the depth on each spot", TRUE),
          checkboxInput("show_track", "Centroid track", TRUE),
          sliderInput("trail_len", "Cell trails (steps)",
                      min = 0, max = 10, value = 3, step = 1)
        ),
        
        # right column: currents
        div(
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
        
      ),
      
      div(
        class = "settings-foot",
        tags$button(
          type = "button", class = "btn btn-primary",
          onclick = "closeSettings()", "Close"
        )
      )
    )
  )
)

############################################################
# 5. SERVER
############################################################

server <- function(input, output, session) {
  
  cur_step  <- reactive(steps[input$step_idx])
  step_data <- reactive(sim[sim$simulation_step == cur_step(), ])
  
  # maximum of the smoothed density over ALL steps, so colours are
  # comparable while the animation plays
  field_max <- reactive({
    sg <- input$heat_sigma
    max(1e-9, vapply(count_list, function(C) max(smooth_grid(C, sg)), numeric(1)))
  })
  
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
      lines <- c(lines, sprintf("Depth: mean %.1f m (range %.1f to %.1f m)",
                                mean(d$z), min(d$z), max(d$z)))
    }
    
    div(class = "info-line", HTML(paste(lines, collapse = "<br>")))
  })
  
  ##########################################################
  # BASE MAP
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
  # FLOW LAYER
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
    
    # arrow direction = drift direction (degrees clockwise from north);
    # at the first step it points towards the next centroid position
    ci <- which(cent$step == last$step)
    if (ci > 1) {
      p0 <- cent[ci - 1, ]; p1 <- cent[ci, ]
    } else if (nrow(cent) > 1) {
      p0 <- cent[1, ]; p1 <- cent[2, ]
    } else {
      p0 <- p1 <- cent[1, ]
    }
    angle <- atan2((p1$x - p0$x) * cos(p1$y * pi / 180), p1$y - p0$y) * 180 / pi
    
    # SVG arrow with the rotation baked in, passed as an image (data URI)
    arrow_svg <- sprintf(
      paste0("<svg xmlns='http://www.w3.org/2000/svg' width='44' height='44' ",
             "viewBox='-22 -22 44 44'>",
             "<g transform='rotate(%.1f)'>",
             "<polygon points='0,-18 12,14 0,7 -12,14' fill='#e63946' ",
             "stroke='#ffffff' stroke-width='2' stroke-linejoin='round'/>",
             "</g></svg>"),
      angle)
    arrow_icon <- makeIcon(
      iconUrl = paste0("data:image/svg+xml;charset=utf-8,",
                       utils::URLencode(arrow_svg, reserved = TRUE)),
      iconWidth = 44, iconHeight = 44,
      iconAnchorX = 22, iconAnchorY = 22
    )
    
    proxy |>
      addCircleMarkers(
        lng = tr$x, lat = tr$y,
        radius = 3, weight = 1, color = "#ffffff",
        fillColor = "#000000", fillOpacity = 1,
        group = "Track"
      ) |>
      addMarkers(
        lng = last$x, lat = last$y,
        icon = arrow_icon,
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
    
    now_ids <- step_data()$id
    d <- d[d$id %in% now_ids, ]
    if (nrow(d) == 0) return()
    
    ids <- unique(d$id)
    if (length(ids) > 300) ids <- sample(ids, 300)
    d <- d[d$id %in% ids, ]
    d <- d[order(d$id, d$simulation_step), ]
    
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
  # CELLS: HEAT MAP + NUMBERS
  #
  # Cells closer than ~MERGE_PX screen pixels are merged into one spot
  # (as before). The heat map is weighted by the number of cells in each
  # spot, so overlapping cells appear as hotter areas. The count is still
  # written on every merged spot, and an invisible marker on top keeps the
  # hover label and click popup.
  ##########################################################
  
  MERGE_PX <- 8
  
  observe({
    proxy <- leafletProxy("map") |>
      clearGroup("Heat") |>
      clearGroup("Cells") |>
      removeControl("heat_legend") |>
      removeControl("depth_legend")
    
    d <- step_data()
    if (nrow(d) == 0) return()
    if (!isTRUE(input$show_heat) && !isTRUE(input$show_cells)) return()
    
    zoom <- input$map_zoom
    if (is.null(zoom)) zoom <- 6
    
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
    
    # ---- repel the spots (numbers / depth labels / hover markers) ----
    if (isTRUE(input$repel_px > 0) && nrow(m) > 1) {
      ky <- 1 / cos(mean(m$y) * pi / 180)     # Mercator: pixels per deg lat
      rp <- repel_xy(m$x / deg_per_px, m$y * ky / deg_per_px, input$repel_px)
      m$x <- rp$x * deg_per_px
      m$y <- rp$y * deg_per_px / ky
    }
    
    # ---- radar-style density field (same colour scale in every step) ----
    if (isTRUE(input$show_heat)) {
      S    <- smooth_grid(count_list[[input$step_idx]], input$heat_sigma)
      fmax <- field_max()
      thr  <- input$heat_min * fmax
      hit  <- which(S > thr, arr.ind = TRUE)
      
      if (nrow(hit) > 0) {
        pal  <- colorNumeric(radar_cols, domain = c(0, fmax), na.color = "transparent")
        vals <- S[hit]
        proxy <- proxy |>
          addRectangles(
            lng1 = rad_xlo + (hit[, 2] - 1) * rdx,
            lng2 = rad_xlo + hit[, 2] * rdx,
            lat1 = rad_ylo + (hit[, 1] - 1) * rdy,
            lat2 = rad_ylo + hit[, 1] * rdy,
            stroke = FALSE, weight = 0,
            fillColor = pal(pmin(vals, fmax)),
            fillOpacity = input$heat_opacity,
            group = "Heat"
          ) |>
          addLegend(
            position = "bottomleft",
            pal = pal, values = c(0, fmax),
            title = "Cell density",
            labFormat = labelFormat(digits = 1),
            opacity = 0.95,
            layerId = "heat_legend"
          )
      }
    }
    
    if (!isTRUE(input$show_cells)) return()
    
    ids_txt <- vapply(split(d$id, grp), function(i) {
      s <- paste(head(i, 8), collapse = ", ")
      if (length(i) > 8) s <- paste0(s, ", ...")
      s
    }, character(1))
    
    popup_txt <- sprintf(
      "<b>%d cell%s here</b><br>ids: %s<br>Mean depth: %.1f m",
      m$n, ifelse(m$n == 1, "", "s"), ids_txt, m$z
    )
    
    hover_txt <- ifelse(
      m$n == 1,
      sprintf("Depth %.1f m", m$z),
      sprintf("%d cells, mean depth %.1f m", m$n, m$z)
    )
    
    # markers coloured by depth on top of the heat map (also carry the
    # hover label and click popup)
    proxy <- proxy |>
      addCircleMarkers(
        lng = m$x, lat = m$y,
        radius = pmin(16, 6 + 2 * log2(m$n)),
        stroke = FALSE,
        fillColor = "#000000",
        fillOpacity = 0,
        popup = popup_txt,
        label = hover_txt,
        group = "Cells"
      )
    
    # count numbers on merged spots
    mm <- m[m$n > 1, , drop = FALSE]
    if (nrow(mm) > 0) {
      proxy <- proxy |>
        addLabelOnlyMarkers(
          lng = mm$x, lat = mm$y,
          label = as.character(mm$n),
          labelOptions = labelOptions(
            noHide = TRUE, direction = "center", textOnly = TRUE,
            style = list("font-weight" = "bold", "font-size" = "12px",
                         "color" = "#111111",
                         "text-shadow" = "0 0 3px #fff, 0 0 3px #fff, 0 0 5px #fff")
          ),
          group = "Cells"
        )
    }
    
    # optional: depth next to every spot
    if (isTRUE(input$show_depth_labels)) {
      proxy |>
        addLabelOnlyMarkers(
          lng = m$x, lat = m$y,
          label = sprintf("%.1f m", m$z),
          labelOptions = labelOptions(
            noHide = TRUE, direction = "bottom", textOnly = TRUE,
            offset = c(0, 10),
            style = list("font-weight" = "600", "font-size" = "10.5px",
                         "color" = "#111111",
                         "text-shadow" = "0 0 3px #fff, 0 0 3px #fff, 0 0 5px #fff")
          ),
          group = "Cells"
        )
    }
  })
  
  ##########################################################
  # RESET VIEW
  ##########################################################
  
  observeEvent(input$reset_view, {
    leafletProxy("map") |>
      fitBounds(r_xlo, r_ylo, r_xhi, r_yhi)
  })
  
  ##########################################################
  # DISTRIBUTIONS
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