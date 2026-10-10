############################################################
# PHYTOPLANKTON PATCH VIEWER (Shiny + Leaflet + leaflet-velocity)
#
#  - Everything comes from phytoplankton_simulation_results.csv
#  - Flow : the u / v columns of the table (current at each cell)
#           interpolated to a grid per step (IDW), because the table
#           only holds u/v where cells are. The flow fades to zero
#           away from the cells.
#  - Patch: centroid track (+ drift speed), fading cell trails, points
#  - Particle set dialog: writes particle_params.csv for the model
#    script; "Reload results" re-reads the csv after a new model run
############################################################

# install.packages(c("shiny", "leaflet", "leaflet.extras2", "jsonlite"))

library(shiny)
library(leaflet)
library(leaflet.extras2)
library(jsonlite)

############################################################
# 1. FILES AND SETTINGS
############################################################

SIM_FILE         <- "phytoplankton_simulation_results.csv"
PARAMS_FILE      <- "particle_params.csv"                      # written by this app
PARAMS_USED_FILE <- "phytoplankton_particle_params_used.csv"   # written by the model

FLOW_PAD   <- 2      # degrees of padding around all positions in the table
D_PAD      <- 2      # degrees of padding for the starting map view
nx         <- 70     # flow grid size
ny         <- 55
MAX_POINTS <- 1500   # cap for speed if the population grows large

############################################################
# 2. FLOW FIELD FROM THE TABLE (u, v of the cells)
############################################################

# The table only has u/v at the cell positions, so a continuous field is
# built per step by inverse-distance weighting (IDW) onto a regular grid.
# The field fades to zero with distance from the nearest cell
# ("influence radius"), so it is only meaningful near the patch.
# (grid_x, grid_y, x_min ... are created by load_results() below.)

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
# 3. PATCH HELPERS
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
# 4. PARTICLE-SET PARAMETERS (shared with the model script)
############################################################

# Same names and defaults as particle_params in phytoplankton_model.R.
# lon_/lat_ min and max = start area; NA = the whole model grid.
PP_DEFAULTS <- list(
  seed         = 101,
  n_particles  = 10,
  start_depth  = 10,
  diameter_min = 3,
  diameter_max = 70,
  rho_min      = 1000,
  rho_max      = 1023,
  age_min      = 0,
  age_max      = 24,
  carbon_min   = 50,
  carbon_max   = 200,
  nitrogen_min = 5,
  nitrogen_max = 30,
  lon_min      = NA_real_,
  lon_max      = NA_real_,
  lat_min      = NA_real_,
  lat_max      = NA_real_
)

# parameters from a name,value csv (missing names keep their defaults)
read_params <- function(file) {
  p <- PP_DEFAULTS
  if (file.exists(file)) {
    df <- tryCatch(read.csv(file, stringsAsFactors = FALSE),
                   error = function(e) NULL)
    if (!is.null(df) && all(c("name", "value") %in% names(df))) {
      for (i in seq_len(nrow(df))) {
        nm <- as.character(df$name[i])
        if (nm %in% names(p)) {
          p[[nm]] <- suppressWarnings(as.numeric(df$value[i]))
        }
      }
    }
  }
  p
}

write_params <- function(p, file) {
  write.csv(
    data.frame(name = names(p), value = as.numeric(unlist(p))),
    file, row.names = FALSE
  )
}

# parameters the last model run used (NULL if the model did not record them)
read_used_params <- function() {
  if (file.exists(PARAMS_USED_FILE)) read_params(PARAMS_USED_FILE) else NULL
}

params_equal <- function(a, b) {
  nm <- names(PP_DEFAULTS)
  isTRUE(all.equal(as.numeric(unlist(a[nm])), as.numeric(unlist(b[nm]))))
}

# returns a character vector of problems (empty = fine)
validate_params <- function(p) {
  msg <- character(0)
  
  required <- setdiff(names(PP_DEFAULTS),
                      c("lon_min", "lon_max", "lat_min", "lat_max"))
  missing_vals <- required[!vapply(required, function(n) {
    length(p[[n]]) == 1 && is.finite(p[[n]])
  }, logical(1))]
  
  if (length(missing_vals) > 0) {
    return(paste("Missing value for:", paste(missing_vals, collapse = ", ")))
  }
  
  if (p$n_particles < 1 || p$n_particles != round(p$n_particles)) {
    msg <- c(msg, "Number of particles must be a whole number, 1 or more.")
  }
  if (p$start_depth < 0) msg <- c(msg, "Start depth cannot be negative.")
  if (p$diameter_min <= 0) msg <- c(msg, "Diameter must be above 0.")
  if (p$rho_min <= 0) msg <- c(msg, "Density must be above 0.")
  if (p$age_min < 0) msg <- c(msg, "Age cannot be negative.")
  
  for (b in c("diameter", "rho", "age", "carbon", "nitrogen")) {
    if (p[[paste0(b, "_min")]] > p[[paste0(b, "_max")]]) {
      msg <- c(msg, sprintf("%s: min is larger than max.", b))
    }
  }
  
  for (a in c("lon", "lat")) {
    lo <- p[[paste0(a, "_min")]]
    hi <- p[[paste0(a, "_max")]]
    if (is.finite(lo) != is.finite(hi)) {
      msg <- c(msg, sprintf("%s: fill both min and max, or leave both empty.", a))
    } else if (is.finite(lo) && lo > hi) {
      msg <- c(msg, sprintf("%s: min is larger than max.", a))
    }
  }
  
  msg
}

############################################################
# 5. LOAD THE RESULTS
#    Called once at start and again by the "Reload results" button.
#    Everything that depends on the table is rebuilt and stored in the
#    app environment, where the UI and the server read it.
############################################################

load_results <- function() {
  
  if (!file.exists(SIM_FILE)) {
    stop("Could not find ", SIM_FILE, " in the app folder.")
  }
  
  sim <- read.csv(SIM_FILE)
  
  needed <- c("simulation_step", "elapsed_hours", "id", "x", "y", "z", "u", "v")
  missing_cols <- setdiff(needed, names(sim))
  if (length(missing_cols) > 0) {
    stop("Table is missing columns: ", paste(missing_cols, collapse = ", "),
         "\nRe-run the model with the u/v edit to create u and v.")
  }
  
  sim <- sim[is.finite(sim$x) & is.finite(sim$y) &
               is.finite(sim$u) & is.finite(sim$v), ]
  if (nrow(sim) == 0) {
    stop("The simulation table has no cells with positions and u/v.")
  }
  
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
  
  # flow grid (covers all positions in the table, padded)
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
  
  # starting map view
  d_xlo <- min(sim$x) - D_PAD
  d_xhi <- max(sim$x) + D_PAD
  d_ylo <- min(sim$y) - D_PAD
  d_yhi <- max(sim$y) + D_PAD
  
  # tight frame around the data (all steps) for the "Reset view" button
  r_padx <- max(0.1, 0.1 * (max(sim$x) - min(sim$x)))
  r_pady <- max(0.1, 0.1 * (max(sim$y) - min(sim$y)))
  r_xlo <- min(sim$x) - r_padx
  r_xhi <- max(sim$x) + r_padx
  r_ylo <- min(sim$y) - r_pady
  r_yhi <- max(sim$y) + r_pady
  
  # depth colours (same scale in every step): shallow = yellow, deep = purple
  z_dom <- range(sim$z[is.finite(sim$z)])
  if (diff(z_dom) == 0) z_dom <- z_dom + c(-0.5, 0.5)
  
  depth_pal <- colorBin(
    palette = "viridis", domain = z_dom, bins = 6, pretty = TRUE,
    reverse = TRUE
  )
  
  # variables for the distribution plots
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
  
  # parameters the model used for this table (NULL if not recorded)
  used_params <- read_used_params()
  
  # store in the environment where this function was defined (the app's)
  keep <- c("sim", "steps", "n_steps", "cent",
            "x_min", "x_max", "y_min", "y_max", "dx", "dy", "grid_x", "grid_y",
            "d_xlo", "d_xhi", "d_ylo", "d_yhi",
            "r_xlo", "r_xhi", "r_ylo", "r_yhi",
            "z_dom", "depth_pal", "var_info", "dist_choices", "used_params")
  list2env(mget(keep, envir = environment()), envir = parent.env(environment()))
  
  invisible(NULL)
}

load_results()

############################################################
# 5b. PARTICLE-SET DIALOG (built per page load, so it shows the saved values)
############################################################

particles_modal <- function() {
  
  p <- read_params(PARAMS_FILE)
  
  num <- function(name, label, ...) {
    numericInput(paste0("pp_", name), label, value = p[[name]],
                 width = "100%", ...)
  }
  pair <- function(base, label_min, label_max, ...) {
    div(
      class = "pp-pair",
      num(paste0(base, "_min"), label_min, ...),
      num(paste0(base, "_max"), label_max, ...)
    )
  }
  
  div(
    id = "particles_modal", class = "settings-modal",
    onclick = "if (event.target === this) closeParticles();",
    
    div(
      class = "settings-dialog",
      
      div(
        class = "settings-head",
        h4("Particle set"),
        tags$button(
          type = "button", class = "settings-close",
          onclick = "closeParticles()", `aria-label` = "Close",
          HTML("&times;")
        )
      ),
      
      div(
        class = "settings-body",
        
        div(
          h5("Population"),
          div(
            class = "pp-pair",
            num("n_particles", "Number of particles", min = 1, step = 1),
            num("seed", "Random seed", step = 1)
          ),
          div(
            class = "pp-pair",
            num("start_depth", "Start depth (m)", min = 0)
          ),
          pair("age", "Age min (h)", "Age max (h)", min = 0),
          
          h5("Start area (empty = whole model grid)"),
          pair("lon", "Longitude min", "Longitude max"),
          pair("lat", "Latitude min", "Latitude max")
        ),
        
        div(
          h5("Cell properties (random between min and max)"),
          pair("diameter", "Diameter min (\u00b5m)", "Diameter max (\u00b5m)", min = 0),
          pair("rho", "Density min (kg/m\u00b3)", "Density max (kg/m\u00b3)", min = 0),
          pair("carbon", "Cell carbon min", "Cell carbon max", min = 0),
          pair("nitrogen", "Cell nitrogen min", "Cell nitrogen max", min = 0)
        )
      ),
      
      div(
        class = "pp-help",
        HTML(paste0(
          "1. <b>Save parameters</b> (writes <code>particle_params.csv</code>)",
          " &rarr; 2. re-run <code>phytoplankton_model.R</code> in the same folder",
          " &rarr; 3. <b>Reload results</b>."
        ))
      ),
      div(class = "pp-status", uiOutput("pp_status")),
      
      div(
        class = "settings-foot split",
        actionButton("pp_defaults", "Defaults"),
        div(
          actionButton("pp_save", "Save parameters", class = "btn-primary"),
          actionButton("reload_results", "Reload results"),
          tags$button(
            type = "button", class = "btn btn-default",
            onclick = "closeParticles()", "Close"
          )
        )
      )
    )
  )
}

############################################################
# 6. UI
############################################################

ui <- function(request) fluidPage(
  
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

    /* simulation-step slider: transparent, medium length, bottom centre */
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
    /* no background, so keep the text readable over the map */
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

    /* the bar is centred, so the map's bottom controls only need to move
       up on narrow screens where they would sit under it */
    .leaflet-bottom { bottom: 0; }
    @media (max-width: 1150px) { .leaflet-bottom { bottom: 84px; } }

    /* ---- settings modal ----
       Always present in the page (so every input has a value from the
       start); only shown/hidden. visibility (not display:none) keeps the
       layout, so sliders and the plot measure their width correctly. */
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
    .panel-box h4 { margin: 0; }
    .panel-buttons { display: flex; gap: 6px; margin-top: 8px; }
    .settings-foot.split {
      display: flex; justify-content: space-between; gap: 8px; flex-wrap: wrap;
    }
    .settings-foot .btn + .btn { margin-left: 4px; }
    .pp-pair { display: grid; grid-template-columns: 1fr 1fr; gap: 0 12px; }
    .pp-help { padding: 8px 18px 0 18px; font-size: 12.5px; color: #555; }
    .pp-status { padding: 6px 18px 10px 18px; font-size: 12.5px; }
  "))),

  tags$script(HTML("
    function openModal(id) {
      document.getElementById(id).classList.add('open');
      // let sliders / the plot re-measure once the dialog is visible
      setTimeout(function() {
        window.dispatchEvent(new Event('resize'));
      }, 60);
    }
    function closeModal(id) {
      document.getElementById(id).classList.remove('open');
    }
    function openSettings()  { openModal('settings_modal'); }
    function closeSettings() { closeModal('settings_modal'); }
    function openParticles()  { openModal('particles_modal'); }
    function closeParticles() { closeModal('particles_modal'); }
    document.addEventListener('keydown', function(e) {
      if (e.key === 'Escape') {
        document.querySelectorAll('.settings-modal.open').forEach(function(m) {
          m.classList.remove('open');
        });
      }
    });
  ")),

  leafletOutput("map", width = "100%", height = "100vh"),

  absolutePanel(
    top = 82, left = 10, width = "auto", fixed = TRUE,
    actionButton("reset_view", "Reset view", icon = icon("crosshairs"))
  ),

  # Simulation step slider + play button: bottom centre, medium length,
  # transparent background
  div(
    class = "bottom-bar",
    sliderInput(
      "step_idx", "Simulation step",
      min = 1, max = n_steps, value = 1, step = 1,
      width = "100%",
      animate = animationOptions(interval = 2000, loop = TRUE)
    )
  ),

  # Compact info panel: title, Settings button and the step read-out
  absolutePanel(
    class = "panel-box",
    top = 12, right = 12, width = 320, draggable = TRUE,
    style = "max-height: calc(100vh - 130px); overflow-y: auto;",

    h4("Phytoplankton patch"),
    div(
      class = "panel-buttons",
      tags$button(
        type = "button", class = "btn btn-default btn-sm",
        onclick = "openSettings()",
        icon("cog"), " Settings"
      ),
      tags$button(
        type = "button", class = "btn btn-default btn-sm",
        onclick = "openParticles()",
        " Particle set"
      )
    ),
    hr(style = "margin: 8px 0;"),
    uiOutput("info"),

    # Distributions: always visible (not part of the settings modal)
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

  # All settings in one modal dialog (click outside, x, Close or Esc to close)
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
          checkboxInput("show_track", "Centroid track", TRUE),
          sliderInput("trail_len", "Cell trails (steps)",
                      min = 0, max = 10, value = 3, step = 1),
          checkboxInput("show_cells", "Show points (individual cells)", TRUE),
          checkboxInput("show_depth_labels", "Write the depth on each point", FALSE)
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
        ),

      ),

      div(
        class = "settings-foot",
        tags$button(
          type = "button", class = "btn btn-primary",
          onclick = "closeSettings()", "Close"
        )
      )
    )
  ),

  # Particle-set dialog: initial particle properties for the model run
  particles_modal()
)

############################################################
# 7. SERVER
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
      lines <- c(lines, sprintf("Depth: mean %.1f m (range %.1f to %.1f m)",
                                mean(d$z), min(d$z), max(d$z)))
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
    proxy <- leafletProxy("map") |>
      clearGroup("Cells") |>
      removeControl("depth_legend")

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

    # hover text: depth of the cell (mean depth for merged markers)
    hover_txt <- ifelse(
      m$n == 1,
      sprintf("Depth %.1f m", m$z),
      sprintf("%d cells, mean depth %.1f m", m$n, m$z)
    )

    # markers coloured by depth
    proxy <- proxy |>
      addCircleMarkers(
        lng = m$x, lat = m$y,
        radius = pmin(22, 5 + 3 * log2(m$n)),
        weight = 1.5,
        color = "#ffffff",
        fillColor = depth_pal(m$z),
        fillOpacity = 0.95,
        popup = popup_txt,
        label = hover_txt,
        group = "Cells"
      )

    # count labels only on merged markers
    mm <- m[m$n > 1, , drop = FALSE]
    if (nrow(mm) > 0) {
      proxy <- proxy |>
        addLabelOnlyMarkers(
          lng = mm$x, lat = mm$y,
          label = as.character(mm$n),
          labelOptions = labelOptions(
            noHide = TRUE, direction = "center", textOnly = TRUE,
            style = list("font-weight" = "bold", "font-size" = "11px",
                         "color" = "#111111",
                         "text-shadow" = "0 0 3px #fff, 0 0 3px #fff")
          ),
          group = "Cells"
        )
    }

    # optional: write the depth next to every marker
    if (isTRUE(input$show_depth_labels)) {
      proxy <- proxy |>
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

    # depth legend (shallow at the top)
    proxy |>
      addLegend(
        position = "bottomleft",
        pal = depth_pal, values = z_dom,
        title = "Depth (m)",
        labFormat = labelFormat(digits = 1),
        opacity = 0.95,
        layerId = "depth_legend"
      )
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

  ##########################################################
  # PARTICLE-SET DIALOG
  ##########################################################

  # parameters saved in particle_params.csv (NULL = no file yet)
  saved_params <- reactiveVal(
    if (file.exists(PARAMS_FILE)) read_params(PARAMS_FILE) else NULL
  )

  # what is typed in the dialog right now
  current_params <- reactive({
    nms <- names(PP_DEFAULTS)
    p <- lapply(nms, function(nm) {
      v <- input[[paste0("pp_", nm)]]
      if (is.null(v) || length(v) == 0) NA_real_ else as.numeric(v)
    })
    names(p) <- nms
    p
  })

  output$pp_status <- renderUI({
    p <- current_params()
    problems <- validate_params(p)

    if (length(problems) > 0) {
      return(div(style = "color: #b00020;",
                 HTML(paste(htmltools::htmlEscape(problems), collapse = "<br>"))))
    }

    if (is.null(used_params)) {
      return(div(style = "color: #666;",
                 "The loaded results have no record of the parameters used."))
    }

    if (params_equal(p, used_params)) {
      return(div(style = "color: #1b7f3b;",
                 "These are the parameters of the loaded results."))
    }

    sv <- saved_params()
    if (!is.null(sv) && params_equal(p, sv)) {
      div(style = "color: #9a6700;",
          "Saved, but different from the loaded results: re-run the model, ",
          "then press Reload results.")
    } else {
      div(style = "color: #9a6700;",
          "Different from the loaded results and not saved yet: ",
          "press Save parameters, re-run the model, then Reload results.")
    }
  })

  observeEvent(input$pp_save, {
    p <- current_params()
    problems <- validate_params(p)

    if (length(problems) > 0) {
      showNotification(paste(problems, collapse = " "),
                       type = "error", duration = 8)
      return()
    }

    ok <- tryCatch({
      write_params(p, PARAMS_FILE)
      TRUE
    }, error = function(e) {
      showNotification(paste("Could not save:", conditionMessage(e)),
                       type = "error", duration = 8)
      FALSE
    })

    if (ok) {
      saved_params(p)
      showNotification(
        paste0("Saved ", PARAMS_FILE, ". Re-run the model, then press Reload results."),
        type = "message", duration = 8
      )
    }
  })

  observeEvent(input$pp_defaults, {
    for (nm in names(PP_DEFAULTS)) {
      updateNumericInput(session, paste0("pp_", nm), value = PP_DEFAULTS[[nm]])
    }
  })

  # re-read the csv written by a new model run, then rebuild the page
  observeEvent(input$reload_results, {
    ok <- tryCatch({
      load_results()
      TRUE
    }, error = function(e) {
      showNotification(paste("Could not reload the results:", conditionMessage(e)),
                       type = "error", duration = 10)
      FALSE
    })

    if (ok) session$reload()
  })
}

shinyApp(ui, server)
