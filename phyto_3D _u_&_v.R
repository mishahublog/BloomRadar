############################################################
# PHYTOPLANKTON 3D LAGRANGIAN MODEL
# HYCOM PHYSICS + COPERNICUS BGC NUTRIENTS
############################################################

library(terra)

############################################################
# 1. LOAD PHYSICAL DATA
############################################################

sal_data  <- rast("sal_data.nc")
temp_data <- rast("temp_data.nc")
u_data    <- rast("u_data.nc")
v_data    <- rast("v_data.nc")

############################################################
# 2. GRID INFORMATION
############################################################

grid_dims <- dim(sal_data)

n_rows <- grid_dims[1]
n_cols <- grid_dims[2]

grid_ext <- ext(sal_data)

x_min <- xmin(grid_ext)
x_max <- xmax(grid_ext)

y_min <- ymin(grid_ext)
y_max <- ymax(grid_ext)

# Raster cell resolution
x_res <- (x_max - x_min) / (n_cols - 1)
y_res <- (y_max - y_min) / (n_rows - 1)

cat("Grid dimensions:", n_rows, "x", n_cols, "\n")
cat("Longitude:", x_min, "-", x_max, "\n")
cat("Latitude :", y_min, "-", y_max, "\n")

############################################################
# 3. DEPTH LEVELS
############################################################

get_depths <- function(r) {
  
  d <- tryCatch(
    unique(r@pntr$depth),
    error = function(e) NULL
  )
  
  if (is.null(d) || length(d) == 0 || all(is.na(d))) {
    stop("Could not determine depth levels.")
  } 
  sort(unique(as.numeric(d)))
}

sal_depth  <- get_depths(sal_data)
temp_depth <- get_depths(temp_data)
u_depth    <- get_depths(u_data)
v_depth    <- get_depths(v_data)

depth_vectors <- list(
  sal  = sal_depth,
  temp = temp_depth,
  u    = u_depth,
  v    = v_depth
)

############################################################
# COMMON DEPTH LEVELS
############################################################

real_depths <- Reduce( intersect,depth_vectors)
real_depths <- sort(real_depths)
n_depths <- length(real_depths)

if (n_depths < 2) {
  stop("At least two common depth levels are required.")
}

cat("Common depth levels:",n_depths,"\n")

cat( "Depth range:",min(real_depths), "-", max(real_depths),"m\n")

############################################################
# 4. DETERMINE COMMON TIME STEPS
############################################################

get_time_layers <- function(r, n_depths) {
  
  n_layers <- nlyr(r)
  
  if (n_layers %% n_depths != 0) {
    
    warning("Layers (",n_layers,") are not exactly divisible by depth levels (",n_depths,").")
  }
  
  floor(n_layers / n_depths)
}

n_time_sal  <- get_time_layers(sal_data, n_depths)
n_time_temp <- get_time_layers(temp_data, n_depths)
n_time_u    <- get_time_layers(u_data, n_depths)
n_time_v    <- get_time_layers(v_data, n_depths)

n_times <- min(
  n_time_sal,
  n_time_temp,
  n_time_u,
  n_time_v
)

cat(
  "Common time steps:",
  n_times,
  "\n"
)

############################################################
# 5. LAYER NAMES
############################################################

sal_names  <- names(sal_data)[1:(n_times * n_depths)]
temp_names <- names(temp_data)[1:(n_times * n_depths)]
u_names    <- names(u_data)[1:(n_times * n_depths)]
v_names    <- names(v_data)[1:(n_times * n_depths)]

############################################################
# 6. CACHE 4D DATA
############################################################

cache_4d_variable <- function(
    stack,
    layer_names,
    n_depths,
    n_times,
    n_rows,
    n_cols
) {
  
  arr_4d <- array(
    NA_real_,
    dim = c(
      n_rows,
      n_cols,
      n_depths,
      n_times
    )
  )
  
  for (t in seq_len(n_times)) {
    
    start_layer <-
      (t - 1) * n_depths + 1
    
    end_layer <-
      t * n_depths
    
    layers <-
      layer_names[start_layer:end_layer]
    
    temp_stack <-
      terra::subset(
        stack,
        layers
      )
    
    arr_4d[, , , t] <-
      as.array(temp_stack)
  }
  
  arr_4d
}

cat("Caching physical data...\n")

master_u <- cache_4d_variable(
  u_data,
  u_names,
  n_depths,
  n_times,
  n_rows,
  n_cols
)

master_v <- cache_4d_variable(
  v_data,
  v_names,
  n_depths,
  n_times,
  n_rows,
  n_cols
)

master_temp <- cache_4d_variable(
  temp_data,
  temp_names,
  n_depths,
  n_times,
  n_rows,
  n_cols
)

master_sal <- cache_4d_variable(
  sal_data,
  sal_names,
  n_depths,
  n_times,
  n_rows,
  n_cols
)

cat("Physical data cached successfully.\n")

############################################################
# 7. LOAD COPERNICUS BIOGEOCHEMICAL DATA
############################################################

# IMPORTANT:
# Replace this with your actual Copernicus raster object or file.
#
# Example:
#
# copernicus_file <- rast(
#   "copernicus_bgc.nc"
# )
#
# Or if you already have an object:
#
# copernicus_file <- my_cms_data

copernicus_file <- my_cms_data

############################################################
# CONVERT TO SPATRASTER
############################################################

bgc_rast <- tryCatch(
  {
    
    if (
      inherits(
        copernicus_file,
        "SpatRaster"
      )
    ) {
      
      copernicus_file
      
    } else {
      
      rast(copernicus_file)
    }
    
  },
  error = function(e) {
    
    warning(
      "Could not load Copernicus BGC data: ",
      e$message
    )
    
    NULL
  }
)

############################################################
# 8. NUTRIENT ARRAYS
############################################################

if (!is.null(bgc_rast)) {
  
  cat(
    "Loading Copernicus BGC nutrients...\n"
  )
  
  ##########################################################
  # PRINT LAYER NAMES
  ##########################################################
  
  bgc_names <- names(bgc_rast)
  
  print(bgc_names)
  
  ##########################################################
  # FIND LAYER INDICES
  ##########################################################
  
  find_layer_index <- function(
    raster,
    pattern
  ) {
    
    idx <- grep(
      pattern,
      names(raster),
      ignore.case = TRUE
    )
    
    if (length(idx) == 0) {
      
      warning(
        "Could not find layer matching: ",
        pattern
      )
      
      return(NA_integer_)
    }
    
    idx[1]
  }
  
  ##########################################################
  # IMPORTANT:
  # terra::subset() requires a layer index or logical vector.
  # Therefore we use numeric indices, not layer names.
  ##########################################################
  
  no3_idx <- find_layer_index(
    bgc_rast,
    "no3"
  )
  
  fe_idx <- find_layer_index(
    bgc_rast,
    "fe"
  )
  
  po4_idx <- find_layer_index(
    bgc_rast,
    "po4"
  )
  
  si_idx <- find_layer_index(
    bgc_rast,
    "si"
  )
  
  ##########################################################
  # CHECK REQUIRED LAYERS
  ##########################################################
  
  if (
    any(
      is.na(
        c(
          no3_idx,
          fe_idx,
          po4_idx,
          si_idx
        )
      )
    )
  ) {
    
    stop(
      "One or more nutrient layers could not be found.\n",
      "Available layers are:\n",
      paste(
        bgc_names,
        collapse = "\n"
      )
    )
  }
  
  ##########################################################
  # EXTRACT LAYERS
  ##########################################################
  
  layer_no3 <- terra::subset(
    bgc_rast,
    no3_idx
  )
  
  layer_fe <- terra::subset(
    bgc_rast,
    fe_idx
  )
  
  layer_po4 <- terra::subset(
    bgc_rast,
    po4_idx
  )
  
  layer_si <- terra::subset(
    bgc_rast,
    si_idx
  )
  
  ##########################################################
  # TARGET GRID
  ##########################################################
  
  target_grid <-
    terra::subset(
      sal_data,
      1
    )
  
  ##########################################################
  # RESAMPLE NUTRIENTS TO HYCOM GRID
  ##########################################################
  
  cat(
    "Resampling nutrients to HYCOM grid...\n"
  )
  
  res_no3 <- resample(
    layer_no3,
    target_grid,
    method = "bilinear"
  )
  
  res_fe <- resample(
    layer_fe,
    target_grid,
    method = "bilinear"
  )
  
  res_po4 <- resample(
    layer_po4,
    target_grid,
    method = "bilinear"
  )
  
  res_si <- resample(
    layer_si,
    target_grid,
    method = "bilinear"
  )
  
  ##########################################################
  # CONVERT TO MATRICES
  ##########################################################
  
  mat_no3_surf <-
    as.matrix(
      res_no3,
      wide = TRUE
    )
  
  mat_fe_surf <-
    as.matrix(
      res_fe,
      wide = TRUE
    )
  
  mat_po4_surf <-
    as.matrix(
      res_po4,
      wide = TRUE
    )
  
  mat_si_surf <-
    as.matrix(
      res_si,
      wide = TRUE
    )
  
  ##########################################################
  # REPLACE NON-FINITE VALUES
  ##########################################################
  
  mat_no3_surf[
    !is.finite(mat_no3_surf)
  ] <- 0
  
  mat_fe_surf[
    !is.finite(mat_fe_surf)
  ] <- 0
  
  mat_po4_surf[
    !is.finite(mat_po4_surf)
  ] <- 0
  
  mat_si_surf[
    !is.finite(mat_si_surf)
  ] <- 0
  
  ##########################################################
  # INITIALIZE 4D NUTRIENT ARRAYS
  ##########################################################
  
  master_no3 <- array(
    NA_real_,
    dim = c(
      n_rows,
      n_cols,
      n_depths,
      n_times
    )
  )
  
  master_nh4 <- array(
    NA_real_,
    dim = c(
      n_rows,
      n_cols,
      n_depths,
      n_times
    )
  )
  
  master_fe <- array(
    NA_real_,
    dim = c(
      n_rows,
      n_cols,
      n_depths,
      n_times
    )
  )
  
  master_po4 <- array(
    NA_real_,
    dim = c(
      n_rows,
      n_cols,
      n_depths,
      n_times
    )
  )
  
  master_si <- array(
    NA_real_,
    dim = c(
      n_rows,
      n_cols,
      n_depths,
      n_times
    )
  )
  
  ##########################################################
  # CREATE DEPTH PROFILES
  ##########################################################
  
  for (
    z_idx in seq_len(n_depths)
  ) {
    
    depth <-
      real_depths[z_idx]
    
    ########################################################
    # NUTRIENT DEPTH ATTENUATION
    ########################################################
    
    decay_factor <-
      exp(
        -0.02 *
          depth
      )
    
    ########################################################
    # NH4 MODEL
    ########################################################
    
    no3_profile <-
      mat_no3_surf *
      decay_factor
    
    nh4_profile <-
      mat_no3_surf *
      0.12 *
      decay_factor
    
    fe_profile <-
      mat_fe_surf *
      decay_factor
    
    po4_profile <-
      mat_po4_surf *
      decay_factor
    
    si_profile <-
      mat_si_surf *
      decay_factor
    
    ########################################################
    # REPEAT FOR ALL TIMES
    ########################################################
    
    for (
      t in seq_len(n_times)
    ) {
      
      master_no3[
        , ,
        z_idx,
        t
      ] <-
        no3_profile
      
      master_nh4[
        , ,
        z_idx,
        t
      ] <-
        nh4_profile
      
      master_fe[
        , ,
        z_idx,
        t
      ] <-
        fe_profile
      
      master_po4[
        , ,
        z_idx,
        t
      ] <-
        po4_profile
      
      master_si[
        , ,
        z_idx,
        t
      ] <-
        si_profile
    }
  }
  
  cat(
    "Nutrient 4D arrays created successfully.\n"
  )
  
} else {
  
  ##########################################################
  # FALLBACK NUTRIENTS
  ##########################################################
  
  warning(
    "Using idealized fallback nutrient concentrations."
  )
  
  master_no3 <-
    master_temp * 0 + 2.5
  
  master_nh4 <-
    master_temp * 0 + 0.3
  
  master_fe <-
    master_temp * 0 + 0.05
  
  master_po4 <-
    master_temp * 0 + 0.2
  
  master_si <-
    master_temp * 0 + 1.5
}

############################################################
# 9. PAR FIELD
############################################################

master_par <-
  master_temp * 0 + 45

############################################################
# 10. 4D INTERPOLATION
############################################################

extract_4d_slice_fast <- function(
    data_arr_4d,
    x,
    y,
    z,
    levels,
    time_idx
) {
  
  if (
    length(x) == 0
  ) {
    
    return(
      numeric(0)
    )
  }
  
  ##########################################################
  # BOUND TIME
  ##########################################################
  
  time_idx <-
    pmax(
      1,
      pmin(
        dim(data_arr_4d)[4],
        time_idx
      )
    )
  
  ##########################################################
  # HORIZONTAL COORDINATES
  ##########################################################
  
  col_f <-
    1 +
    (x - x_min) /
    x_res
  
  row_f <-
    1 +
    (y_max - y) /
    y_res
  
  col_f <-
    pmax(
      1,
      pmin(
        n_cols,
        col_f
      )
    )
  
  row_f <-
    pmax(
      1,
      pmin(
        n_rows,
        row_f
      )
    )
  
  c1 <-
    floor(col_f)
  
  r1 <-
    floor(row_f)
  
  c2 <-
    pmin(
      n_cols,
      c1 + 1
    )
  
  r2 <-
    pmin(
      n_rows,
      r1 + 1
    )
  
  wc2 <-
    col_f - c1
  
  wc1 <-
    1 - wc2
  
  wr2 <-
    row_f - r1
  
  wr1 <-
    1 - wr2
  
  ##########################################################
  # VERTICAL COORDINATES
  ##########################################################
  
  z <-
    pmax(
      min(levels),
      pmin(
        max(levels),
        z
      )
    )
  
  lower_idx <-
    findInterval(
      z,
      levels
    )
  
  lower_idx <-
    pmax(
      1,
      pmin(
        length(levels) - 1,
        lower_idx
      )
    )
  
  upper_idx <-
    lower_idx + 1
  
  deepest <-
    z >= max(levels)
  
  lower_idx[deepest] <-
    length(levels) - 1
  
  upper_idx[deepest] <-
    length(levels)
  
  z_low <-
    levels[lower_idx]
  
  z_high <-
    levels[upper_idx]
  
  dz <-
    z_high - z_low
  
  wz2 <-
    ifelse(
      dz == 0,
      0,
      (z - z_low) /
        dz
    )
  
  wz1 <-
    1 - wz2
  
  ##########################################################
  # DATA EXTRACTION
  ##########################################################
  
  get_value <- function(
    r,
    c,
    d
  ) {
    
    idx <-
      cbind(
        r,
        c,
        d,
        time_idx
      )
    
    data_arr_4d[idx]
  }
  
  v111 <-
    get_value(
      r1,
      c1,
      lower_idx
    )
  
  v121 <-
    get_value(
      r1,
      c2,
      lower_idx
    )
  
  v211 <-
    get_value(
      r2,
      c1,
      lower_idx
    )
  
  v221 <-
    get_value(
      r2,
      c2,
      lower_idx
    )
  
  v112 <-
    get_value(
      r1,
      c1,
      upper_idx
    )
  
  v122 <-
    get_value(
      r1,
      c2,
      upper_idx
    )
  
  v212 <-
    get_value(
      r2,
      c1,
      upper_idx
    )
  
  v222 <-
    get_value(
      r2,
      c2,
      upper_idx
    )
  
  ##########################################################
  # HORIZONTAL INTERPOLATION
  ##########################################################
  
  lower <-
    (
      v111 * wc1 +
        v121 * wc2
    ) *
    wr1 +
    
    (
      v211 * wc1 +
        v221 * wc2
    ) *
    wr2
  
  upper <-
    (
      v112 * wc1 +
        v122 * wc2
    ) *
    wr1 +
    
    (
      v212 * wc1 +
        v222 * wc2
    ) *
    wr2
  
  ##########################################################
  # VERTICAL INTERPOLATION
  ##########################################################
  
  result <-
    lower * wz1 +
    upper * wz2
  
  result[
    !is.finite(result)
  ] <- 0
  
  result
}

############################################################
# 11. ADVECTION
############################################################

kernel_advection_rk4_3d_fast <- function(
    particles,
    arr_u,
    arr_v,
    levels,
    dt,
    time_idx
) {
  
  if (
    nrow(particles) == 0
  ) {
    
    return(
      particles
    )
  }
  
  deg_lat_scale <-
    1 / 111000
  
  deg_lon_scale <-
    1 /
    (
      111000 *
        pmax(
          0.1,
          cos(
            particles$y *
              pi /
              180
          )
        )
    )
  
  ##########################################################
  # K1
  ##########################################################
  
  u1 <-
    extract_4d_slice_fast(
      arr_u,
      particles$x,
      particles$y,
      particles$z,
      levels,
      time_idx
    )
  
  v1 <-
    extract_4d_slice_fast(
      arr_v,
      particles$x,
      particles$y,
      particles$z,
      levels,
      time_idx
    )
  
  ##########################################################
  # K2
  ##########################################################
  
  x2 <-
    particles$x +
    u1 *
    dt /
    2 *
    deg_lon_scale
  
  y2 <-
    particles$y +
    v1 *
    dt /
    2 *
    deg_lat_scale
  
  u2 <-
    extract_4d_slice_fast(
      arr_u,
      x2,
      y2,
      particles$z,
      levels,
      time_idx
    )
  
  v2 <-
    extract_4d_slice_fast(
      arr_v,
      x2,
      y2,
      particles$z,
      levels,
      time_idx
    )
  
  ##########################################################
  # K3
  ##########################################################
  
  x3 <-
    particles$x +
    u2 *
    dt /
    2 *
    deg_lon_scale
  
  y3 <-
    particles$y +
    v2 *
    dt /
    2 *
    deg_lat_scale
  
  u3 <-
    extract_4d_slice_fast(
      arr_u,
      x3,
      y3,
      particles$z,
      levels,
      time_idx
    )
  
  v3 <-
    extract_4d_slice_fast(
      arr_v,
      x3,
      y3,
      particles$z,
      levels,
      time_idx
    )
  
  ##########################################################
  # K4
  ##########################################################
  
  x4 <-
    particles$x +
    u3 *
    dt *
    deg_lon_scale
  
  y4 <-
    particles$y +
    v3 *
    dt *
    deg_lat_scale
  
  u4 <-
    extract_4d_slice_fast(
      arr_u,
      x4,
      y4,
      particles$z,
      levels,
      time_idx
    )
  
  v4 <-
    extract_4d_slice_fast(
      arr_v,
      x4,
      y4,
      particles$z,
      levels,
      time_idx
    )
  
  ##########################################################
  # RK4 UPDATE
  ##########################################################
  
  particles$x <-
    particles$x +
    
    (
      u1 +
        2 * u2 +
        2 * u3 +
        u4
    ) *
    dt /
    6 *
    deg_lon_scale
  
  particles$y <-
    particles$y +
    
    (
      v1 +
        2 * v2 +
        2 * v3 +
        v4
    ) *
    dt /
    6 *
    deg_lat_scale
  
  ##########################################################
  # HORIZONTAL DIFFUSION
  ##########################################################
  
  Kh <-
    2
  
  sigma <-
    sqrt(
      2 *
        Kh *
        dt
    )
  
  particles$x <-
    particles$x +
    
    rnorm(
      nrow(particles),
      0,
      sigma
    ) *
    deg_lon_scale
  
  particles$y <-
    particles$y +
    
    rnorm(
      nrow(particles),
      0,
      sigma
    ) *
    deg_lat_scale
  
  ##########################################################
  # BOUNDARIES
  ##########################################################
  
  particles$x <-
    pmax(
      x_min,
      pmin(
        x_max,
        particles$x
      )
    )
  
  particles$y <-
    pmax(
      y_min,
      pmin(
        y_max,
        particles$y
      )
    )
  
  particles
}

############################################################
# 12. STOKES SINKING
############################################################

kernel_stokes_sinking_3d_fast <- function(
    particles,
    arr_temp,
    arr_sal,
    levels,
    dt_sec,
    time_idx
) {
  
  if (
    nrow(particles) == 0
  ) {
    
    return(
      particles
    )
  }
  
  env_t <-
    extract_4d_slice_fast(
      arr_temp,
      particles$x,
      particles$y,
      particles$z,
      levels,
      time_idx
    )
  
  env_s <-
    extract_4d_slice_fast(
      arr_sal,
      particles$x,
      particles$y,
      particles$z,
      levels,
      time_idx
    )
  
  ##########################################################
  # SEAWATER DENSITY
  ##########################################################
  
  rho_f <-
    1000 +
    0.78 *
    env_s -
    0.07 *
    env_t -
    0.0045 *
    env_t^2
  
  ##########################################################
  # VISCOSITY
  ##########################################################
  
  mu <-
    0.001779 /
    (
      1 +
        0.03368 *
        env_t +
        0.000221 *
        env_t^2
    )
  
  ##########################################################
  # DIAMETER IS IN MICROMETERS
  ##########################################################
  
  radius_m <-
    particles$diameter *
    1e-6 /
    2
  
  ##########################################################
  # STOKES VELOCITY
  ##########################################################
  
  w_s <-
    (
      2 *
        9.81 *
        radius_m^2 *
        (
          particles$rho_p -
            rho_f
        )
    ) /
    (
      9 *
        mu
    )
  
  particles$settling_velocity_ms <-
    w_s
  
  particles$z <-
    particles$z +
    w_s *
    dt_sec
  
  particles$z <-
    pmax(
      min(levels),
      pmin(
        max(levels),
        particles$z
      )
    )
  
  particles
}

############################################################
# 13. PHYTOPLANKTON DYNAMICS
############################################################

kernel_phytoplankton_dynamics <- function(
    particles,
    dt,
    env_no3,
    env_nh4,
    env_fe,
    env_po4,
    env_si,
    env_par,
    env_temp
) {
  
  if (
    nrow(particles) == 0
  ) {
    
    return(
      particles
    )
  }
  
  ##########################################################
  # AGE
  ##########################################################
  
  dt_hours <-
    dt /
    3600
  
  particles$age_hours <-
    particles$age_hours +
    dt_hours
  
  ##########################################################
  # MORTALITY
  ##########################################################
  
  m_early <-
    0.05 /86400
  
  prob_die_early <-
    m_early *
    exp(
      -0.15 *
        particles$age_hours
    ) *
    dt
  
  viability <-
    1 /
    (
      1 +
        exp(
          0.12 *
            (
              particles$age_hours -
                72
            )
        )
    )
  
  prob_die_late <-
    (
      1 -
        viability
    ) *
    (
      0.01 /
        86400
    ) *
    dt
  
  death_probability <-
    pmin(
      0.95,
      prob_die_early +
        prob_die_late
    )
  
  surviving_mask <-
    runif(
      nrow(particles)
    ) >
    death_probability
  
  ##########################################################
  # KEEP SURVIVORS
  ##########################################################
  
  particles <-
    particles[
      surviving_mask,
      ,
      drop = FALSE
    ]
  
  if (
    nrow(particles) == 0
  ) {
    
    return(
      particles
    )
  }
  
  ##########################################################
  # MATCH ENVIRONMENT TO SURVIVORS
  ##########################################################
  
  viability <-
    viability[
      surviving_mask
    ]
  
  e_no3 <-
    env_no3[
      surviving_mask
    ]
  
  e_nh4 <-
    env_nh4[
      surviving_mask
    ]
  
  e_fe <-
    env_fe[
      surviving_mask
    ]
  
  e_po4 <-
    env_po4[
      surviving_mask
    ]
  
  e_si <-
    env_si[
      surviving_mask
    ]
  
  e_par <-
    env_par[
      surviving_mask
    ]
  
  e_temp <-
    env_temp[
      surviving_mask
    ]
  
  ##########################################################
  # NUTRIENT CONSTANTS
  ##########################################################
  
  K_NO3 <-
    0.5
  
  K_NH4 <-
    0.1
  
  K_Fe <-
    0.2
  
  K_PO4 <-
    0.05
  
  K_Si <-
    1.0
  
  Psi <-
    1.5
  
  ##########################################################
  # NITROGEN LIMITATION
  ##########################################################
  
  lim_NH4 <-
    e_nh4 /
    (
      e_nh4 +
        K_NH4
    )
  
  lim_NO3 <-
    (
      e_no3 /
        (
          e_no3 +
            K_NO3
        )
    ) *
    exp(
      -Psi *
        e_nh4
    )
  
  lim_N <-
    pmin(
      1,
      lim_NH4 +
        lim_NO3
    )
  
  ##########################################################
  # OTHER NUTRIENTS
  ##########################################################
  
  lim_Fe <-
    e_fe /
    (
      e_fe +
        K_Fe
    )
  
  lim_PO4 <-
    e_po4 /
    (
      e_po4 +
        K_PO4
    )
  
  lim_Si <-
    e_si /
    (
      e_si +
        K_Si
    )
  
  ##########################################################
  # LIEBIG MINIMUM
  ##########################################################
  
  lim_Nutrient <-
    pmin(
      lim_N,
      lim_Fe,
      lim_PO4,
      lim_Si
    )
  
  ##########################################################
  # LIGHT
  ##########################################################
  
  I_opt <-
    50
  
  beta_I <-
    0.4
  
  I_rel <-
    e_par /
    I_opt
  
  lim_I <-
    (
      2 *
        (
          1 +
            beta_I
        ) *
        I_rel
    ) /
    (
      I_rel^2 +
        2 *
        beta_I *
        I_rel +
        1
    )
  
  lim_I <-
    pmax(
      0,
      pmin(
        1,
        lim_I
      )
    )
  
  ##########################################################
  # TEMPERATURE
  ##########################################################
  
  T_opt <-
    24
  
  T_let <-
    12
  
  beta_T <-
    0.5
  
  theta <-
    pmax(
      0,
      (
        e_temp -
          T_let
      ) /
        (
          T_opt -
            T_let
        )
    )
  
  lim_T <-
    (
      2 *
        (
          1 +
            beta_T
        ) *
        theta
    ) /
    (
      theta^2 +
        2 *
        beta_T *
        theta +
        1
    )
  
  lim_T <-
    pmax(
      0,
      pmin(
        1,
        lim_T
      )
    )
  
  ##########################################################
  # GROWTH RATE
  ##########################################################
  
  mu_max <-
    1.5 /
    86400
  
  mu <-
    mu_max *
    lim_I *
    lim_T *
    lim_Nutrient *
    viability
  
  ##########################################################
  # DIVISION PROBABILITY
  ##########################################################
  
  division_probability <-
    pmin(
      1,
      mu *
        dt
    )
  
  dividing_mask <-
    runif(
      nrow(particles)
    ) <
    division_probability
  
  ##########################################################
  # CELL DIVISION
  ##########################################################
  
  if (
    any(dividing_mask)
  ) {
    
    new_cells <-
      particles[
        dividing_mask,
        ,
        drop = FALSE
      ]
    
    n_new <-
      nrow(new_cells)
    
    max_id <-
      ifelse(
        nrow(particles) > 0,
        max(particles$id),
        0
      )
    
    new_cells$id <-
      seq(
        max_id + 1,
        max_id + n_new
      )
    
    new_cells$age_hours <-
      0
    
    ########################################################
    # SMALL DAUGHTER CELL SEPARATION
    ########################################################
    
    nudge <-
      0.0009
    
    new_cells$x <-
      new_cells$x +
      rnorm(
        n_new,
        0,
        nudge
      )
    
    new_cells$y <-
      new_cells$y +
      rnorm(
        n_new,
        0,
        nudge
      )
    
    ########################################################
    # CRITICAL FIX:
    # FORCE IDENTICAL COLUMNS BEFORE RBIND
    ########################################################
    
    common_cols <-
      intersect(
        names(particles),
        names(new_cells)
      )
    
    new_cells <-
      new_cells[
        ,
        common_cols,
        drop = FALSE
      ]
    
    particles <-
      particles[
        ,
        common_cols,
        drop = FALSE
      ]
    
    particles <-
      rbind(
        particles,
        new_cells
      )
  }
  
  particles
}

############################################################
# 14. PATCH COHESION
############################################################

apply_patch_movement_cohesion_fast <- function(
    particles,
    patch_radius = 0.05,
    cohesion_factor = 0.1
) {
  
  if (
    nrow(particles) <= 1
  ) {
    
    return(
      particles
    )
  }
  
  grid_x <-
    round(
      particles$x /
        patch_radius
    ) *
    patch_radius
  
  grid_y <-
    round(
      particles$y /
        patch_radius
    ) *
    patch_radius
  
  patch_key <-
    paste(
      grid_x,
      grid_y,
      sep = "_"
    )
  
  particles$patch_id <-
    match(
      patch_key,
      unique(patch_key)
    )
  
  counts <-
    tabulate(
      particles$patch_id
    )
  
  mean_x <-
    rowsum(
      particles$x,
      particles$patch_id
    ) /
    counts
  
  mean_y <-
    rowsum(
      particles$y,
      particles$patch_id
    ) /
    counts
  
  particles$x <-
    particles$x +
    
    (
      mean_x[
        particles$patch_id,
        1
      ] -
        particles$x
    ) *
    cohesion_factor
  
  particles$y <-
    particles$y +
    
    (
      mean_y[
        particles$patch_id,
        1
      ] -
        particles$y
    ) *
    cohesion_factor
  
  particles$patch_id <-
    NULL
  
  particles
}

############################################################
# 15. INITIAL PARTICLES
############################################################

set.seed(101)

n_particles <-
  10

particle_set <-
  data.frame(
    
    id =
      seq_len(
        n_particles
      ),
    
    x =
      runif(
        n_particles,
        x_min,
        x_max
      ),
    
    y =
      runif(
        n_particles,
        y_min,
        y_max
      ),
    
    z =
      rep(
        10,
        n_particles
      ),
    
    ########################################################
    # CELL DIAMETER IN MICROMETERS
    ########################################################
    
    diameter =
      runif(
        n_particles,
        3,
        70
      ),
    
    ########################################################
    # CELL DENSITY
    ########################################################
    
    rho_p =
      runif(
        n_particles,
        1000,
        1023
      ),
    
    settling_velocity_ms =
      rep(
        0,
        n_particles
      ),
    
    age_hours =
      runif(
        n_particles,
        0,
        24
      ),
    
    cell_carbon =
      runif(
        n_particles,
        50,
        200
      ),
    
    cell_nitrogen =
      runif(
        n_particles,
        5,
        30
      )
  )

############################################################
# 16. SIMULATION
############################################################

dt <-
  10800

full_simulation_archive <-
  vector(
    "list",
    n_times
  )

cat(
  "Running simulation...\n"
)

for (
  step in seq_len(n_times)
) {
  
  cat(
    "Step",
    step,
    "of",
    n_times,
    "\n"
  )
  
  current_time_idx <-
    step
  
  ##########################################################
  # PHYSICAL MOVEMENT
  ##########################################################
  
  if (
    step > 1 &&
    nrow(particle_set) > 0
  ) {
    
    particle_set <-
      kernel_advection_rk4_3d_fast(
        particle_set,
        master_u,
        master_v,
        real_depths,
        dt,
        current_time_idx
      )
    
    particle_set <-
      kernel_stokes_sinking_3d_fast(
        particle_set,
        master_temp,
        master_sal,
        real_depths,
        dt,
        current_time_idx
      )
    
    ########################################################
    # ENVIRONMENTAL EXTRACTION
    ########################################################
    
    p_no3 <-
      extract_4d_slice_fast(
        master_no3,
        particle_set$x,
        particle_set$y,
        particle_set$z,
        real_depths,
        current_time_idx
      )
    
    p_nh4 <-
      extract_4d_slice_fast(
        master_nh4,
        particle_set$x,
        particle_set$y,
        particle_set$z,
        real_depths,
        current_time_idx
      )
    
    p_fe <-
      extract_4d_slice_fast(
        master_fe,
        particle_set$x,
        particle_set$y,
        particle_set$z,
        real_depths,
        current_time_idx
      )
    
    p_po4 <-
      extract_4d_slice_fast(
        master_po4,
        particle_set$x,
        particle_set$y,
        particle_set$z,
        real_depths,
        current_time_idx
      )
    
    p_si <-
      extract_4d_slice_fast(
        master_si,
        particle_set$x,
        particle_set$y,
        particle_set$z,
        real_depths,
        current_time_idx
      )
    
    p_par <-
      extract_4d_slice_fast(
        master_par,
        particle_set$x,
        particle_set$y,
        particle_set$z,
        real_depths,
        current_time_idx
      )
    
    p_temp <-
      extract_4d_slice_fast(
        master_temp,
        particle_set$x,
        particle_set$y,
        particle_set$z,
        real_depths,
        current_time_idx
      )
    
    ########################################################
    # PHYTOPLANKTON GROWTH AND DIVISION
    ########################################################
    
    particle_set <-
      kernel_phytoplankton_dynamics(
        particle_set,
        dt,
        p_no3,
        p_nh4,
        p_fe,
        p_po4,
        p_si,
        p_par,
        p_temp
      )
    
    ########################################################
    # PATCH COHESION
    ########################################################
    
    if (
      nrow(particle_set) > 1
    ) {
      
      particle_set <-
        apply_patch_movement_cohesion_fast(
          particle_set,
          patch_radius = 0.005,
          cohesion_factor = 0.7
        )
    }
  }
  
  ##########################################################
  # ARCHIVE
  ##########################################################
  
  if (
    nrow(particle_set) > 0
  ) {
    
    snapshot <-
      particle_set
    
    ########################################################
    # ADD U AND V (INPUT CURRENTS AT PARTICLE POSITION)
    ########################################################
    
    snapshot$u <-
      extract_4d_slice_fast(
        master_u,
        snapshot$x,
        snapshot$y,
        snapshot$z,
        real_depths,
        current_time_idx
      )
    
    snapshot$v <-
      extract_4d_slice_fast(
        master_v,
        snapshot$x,
        snapshot$y,
        snapshot$z,
        real_depths,
        current_time_idx
      )
    
    snapshot$simulation_step <-
      step
    
    snapshot$elapsed_hours <-
      (
        step - 1
      ) *
      3
    
  } else {
    
    snapshot <-
      data.frame(
        simulation_step =
          step,
        
        elapsed_hours =
          (
            step - 1
          ) *
          3
      )
  }
  
  full_simulation_archive[[step]] <-
    snapshot
}

cat(
  "Simulation completed successfully.\n"
)

############################################################
# 17. COMBINE ARCHIVE
############################################################

simulation_results <-
  do.call(
    rbind,
    full_simulation_archive
  )

############################################################
# 18. BASIC SUMMARY
############################################################

cat(
  "Final number of cells:",
  nrow(particle_set),
  "\n"
)

cat(
  "Maximum population:",
  max(
    sapply(
      full_simulation_archive,
      nrow
    )
  ),
  "\n"
)

############################################################
# 19. SAVE RESULTS
############################################################

write.csv(
  simulation_results,
  "phytoplankton_simulation_results.csv",
  row.names = FALSE
)

saveRDS(
  full_simulation_archive,
  "phytoplankton_simulation_archive.rds"
)

cat(
  "Results saved successfully.\n"
)