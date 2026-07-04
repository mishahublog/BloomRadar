library(terra)

# ==============================================================================
# 1. READ & ULTRA-FAST 4D MEMORY CACHING
# ==============================================================================
nc_file <- "~/Downloads/runs_FMRC_ESPC-D-V02_all_RUN_2026-06-27T12_00_00Z.nc4"
master_stack <- rast(nc_file)

real_depths <- c(0, 2, 4, 6, 8, 10, 15, 20, 25, 30, 40, 50, 60, 75, 100, 125, 150, 
                 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000, 1250, 1500, 
                 1750, 2000, 2500, 3000, 3500, 4000, 4250, 4500, 4750, 4900, 5000)

n_depths  <- length(real_depths)
grid_ext  <- ext(master_stack)
grid_dims <- dim(master_stack) 
n_steps   <- 12                # 12 steps * 3 hours = 36 total forecast hours covered

x_min <- grid_ext$xmin; x_max <- grid_ext$xmax; x_res <- (x_max - x_min) / grid_dims[2]
y_min <- grid_ext$ymin; y_max <- grid_ext$ymax; y_res <- (y_max - y_min) / grid_dims[1]

# Filter layers once
sal_names  <- grep("sal", names(master_stack), value = TRUE)
temp_names <- grep("temp", names(master_stack), value = TRUE)
u_names    <- grep("water_u", names(master_stack), value = TRUE)
v_names    <- grep("water_v", names(master_stack), value = TRUE)

# --- Blazing Fast 4D Array Caching Structure [Rows, Cols, Depths, Time] ---
cache_4d_variable <- function(stack, layer_names, n_depths) {
  arr_4d <- array(0, dim = c(grid_dims[1], grid_dims[2], n_depths, 37))
  for (t in 1:37) {
    start_layer <- (t - 1) * n_depths + 1
    end_layer   <- t * n_depths
    if (end_layer <= length(layer_names)) {
      sub_stack <- subset(stack, layer_names[start_layer:end_layer])
      arr_4d[, , , t] <- as.array(sub_stack)
    }
  }
  return(arr_4d)
}

cat("Caching NetCDF data into RAM...\n")
master_u    <- cache_4d_variable(master_stack, u_names, n_depths)
master_v    <- cache_4d_variable(master_stack, v_names, n_depths)
master_temp <- cache_4d_variable(master_stack, temp_names, n_depths)
master_sal  <- cache_4d_variable(master_stack, sal_names, n_depths)

# Mock baseline data structures for nutrient & light limits if missing from physical file
master_no3  <- master_temp * 0 + 2.5   # Standardized ambient concentrations
master_nh4  <- master_temp * 0 + 0.3
master_par  <- master_temp * 0 + 45.0  # Steady state PAR for daylight calculations

# ==============================================================================
# 2. OPTIMIZED 4D MEMORY TRILINEAR INTERPOLATION ENGINE
# ==============================================================================
extract_4d_slice_fast <- function(data_arr_4d, x, y, z, levels, time_idx) {
  col_f <- 1 + (x - x_min) / x_res
  row_f <- 1 + (y_max - y) / y_res
  
  col_f <- pmax(1, pmin(grid_dims[2], col_f))
  row_f <- pmax(1, pmin(grid_dims[1], row_f))
  
  c1 <- floor(col_f); c2 <- pmin(grid_dims[2], c1 + 1)
  r1 <- floor(row_f); r2 <- pmin(grid_dims[1], r1 + 1)
  
  w_c2 <- col_f - c1; w_c1 <- 1 - w_c2
  w_r2 <- row_f - r1; w_r1 <- 1 - w_r2
  
  z_bounded <- pmax(min(levels), pmin(max(levels), z))
  l_idx <- findInterval(z_bounded, levels)
  l_idx <- pmax(1, pmin(length(levels) - 1, l_idx))
  u_idx <- l_idx + 1
  
  exact_max <- (z_bounded == max(levels))
  l_idx[exact_max] <- length(levels) - 1
  u_idx[exact_max] <- length(levels)
  
  z_low <- levels[l_idx]
  z_upp <- levels[u_idx]
  w_z2  <- (z_bounded - z_low) / (z_upp - z_low)
  w_z1  <- 1 - w_z2
  
  idx_r1_c1_l <- cbind(r1, c1, l_idx, time_idx)
  idx_r1_c1_u <- cbind(r1, c1, u_idx, time_idx)
  idx_r1_c2_l <- cbind(r1, c2, l_idx, time_idx)
  idx_r1_c2_u <- cbind(r1, c2, u_idx, time_idx)
  idx_r2_c1_l <- cbind(r2, c1, l_idx, time_idx)
  idx_r2_c1_u <- cbind(r2, c1, u_idx, time_idx)
  idx_r2_c2_l <- cbind(r2, c2, l_idx, time_idx)
  idx_r2_c2_u <- cbind(r2, c2, u_idx, time_idx)
  
  v_l <- (data_arr_4d[idx_r1_c1_l] * w_c1 + data_arr_4d[idx_r1_c2_l] * w_c2) * w_r1 +
    (data_arr_4d[idx_r2_c1_l] * w_c1 + data_arr_4d[idx_r2_c2_l] * w_c2) * w_r2
  
  # FIX: Corrected final weight from w_r1 to w_r2 to complete vertical matrix map
  v_u <- (data_arr_4d[idx_r1_c1_u] * w_c1 + data_arr_4d[idx_r1_c2_u] * w_c2) * w_r1 +
    (data_arr_4d[idx_r2_c1_u] * w_c1 + data_arr_4d[idx_r2_c2_u] * w_c2) * w_r2
  
  output_values <- v_l * w_z1 + v_u * w_z2
  output_values[is.na(output_values)] <- 0
  
  return(output_values)
}

# ==============================================================================
# 3. HIGH-PERFORMANCE 3D OCEAN KERNELS
# ==============================================================================
kernel_advection_rk4_3d_fast <- function(particles, arr_u, arr_v, levels, dt, time_idx) {
  deg_lat_scale <- 1 / 111000 
  deg_lon_scale <- 1 / (111000 * cos(particles$y * pi / 180))
  
  u1 <- extract_4d_slice_fast(arr_u, particles$x, particles$y, particles$z, levels, time_idx)
  v1 <- extract_4d_slice_fast(arr_v, particles$x, particles$y, particles$z, levels, time_idx)
  
  pos_k2_x <- pmax(59.96, pmin(75.0, particles$x + u1 * (dt / 2) * deg_lon_scale))
  pos_k2_y <- pmax(14.98, pmin(22.02, particles$y + v1 * (dt / 2) * deg_lat_scale))
  u2 <- extract_4d_slice_fast(arr_u, pos_k2_x, pos_k2_y, particles$z, levels, time_idx)
  v2 <- extract_4d_slice_fast(arr_v, pos_k2_x, pos_k2_y, particles$z, levels, time_idx)
  
  pos_k3_x <- pmax(59.96, pmin(75.0, particles$x + u2 * (dt / 2) * deg_lon_scale))
  pos_k3_y <- pmax(14.98, pmin(22.02, particles$y + v2 * (dt / 2) * deg_lat_scale))
  u3 <- extract_4d_slice_fast(arr_u, pos_k3_x, pos_k3_y, particles$z, levels, time_idx)
  v3 <- extract_4d_slice_fast(arr_v, pos_k3_x, pos_k3_y, particles$z, levels, time_idx)
  
  pos_k4_x <- pmax(59.96, pmin(75.0, particles$x + u3 * dt * deg_lon_scale))
  pos_k4_y <- pmax(14.98, pmin(22.02, particles$y + v3 * dt * deg_lat_scale))
  u4 <- extract_4d_slice_fast(arr_u, pos_k4_x, pos_k4_y, particles$z, levels, time_idx)
  v4 <- extract_4d_slice_fast(arr_v, pos_k4_x, pos_k4_y, particles$z, levels, time_idx)
  
  particles$x <- particles$x + (u1 + 2*u2 + 2*u3 + u4) * (dt / 6) * deg_lon_scale
  particles$y <- particles$y + (v1 + 2*v2 + 2*v3 + v4) * (dt / 6) * deg_lat_scale
  
  Kh <- 2.0 
  sigma <- sqrt(2 * Kh * dt) 
  
  particles$x <- particles$x + rnorm(nrow(particles), 0, sigma) * deg_lon_scale
  particles$y <- particles$y + rnorm(nrow(particles), 0, sigma) * deg_lat_scale
  
  particles$x <- pmax(59.96, pmin(74.99, particles$x))
  particles$y <- pmax(14.98, pmin(22.01, particles$y))
  
  return(particles)
}

kernel_stokes_sinking_3d_fast <- function(particles, arr_temp, arr_sal, levels, dt_sec, time_idx) {
  env_t <- extract_4d_slice_fast(arr_temp, particles$x, particles$y, particles$z, levels, time_idx)
  env_s <- extract_4d_slice_fast(arr_sal, particles$x, particles$y, particles$z, levels, time_idx)
  
  g <- 9.81
  rho_f <- 1000 + (0.78 * env_s) - (0.07 * env_t) - (0.0045 * env_t^2)
  mu    <- 0.001779 / (1 + 0.03368 * env_t + 0.000221 * env_t^2)
  
  w_s <- (g * (particles$rho_p - rho_f) * (particles$diameter)^2) / (18 * mu)
  particles$settling_velocity_ms <- w_s
  
  particles$z <- particles$z + (w_s * dt_sec)
  particles$z <- pmax(0, pmin(max(levels), particles$z)) 
  
  return(particles)
}

kernel_phytoplankton_dynamics <- function(particles, dt, env_no3, env_nh4, env_par, env_temp) {
  if (nrow(particles) == 0) return(particles)
  
  # ----------------------------------------------------------------------------
  # A. AGE PROGRESSION
  # ----------------------------------------------------------------------------
  dt_hours <- dt / 3600
  particles$age_hours <- particles$age_hours + dt_hours
  
  # ----------------------------------------------------------------------------
  # B. STOCHASTIC MORTALITY & VIABILITY FILTER
  # ----------------------------------------------------------------------------
  # 1. Early-stage mortality: High initial risk that decays exponentially
  # Drops significantly after the first 6 hours of cell division
  m_early_base <- 0.05 / 3600  # Base hourly risk converted to per-second
  prob_die_early <- (m_early_base * exp(-0.15 * particles$age_hours)) * dt
  
  # 2. Late-stage viability: Logistic drop-off (inflection point at 72 hours)
  # Viability stays close to 1.0 early on, then collapses towards 0.0 late-stage
  viability <- 1 / (1 + exp(0.12 * (particles$age_hours - 72)))
  prob_die_late <- (1 - viability) * (0.04 / 3600) * dt # Scaled death risk
  
  # Combine background and age-dependent risks
  total_death_prob <- pmin(0.95, prob_die_early + prob_die_late)
  surviving_mask   <- runif(nrow(particles)) > total_death_prob
  
  # Apply immediate mortality pruning
  particles <- particles[surviving_mask, ]
  if (nrow(particles) == 0) return(particles)
  
  # ----------------------------------------------------------------------------
  # C. REPRODUCTIVE DIVISION ENGINE (MODECOGeL Based)
  # ----------------------------------------------------------------------------
  mu_max <- 1.5 / 86400  
  K_NO3  <- 0.5; K_NH4 <- 0.1; Psi <- 1.5          
  I_opt  <- 50;  beta_I <- 0.4          
  T_opt  <- 24;  T_let <- 12;  beta_T <- 0.5          
  
  lim_NH4 <- env_nh4 / (env_nh4 + K_NH4)
  lim_NO3 <- (env_no3 / (env_no3 + K_NO3)) * exp(-Psi * env_nh4)
  lim_N   <- lim_NO3 + lim_NH4
  
  I_rel   <- env_par / I_opt
  lim_I   <- (2 * (1 + beta_I) * I_rel) / (I_rel^2 + 2 * beta_I * I_rel + 1)
  lim_I   <- pmax(0, pmin(1, lim_I))
  
  theta   <- pmax(0, (env_temp - T_let) / (T_opt - T_let)) 
  lim_T   <- (2 * (1 + beta_T) * theta) / (theta^2 + 2 * beta_T * theta + 1)
  lim_T   <- pmax(0, lim_T)
  
  # Only viable/healthy cells execute high-efficiency division
  mu <- mu_max * lim_I * lim_T * lim_N * viability[surviving_mask]
  
  division_prob <- mu * dt
  dividing_mask <- runif(nrow(particles)) < division_prob
  
  if (any(dividing_mask)) {
    new_cells <- particles[dividing_mask, ]
    max_id    <- max(particles$id)
    new_cells$id <- (max_id + 1):(max_id + nrow(new_cells))
    
    # RESET AGE FOR NEWBORN DAUGHTER CELLS
    new_cells$age_hours <- 0
    
    # Small spatial displacement at birth
    nudge_scale <- 0.1
    new_cells$x <- new_cells$x + rnorm(nrow(new_cells), 0, nudge_scale)
    new_cells$y <- new_cells$y + rnorm(nrow(new_cells), 0, nudge_scale)
    
    particles <- rbind(particles, new_cells)
  }
  
  return(particles)
}


#--------
# this is a test
#----------

# apply_patch_movement_cohesion <- function(particles, patch_radius = 0.1, cohesion_factor = 0.2) {
#   # Cluster nearby cells into dynamic irregular patches
#   hc <- hclust(dist(particles[, c("x", "y")]), method = "quick")
#   # Group into discrete irregular patches based on your spatial threshold
#   particles$patch_id <- cutree(hc, h = patch_radius)
#   
#   # Pull individual cells slightly toward their specific irregular patch center
#   for(p in unique(particles$patch_id)) {
#     idx <- which(particles$patch_id == p)
#     if(length(idx) > 1) {
#       mean_x <- mean(particles$x[idx])
#       mean_y <- mean(particles$y[idx])
#       
#       # Nudge cells toward the patch gravity center to maintain patch boundaries
#       particles$x[idx] <- particles$x[idx] + (mean_x - particles$x[idx]) * cohesion_factor
#       particles$y[idx] <- particles$y[idx] + (mean_y - particles$y[idx]) * cohesion_factor
#     }
#   }
#   return(particles)
# }


# ==============================================================================
# 4. PARTICLESET INITIALIZATION 
# ==============================================================================
set.seed(101)
n_particles <- 10

particle_set <- data.frame(
  id                   = 1:n_particles,
  x                    = runif(n_particles, 65.0, 70.0),
  y                    = runif(n_particles, 18.0, 20.0),
  z                    = rep(0.5, n_particles),
  diameter             = runif(n_particles, 30e-6, 90e-6), 
  rho_p                = runif(n_particles, 1035, 1045),  
  settling_velocity_ms = 0,
  age_hours            = runif(n_particles, 0, 24) # Seed cells have mixed initial ages
)
# ==============================================================================
# 5. SIMULATION LOGIC CORE & DATA ARCHIVING
# ==============================================================================
dt <- 10800  # 3 hours per step

full_simulation_archive <- list()

cat("Running physics engine loops...\n")
for (step in 0:n_steps) {
  current_time_idx <- (step * 3) + 1 
  if(current_time_idx > 37) current_time_idx <- 37
  
  if (step > 0) {
    # Extract structural 3D matrix localized vectors right before calculating physiological growth
    p_no3  <- extract_4d_slice_fast(master_no3, particle_set$x, particle_set$y, particle_set$z, real_depths, current_time_idx)
    p_nh4  <- extract_4d_slice_fast(master_nh4, particle_set$x, particle_set$y, particle_set$z, real_depths, current_time_idx)
    p_par  <- extract_4d_slice_fast(master_par, particle_set$x, particle_set$y, particle_set$z, real_depths, current_time_idx)
    p_temp <- extract_4d_slice_fast(master_temp, particle_set$x, particle_set$y, particle_set$z, real_depths, current_time_idx)
    
    # Execution pipe
    particle_set <- particle_set |> 
      kernel_advection_rk4_3d_fast(master_u, master_v, real_depths, dt, current_time_idx) |> 
      kernel_stokes_sinking_3d_fast(master_temp, master_sal, real_depths, dt, current_time_idx) |> 
      kernel_phytoplankton_dynamics(dt, p_no3, p_nh4, p_par, p_temp) # FIX: Added required extracted environmental variables
  }
  
  snapshot <- particle_set
  snapshot$simulation_step <- step
  snapshot$elapsed_hours    <- step * 3
  full_simulation_archive[[step + 1]] <- snapshot
}

# ==============================================================================
# 6. EXPORT UNIFIED MASTER TRAJECTORY DATASET
# ==============================================================================
cat("Assembling final dataset...\n")
master_export_df <- do.call(rbind, full_simulation_archive)

master_export_df <- master_export_df[, c("simulation_step", "elapsed_hours", "id", 
                                         "x", "y", "z", "diameter", "rho_p", "age_hours",
                                         "settling_velocity_ms")]

write.csv(master_export_df, "simulated_phytoplankton_trajectories.csv", row.names = FALSE)
saveRDS(master_export_df, "simulated_phytoplankton_trajectories.rds")

cat("Done! 'simulated_phytoplankton_trajectories.csv' successfully generated with growth data.\n")

growth_summary <- aggregate(id ~ elapsed_hours + simulation_step, 
                            data = master_export_df, 
                            FUN = length)

names(growth_summary)[3] <- "total_cells"
print(growth_summary)