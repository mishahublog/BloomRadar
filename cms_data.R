library(CopernicusMarine)

# Download the subset
my_cms_data <- cms_download_subset(
  product = "GLOBAL_ANALYSISFORECAST_BGC_001_028", # Ensure this is the correct Product ID
  layer = "cmems_mod_glo_bgc-nut_anfc_0.25deg_P1D-m",
  variable = c("fe", "no3", "po4", "si"),
  region = c(60, 65, 15, 20),           # c(xmin, ymin, xmax, ymax)
  timerange = c("2026-07-27", "2026-07-27"),
  verticalrange = c(0.494, -40.34)
)

my_phy_data <- cms_download_subset(
  product = "GLOBAL_ANALYSISFORECAST_PHY_001_024", # Ensure this is the correct Product ID
  layer = "cmems_mod_glo_phy_anfc_0.083deg_PT1H-m",
  variable = c("so", "thetao", "uo", "vo", "zos"),
  region = c(60, 65, 15, 20),           # c(xmin, ymin, xmax, ymax)
  timerange = c("2026-10-25", "2026-11-10"),
  verticalrange = c(0.494, 0.494)
)

