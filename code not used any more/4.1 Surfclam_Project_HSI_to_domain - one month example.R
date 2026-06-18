library(tidyverse)
library(terra)


# 1. Load aligned environmental rasters ---- ----

sediment_phi <- rast(
  "../data/sediment data/kriged.sediment/sediment_phi_kriged_1km.tif"
)
depth_rast <- rast("../data/bathymetry data/depth_gebco_1km_domain.tif")
temp_rast <- rast("../data/temperature data/bottomT_insitu_1km_domain.tif")

# quick geometry check: all should match exactly
compareGeom(sediment_phi, depth_rast, stopOnError = FALSE)
compareGeom(sediment_phi, temp_rast, stopOnError = FALSE)

# ---------------------------------------------------------------------------------------------------------- #

# 2. Load SI tables ----

sediment_SI <- readRDS("../results/HSI/surfclam/surfclam_sediment_SI.rds")
depth_SI <- readRDS("../results/HSI/surfclam/surfclam_depth_SI.rds")
temp_SI <- readRDS("../results/HSI/surfclam/surfclam_temp_SI.rds")

# ------------------------------------------------------------ #

# 3. Load BRT weights ----

brt_weights <- read.csv("../results/HSI/surfclam/surfclam_BRT_weights.csv")

w_sed <- brt_weights$weight[brt_weights$var == "SED_PHI_MEAN"]
w_dep <- brt_weights$weight[brt_weights$var == "DEPTH"]
w_tmp <- brt_weights$weight[brt_weights$var == "TEMP"]

# ------------------------------------------------------------ #

# 4. Helper function: convert raster env values to SI values ----

# Uses linear interpolation from the saved SI curve for between the bin values
# rule = 2 extends edge values flat beyond the observed range

make_si_raster <- function(r, si_df, x_col, y_col) {
  # keep only needed columns and sort by x
  lut <- si_df |>
    dplyr::select(x = all_of(x_col), y = all_of(y_col)) |>
    distinct() |>
    arrange(x)

  terra::app(
    r,
    fun = function(x) {
      approx(
        x = lut$x,
        y = lut$y,
        xout = x,
        method = "linear",
        rule = 2 # extend edge values flat beyond observed range
      )$y
    }
  )
}

# ------------------------------------------------------------ #

# 5. Project each environmental raster to SI raster ----

SI_sed_rast <- make_si_raster(
  r = sediment_phi,
  si_df = sediment_SI,
  x_col = "mid",
  y_col = "SI_sediment"
)

SI_depth_rast <- make_si_raster(
  r = depth_rast,
  si_df = depth_SI,
  x_col = "mid",
  y_col = "SI_depth"
)

SI_temp_rast <- make_si_raster(
  r = temp_rast,
  si_df = temp_SI,
  x_col = "mid",
  y_col = "SI_temp"
)

# optional visual checks
plot(SI_sed_rast)
plot(SI_depth_rast)
plot(SI_temp_rast)

# ------------------------------------------------------------ #

# 6. Compute weighted HSI raster ----

HSI_weighted <- w_sed *
  SI_sed_rast +
  w_dep * SI_depth_rast +
  w_tmp * SI_temp_rast

names(HSI_weighted) <- "HSI_weighted"

# quick check
summary(values(HSI_weighted))
plot(HSI_weighted)

# ------------------------------------------------------------ #

# 7. Save final weighted HSI raster ----

writeRaster(
  HSI_weighted,
  "../results/HSI/surfclam/surfclam_HSI_weighted_1km.tif",
  overwrite = TRUE
)

# ------------------------------------------------------------ #
