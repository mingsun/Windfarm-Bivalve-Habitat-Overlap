library(tidyverse)
library(terra)

source("../4.bivalve_habitat_overlap_R/functions/check_temp_extrapolation.R")

# 1. Load aligned environmental rasters ---- ----

sediment_phi <- rast(
  "../data/sediment data/kriged.sediment/sediment_phi_kriged_1km.tif"
)
depth_rast <- rast("../data/bathymetry data/depth_gebco_1km_domain.tif")

# quick geometry check: all should match exactly
compareGeom(sediment_phi, depth_rast, stopOnError = FALSE)

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

#  Get training temperature range
temp_min_train <- min(temp_SI$mid, na.rm = TRUE)
temp_max_train <- max(temp_SI$mid, na.rm = TRUE)

cat("Training temp range:", temp_min_train, "to", temp_max_train, "\n")


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

# optional visual checks
plot(SI_sed_rast)
plot(SI_depth_rast)

# ------------------------------------------------------------ #

# 6. Monthly temperature files ----

months <- sprintf("%02d", 1:12)

temp_files <- paste0(
  "../data/temperature data/bottomT_insitu_2025",
  months,
  "_1km_domain.tif"
)

# ------------------------------------------------------------ #

# 7. Loop through months and generate HSI ----

# create a results container
extrap_list <- list()

for (i in seq_along(temp_files)) {
  month_id <- months[i]

  cat("Processing month:", month_id, "\n")

  temp_rast <- rast(temp_files[i])

  extrap_list[[month_id]] <- calc_temp_extrapolation(
    temp_rast = temp_rast,
    temp_SI = temp_SI,
    month_id = month_id
  )

  # ensure geometry matches
  compareGeom(sediment_phi, temp_rast)

  # convert temperature to SI
  SI_temp_rast <- make_si_raster(
    r = temp_rast,
    si_df = temp_SI,
    x_col = "mid",
    y_col = "SI_temp"
  )

  # compute weighted HSI
  HSI_weighted <- w_sed *
    SI_sed_rast +
    w_dep * SI_depth_rast +
    w_tmp * SI_temp_rast

  names(HSI_weighted) <- "HSI_weighted"

  # quick check
  print(summary(values(HSI_weighted)))

  # save monthly HSI
  out_file <- paste0(
    "../results/HSI/surfclam/HSI_weighted/surfclam_HSI_weighted_2025",
    month_id,
    ".tif"
  )

  writeRaster(
    HSI_weighted,
    out_file,
    overwrite = TRUE
  )

  cat("Saved:", out_file, "\n")
}

extrap_df <- bind_rows(extrap_list)

write.csv(
  extrap_df,
  "../results/HSI/surfclam/HSI_weighted/temperature_extrapolation_summary_2025.csv",
  row.names = FALSE
)

# ------------------------------------------------------------ #
