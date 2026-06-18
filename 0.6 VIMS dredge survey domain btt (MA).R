library(terra)
library(gsw)
library(tidyverse)


# 1. Load static template rasters ----

sediment_phi <- rast(
  "../data/sediment data/kriged.sediment/VIMS_sediment_phi_kriged_1km.tif"
)
depth_resampled <- rast(
  "../data/bathymetry data/VIMS_depth_gebco_1km_domain.tif"
)

# make sure depth matches template exactly
depth_resampled <- project(depth_resampled, sediment_phi, method = "bilinear")

# ------------------------------------------------------------------------------------------- #

# 2. Precompute static grid information once ----
# Only keep cells inside the analysis domain (non-NA sediment and non-NA depth)

template_df <- as.data.frame(
  sediment_phi,
  xy = TRUE,
  cells = TRUE,
  na.rm = FALSE
)
names(template_df)[4] <- "sediment_val"

depth_vals <- values(depth_resampled)

template_df$depth_m <- depth_vals[template_df$cell]

# valid analysis cells
template_df <- template_df %>%
  filter(!is.na(sediment_val) & !is.na(depth_m))

# convert x/y to lon/lat once
pts <- vect(template_df[, c("x", "y")], crs = crs(sediment_phi))
pts_ll <- project(pts, "EPSG:4326")
ll <- crds(pts_ll)

template_df$lon <- ll[, 1]
template_df$lat <- ll[, 2]

# pressure depends only on depth and latitude, so compute once
template_df$p <- gsw_p_from_z(
  z = -template_df$depth_m,
  latitude = template_df$lat
)

# keep only columns needed later
template_df <- template_df |>
  dplyr::select(cell, x, y, lon, lat, depth_m, p)

# ------------------------------------------------------------------------------------------- #

# 3. Helper function for deepest non-NA salinity ----

deepest_non_na <- function(x) {
  idx <- max(which(!is.na(x)))
  if (is.finite(idx)) {
    return(x[idx])
  } else {
    return(NA_real_)
  }
}

# ------------------------------------------------------------------------------------------- #

# 4. Monthly files ----

months <- sprintf("%02d", 1:12)

glorys_files <- paste0(
  "../data/temperature data/mercatorglorys12v1_gl12_mean_2025",
  months,
  ".nc"
)

# ------------------------------------------------------------------------------------------- #

# 5. Process each month ----

for (i in seq_along(glorys_files)) {
  glorys_file <- glorys_files[i]
  month_id <- months[i]

  cat("Processing month:", month_id, "\n")

  # --------------------
  # 5.1 Load monthly bottom potential temperature
  # --------------------

  bottomT <- rast(glorys_file, subds = "bottomT")
  values(bottomT)[!is.finite(values(bottomT))] <- NA

  bottomT_resampled <- project(
    bottomT,
    sediment_phi,
    method = "bilinear"
  )

  # --------------------
  # 5.2 Load monthly salinity and derive bottom salinity
  # --------------------

  so <- rast(glorys_file, subds = "so")

  so_bottom <- app(so, deepest_non_na)

  so_bottom_resampled <- project(
    so_bottom,
    sediment_phi,
    method = "bilinear"
  )

  # --------------------
  # 5.3 Extract monthly values only at precomputed valid cells
  # --------------------

  pt0_vals <- values(bottomT_resampled)[template_df$cell]
  sp_vals <- values(so_bottom_resampled)[template_df$cell]

  month_df <- template_df %>%
    mutate(
      pt0 = pt0_vals,
      SP = sp_vals
    ) %>%
    filter(!is.na(pt0) & !is.na(SP))

  # --------------------
  # 5.4 TEOS-10 conversion for this month
  # --------------------

  month_df$SA <- gsw_SA_from_SP(
    SP = month_df$SP,
    p = month_df$p,
    longitude = month_df$lon,
    latitude = month_df$lat
  )

  month_df$CT <- gsw_CT_from_pt(
    SA = month_df$SA,
    pt = month_df$pt0
  )

  month_df$t_insitu <- gsw_t_from_CT(
    SA = month_df$SA,
    CT = month_df$CT,
    p = month_df$p
  )

  # --------------------
  # 5.5 Write back to raster
  # --------------------

  bottomT_insitu <- sediment_phi
  values(bottomT_insitu) <- NA_real_

  vals <- values(bottomT_insitu)
  vals[month_df$cell] <- month_df$t_insitu
  values(bottomT_insitu) <- vals

  # --------------------
  # 5.6 Save monthly raster
  # --------------------

  out_file <- paste0(
    "../data/temperature data/VIMS_bottomT_insitu_2025",
    month_id,
    "_1km_domain.tif"
  )

  writeRaster(
    bottomT_insitu,
    out_file,
    overwrite = TRUE
  )

  # optional quick summary
  print(summary(month_df$t_insitu))
}

# ------------------------------------------------------------------------------------------- #

# -----------------------
# 6) Optional check: load one result
# -----------------------

# test_rast <- rast("../data/temperature data/bottomT_insitu_202508_1km_domain.tif")
# plot(test_rast)
# summary(values(test_rast), na.rm = TRUE)
