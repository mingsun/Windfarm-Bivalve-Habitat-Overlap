library(tidyverse)
library(sf)
library(terra)

# 1. Load weighted HSI raster ----

# use January HSI as the template; all 12 monthly HSI rasters should match
HSI_template <- rast(
  "../results/HSI/scallop/HSI_weighted/scallop_HSI_weighted_202501.tif"
)
plot(HSI_template)

# ----------------------------------------------------------------- #

# 2. WEA shapefile ----

# 2.1 load WEA shapefile and check CRS
wea_sf <- st_read(
  "../data/boem-renewable-energy-shapefiles_2.5.2025/Offshore_Wind_Leases_outlines.shp",
  quiet = TRUE
)

# fix invalid geometries
wea_sf <- st_make_valid(wea_sf)

plot(st_geometry(wea_sf))

# ----------------------------------------------------------------- #

# 3. crop WEA to HSI extent ----

# convert HSI raster extent to polygon
bbox_vect <- terra::as.polygons(
  terra::ext(HSI_template),
  crs = terra::crs(HSI_template)
)
bbox_poly <- sf::st_as_sf(bbox_vect)

# project crop box to WEA CRS
bbox_poly_wea <- sf::st_transform(bbox_poly, sf::st_crs(wea_sf))

# crop WEA polygons to study-region rectangle
wea_crop <- sf::st_crop(wea_sf, sf::st_bbox(bbox_poly_wea))

plot(st_geometry(wea_crop))

# ----------------------------------------------------------------- #

# 4. Project cropped WEA to the HSI raster CRS ----

wea_hsi <- st_transform(wea_crop, crs(HSI_template))

# quick check
plot(st_geometry(wea_hsi))
plot(as.polygons(HSI_template), add = TRUE)

# ----------------------------------------------------------------- #

# 5. Rasterize WEA to the exact HSI grid ----

# cells inside WEA = 1
# cells outside WEA = NA initially, then convert to 0

wea_vect <- vect(wea_hsi)

WEA_raster <- rasterize(
  wea_vect,
  HSI_template,
  field = 1,
  background = NA
)

# convert outside-domain NA values to 0 only where HSI exists
WEA_raster <- ifel(!is.na(HSI_template) & is.na(WEA_raster), 0, WEA_raster)

names(WEA_raster) <- "in_WEA"

plot(WEA_raster)

# save
writeRaster(
  WEA_raster,
  "../results/WEA_overlay/scallop_WEA_on_HSI_grid.tif",
  overwrite = TRUE
)

# ----------------------------------------------------------------- #

# 6. Monthly HSI files ----

months <- sprintf("%02d", 1:12)

hsi_files <- paste0(
  "../results/HSI/scallop/HSI_weighted/scallop_HSI_weighted_2025",
  months,
  ".tif"
)
# ----------------------------------------------------------------- #

# 7. Loop through months and build overlay_df ----

overlay_list <- vector("list", length(months))

for (i in seq_along(months)) {
  month_id <- months[i]
  cat("Processing month:", month_id, "\n")

  HSI_weighted <- rast(hsi_files[i])

  # make sure geometry matches
  compareGeom(HSI_template, HSI_weighted)

  # build monthly overlay stack
  overlay_stack <- c(HSI_weighted, WEA_raster)
  names(overlay_stack) <- c("HSI", "in_WEA")

  overlay_df <- as.data.frame(overlay_stack, xy = TRUE, na.rm = FALSE) |>
    filter(!is.na(HSI)) |>
    mutate(
      in_WEA = ifelse(is.na(in_WEA), 0, in_WEA),
      in_WEA = factor(
        in_WEA,
        levels = c(0, 1),
        labels = c("Outside WEA", "Inside WEA")
      ),
      month = month_id
    )

  # save monthly cell level overlay
  write.csv(
    overlay_df,
    paste0(
      "../results/WEA_overlay/scallop_HSI_WEA_cell_overlay_2025",
      month_id,
      ".csv"
    ),
    row.names = FALSE
  )

  overlay_list[[i]] <- overlay_df
}

# ----------------------------------------------------------------- #

# 8. Combine all months ----

overlay_df_all <- bind_rows(overlay_list)

summary(overlay_df_all$HSI)
table(overlay_df_all$in_WEA)
table(overlay_df_all$month)

# save combined table
saveRDS(
  overlay_df_all,
  "../results/WEA_overlay/scallop_HSI_WEA_cell_overlay_2025_all_months.rds"
)

# ----------------------------------------------------------------- #

# 9. Optional quick summaries ----

# mean HSI by month and WEA status
overlay_df_all |>
  group_by(month, in_WEA) |>
  summarise(
    mean_HSI = mean(HSI, na.rm = TRUE),
    median_HSI = median(HSI, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )

# total number of cells inside vs outside WEA by month
overlay_df_all |>
  count(month, in_WEA)

# ----------------------------------------------------------------- #
