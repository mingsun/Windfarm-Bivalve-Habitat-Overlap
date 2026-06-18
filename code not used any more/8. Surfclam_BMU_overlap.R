library(tidyverse)
library(sf)
library(terra)
library(kohonen)

# ============================================================ #
# 1. Settings and input files ----
# ============================================================ #

## file paths and options ----

species_name <- "surfclam"
year_tag <- "2025"
months <- sprintf("%02d", 1:12)

tow_file <- "../results/HSI/surfclam/SC_tows_with_composite_SI.csv"

sediment_file <- "../data/sediment data/kriged.sediment/sediment_phi_kriged_1km.tif"
depth_file <- "../data/bathymetry data/depth_gebco_1km_domain.tif"
temp_dir <- "../data/temperature data"
hsi_dir <- "../results/HSI/surfclam/HSI_weighted"

wea_shp <- "../data/boem-renewable-energy-shapefiles_2.5.2025/Offshore_Wind_Leases_outlines.shp"
wea_rast_file <- "../results/WEA_overlay/surfclam_WEA_on_HSI_grid.tif"

result_dir <- file.path("../results/BMU", species_name)
plot_dir <- "../plot/BMU"

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(
  file.path(result_dir, "model"),
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  file.path(result_dir, "rasters"),
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  file.path(result_dir, "cell_tables"),
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  file.path(result_dir, "monthly"),
  recursive = TRUE,
  showWarnings = FALSE
)

som_vars <- c("TEMP", "DEPTH", "SED_PHI_MEAN")

# 3 x 3 grid = 9 BMUs
# chosen to balance habitat heterogeneity and interpretability
som_grid <- somgrid(xdim = 3, ydim = 3, topo = "hexagonal")

set.seed(123)

## helper functions ----

nearest_bmu <- function(x, codebook) {
  d2 <- outer(
    rowSums(x^2),
    rowSums(codebook^2),
    FUN = "+"
  ) -
    2 * x %*% t(codebook)

  list(
    unit = max.col(-d2, ties.method = "first"),
    dist = sqrt(apply(d2, 1, min))
  )
}

class_from_rank <- function(rank_vec, low_lab, mid_lab, high_lab) {
  case_when(
    rank_vec == 1 ~ low_lab,
    rank_vec == 2 ~ mid_lab,
    rank_vec == 3 ~ high_lab,
    TRUE ~ NA_character_
  )
}

# ============================================================ #
# 2. Fit SOM and build BMU interpretation table ----
# ============================================================ #

## load tow data and fit SOM ----

tow <- read.csv(tow_file) |>
  filter(
    !is.na(LON),
    !is.na(LAT),
    !is.na(TEMP),
    !is.na(DEPTH),
    !is.na(SED_PHI_MEAN)
  ) |>
  mutate(row_id = row_number(), .before = 1)

x <- tow |>
  dplyr::select(all_of(som_vars)) |>
  as.matrix()

x_scaled <- scale(x)

scale_info <- list(
  vars = som_vars,
  center = attr(x_scaled, "scaled:center"),
  scale = attr(x_scaled, "scaled:scale")
)

som_model <- kohonen::som(
  X = x_scaled,
  grid = som_grid,
  rlen = 200,
  alpha = c(0.05, 0.01),
  keep.data = TRUE
)

saveRDS(scale_info, file.path(result_dir, "model", "SOM_scale_info.rds"))
saveRDS(som_model, file.path(result_dir, "model", "SOM_model.rds"))

## assign BMUs to tows and build interpretation table ----

tow$BMU <- som_model$unit.classif

codebook <- as_tibble(som_model$codes[[1]]) |>
  setNames(paste0("scaled_", som_vars)) |>
  mutate(BMU = row_number())

for (v in som_vars) {
  codebook[[v]] <- codebook[[paste0("scaled_", v)]] *
    scale_info$scale[v] +
    scale_info$center[v]
}

bmu_table <- tow |>
  group_by(BMU) |>
  summarise(
    n_tows = n(),
    mean_HSI = mean(HSI_weighted, na.rm = TRUE),
    .groups = "drop"
  ) |>
  right_join(codebook, by = "BMU") |>
  mutate(n_tows = replace_na(n_tows, 0))

# qualitative classes based on ranking among BMUs
bmu_table <- bmu_table |>
  mutate(
    temp_rank = ntile(TEMP, 3),
    depth_rank = ntile(DEPTH, 3),
    sed_rank = ntile(SED_PHI_MEAN, 3),
    temp_class = class_from_rank(temp_rank, "cold", "moderate", "warm"),
    depth_class = class_from_rank(depth_rank, "shallow", "mid-depth", "deep"),
    sed_class = class_from_rank(sed_rank, "coarse", "intermediate", "muddy"),
    BMU_label = paste(temp_class, depth_class, sed_class, sep = "-")
  ) |>
  arrange(temp_rank, depth_rank, sed_rank, BMU) |>
  mutate(BMU_order = row_number()) |>
  dplyr::select(
    BMU,
    BMU_order,
    BMU_label,
    temp_class,
    depth_class,
    sed_class,
    TEMP,
    DEPTH,
    SED_PHI_MEAN,
    n_tows,
    mean_HSI,
    starts_with("scaled_")
  )

tow <- tow |>
  left_join(
    bmu_table |>
      dplyr::select(
        BMU,
        BMU_order,
        BMU_label,
        temp_class,
        depth_class,
        sed_class
      ),
    by = "BMU"
  )

write.csv(tow, file.path(result_dir, "SC_tows_with_BMU.csv"), row.names = FALSE)
write.csv(
  bmu_table,
  file.path(result_dir, "SOM_BMU_interpretation_table.csv"),
  row.names = FALSE
)

# optional diagnostics saved for record
png(
  file.path(plot_dir, paste0(species_name, "_SOM_training_changes.png")),
  width = 7,
  height = 5,
  units = "in",
  res = 300
)
plot(
  som_model,
  type = "changes",
  main = paste(species_name, "SOM training progress")
)
dev.off()

png(
  file.path(plot_dir, paste0(species_name, "_SOM_counts.png")),
  width = 7,
  height = 5,
  units = "in",
  res = 300
)
plot(
  som_model,
  type = "counts",
  main = paste(species_name, "Counts per SOM unit")
)
dev.off()

png(
  file.path(plot_dir, paste0(species_name, "_SOM_quality.png")),
  width = 7,
  height = 5,
  units = "in",
  res = 300
)
plot(
  som_model,
  type = "quality",
  main = paste(species_name, "Mean distance to BMU")
)
dev.off()

# ============================================================ #
# 3. Prepare static rasters and WEA layer ----
# ============================================================ #

## load static rasters ----

sediment_rast <- rast(sediment_file)
depth_rast <- rast(depth_file)
wea_rast <- rast(wea_rast_file)

compareGeom(sediment_rast, depth_rast, wea_rast)

## prepare WEA outline for plotting ----

wea_plot <- st_read(wea_shp, quiet = TRUE) |>
  st_make_valid()

bbox_vect <- terra::as.polygons(
  terra::ext(sediment_rast),
  crs = terra::crs(sediment_rast)
)
bbox_sf <- st_as_sf(bbox_vect)

wea_plot <- wea_plot |>
  st_transform(st_crs(bbox_sf)) |>
  st_crop(st_bbox(bbox_sf)) |>
  st_transform(crs(sediment_rast))

# ============================================================ #
# 4. Project BMUs month by month and save raster outputs ----
# ============================================================ #

for (mm in months) {
  message("Processing month: ", mm)

  temp_file <- file.path(
    temp_dir,
    paste0("bottomT_insitu_", year_tag, mm, "_1km_domain.tif")
  )

  hsi_file <- file.path(
    hsi_dir,
    paste0(species_name, "_HSI_weighted_", year_tag, mm, ".tif")
  )

  temp_rast <- rast(temp_file)
  hsi_rast <- rast(hsi_file)

  compareGeom(sediment_rast, depth_rast, temp_rast, hsi_rast, wea_rast)

  ## classify raster cells into BMUs ----

  env_df <- c(sediment_rast, depth_rast, temp_rast) |>
    `names<-`(c("SED_PHI_MEAN", "DEPTH", "TEMP")) |>
    as.data.frame(xy = TRUE, cells = TRUE, na.rm = FALSE) |>
    filter(
      !is.na(TEMP),
      !is.na(DEPTH),
      !is.na(SED_PHI_MEAN)
    )

  env_scaled <- env_df |>
    dplyr::select(all_of(som_vars)) |>
    as.matrix()

  for (v in som_vars) {
    env_scaled[, v] <- (env_scaled[, v] - scale_info$center[v]) /
      scale_info$scale[v]
  }

  bmu_out <- nearest_bmu(env_scaled, som_model$codes[[1]])

  env_df$BMU <- bmu_out$unit
  env_df$BMU_dist <- bmu_out$dist

  env_df <- env_df |>
    left_join(
      bmu_table |>
        dplyr::select(
          BMU,
          BMU_order,
          BMU_label,
          temp_class,
          depth_class,
          sed_class
        ),
      by = "BMU"
    )

  ## build and save BMU rasters ----

  bmu_rast <- rast(hsi_rast)
  values(bmu_rast) <- NA
  values(bmu_rast)[env_df$cell] <- env_df$BMU
  names(bmu_rast) <- "BMU"

  bmu_order_rast <- rast(hsi_rast)
  values(bmu_order_rast) <- NA
  values(bmu_order_rast)[env_df$cell] <- env_df$BMU_order
  names(bmu_order_rast) <- "BMU_order"

  bmu_dist_rast <- rast(hsi_rast)
  values(bmu_dist_rast) <- NA
  values(bmu_dist_rast)[env_df$cell] <- env_df$BMU_dist
  names(bmu_dist_rast) <- "BMU_dist"

  writeRaster(
    bmu_rast,
    file.path(result_dir, "rasters", paste0("SOM_BMU_", year_tag, mm, ".tif")),
    overwrite = TRUE
  )

  writeRaster(
    bmu_order_rast,
    file.path(
      result_dir,
      "rasters",
      paste0("SOM_BMU_order_", year_tag, mm, ".tif")
    ),
    overwrite = TRUE
  )

  writeRaster(
    bmu_dist_rast,
    file.path(
      result_dir,
      "rasters",
      paste0("SOM_BMU_distance_", year_tag, mm, ".tif")
    ),
    overwrite = TRUE
  )

  write.csv(
    env_df,
    file.path(
      result_dir,
      "cell_tables",
      paste0("BMU_cells_", year_tag, mm, ".csv")
    ),
    row.names = FALSE
  )

  ## monthly BMU map for appendix ----

  png(
    file.path(
      plot_dir,
      paste0(species_name, "_BMU_map_", year_tag, mm, ".png")
    ),
    width = 8,
    height = 6,
    units = "in",
    res = 300
  )
  plot(
    bmu_order_rast,
    main = paste(species_name, "BMU map -", year_tag, mm),
    type = "classes"
  )
  plot(st_geometry(wea_plot), add = TRUE, border = "black", lwd = 0.8)
  dev.off()

  rm(
    temp_rast,
    hsi_rast,
    env_df,
    env_scaled,
    bmu_out,
    bmu_rast,
    bmu_order_rast,
    bmu_dist_rast
  )
}

# ============================================================ #
# 5. Overlay BMU with HSI and WEA and calculate enrichment ----
# ============================================================ #

overlay_all <- vector("list", length(months))
enrich_all <- vector("list", length(months))

for (i in seq_along(months)) {
  mm <- months[i]

  bmu_rast <- rast(file.path(
    result_dir,
    "rasters",
    paste0("SOM_BMU_", year_tag, mm, ".tif")
  ))
  hsi_rast <- rast(file.path(
    hsi_dir,
    paste0(species_name, "_HSI_weighted_", year_tag, mm, ".tif")
  ))

  compareGeom(bmu_rast, hsi_rast, wea_rast)

  overlay_df <- c(bmu_rast, hsi_rast, wea_rast) |>
    `names<-`(c("BMU", "HSI", "in_WEA")) |>
    as.data.frame(xy = TRUE, cells = TRUE, na.rm = FALSE) |>
    filter(!is.na(BMU), !is.na(HSI)) |>
    mutate(
      in_WEA = ifelse(is.na(in_WEA), 0, in_WEA),
      in_WEA = factor(
        in_WEA,
        levels = c(0, 1),
        labels = c("Outside WEA", "Inside WEA")
      ),
      month = mm
    ) |>
    left_join(
      bmu_table |>
        dplyr::select(
          BMU,
          BMU_order,
          BMU_label,
          temp_class,
          depth_class,
          sed_class
        ),
      by = "BMU"
    )

  write.csv(
    overlay_df,
    file.path(
      result_dir,
      "monthly",
      paste0("BMU_HSI_WEA_overlay_", year_tag, mm, ".csv")
    ),
    row.names = FALSE
  )

  write.csv(
    overlay_df |>
      count(month, BMU, BMU_order, BMU_label, in_WEA, name = "n_cells") |>
      arrange(month, BMU_order, in_WEA),
    file.path(
      result_dir,
      "monthly",
      paste0("BMU_counts_", year_tag, mm, ".csv")
    ),
    row.names = FALSE
  )

  write.csv(
    overlay_df |>
      group_by(month, BMU, BMU_order, BMU_label, in_WEA) |>
      summarise(
        mean_HSI = mean(HSI, na.rm = TRUE),
        median_HSI = median(HSI, na.rm = TRUE),
        n_cells = n(),
        .groups = "drop"
      ) |>
      arrange(month, BMU_order, in_WEA),
    file.path(
      result_dir,
      "monthly",
      paste0("BMU_HSI_summary_", year_tag, mm, ".csv")
    ),
    row.names = FALSE
  )

  enrich_df <- overlay_df |>
    count(
      BMU,
      BMU_order,
      BMU_label,
      temp_class,
      depth_class,
      sed_class,
      in_WEA,
      name = "n_cells"
    ) |>
    group_by(in_WEA) |>
    mutate(prop = n_cells / sum(n_cells)) |>
    ungroup() |>
    mutate(
      in_WEA = ifelse(in_WEA == "Inside WEA", "Inside_WEA", "Outside_WEA")
    ) |>
    pivot_wider(
      names_from = in_WEA,
      values_from = c(n_cells, prop),
      values_fill = 0
    ) |>
    mutate(
      enrichment_ratio = ifelse(
        prop_Outside_WEA > 0,
        prop_Inside_WEA / prop_Outside_WEA,
        NA_real_
      ),
      log2_enrichment = log2(enrichment_ratio),
      month = mm
    ) |>
    arrange(BMU_order)

  overlay_all[[i]] <- overlay_df
  enrich_all[[i]] <- enrich_df

  rm(bmu_rast, hsi_rast, overlay_df, enrich_df)
}

## combined final outputs ----

enrich_all_df <- bind_rows(enrich_all)

write.csv(
  enrich_all_df,
  file.path(result_dir, paste0("BMU_enrichment_", year_tag, "_all_months.csv")),
  row.names = FALSE
)

enrich_mean_df <- enrich_all_df |>
  group_by(BMU, BMU_order, BMU_label, temp_class, depth_class, sed_class) |>
  summarise(
    mean_enrichment_ratio = mean(enrichment_ratio, na.rm = TRUE),
    mean_log2_enrichment = mean(log2_enrichment, na.rm = TRUE),
    mean_prop_inside = mean(prop_Inside_WEA, na.rm = TRUE),
    mean_prop_outside = mean(prop_Outside_WEA, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(BMU_order)

write.csv(
  enrich_mean_df,
  file.path(
    result_dir,
    paste0("BMU_enrichment_", year_tag, "_mean_across_months.csv")
  ),
  row.names = FALSE
)
