library(tidyverse)
library(sf)
library(terra)
library(kohonen)

# ============================================================ #
# 1. Settings and input files ----
# ============================================================ #

## file paths and options ----

year_tag <- "2025"
months <- sprintf("%02d", 1:12)

# --- aggregated tow files for all species ---
surfclam_tow <- read.csv(
  "../results/HSI/surfclam/SC_tows_with_composite_SI.csv"
) |>
  mutate(species = "surfclam")
quahog_tow <- read.csv("../results/HSI/quahog/OQ_tows_with_composite_SI.csv") |>
  mutate(species = "quahog")
scallop_tow <- read.csv(
  "../results/HSI/scallop/SS_tows_with_composite_SI.csv"
) |>
  mutate(species = "scallop")

scallop_tow <- scallop_tow |>
  rename(
    CRUISE6 = CruiseID,
    YR = YEAR,
    STRATUM = Month,
    STATION = StationID,
    Dist_m = TowDist
  )


# --- raster inputs ---
sediment_file <- "../data/sediment data/kriged.sediment/sediment_phi_kriged_1km.tif"
depth_file <- "../data/bathymetry data/depth_gebco_1km_domain.tif"
temp_dir <- "../data/temperature data"

# choose one common raster domain for the aggregated BMU analysis
# this should be the same domain to map BMUs and evaluate WEA overlap
# *scallop have different HSI domains, here we use a common environmental domain
template_hsi_dir <- "../results/HSI/surfclam/HSI_weighted"
template_hsi_species <- "surfclam"

wea_shp <- "../data/boem-renewable-energy-shapefiles_2.5.2025/Offshore_Wind_Leases_outlines.shp"
wea_rast_file <- "../results/WEA_overlay/surfclam_WEA_on_HSI_grid.tif"

result_dir <- "../results/BMU/aggregated_species"
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
som_grid <- kohonen::somgrid(xdim = 3, ydim = 3, topo = "hexagonal")

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
# 2. Build aggregated tow dataset and fit SOM ----
# ============================================================ #

## load and combine tow data ----
## required columns in each file: LON, LAT, TEMP, DEPTH, SED_PHI_MEAN

tow <- bind_rows(surfclam_tow, quahog_tow, scallop_tow) |>
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
    n_tows_total = n(),
    n_surfclam = sum(species == "surfclam"),
    n_ocean_quahog = sum(species == "ocean_quahog"),
    n_scallop = sum(species == "scallop"),
    .groups = "drop"
  ) |>
  right_join(codebook, by = "BMU") |>
  mutate(
    n_tows_total = replace_na(n_tows_total, 0),
    n_surfclam = replace_na(n_surfclam, 0),
    n_ocean_quahog = replace_na(n_ocean_quahog, 0),
    n_scallop = replace_na(n_scallop, 0)
  )

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
  arrange(TEMP, BMU) |>
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
    n_tows_total,
    n_surfclam,
    n_ocean_quahog,
    n_scallop,
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

write.csv(
  tow,
  file.path(result_dir, "aggregated_tows_with_BMU.csv"),
  row.names = FALSE
)
write.csv(
  bmu_table,
  file.path(result_dir, "BMU_interpretation_table.csv"),
  row.names = FALSE
)

# optional diagnostics
png(
  file.path(plot_dir, "aggregated_species_SOM_training_changes.png"),
  width = 7,
  height = 5,
  units = "in",
  res = 300
)
plot(
  som_model,
  type = "changes",
  main = "Aggregated species SOM training progress"
)
dev.off()

png(
  file.path(plot_dir, "aggregated_species_SOM_counts.png"),
  width = 7,
  height = 5,
  units = "in",
  res = 300
)
plot(
  som_model,
  type = "counts",
  main = "Aggregated species counts per SOM unit"
)
dev.off()

png(
  file.path(plot_dir, "aggregated_species_SOM_quality.png"),
  width = 7,
  height = 5,
  units = "in",
  res = 300
)
plot(
  som_model,
  type = "quality",
  main = "Aggregated species mean distance to BMU"
)
dev.off()

# ============================================================ #
# 3. Prepare common static rasters and WEA layer ----
# ============================================================ #

sediment_rast <- rast(sediment_file)
depth_rast <- rast(depth_file)
wea_rast <- rast(wea_rast_file)

compareGeom(sediment_rast, depth_rast, wea_rast)

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
# 4. Project BMUs to monthly rasters ----
# ============================================================ #

# to save plot
map_df_list <- vector("list", length(months))

for (mm in months) {
  message("Processing month: ", mm)

  temp_file <- file.path(
    temp_dir,
    paste0("bottomT_insitu_", year_tag, mm, "_1km_domain.tif")
  )

  # only used as a template raster for geometry / plotting
  template_hsi_file <- file.path(
    template_hsi_dir,
    paste0(template_hsi_species, "_HSI_weighted_", year_tag, mm, ".tif")
  )

  temp_rast <- rast(temp_file)
  template_rast <- rast(template_hsi_file)

  compareGeom(sediment_rast, depth_rast, temp_rast, template_rast, wea_rast)

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

  bmu_rast <- rast(template_rast)
  values(bmu_rast) <- NA
  values(bmu_rast)[env_df$cell] <- env_df$BMU
  names(bmu_rast) <- "BMU"

  bmu_order_rast <- rast(template_rast)
  values(bmu_order_rast) <- NA
  values(bmu_order_rast)[env_df$cell] <- env_df$BMU_order
  names(bmu_order_rast) <- "BMU_order"

  bmu_dist_rast <- rast(template_rast)
  values(bmu_dist_rast) <- NA
  values(bmu_dist_rast)[env_df$cell] <- env_df$BMU_dist
  names(bmu_dist_rast) <- "BMU_dist"

  writeRaster(
    bmu_rast,
    file.path(
      result_dir,
      "rasters",
      paste0("aggregated_BMU_", year_tag, mm, ".tif")
    ),
    overwrite = TRUE
  )

  writeRaster(
    bmu_order_rast,
    file.path(
      result_dir,
      "rasters",
      paste0("aggregated_BMU_order_", year_tag, mm, ".tif")
    ),
    overwrite = TRUE
  )

  writeRaster(
    bmu_dist_rast,
    file.path(
      result_dir,
      "rasters",
      paste0("aggregated_BMU_distance_", year_tag, mm, ".tif")
    ),
    overwrite = TRUE
  )

  write.csv(
    env_df,
    file.path(
      result_dir,
      "cell_tables",
      paste0("aggregated_BMU_cells_", year_tag, mm, ".csv")
    ),
    row.names = FALSE
  )

  map_df_list[[which(months == mm)]] <- env_df |>
    mutate(month = mm)

  png(
    file.path(
      plot_dir,
      paste0("aggregated_species_BMU_map_", year_tag, mm, ".png")
    ),
    width = 8,
    height = 6,
    units = "in",
    res = 300
  )

  plot(
    bmu_order_rast,
    main = paste("Aggregated species BMU map -", year_tag, mm),
    type = "classes"
  )
  plot(st_geometry(wea_plot), add = TRUE, border = "black", lwd = 0.8)
  dev.off()

  rm(
    temp_rast,
    template_rast,
    env_df,
    env_scaled,
    bmu_out,
    bmu_rast,
    bmu_order_rast,
    bmu_dist_rast
  )
}

# ============================================================ #
# 5. Monthly and yearly BMU-WEA overlap / enrichment tables ----
# ============================================================ #

monthly_table_list <- vector("list", length(months))

for (i in seq_along(months)) {
  mm <- months[i]

  bmu_rast <- rast(file.path(
    result_dir,
    "rasters",
    paste0("aggregated_BMU_", year_tag, mm, ".tif")
  ))

  compareGeom(bmu_rast, wea_rast)

  overlay_df <- c(bmu_rast, wea_rast) |>
    `names<-`(c("BMU", "in_WEA")) |>
    as.data.frame(xy = TRUE, cells = TRUE, na.rm = FALSE) |>
    filter(!is.na(BMU)) |>
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
          sed_class,
          TEMP,
          DEPTH,
          SED_PHI_MEAN
        ),
      by = "BMU"
    )

  monthly_counts <- overlay_df |>
    count(
      BMU,
      BMU_order,
      BMU_label,
      temp_class,
      depth_class,
      sed_class,
      TEMP,
      DEPTH,
      SED_PHI_MEAN,
      in_WEA,
      name = "n_cells"
    ) |>
    pivot_wider(
      names_from = in_WEA,
      values_from = n_cells,
      values_fill = 0
    ) |>
    rename(
      inside_wea_cells = `Inside WEA`,
      outside_wea_cells = `Outside WEA`
    ) |>
    mutate(
      total_cells = inside_wea_cells + outside_wea_cells,
      month = mm
    ) |>
    arrange(BMU_order)

  total_inside <- sum(monthly_counts$inside_wea_cells, na.rm = TRUE)
  total_outside <- sum(monthly_counts$outside_wea_cells, na.rm = TRUE)

  monthly_counts <- monthly_counts |>
    mutate(
      prop_inside = inside_wea_cells / total_inside,
      prop_outside = outside_wea_cells / total_outside,
      enrichment_ratio = prop_inside / prop_outside
    ) |>
    dplyr::select(
      month,
      BMU,
      BMU_order,
      BMU_label,
      temp_class,
      depth_class,
      sed_class,
      TEMP,
      DEPTH,
      SED_PHI_MEAN,
      total_cells,
      inside_wea_cells,
      outside_wea_cells,
      prop_inside,
      prop_outside,
      enrichment_ratio
    )

  write.csv(
    monthly_counts,
    file.path(
      result_dir,
      "monthly",
      paste0("aggregated_BMU_overlap_enrichment_", year_tag, mm, ".csv")
    ),
    row.names = FALSE
  )

  monthly_table_list[[i]] <- monthly_counts

  rm(bmu_rast, overlay_df, monthly_counts)
}

monthly_table_all <- bind_rows(monthly_table_list)

write.csv(
  monthly_table_all,
  file.path(
    result_dir,
    paste0("aggregated_BMU_overlap_enrichment_", year_tag, "_all_months.csv")
  ),
  row.names = FALSE
)

yearly_summary_table <- monthly_table_all |>
  group_by(
    BMU,
    BMU_order,
    BMU_label,
    temp_class,
    depth_class,
    sed_class,
    TEMP,
    DEPTH,
    SED_PHI_MEAN
  ) |>
  summarise(
    mean_total_cells = mean(total_cells, na.rm = TRUE),
    mean_inside_wea_cells = mean(inside_wea_cells, na.rm = TRUE),
    mean_outside_wea_cells = mean(outside_wea_cells, na.rm = TRUE),
    mean_prop_inside = mean(prop_inside, na.rm = TRUE),
    mean_prop_outside = mean(prop_outside, na.rm = TRUE),
    mean_enrichment_ratio = mean(enrichment_ratio, na.rm = TRUE),
    sd_enrichment_ratio = sd(enrichment_ratio, na.rm = TRUE),
    log2_mean_enrichment_ratio = log2(mean(enrichment_ratio, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  arrange(BMU_order)

write.csv(
  yearly_summary_table,
  file.path(
    result_dir,
    paste0(
      "aggregated_BMU_overlap_enrichment_",
      year_tag,
      "_yearly_summary.csv"
    )
  ),
  row.names = FALSE
)


# ============================================================ #
# 6. Combined 12-month BMU map figure ----
# ============================================================ #

map_df_all <- bind_rows(map_df_list)

# build legend labels: "BMU 1: cold-shallow-coarse"
bmu_legend_df <- bmu_table |>
  arrange(BMU_order) |>
  mutate(
    legend_label = paste0("BMU ", BMU_order, ": ", BMU_label)
  )

map_df_all <- map_df_all |>
  left_join(
    bmu_legend_df |>
      dplyr::select(BMU, BMU_order, legend_label),
    by = c("BMU", "BMU_order")
  ) |>
  mutate(
    month = factor(month, levels = months),
    legend_label = factor(
      legend_label,
      levels = bmu_legend_df$legend_label
    )
  )

wea_plot_df <- sf::st_as_sf(wea_plot)

# define palette
cb_palette <- c(
  "#0072B2",
  "#E69F00",
  "#009E73",
  "#D55E00",
  "#CC79A7",
  "#F0E442",
  "#56B4E9",
  "#8C8C8C",
  "#A6CEE3"
)

cb_palette <- viridisLite::viridis(9)
names(cb_palette) <- levels(map_df_all$legend_label)

p_bmu_12mo <- ggplot() +
  geom_raster(
    data = map_df_all,
    aes(x = x, y = y, fill = legend_label)
  ) +
  geom_sf(
    data = wea_plot_df,
    fill = NA,
    color = "black",
    linewidth = 0.25
  ) +
  scale_fill_manual(values = cb_palette) +
  facet_wrap(~month, ncol = 4) +
  coord_sf() +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7),
    axis.ticks = element_line(linewidth = 0.3),
    strip.background = element_rect(fill = "grey95"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 8)
  )

ggsave(
  filename = file.path(
    plot_dir,
    paste0("7. aggregated_species_BMU_map_", year_tag, "_12months.png")
  ),
  plot = p_bmu_12mo,
  width = 14,
  height = 10,
  dpi = 300
)
