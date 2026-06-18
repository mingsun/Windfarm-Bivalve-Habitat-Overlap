library(tidyverse)
library(sf)
library(terra)

# 1. Load weighted HSI raster ----

HSI_weighted <- rast("../results/HSI/surfclam/surfclam_HSI_weighted_1km.tif")
plot(HSI_weighted)

# ----------------------------------------------------------------- #

# 2. WEA shapefile ----

# 2.1 load WEA shapefile and check CRS
wea_sf <- st_read(
  "../data/boem-renewable-energy-shapefiles_2.5.2025/Offshore_Wind_Leases_outlines.shp",
  quiet = TRUE
)

# wea_sf <- st_read(
#   "../data/boem-renewable-energy-shapefiles_2.5.2025/Offshore_Wind_Leases.shp",
#   quiet = TRUE
# )

plot(st_geometry(wea_sf))
wea_sf <- st_make_valid(wea_sf) # fix invalid geometries first

# ----------------------------------------------------------------- #

# 3. crop WEA to HSI extent ----

# convert raster extent to polygon using terra
bbox_vect <- as.polygons(ext(HSI_weighted), crs = crs(HSI_weighted))

# convert to sf
bbox_poly <- st_as_sf(bbox_vect)
plot(bbox_poly)

# project crop box to WEA CRS for intersection
bbox_poly_wea <- st_transform(bbox_poly, st_crs(wea_sf))

# crop WEA polygons
wea_crop <- st_intersection(wea_sf, bbox_poly_wea)

plot(st_geometry(wea_crop))

# ----------------------------------------------------------------- #

# 4. Project cropped WEA to the HSI raster CRS ----

wea_hsi <- st_transform(wea_crop, crs(HSI_weighted))

# quick check
plot(st_geometry(wea_hsi))
plot(as.polygons(HSI_weighted), add = TRUE)

# ----------------------------------------------------------------- #

# 5. Rasterize WEA to the exact HSI grid ----

# cells inside WEA = 1
# cells outside WEA = NA initially, then convert to 0

wea_vect <- vect(wea_hsi)

WEA_raster <- rasterize(
  wea_vect,
  HSI_weighted,
  field = 1,
  background = NA
)

# convert outside-domain NA values to 0 only where HSI exists
WEA_raster <- ifel(!is.na(HSI_weighted) & is.na(WEA_raster), 0, WEA_raster)

names(WEA_raster) <- "in_WEA"

plot(WEA_raster)

# save
writeRaster(
  WEA_raster,
  "../results/WEA_overlay/surfclam_WEA_on_HSI_grid.tif",
  overwrite = TRUE
)

# ----------------------------------------------------------------- #

# 6. Create cell-by-cell dataframe of HSI and WEA overlap ----

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
    )
  )

summary(overlay_df$HSI)
table(overlay_df$in_WEA)

write.csv(
  overlay_df,
  "../results/WEA_overlay/surfclam_HSI_WEA_cell_overlay.csv",
  row.names = FALSE
)

# ----------------------------------------------------------------- #

remove(list = ls())

# 7. Bin HSI values and summarize frequency ----

overlay_df <- read.csv(
  "../results/WEA_overlay/surfclam_HSI_WEA_cell_overlay.csv"
)

# choose number of bins here
n_bins <- 20

hsi_breaks <- seq(0, 1, length.out = n_bins + 1)

overlay_df <- overlay_df |>
  mutate(
    HSI_bin = cut(HSI, breaks = hsi_breaks, include.lowest = TRUE, right = TRUE)
  )

hsi_freq <- overlay_df |>
  group_by(HSI_bin, in_WEA) |>
  summarise(freq = n(), .groups = "drop")

# total frequency per HSI bin
hsi_total <- overlay_df |>
  group_by(HSI_bin) |>
  summarise(total_freq = n(), .groups = "drop")

# merge if you want proportions later
hsi_freq <- left_join(hsi_freq, hsi_total, by = "HSI_bin") |>
  mutate(prop_within_bin = freq / total_freq)

hsi_freq

# ----------------------------------------------------------------- #

# 8. Plot frequency histogram with WEA portion highlighted ----

# stacked bars:
# total cells by HSI bin, with inside-WEA portion highlighted

ggplot(hsi_freq, aes(x = HSI_bin, y = freq, fill = in_WEA)) +
  geom_col(position = "stack") +
  labs(
    x = "Weighted HSI bin",
    y = "Number of raster cells",
    fill = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# ----------------------------------------------------------------- #

# 9. plot proportions instead of counts ----

# this shows the fraction of each HSI bin that falls inside WEA

ggplot(
  hsi_freq,
  aes(x = HSI_bin, y = prop_within_bin, fill = in_WEA)
) +
  geom_col(position = "fill") +
  labs(
    x = "Weighted HSI bin",
    y = "Proportion of cells",
    fill = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# ----------------------------------------------------------------- #

# 10. Mean HSI inside vs outside WEA ----

# summarize mean HSI and 95% CI by WEA status
hsi_summary <- overlay_df %>%
  group_by(in_WEA) %>%
  summarise(
    n = n(),
    mean_HSI = mean(HSI, na.rm = TRUE),
    sd_HSI = sd(HSI, na.rm = TRUE),
    se_HSI = sd_HSI / sqrt(n),
    ci_low = mean_HSI - 1.96 * se_HSI,
    ci_high = mean_HSI + 1.96 * se_HSI,
    .groups = "drop"
  )

hsi_summary

# bar plot with 95% CI
p_mean_ci <- ggplot(hsi_summary, aes(x = in_WEA, y = mean_HSI, fill = in_WEA)) +
  geom_col(width = 0.6) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    width = 0.15,
    linewidth = 0.7
  ) +
  labs(
    x = NULL,
    y = "Mean weighted HSI"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none"
  )

p_mean_ci

p_box <- ggplot(overlay_df, aes(x = in_WEA, y = HSI, fill = in_WEA)) +
  geom_boxplot(outlier.size = 0.3) +
  labs(
    x = NULL,
    y = "Weighted HSI"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none"
  )

p_box

p_violin <- ggplot(overlay_df, aes(x = in_WEA, y = HSI, fill = in_WEA)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_boxplot(width = 0.12, outlier.size = 0.2, fill = "white") +
  labs(
    x = NULL,
    y = "Weighted HSI"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none"
  )

p_violin
