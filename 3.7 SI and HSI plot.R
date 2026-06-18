library(tidyverse)
library(terra)
library(viridis)

# 1. set up loop variables ----

species_names <- c("surfclam", "quahog", "scallop")
species_titles <- c("Atlantic Surfclam", "Ocean Quahog", "Atlantic Sea Scallop")


# -------------------------------------------------------------------------- #

# 2. loop to plot SI ----
i = 1

for (i in seq_along(species_names)) {
  sp <- species_names[i]

  # load RDS (these preserve factor levels!)
  sed_plot <- readRDS(paste0(
    "../results/HSI/",
    sp,
    "/",
    sp,
    "_sediment_SI.rds"
  )) %>%
    rename(bin = sediment_bins, SI = SI_sediment) %>%
    mutate(var = "Sediment (phi)")

  depth_plot <- readRDS(paste0(
    "../results/HSI/",
    sp,
    "/",
    sp,
    "_depth_SI.rds"
  )) %>%
    rename(bin = depth_bins, SI = SI_depth) %>%
    mutate(var = "Depth (m)")

  temp_plot <- readRDS(paste0(
    "../results/HSI/",
    sp,
    "/",
    sp,
    "_temp_SI.rds"
  )) %>%
    rename(bin = temp_bins, SI = SI_temp) %>%
    mutate(var = "Temperature (°C)")

  plot_df <- bind_rows(sed_plot, depth_plot, temp_plot)

  #  # plot by bin
  #     p <- ggplot(plot_df, aes(x = bin, y = SI, group = 1)) +
  #     geom_line() +
  #     geom_point(size = 1.5, alpha = 0.8) +
  #     facet_wrap(~var, scales = "free_x", ncol = 1) +
  #     labs(
  #       title = species_titles[i],
  #       x = NULL,
  #       y = "Suitability index (0–1)"
  #     ) +
  #     theme_minimal() +
  #     theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # plot by mid point
  p <- ggplot(plot_df, aes(x = mid, y = SI, group = 1)) +
    geom_line() +
    geom_point(size = 1.5, alpha = 0.8) +
    scale_x_continuous(n.breaks = 10) +
    facet_wrap(~var, scales = "free_x", ncol = 1) +
    labs(
      title = species_titles[i],
      x = NULL,
      y = "Suitability index (0–1)"
    ) +
    theme_minimal()

  ggsave(
    filename = paste0("../plot/2. SI curve/", species_names[i], "_SI_plot.png"),
    plot = p,
    width = 3,
    height = 5,
    dpi = 300
  )
}

# -------------------------------------------------------------------------- #

# 3. plot HSI ----

## 3.1 All month domain  ----

months_all <- sprintf("%02d", 1:12)
month_labels <- c(
  "Jan",
  "Feb",
  "Mar",
  "Apr",
  "May",
  "Jun",
  "Jul",
  "Aug",
  "Sep",
  "Oct",
  "Nov",
  "Dec"
)

sp <- "surfclam"

for (sp in species_names) {
  files <- tibble(
    month = months_all,
    file = paste0(
      "../results/HSI/",
      sp,
      "/HSI_weighted/",
      sp,
      "_HSI_weighted_2025",
      months_all,
      ".tif"
    )
  )

  hsi_df <- map2_dfr(files$file, files$month, function(f, m) {
    r <- rast(f)
    df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
    names(df)[3] <- "HSI"
    df$month <- m
    df
  })

  extrap_df <- read.csv(
    paste0(
      "../results/HSI/",
      sp,
      "/HSI_weighted/temperature_extrapolation_summary_2025.csv"
    )
  ) %>%
    mutate(
      month = sprintf("%02d", month),
      extrap_total = pct_below + pct_above,
      label = paste0("Extrapolation: ", sprintf("%.2f", extrap_total), "%")
    )

  hsi_df <- hsi_df %>%
    left_join(extrap_df %>% select(month, label), by = "month") %>%
    mutate(month = factor(month, levels = months_all, labels = month_labels))

  label_df <- hsi_df %>%
    group_by(month) %>%
    summarize(
      x = min(x, na.rm = TRUE) +
        0.03 * (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)),
      y = max(y, na.rm = TRUE) -
        0.03 * (max(y, na.rm = TRUE) - min(y, na.rm = TRUE)),
      label = first(label),
      .groups = "drop"
    )

  p_appendix <- ggplot(hsi_df, aes(x = x, y = y, fill = HSI)) +
    geom_raster() +
    geom_text(
      data = label_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 3
    ) +
    facet_wrap(~month, ncol = 4) +
    coord_equal() +
    scale_fill_viridis_c(
      option = "C",
      limits = c(0, 1),
      na.value = "white",
      name = "HSI"
    ) +
    labs(
      title = species_titles[species_names == sp],
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right"
    )

  ggsave(
    paste0("../plot/appendix/", sp, "_HSI_2025_all_months.png"),
    p_appendix,
    width = 12,
    height = 9,
    dpi = 500
  )
}

## ---------------------------------------------------------------- ##

remove(list = ls())

## 3.2 surfclam example  ----

# august has the best data

data <- rast(
  "../results/HSI/surfclam/HSI_weighted/surfclam_HSI_weighted_202507.tif"
)

r_lonlat <- project(data, "EPSG:4326") # Reproject to lon/lat (EPSG:4326)

data <- as.data.frame(r_lonlat, xy = TRUE, na.rm = FALSE)
names(data)[3] <- "HSI"

p_surfclam <- ggplot(data, aes(x = x, y = y, fill = HSI)) +
  geom_raster() +
  coord_equal() +
  scale_fill_viridis_c(
    option = "C",
    limits = c(0, 1),
    na.value = NA,
    name = "HSI"
  ) +
  labs(
    title = "Atlantic Surfclam [August]",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    # panel.grid = element_blank(),
    # axis.text = element_blank(),
    # axis.ticks = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none"
  )

ggsave(
  paste0("../plot/3. HSI month example/surfclam_HSI_2025_August.png"),
  p_surfclam,
  width = 4,
  height = 4,
  dpi = 500
)

## ---------------------------------------------------------------- ##

remove(list = ls())

## 3.3 quahog example  ----

# august has the best data

data <- rast(
  "../results/HSI/quahog/HSI_weighted/quahog_HSI_weighted_202508.tif"
)

r_lonlat <- project(data, "EPSG:4326") # Reproject to lon/lat (EPSG:4326)

data <- as.data.frame(r_lonlat, xy = TRUE, na.rm = FALSE)
names(data)[3] <- "HSI"

p_surfclam <- ggplot(data, aes(x = x, y = y, fill = HSI)) +
  geom_raster() +
  coord_equal() +
  scale_fill_viridis_c(
    option = "C",
    limits = c(0, 1),
    na.value = NA,
    name = "HSI"
  ) +
  labs(
    title = "Ocean Quahog [August]",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    # panel.grid = element_blank(),
    # axis.text = element_blank(),
    # axis.ticks = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none"
  )

ggsave(
  paste0("../plot/3. HSI month example/quahog_HSI_2025_August.png"),
  p_surfclam,
  width = 4,
  height = 4,
  dpi = 500
)

## ---------------------------------------------------------------- ##

remove(list = ls())

## 3.4 scallop example  ----

# May has the best data

data <- rast(
  "../results/HSI/scallop/HSI_weighted/scallop_HSI_weighted_202505.tif"
)

r_lonlat <- project(data, "EPSG:4326") # Reproject to lon/lat (EPSG:4326)

data <- as.data.frame(r_lonlat, xy = TRUE, na.rm = FALSE)
names(data)[3] <- "HSI"

p_surfclam <- ggplot(data, aes(x = x, y = y, fill = HSI)) +
  geom_raster() +
  coord_equal() +
  scale_fill_viridis_c(
    option = "C",
    limits = c(0, 1),
    na.value = NA,
    name = "HSI"
  ) +
  labs(
    title = "Atlantic Sea Scallop [May]",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    # panel.grid = element_blank(),
    # axis.text = element_blank(),
    # axis.ticks = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none"
  )

ggsave(
  paste0("../plot/3. HSI month example/scallop_HSI_2025_May.png"),
  p_surfclam,
  width = 4,
  height = 4,
  dpi = 500
)

## ---------------------------------------------------------------- ##

# -------------------------------------------------------------------------- #
