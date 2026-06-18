library(tidyverse)

# 1. load data ----

surfclam_overlay_df <- readRDS(
  "../results/WEA_overlay/surfclam_HSI_WEA_cell_overlay_2025_all_months.rds"
) |>
  add_column(Species = "Atlantic Surfclam", .before = 1)

quahog_overlay_df <- readRDS(
  "../results/WEA_overlay/quahog_HSI_WEA_cell_overlay_2025_all_months.rds"
) |>
  add_column(Species = "Ocean Quahog", .before = 1)

scallop_overlay_df <- readRDS(
  "../results/WEA_overlay/scallop_HSI_WEA_cell_overlay_2025_all_months.rds"
) |>
  add_column(Species = "Atlantic Sea Scallop", .before = 1)

overlay_df <- bind_rows(
  surfclam_overlay_df,
  quahog_overlay_df,
  scallop_overlay_df
)
remove(surfclam_overlay_df, quahog_overlay_df, scallop_overlay_df)

# --------------------------------------------------------- #

# 2. bin and frequency ----

n_bins <- 20

hsi_breaks <- seq(0, 1, length.out = n_bins + 1)

overlay_df <- overlay_df |>
  mutate(
    HSI_bin = cut(HSI, breaks = hsi_breaks, include.lowest = TRUE, right = TRUE)
  )

hsi_freq <- overlay_df |>
  group_by(Species, HSI_bin, in_WEA) |>
  summarise(freq = n(), .groups = "drop")

# total frequency per HSI bin
hsi_total <- overlay_df |>
  group_by(Species, HSI_bin) |>
  summarise(total_freq = n(), .groups = "drop")

# merge for proportions later
hsi_freq <- left_join(hsi_freq, hsi_total, by = c("Species", "HSI_bin")) |>
  mutate(prop_within_bin = freq / total_freq)

# ----------------------------------------------------------------- #

# 3. Plot frequency histogram with WEA portion highlighted ----

# first summarize the proportion of cells by category (poor, fair, good)

h.status_freq_summary <- overlay_df |>
  mutate(
    H_Status = case_when(
      HSI <= 0.3 ~ "Poor",
      HSI > 0.3 & HSI <= 0.7 ~ "Fair",
      HSI > 0.7 ~ "Good",
      TRUE ~ NA_character_
    )
  ) |>
  group_by(Species, H_Status) |>
  summarize(
    n = n(),
    prop_in_WEA = sum(in_WEA == "Inside WEA", na.rm = TRUE) / n
  )

write.csv(
  h.status_freq_summary,
  "../plot/4. h.status_freq_summary.csv",
  row.names = FALSE
)

# stacked bars:
# total cells by HSI bin, with inside-WEA portion highlighted

HSI_WEA_p <- ggplot(hsi_freq, aes(x = HSI_bin, y = freq, fill = in_WEA)) +
  # geom_vline(xintercept = c("(0.25,0.3]", "(0.65,0.7]"), linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = c(5.5, 13.5), linetype = "dashed", color = "gray40") +
  geom_col(position = "stack") +
  labs(
    x = "Weighted HSI bin",
    y = "Number of raster cells",
    fill = NULL
  ) +
  facet_wrap(~Species, ncol = 1, scales = "free_y") +
  scale_fill_manual(
    values = c(
      "Inside WEA" = "indianred2",
      "Outside WEA" = "steelblue3"
    ),
    name = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave(
  "../plot/4. HSI_WEA_histogram.png",
  HSI_WEA_p,
  width = 6,
  height = 8,
  dpi = 300
)

# ----------------------------------------------------------------- #

# 4. plot proportions ----

# this shows the fraction of each HSI bin that falls inside WEA

HSI_WEA_proportion_p <- ggplot(
  hsi_freq,
  aes(x = HSI_bin, y = prop_within_bin, fill = in_WEA)
) +
  # geom_vline(xintercept = c("(0.25,0.3]", "(0.65,0.7]"), linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = c(5.5, 13.5), linetype = "dashed", color = "gray40") +
  geom_col(position = "fill") +
  labs(
    x = "Weighted HSI bin",
    y = "Number of raster cells",
    fill = NULL
  ) +
  facet_wrap(~Species, ncol = 1) +
  scale_fill_manual(
    values = c(
      "Inside WEA" = "indianred2",
      "Outside WEA" = "steelblue3"
    ),
    name = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

ggsave(
  "../plot/appendix/HSI_WEA_proportion.png",
  HSI_WEA_proportion_p,
  width = 6,
  height = 8,
  dpi = 300
)

# ----------------------------------------------------------------- #

# 5. Pie chart for habitat quality for cells inside WEA ----

pie_df <- overlay_df |>
  filter(in_WEA == "Inside WEA") |>
  mutate(
    H_Status = case_when(
      HSI <= 0.3 ~ "Poor",
      HSI > 0.3 & HSI <= 0.7 ~ "Fair",
      HSI > 0.7 ~ "Good",
      TRUE ~ NA_character_
    )
  ) |>
  group_by(Species, H_Status) |>
  summarize(n = n(), .groups = "drop") |>
  group_by(Species) |>
  mutate(prop = n / sum(n))

pie_df$H_Status <- factor(pie_df$H_Status, levels = c("Good", "Fair", "Poor"))

pie_p <- ggplot(pie_df, aes(x = "", y = prop, fill = H_Status)) +
  coord_polar(theta = "y") +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(
    values = c(
      "#DDEEE5", # Good (very light sage)
      "#7FAF9B", # Fair (muted green)
      "#2F5D50" # Poor (deep green)
    )
  ) +
  geom_text(
    aes(label = scales::percent(prop)),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  facet_wrap(~Species, ncol = 2) +
  theme_void() +
  labs(fill = "H Status")

ggsave(
  "../plot/appendix/HSI_WEA_PIE.png",
  pie_p,
  width = 6,
  height = 8,
  dpi = 300
)

# ----------------------------------------------------------------- #

# 6. Mean HSI inside vs outside WEA ----

# summarize mean HSI and 95% CI by WEA status
hsi_summary <- overlay_df %>%
  group_by(Species, in_WEA) %>%
  summarise(
    n = n(),
    mean_HSI = mean(HSI, na.rm = TRUE),
    sd_HSI = sd(HSI, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(
      # "n = ", n, "\n",
      "Mean = ",
      round(mean_HSI, 2),
      "\n",
      "SD = ",
      round(sd_HSI, 2)
    )
  )

write.csv(
  hsi_summary,
  "../plot/5. HSI_WEA_summary.csv",
  row.names = FALSE
)

# violin plot to show the distribution of HSI values inside vs outside WEA for each species
p_violin <- ggplot(overlay_df, aes(x = in_WEA, y = HSI, fill = in_WEA)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_boxplot(width = 0.12, outlier.size = 0.2, fill = "white") +
  geom_text(
    data = hsi_summary,
    aes(x = in_WEA, y = Inf, label = label),
    vjust = 1.2,
    size = 2,
    color = "black"
  ) +
  facet_wrap(~Species, ncol = 3) +
  scale_fill_manual(
    values = c(
      "Inside WEA" = "indianred2",
      "Outside WEA" = "steelblue3"
    ),
    name = NULL
  ) +
  labs(
    x = NULL,
    y = "Weighted HSI"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none"
  ) +
  coord_cartesian(ylim = c(0, 1))

ggsave(
  "../plot/5. HSI_WEA_violin.png",
  p_violin,
  width = 6,
  height = 8,
  dpi = 300
)

# wilcox test all siginificant
wilcox.test(
  HSI ~ in_WEA,
  data = subset(overlay_df, Species == "Atlantic Sea Scallop")
)

wilcox.test(
  HSI ~ in_WEA,
  data = subset(overlay_df, Species == "Atlantic Surfclam")
)

wilcox.test(
  HSI ~ in_WEA,
  data = subset(overlay_df, Species == "Ocean Quahog")
)

# # bar plot to show mean and sd
# p_mean_HSI <- ggplot(
#   hsi_summary,
#   aes(x = in_WEA, y = mean_HSI, fill = in_WEA)
# ) +
#   geom_col(width = 0.6) +
#   geom_errorbar(
#     aes(
#       ymin = mean_HSI - sd_HSI,
#       ymax = mean_HSI + sd_HSI
#     ),
#     width = 0.2
#   ) +
#   facet_wrap(~Species, ncol = 3) +
#   scale_fill_manual(
#     values = c(
#       "Inside WEA" = "steelblue3",
#       "Outside WEA" = "indianred2"
#     ),
#     name = NULL
#   ) +
#   labs(
#     x = NULL,
#     y = "Mean weighted HSI"
#   ) +
#   theme_minimal() +
#   theme(
#     legend.position = "none"
#   )
