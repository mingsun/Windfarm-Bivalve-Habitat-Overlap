library(tidyverse)

# 1. settings and load data ----

species_cfg <- tribble(
  ~species   , ~overlay_rds                                                               , ~survey_file                                 ,
  "surfclam" , "../results/WEA_overlay/surfclam_HSI_WEA_cell_overlay_2025_all_months.rds" , "../data/surfclam data/SC_tows_MA_final.csv" ,
  "quahog"   , "../results/WEA_overlay/quahog_HSI_WEA_cell_overlay_2025_all_months.rds"   , "../data/quahog data/OQ_tows_MA_final.csv"   ,
  "scallop"  , "../results/WEA_overlay/scallop_HSI_WEA_cell_overlay_2025_all_months.rds"  , "../data/scallop data/SS_tows_MA_final.csv"
)

# bootstrap settings
n_boot <- 2000
set.seed(100)

# output directory
out_dir <- "../results/WEA_effect_size"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 2. build functions ----

## ---- 2.1 derive survey year from common column names ----
get_survey_year <- function(df) {
  nms <- names(df)

  if ("CRUISE6" %in% nms) {
    return(substr(as.character(df$CRUISE6), 1, 4))
  } else if ("YEAR" %in% nms) {
    return(as.character(df$YEAR))
  } else if ("YR" %in% nms) {
    return(as.character(df$YR))
  } else {
    stop("Could not find a year column. Need CRUISE6, YEAR, or YR.")
  }
}

## -------------------------------------------- ##

## ---- 2.2 average annual number of survey tows ----
get_avg_annual_tows <- function(survey_file) {
  survey_df <- read.csv(survey_file)

  survey_df$survey_year <- get_survey_year(survey_df)

  annual_counts <- survey_df %>%
    filter(!is.na(survey_year)) %>%
    count(survey_year, name = "n_tows")

  mean(annual_counts$n_tows, na.rm = TRUE)
}

## -------------------------------------------- ##

## ---- 2.3 Hedges' g ----
hedges_g <- function(x_in, x_out) {
  x_in <- x_in[is.finite(x_in)]
  x_out <- x_out[is.finite(x_out)]

  n1 <- length(x_in)
  n2 <- length(x_out)

  if (n1 < 2 || n2 < 2) {
    return(NA_real_)
  }

  m1 <- mean(x_in)
  m2 <- mean(x_out)
  s1 <- stats::sd(x_in)
  s2 <- stats::sd(x_out)

  sp <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))

  if (!is.finite(sp) || sp == 0) {
    return(NA_real_)
  }

  d <- (m1 - m2) / sp

  # small-sample correction
  J <- 1 - 3 / (4 * (n1 + n2) - 9)

  g <- J * d
  g
}

## -------------------------------------------- ##

## ---- 2.4 bootstrap Hedges' g using simulated survey sample size ----
bootstrap_hedges_g <- function(x_in, x_out, n_in, n_out, n_boot = 2000) {
  x_in <- x_in[is.finite(x_in)]
  x_out <- x_out[is.finite(x_out)]

  if (length(x_in) < 2 || length(x_out) < 2) {
    return(tibble(
      g_boot = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_
    ))
  }

  g_vals <- replicate(n_boot, {
    boot_in <- sample(x_in, size = n_in, replace = TRUE)
    boot_out <- sample(x_out, size = n_out, replace = TRUE)
    hedges_g(boot_in, boot_out)
  })

  tibble(
    g_boot = mean(g_vals, na.rm = TRUE),
    ci_low = quantile(g_vals, probs = 0.025, na.rm = TRUE),
    ci_high = quantile(g_vals, probs = 0.975, na.rm = TRUE)
  )
}

## -------------------------------------------- ##

## ---- 2.5 summarize full-cell means and SDs by month + WEA group ----
summarise_group_stats <- function(df) {
  df %>%
    group_by(month, in_WEA) %>%
    summarise(
      n_cells = n(),
      mean_HSI = mean(HSI, na.rm = TRUE),
      sd_HSI = sd(HSI, na.rm = TRUE),
      .groups = "drop"
    )
}

## -------------------------------------------- ##

## ---- 2.6 monthly effect size calculation for one species ----
calc_monthly_effect_size <- function(
  overlay_df_all,
  avg_annual_tows,
  species_name,
  n_boot = 2000
) {
  # standardize month format
  overlay_df_all <- overlay_df_all %>%
    mutate(
      month = sprintf("%02d", as.integer(month))
    )

  # static inside/outside proportion based on unique cells
  # this is the allocation ratio used for the simulated survey
  cell_df <- overlay_df_all %>%
    distinct(x, y, in_WEA)

  prop_inside <- mean(cell_df$in_WEA == "Inside WEA", na.rm = TRUE)
  prop_outside <- mean(cell_df$in_WEA == "Outside WEA", na.rm = TRUE)

  n_total <- round(avg_annual_tows)
  n_in <- max(1, round(n_total * prop_inside))
  n_out <- max(1, n_total - n_in)

  # group stats from all cells
  group_stats <- summarise_group_stats(overlay_df_all) %>%
    mutate(species = species_name)

  # bootstrap effect size by month
  effect_df <- overlay_df_all %>%
    group_by(month) %>%
    group_modify(
      ~ {
        x_in <- .x %>% filter(in_WEA == "Inside WEA") %>% pull(HSI)
        x_out <- .x %>% filter(in_WEA == "Outside WEA") %>% pull(HSI)

        boot_res <- bootstrap_hedges_g(
          x_in = x_in,
          x_out = x_out,
          n_in = n_in,
          n_out = n_out,
          n_boot = n_boot
        )

        tibble(
          n_survey_total = n_total,
          n_survey_inside = n_in,
          n_survey_outside = n_out,
          prop_cells_inside = prop_inside,
          prop_cells_outside = prop_outside,
          hedges_g = boot_res$g_boot,
          ci_low = boot_res$ci_low,
          ci_high = boot_res$ci_high
        )
      }
    ) %>%
    ungroup() %>%
    mutate(species = species_name)

  list(
    effect_df = effect_df,
    group_stats = group_stats
  )
}

## -------------------------------------------- ##

# ----------------------------------------------------------------- #

# 3. run for all species ----

all_effect_list <- list()
all_groupstats_list <- list()

for (i in seq_len(nrow(species_cfg))) {
  sp_name <- species_cfg$species[i]
  overlay_file <- species_cfg$overlay_rds[i]
  survey_file <- species_cfg$survey_file[i]

  cat("Processing species:", sp_name, "\n")

  overlay_df_all <- readRDS(overlay_file)

  avg_annual_tows <- get_avg_annual_tows(survey_file)

  cat("Average annual tows:", avg_annual_tows, "\n")

  res <- calc_monthly_effect_size(
    overlay_df_all = overlay_df_all,
    avg_annual_tows = avg_annual_tows,
    species_name = sp_name,
    n_boot = n_boot
  )

  all_effect_list[[sp_name]] <- res$effect_df
  all_groupstats_list[[sp_name]] <- res$group_stats
}

effect_size_df <- bind_rows(all_effect_list)
group_stats_df <- bind_rows(all_groupstats_list)

# ----------------------------------------------------------------- #

# 4. Reshape group stats to wide format (inside vs outside) ----

group_stats_wide <- group_stats_df %>%
  mutate(in_WEA = ifelse(in_WEA == "Inside WEA", "inside", "outside")) %>%
  pivot_wider(
    names_from = in_WEA,
    values_from = c(n_cells, mean_HSI, sd_HSI)
  ) %>%
  arrange(species, month)

# combine effect size + group summaries
final_effect_table <- effect_size_df %>%
  left_join(group_stats_wide, by = c("species", "month")) %>%
  arrange(species, month)

# save
write.csv(
  final_effect_table,
  file.path(out_dir, "monthly_HSI_effect_size_hedges_g.csv"),
  row.names = FALSE
)

write.csv(
  group_stats_df,
  file.path(out_dir, "monthly_HSI_group_stats_long.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------------- #

remove(list = ls())

# 5. plot ----

out_dir <- "../results/WEA_effect_size"
final_effect_table <- read.csv(file.path(
  out_dir,
  "monthly_HSI_effect_size_hedges_g.csv"
))

group_stats_df <- read.csv(file.path(
  out_dir,
  "monthly_HSI_group_stats_long.csv"
))

## 5.1 Hedges' g by month, faceted by species ----

plot_df <- final_effect_table %>%
  mutate(
    month_num = as.integer(month),
    month_lab = factor(
      month_num,
      levels = 1:12,
      labels = c(
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
    ),
    # reverse order so Jan is at top and Dec at bottom
    month_lab = factor(month_lab, levels = rev(levels(month_lab))),
    effect_dir = ifelse(
      hedges_g >= 0,
      "Higher HSI inside WEA",
      "Higher HSI outside WEA"
    ),
    species = recode(
      species,
      "surfclam" = "Atlantic Surfclam",
      "quahog" = "Ocean Quahog",
      "scallop" = "Atlantic Sea Scallop"
    )
  )

p_g_month <- ggplot(plot_df, aes(x = hedges_g, y = month_lab)) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey40"
  ) +
  geom_errorbarh(
    aes(xmin = ci_low, xmax = ci_high, color = effect_dir),
    height = 0.18,
    linewidth = 0.8
  ) +
  geom_point(aes(color = effect_dir), size = 3.8) +
  facet_wrap(~species, ncol = 3) +
  scale_color_manual(
    values = c(
      "Higher HSI inside WEA" = "indianred2",
      "Higher HSI outside WEA" = "steelblue3"
    ),
    name = NULL
  ) +
  labs(x = "Hedges' g", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold", size = 13),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(xlim = c(-1.2, 1.2))

ggsave(
  "../plot/6. monthly_HSI_hedges_g_by_species.png",
  p_g_month,
  width = 8,
  height = 7,
  dpi = 300
)

## -------------------------------------------- ##

# 5.2  Mean HSI with SD by month ----

plot_mean_sd <- group_stats_df %>%
  mutate(
    month_num = as.integer(month),
    month_lab = factor(
      month_num,
      levels = 1:12,
      labels = c(
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
    ),
    month_lab = factor(month_lab, levels = rev(levels(month_lab))),
    species = recode(
      species,
      "surfclam" = "Atlantic Surfclam",
      "quahog" = "Ocean Quahog",
      "scallop" = "Atlantic Sea Scallop"
    ),
    sd_low = mean_HSI - sd_HSI,
    sd_high = mean_HSI + sd_HSI
  )

pd <- position_dodge(width = 0.5)

p_mean_sd <- ggplot(
  plot_mean_sd,
  aes(x = mean_HSI, y = month_lab, color = in_WEA)
) +
  geom_errorbarh(
    aes(xmin = sd_low, xmax = sd_high),
    height = 0.18,
    linewidth = 0.8,
    position = pd
  ) +
  geom_point(
    size = 3.5,
    position = pd
  ) +
  facet_wrap(~species, ncol = 3) +
  scale_color_manual(
    values = c(
      "Inside WEA" = "indianred2",
      "Outside WEA" = "steelblue3"
    ),
    name = NULL
  ) +
  labs(x = "Mean HSI ± 1 SD", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold", size = 13),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(xlim = c(0, 1))

ggsave(
  "../plot/appendix/monthly_HSI_mean_sd_by_species.png",
  p_mean_sd,
  width = 8,
  height = 8,
  dpi = 300
)

# ----------------------------------------------------------------- #
