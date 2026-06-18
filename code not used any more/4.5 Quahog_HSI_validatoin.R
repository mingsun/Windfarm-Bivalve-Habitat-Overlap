library(tidyverse)
library(terra)
library(sf)
library(quantreg)

# 1. Case settings ----

species_name <- "quahog"

survey_file <- "../data/quahog data/OQ_tows_MA_final.csv"

hsi_dir <- "../results/HSI/quahog/HSI_weighted"
hsi_pattern <- "quahog_HSI_weighted_2025"

# survey columns
lon_col <- "LON"
lat_col <- "LAT"
abund_col <- "NperMsq"

# survey data
survey_df <- read.csv(survey_file) |>
  filter(YR >= 2013)


# number of HSI bins for validation
n_bins <- 10

# months used for validation
months <- sprintf("%02d", 7:9)

# ----------------------------------------------------------------- #

# 2. Load monthly HSI rasters ----

hsi_files <- paste0(
  hsi_dir,
  "/",
  hsi_pattern,
  months,
  ".tif"
)

hsi_list <- lapply(hsi_files, rast)
names(hsi_list) <- months

# template raster
hsi_template <- hsi_list[[1]]

# convert survey to sf and then terra vect for spatial operations
survey_sf <- st_as_sf(
  survey_df,
  coords = c(lon_col, lat_col),
  crs = 4326,
  remove = FALSE
)

survey_sf <- st_transform(survey_sf, crs(hsi_template))
survey_vect <- vect(survey_sf)

# ----------------------------------------------------------------- #

# 3. Helper function: calculate validation metrics ----

calc_validation_metrics <- function(
  df,
  hsi_col = "HSI",
  abund_col = "NperMsq",
  n_bins = 10
) {
  df <- df |>
    filter(!is.na(.data[[hsi_col]]), !is.na(.data[[abund_col]]))

  # if too few rows, return NA
  if (nrow(df) < 10) {
    return(
      tibble(
        n = nrow(df),
        binned_spearman = NA_real_,
        binned_spearman_p = NA_real_,
        quantile_slope_tau_0.9 = NA_real_,
        quantile_slope_p = NA_real_,
        variance_explained_bins = NA_real_
      )
    )
  }

  ## 3.1 Binned monotonic correlation ----

  # use equal-width HSI bins from 0 to 1
  hsi_breaks <- seq(0, 1, length.out = n_bins + 1)

  df_bin <- df |>
    mutate(
      HSI_bin = cut(
        .data[[hsi_col]],
        breaks = hsi_breaks,
        include.lowest = TRUE,
        right = TRUE
      )
    ) |>
    filter(!is.na(HSI_bin))

  bin_summary <- df_bin |>
    group_by(HSI_bin) |>
    summarise(
      mean_abundance = mean(.data[[abund_col]], na.rm = TRUE),
      mid_HSI = mean(.data[[hsi_col]], na.rm = TRUE),
      n_bin = n(),
      .groups = "drop"
    ) |>
    filter(!is.na(mean_abundance), !is.na(mid_HSI))

  if (nrow(bin_summary) >= 3) {
    cor_out <- suppressWarnings(
      cor.test(
        bin_summary$mid_HSI,
        bin_summary$mean_abundance,
        method = "spearman",
        exact = FALSE
      )
    )
    binned_spearman <- unname(cor_out$estimate)
    binned_spearman_p <- cor_out$p.value
  } else {
    binned_spearman <- NA_real_
    binned_spearman_p <- NA_real_
  }

  ## -------------------------------------------- ##

  ## 3.2 Quantile regression slope (tau = 0.9) ----

  rq_fit <- tryCatch(
    quantreg::rq(
      formula = as.formula(paste0(abund_col, " ~ ", hsi_col)),
      tau = 0.9,
      data = df
    ),
    error = function(e) NULL
  )

  if (!is.null(rq_fit)) {
    rq_sum <- tryCatch(summary(rq_fit, se = "nid"), error = function(e) NULL)

    if (!is.null(rq_sum)) {
      quantile_slope <- coef(rq_fit)[2]
      quantile_slope_p <- rq_sum$coefficients[2, 4]
    } else {
      quantile_slope <- coef(rq_fit)[2]
      quantile_slope_p <- NA_real_
    }
  } else {
    quantile_slope <- NA_real_
    quantile_slope_p <- NA_real_
  }

  ## -------------------------------------------- ##

  ## 3.3 Variance explained in bins ----
  # between-bin SS / total SS

  grand_mean <- mean(df_bin[[abund_col]], na.rm = TRUE)

  bin_ss <- df_bin |>
    group_by(HSI_bin) |>
    summarise(
      n_bin = n(),
      mean_abundance = mean(.data[[abund_col]], na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      ss_between = n_bin * (mean_abundance - grand_mean)^2
    )

  ss_between <- sum(bin_ss$ss_between, na.rm = TRUE)
  ss_total <- sum((df_bin[[abund_col]] - grand_mean)^2, na.rm = TRUE)

  variance_explained_bins <- ifelse(
    ss_total > 0,
    ss_between / ss_total,
    NA_real_
  )

  ## -------------------------------------------- ##

  ## 3.4  return ----

  tibble(
    n = nrow(df),
    binned_spearman = binned_spearman,
    binned_spearman_p = binned_spearman_p,
    quantile_slope_tau_0.9 = quantile_slope,
    quantile_slope_p = quantile_slope_p,
    variance_explained_bins = variance_explained_bins
  )
}

## -------------------------------------------- ##

# ----------------------------------------------------------------- #

# 4. Extract monthly HSI to survey locations ----

validation_list <- vector("list", length(months))

for (i in seq_along(months)) {
  month_id <- months[i]
  hsi_rast <- hsi_list[[i]]

  cat("Processing month:", month_id, "\n")

  # extract monthly HSI at the same survey locations
  ext_df <- terra::extract(hsi_rast, survey_vect)

  validation_list[[i]] <- survey_df |>
    mutate(
      month = month_id,
      HSI = ext_df[, 2]
    ) |>
    filter(!is.na(HSI))
}

validation_df <- bind_rows(validation_list)

# quick checks
validation_df |>
  count(month)

summary(validation_df$HSI)
summary(validation_df[[abund_col]])

# ----------------------------------------------------------------- #

# 5. Final validation metrics monthly ----

validation_list <- vector("list", length(months))

for (i in seq_along(months)) {
  month_id <- months[i]
  hsi_rast <- hsi_list[[i]]

  cat("Processing month:", month_id, "\n")

  # extract monthly HSI at the same survey locations
  ext_df <- terra::extract(hsi_rast, survey_vect)

  validation_list[[i]] <- survey_df |>
    mutate(
      month = month_id,
      HSI = ext_df[, 2]
    ) |>
    filter(!is.na(HSI))
}

validation_df <- bind_rows(validation_list)

# quick checks
validation_df |>
  count(month)

summary(validation_df$HSI)
summary(validation_df[[abund_col]])

# ----------------------------------------------------------------- #

# 7. Monthly validation metrics ----

validation_metrics_monthly <- validation_df |>
  group_by(month) |>
  group_modify(
    ~ calc_validation_metrics(
      .x,
      hsi_col = "HSI",
      abund_col = abund_col,
      n_bins = n_bins
    )
  ) |>
  ungroup() |>
  mutate(
    species = species_name,
    season = "summer monthly"
  ) |>
  dplyr::select(species, season, month, everything())

validation_metrics_monthly

# ----------------------------------------------------------------- #

# 6. save ----

write.csv(
  validation_metrics_monthly,
  "../results/HSI/quahog/quahog_validation_metrics_summer_by_month.csv",
  row.names = FALSE
)


saveRDS(
  validation_df,
  "../results/HSI/quahog/quahog_validation_df_summer_by_month.rds"
)
# ----------------------------------------------------------------- #
