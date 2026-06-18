library(tidyverse)
library(quantreg)

# 1. Case settings ----

species_name <- "quahog"

tow_file <- "../results/HSI/quahog/OQ_tows_with_composite_SI.csv"

# columns
abund_col <- "NperMsq"
hsi_col_name <- "HSI_weighted"

# tow-level data with fitted HSI already included
tow_df <- read.csv(tow_file) |>
  filter(YR >= 2013)

# ----------------------------------------------------------------- #

# 2. Helper function: calculate validation metrics using raw tow-level values ----

calc_validation_metrics <- function(
  df,
  hsi_col = "HSI_weighted",
  abund_col = "NperMsq"
) {
  df <- df |>
    filter(!is.na(.data[[hsi_col]]), !is.na(.data[[abund_col]]))

  # if too few rows, return NA
  if (nrow(df) < 10) {
    return(
      tibble(
        n = nrow(df),
        spearman = NA_real_,
        spearman_p = NA_real_,
        quantile_slope_tau_0.9 = NA_real_,
        quantile_slope_p = NA_real_,
        lm_slope = NA_real_,
        lm_slope_p = NA_real_,
        lm_r2 = NA_real_
      )
    )
  }

  ## 2.1 Raw Spearman correlation ----
  cor_out <- suppressWarnings(
    cor.test(
      df[[hsi_col]],
      df[[abund_col]],
      method = "spearman",
      exact = FALSE
    )
  )

  spearman <- unname(cor_out$estimate)
  spearman_p <- cor_out$p.value

  ## 2.2 Quantile regression slope (tau = 0.9) ----
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

  ## 2.3 Simple linear regression ----
  lm_fit <- tryCatch(
    lm(
      formula = as.formula(paste0(hsi_col, " ~ ", abund_col)),
      data = df
    ),
    error = function(e) NULL
  )

  if (!is.null(lm_fit)) {
    lm_sum <- summary(lm_fit)
    lm_slope <- coef(lm_fit)[2]
    lm_slope_p <- lm_sum$coefficients[2, 4]
    lm_r2 <- lm_sum$r.squared
  } else {
    lm_slope <- NA_real_
    lm_slope_p <- NA_real_
    lm_r2 <- NA_real_
  }

  ## 2.4 return ----
  tibble(
    n = nrow(df),
    spearman = spearman,
    spearman_p = spearman_p,
    quantile_slope_tau_0.9 = quantile_slope,
    quantile_slope_p = quantile_slope_p,
    lm_slope = lm_slope,
    lm_slope_p = lm_slope_p,
    lm_r2 = lm_r2
  )
}

# ----------------------------------------------------------------- #

# 3. Build validation data directly from tow-level fitted HSI ----

validation_df <- tow_df |>
  filter(!is.na(.data[[hsi_col_name]]), !is.na(.data[[abund_col]]))

# quick checks
summary(validation_df[[hsi_col_name]])
summary(validation_df[[abund_col]])

# ----------------------------------------------------------------- #

# 4. Validation metrics for all tows ----

validation_metrics <- calc_validation_metrics(
  validation_df,
  hsi_col = hsi_col_name,
  abund_col = abund_col
) |>
  mutate(
    species = species_name,
    season = "all tows"
  ) |>
  dplyr::select(species, season, everything())

validation_metrics

# ----------------------------------------------------------------- #

# 5. save ----

write.csv(
  validation_metrics,
  "../results/HSI/quahog/quahog_validation_metrics.csv",
  row.names = FALSE
)

# ----------------------------------------------------------------- #
