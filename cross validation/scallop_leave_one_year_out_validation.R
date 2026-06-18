library(tidyverse)
library(classInt)
library(quantreg)
library(gbm)
library(dismo)

source("functions/fit_bell_SI.R")

# 1. Case settings ----

species_name <- "scallop"

tow_file <- "../data/scallop data/SS_tows_MA_final.csv"
out_dir <- "../results/HSI/scallop"

abund_col <- "NperMsq"
year_col <- "YEAR"

env_cols <- c("SED_PHI_MEAN", "DEPTH", "TEMP")

n_bins <- 20

# BRT settings: match final HSI fitting script
brt_tree_complexity <- 3
brt_learning_rate <- 0.008
brt_bag_fraction <- 0.75

# ----------------------------------------------------------------- #

# 2. Load data ----

data_all <- read.csv(tow_file) |>
  filter(YEAR >= 2013) |>
  filter(
    !is.na(.data[[abund_col]]),
    !is.na(TEMP),
    !is.na(DEPTH),
    !is.na(SED_PHI_MEAN),
    !is.na(.data[[year_col]])
  )

years <- sort(unique(data_all[[year_col]]))

# ----------------------------------------------------------------- #

# 3. Function: build SI curve from training data ----

build_si_curve <- function(train_df, env_col, abund_col, si_col_name) {
  x <- as.numeric(train_df[[env_col]])

  # protect against rare cases with too few unique values
  n_bins_use <- min(n_bins, length(unique(x)) - 1)

  brks <- classIntervals(x, n_bins_use, style = "fisher")$brks
  brks[1] <- brks[1] - .Machine$double.eps

  train_tmp <- train_df |>
    mutate(
      env_bin = cut(.data[[env_col]], breaks = brks, include.lowest = TRUE)
    )

  si_df <- train_tmp |>
    group_by(env_bin) |>
    summarize(
      abundance = mean(.data[[abund_col]], na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      !!si_col_name := {
        mn <- min(abundance, na.rm = TRUE)
        mx <- max(abundance, na.rm = TRUE)
        if (mx == mn) {
          rep(0, length(abundance))
        } else {
          (abundance - mn) / (mx - mn)
        }
      },
      mid = (brks[-1] + brks[-length(brks)])[as.integer(env_bin)] / 2
    )

  # fit bell-shaped SI curve
  bell_fit <- fit_bell_SI(
    SI_df = si_df,
    SI_col = !!rlang::sym(si_col_name)
  )

  # predict smoothed SI and rescale to 0-1
  si_df <- si_df |>
    mutate(
      !!si_col_name := predict(bell_fit, newdata = si_df)
    ) |>
    mutate(
      !!si_col_name := {
        z <- .data[[si_col_name]]
        mn <- min(z, na.rm = TRUE)
        mx <- max(z, na.rm = TRUE)
        if (mx == mn) {
          rep(0, length(z))
        } else {
          (z - mn) / (mx - mn)
        }
      }
    ) |>
    dplyr::select(mid, all_of(si_col_name)) |>
    arrange(mid)

  si_df
}

# ----------------------------------------------------------------- #

# 4. Function: predict SI from SI curve ----

predict_si <- function(x, si_df, si_col) {
  approx(
    x = si_df$mid,
    y = si_df[[si_col]],
    xout = x,
    method = "linear",
    rule = 2
  )$y
}

# ----------------------------------------------------------------- #

# 5. Function: calculate BRT weights from training data ----

calc_brt_weights <- function(train_df, abund_col) {
  set.seed(123)

  brt_fit <- gbm.step(
    data = train_df,
    gbm.x = env_cols,
    gbm.y = abund_col,
    family = "gaussian",
    tree.complexity = brt_tree_complexity,
    learning.rate = brt_learning_rate,
    bag.fraction = brt_bag_fraction
  )

  brt_fit$contributions |>
    as.data.frame() |>
    transmute(
      var = var,
      weight = rel.inf / sum(rel.inf)
    )
}

# ----------------------------------------------------------------- #

# 6. Function: calculate validation metrics ----

calc_validation_metrics <- function(
  df,
  hsi_col = "HSI_leave_one_year",
  abund_col = "NperMsq"
) {
  df <- df |>
    filter(!is.na(.data[[hsi_col]]), !is.na(.data[[abund_col]]))

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

  # raw Spearman correlation
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

  # quantile regression slope, tau = 0.9
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

  # simple linear regression
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

# 7. Leave-one-year-out HSI prediction ----

cv_list <- vector("list", length(years))

for (i in seq_along(years)) {
  hold_year <- years[i]

  message("Holding out year: ", hold_year)

  train_df <- data_all |>
    filter(.data[[year_col]] != hold_year)

  test_df <- data_all |>
    filter(.data[[year_col]] == hold_year)

  # build SI curves using training years only
  sediment_SI <- build_si_curve(
    train_df = train_df,
    env_col = "SED_PHI_MEAN",
    abund_col = abund_col,
    si_col_name = "SI_sediment"
  )

  depth_SI <- build_si_curve(
    train_df = train_df,
    env_col = "DEPTH",
    abund_col = abund_col,
    si_col_name = "SI_depth"
  )

  temp_SI <- build_si_curve(
    train_df = train_df,
    env_col = "TEMP",
    abund_col = abund_col,
    si_col_name = "SI_temp"
  )

  # refit BRT weights using training years only
  brt_weights <- calc_brt_weights(
    train_df = train_df,
    abund_col = abund_col
  )

  w_sed <- brt_weights$weight[brt_weights$var == "SED_PHI_MEAN"]
  w_dep <- brt_weights$weight[brt_weights$var == "DEPTH"]
  w_tmp <- brt_weights$weight[brt_weights$var == "TEMP"]

  # calculate HSI for held-out year
  cv_list[[i]] <- test_df |>
    mutate(
      SI_sediment_cv = predict_si(
        SED_PHI_MEAN,
        sediment_SI,
        "SI_sediment"
      ),
      SI_depth_cv = predict_si(
        DEPTH,
        depth_SI,
        "SI_depth"
      ),
      SI_temp_cv = predict_si(
        TEMP,
        temp_SI,
        "SI_temp"
      ),
      HSI_leave_one_year = w_sed *
        SI_sediment_cv +
        w_dep * SI_depth_cv +
        w_tmp * SI_temp_cv,
      held_out_year = hold_year
    )
}

cv_df <- bind_rows(cv_list)

# ----------------------------------------------------------------- #

# 8. Final leave-one-year-out validation metrics ----

leave_one_year_validation_metrics <- calc_validation_metrics(
  cv_df,
  hsi_col = "HSI_leave_one_year",
  abund_col = abund_col
) |>
  mutate(
    species = species_name,
    validation = "leave-one-year-out"
  ) |>
  dplyr::select(species, validation, everything())

leave_one_year_validation_metrics

write.csv(
  leave_one_year_validation_metrics,
  file.path(out_dir, "scallop_leave_one_year_validation_metrics.csv"),
  row.names = FALSE
)
