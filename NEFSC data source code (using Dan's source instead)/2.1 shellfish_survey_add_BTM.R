library(tidyverse)
# library(ncdf4)
library(terra)

# variabls needed for habitat modeling include: latitude, bottom temperature, and depth
# bottom temperature need to be interpolated to the station locations, so we will use kriging for that

cat.df <- read.csv("../results/tow data/Shellfis_catch_by_tow.csv")

# mannually add depth for missing entries based on geographic location
cat.df <- cat.df |>
  mutate(AVGDEPTH = ifelse(!is.na(AVGDEPTH), AVGDEPTH, 45)) |>
  rename(LAT = DECDEG_BEGLAT, LON = DECDEG_BEGLON)


# no way to krige the BTM data because of the large number of missing values
na_by_year <- cat.df |>
  group_by(YEAR) |>
  summarise(
    n_na = sum(is.na(BOTTEMP)),
    n_available = sum(!is.na(BOTTEMP)),
    .groups = "drop"
  )

remove(na_by_year)

# test <- nc_open(
#   "../data/temperature data/mercatorglorys12v1_gl12_mean_202504.nc"
# )
# test <- rast(
#   "../data/temperature data/mercatorglorys12v1_gl12_mean_202504.nc",
#   subds = "thetao"
# )

## =========================
## 0) User settings (edit)
## =========================
data_dir <- "E:/Copernicus/" # folder where .nc files are
file_template <- "GLORYS_thetao_%s.nc" # %s will be CRUISE6, e.g., 198204
subds_name <- "thetao" # Copernicus temp variable
lon_col <- "LON"
lat_col <- "LAT"
cruise_col <- "CRUISE6"
out_col <- "btm_temp_C"

test_year <- 1982 # test with one year first

## =========================
## 1) Preparation
## =========================

# Ensure CRUISE6 is character (so sprintf works cleanly)
cat.df[[cruise_col]] <- as.character(cat.df[[cruise_col]])

# Create output column
cat.df[[out_col]] <- NA_real_

# Build spatial points once (assumes LON/LAT in EPSG:4326)
pts_all <- vect(cat.df, geom = c(lon_col, lat_col), crs = "EPSG:4326")

# Filter CRUISE6 values for the year you want to test
# (CRUISE6 is like "198204")
cruise_vals_all <- sort(unique(cat.df[[cruise_col]]))
cruise_vals <- cruise_vals_all[
  substr(cruise_vals_all, 1, 4) == as.character(test_year)
]

cruise_vals
# You should see something like: "198201" "198202" ... for months that exist in cat.df

## =========================
## 2) Loop over CRUISE6 (YYYYMM)
## =========================
for (cc in cruise_vals) {
  # Build file path for this year-month
  f <- file.path(data_dir, sprintf(file_template, cc))

  if (!file.exists(f)) {
    warning(sprintf("File not found for CRUISE6=%s: %s (skipping)", cc, f))
    next
  }

  message("Processing CRUISE6 = ", cc, " | file = ", f)

  # Load temperature cube for that month (usually 50 depth layers)
  r <- rast(f, subds = subds_name)

  # Compute bottom temperature raster:
  # per grid cell, take deepest non-NA across depth layers
  bottom_r <- app(r, fun = function(v) {
    i <- which(!is.na(v))
    if (length(i) == 0) {
      return(NA_real_)
    }
    v[max(i)]
  })

  # Extract for rows in cat.df that match this CRUISE6
  idx <- which(cat.df[[cruise_col]] == cc)
  bt <- extract(bottom_r, pts_all[idx])

  vals <- bt[, 2]

  # Auto-convert Kelvin -> Celsius if needed
  medv <- suppressWarnings(median(vals, na.rm = TRUE))
  if (is.finite(medv) && medv > 100) {
    vals <- vals - 273.15
  }

  # Fill back into cat.df
  cat.df[[out_col]][idx] <- vals

  # Cleanup
  rm(r, bottom_r, bt, vals)
  gc()
}

## =========================
## 3) Quick checks
## =========================
summary(cat.df[[out_col]])
sum(is.na(cat.df[[out_col]]))
