library(terra)
library(ncdf4)
library(gsw)
library(tidyverse)

# 1. Load GLORYS NetCDF and inspect variables ----

glorys_file <- "../data/temperature data/mercatorglorys12v1_gl12_mean_202508.nc "

# open with ncdf4, so you can inspect dimensions and variable names more clearly
nc <- nc_open(glorys_file)

# print summary of the file
print(nc)

# list variable names only
names(nc$var)

# close connection when done inspecting
nc_close(nc)

# ---------------------------------------------------------------------------------------------------------- #

# 2. Extract bottom potential temperature from GLORYS

# load only bottomT
bottomT <- rast(glorys_file, subds = "bottomT")

# basic checks
bottomT
names(bottomT)
crs(bottomT)
ext(bottomT)
res(bottomT)
minmax(bottomT)

# replace weird values with NA
values(bottomT)[!is.finite(values(bottomT))] <- NA

# quick plot
plot(bottomT)


# ---------------------------------------------------------------------------------------------------------- #

# 3. Project bottom potential temperature to the survey domain ----

# load the sediment raster you already created
sediment_phi <- rast(
  "../data/sediment data/kriged.sediment/sediment_phi_kriged_1km.tif"
)

# project GLORYS bottom potential temperature to the same
# CRS / extent / resolution as sediment raster
bottomT_match <- project(
  bottomT,
  sediment_phi,
  method = "bilinear"
)

# project GLORYS bottom potential temperature to the same
# CRS / extent / resolution as sediment raster
bottomT_resampled <- resample(
  bottomT_match,
  sediment_phi,
  method = "bilinear"
)

# quick checks: these should match the sediment raster
crs(bottomT_resampled)
res(bottomT_resampled)
ext(bottomT_resampled)

crs(sediment_phi)
res(sediment_phi)
ext(sediment_phi)

# plot for visual check
plot(bottomT_match)
plot(bottomT_resampled)

# ---------------------------------------------------------------------------------------------------------- #

# 4. Extract bottom salinity for temperature conversion ----

# load the full 3D salinity variable
so <- rast(glorys_file, subds = "so")

# inspect layer names / count
so
names(so)
nlyr(so)

# function: return deepest non-NA value in the vertical profile
deepest_non_na <- function(x) {
  idx <- max(which(!is.na(x)))
  if (is.finite(idx)) {
    return(x[idx])
  } else {
    return(NA_real_)
  }
}

# apply across layers to get bottom salinity
so_bottom <- app(so, deepest_non_na)

# check result
plot(so_bottom)
summary(values(so_bottom))

# ---------------------------------------------------------------------------------------------------------- #

# 5. Project bottom salinity to the survey domain ----

# project bottom salinity to the same CRS / extent / resolution
# as the sediment raster
so_bottom_match <- project(
  so_bottom,
  sediment_phi,
  method = "bilinear"
)

# optional extra resample step to mirror the earlier workflow
so_bottom_resampled <- resample(
  so_bottom_match,
  sediment_phi,
  method = "bilinear"
)

# quick checks: these should match sediment raster
crs(so_bottom_resampled)
res(so_bottom_resampled)
ext(so_bottom_resampled)

crs(sediment_phi)
res(sediment_phi)
ext(sediment_phi)

# visual check
plot(so_bottom_resampled)

# summary check
summary(values(so_bottom_resampled))

# ---------------------------------------------------------------------------------------------------------- #

# 6. Convert GLORYS bottom potential temperature to bottom in-situ temperature using TEOS-10 ----

# load final processed depth raster from the depth script
depth_resampled <- rast("../data/bathymetry data/depth_gebco_1km_domain.tif")
plot(depth_resampled)

# convert the aligned rasters to a dataframe with x/y coordinates
# x = longitude/easting, y = latitude/northing in the current raster CRS
temp_df <- as.data.frame(bottomT_resampled, xy = TRUE, na.rm = FALSE)
salt_df <- as.data.frame(so_bottom_resampled, xy = TRUE, na.rm = FALSE)
depth_df <- as.data.frame(depth_resampled, xy = TRUE, na.rm = FALSE)

# combine into one dataframe
# the third column in each dataframe is the raster value column
tsd_df <- data.frame(
  x = temp_df$x,
  y = temp_df$y,
  pt0 = temp_df[, 3], # bottom potential temperature
  SP = salt_df[, 3], # bottom salinity (Practical Salinity)
  depth_m = depth_df[, 3] # bottom depth in meters (positive)
)

# keep only complete rows
tsd_df <- tsd_df %>%
  filter(!is.na(pt0) & !is.na(SP) & !is.na(depth_m))

# convert raster coordinates to lon/lat (required for TEOS-10)
pts <- vect(tsd_df[, c("x", "y")], crs = crs(bottomT_resampled))
pts_ll <- project(pts, "EPSG:4326")
ll <- crds(pts_ll)

tsd_df$lon <- ll[, 1]
tsd_df$lat <- ll[, 2]

# compute pressure from depth (TEOS-10 uses negative z)
tsd_df$p <- gsw_p_from_z(
  z = -tsd_df$depth_m,
  latitude = tsd_df$lat
)

# compute pressure from depth (TEOS-10 uses negative z)
tsd_df$SA <- gsw_SA_from_SP(
  SP = tsd_df$SP,
  p = tsd_df$p,
  longitude = tsd_df$lon,
  latitude = tsd_df$lat
)

# convert potential temperature to Conservative Temperature
tsd_df$CT <- gsw_CT_from_pt(
  SA = tsd_df$SA,
  pt = tsd_df$pt0
)

# convert Conservative Temperature to in-situ temperature
tsd_df$t_insitu <- gsw_t_from_CT(
  SA = tsd_df$SA,
  CT = tsd_df$CT,
  p = tsd_df$p
)

# quick check
summary(tsd_df$t_insitu)

# write in-situ temperature values back to raster
bottomT_insitu <- sediment_phi
values(bottomT_insitu) <- NA_real_

cell_id <- cellFromXY(bottomT_insitu, tsd_df[, c("x", "y")])
values(bottomT_insitu)[cell_id] <- tsd_df$t_insitu

# visual check
plot(bottomT_insitu)
res(bottomT_insitu)

terra::writeRaster(
  bottomT_insitu,
  "../data/temperature data/bottomT_insitu_1km_domain.tif",
  overwrite = TRUE
)

# ---------------------------------------------------------------------------------------------------------- #

# a quick comparison

catch_by_tow_df <- read.csv("../data/surfclam data/SC_tows_MA_final.csv") |>
  filter(YR >= 2013) |> # filter out
  select(YR, LON, LAT, TEMP) |>
  filter(!is.na(TEMP))

summary(catch_by_tow_df$TEMP)
summary(values(bottomT_insitu), na.rm = TRUE)

plot(LAT ~ LON, catch_by_tow_df, color = TEMP)
