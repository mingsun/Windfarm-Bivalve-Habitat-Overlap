library(tidyverse)
library(sf)
library(gstat)
library(terra)
library(units)

# ============================================================
# Tow-level sediment (phi) via kriging from USGS ecstdb2014
# - Builds an irregular "survey domain" polygon from tow points
# - Kriges MEAN (phi) onto 1 km grid inside domain
# - Extracts kriged values back to each tow
# ============================================================

# 1) Read tow data ----

tows <- read.csv("../data/surfclam data/SC_tows_MA.csv") |>
  filter(!is.na(LAT) & !is.na(LON)) # drop rows with missing coordinates

# Convert to sf points (WGS84)
tows_sf <- st_as_sf(
  tows,
  coords = c("LON", "LAT"),
  crs = 4326,
  remove = FALSE
)

# ------------------------------------------------------------------------------------------ #

# 2) Choose a projected CRS (for 1 km grid + variogram/kriging) ----
#    Use UTM zone based on tow centroid

cent <- st_coordinates(st_centroid(st_union(tows_sf)))
lon0 <- cent[1]
lat0 <- cent[2]
utm_zone <- floor((lon0 + 180) / 6) + 1
epsg_utm <- if (lat0 >= 0) 32600 + utm_zone else 32700 + utm_zone

message(sprintf(
  "Using projected CRS EPSG:%s (UTM zone %s).",
  epsg_utm,
  utm_zone
))

remove(cent, lon0, lat0, utm_zone)

tows_utm <- st_transform(tows_sf, crs = epsg_utm)

# ------------------------------------------------------------------------------------------ #

# 3. assign sediment to tow ----

# load kriged raster for the survey domain
sediment_phi <- rast(
  "../data/sediment data/kriged.sediment/sediment_phi_kriged_1km.tif"
)

# Extract to tow points
tows_vect <- terra::vect(tows_utm)
phi_at_tows <- terra::extract(sediment_phi, tows_vect)[, 2] # second column is extracted values

# Attach and write output CSV
tows_out <- tows %>%
  mutate(SED_PHI_MEAN = phi_at_tows) |>
  filter(!is.na(SED_PHI_MEAN) & !is.na(TEMP))

write.csv(
  tows_out,
  "../data/surfclam data/SC_tows_MA_final.csv",
  row.names = FALSE
)

# ------------------------------------------------------------------------------------------ #
