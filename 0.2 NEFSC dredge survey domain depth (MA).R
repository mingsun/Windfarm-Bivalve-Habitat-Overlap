library(tidyverse)
library(sf)
library(terra)
# use dredge survey full station to define a spatial domain for interpolation of sediment data

# 1. Load survey domain polygon ----

dom_buf <- st_read(
  "../data/sediment data/NEFSC_survey_domain_polygon.gpkg",
  quiet = TRUE
)

st_crs(dom_buf)$epsg

# ---------------------------------------------------------------------------------------------------------- #

# 2. Load sediment raster as master template ----

sediment_phi <- rast(
  "../data/sediment data/kriged.sediment/sediment_phi_kriged_1km.tif"
)

# check template properties
crs(sediment_phi)
res(sediment_phi)
ext(sediment_phi)

# ---------------------------------------------------------------------------------------------------------- #

# 3. Load GEBCO bathymetry raster ----

gebco <- rast(
  "../data/bathymetry data/GEBCO_06_Mar_2026_81bed68aa2ac/gebco_2025_n42.5061_s33.9807_w-77.2119_e-68.9062.nc"
)

# GEBCO uses elevation:
# negative = below sea level (water depth)
# positive = land elevation
names(gebco)
plot(gebco)

# convert to positive depth
# keep only ocean cells; land becomes NA
depth <- ifel(gebco < 0, -gebco, NA)

# ---------------------------------------------------------------------------------------------------------- #

# 4. Project bathymetry to the same CRS as the survey domain ----

depth_match <- project(
  depth,
  sediment_phi,
  method = "bilinear"
)


# ---------------------------------------------------------------------------------------------------------- #

# 5. Mask depth raster to the irregular survey domain ----

dom_vect <- vect(dom_buf)

depth_mask <- mask(
  depth_match,
  dom_vect
)

# ---------------------------------------------------------------------------------------------------------- #

# 6. Check that depth and sediment match ----

compareGeom(sediment_phi, depth_mask, stopOnError = FALSE)
crs(depth_mask)
res(depth_mask)
ext(depth_mask)

# visual check
plot(depth_mask)
plot(st_geometry(dom_buf), add = TRUE, border = "black")

writeRaster(
  depth_mask,
  "../data/bathymetry data/depth_gebco_1km_domain.tif",
  overwrite = TRUE
)

# ---------------------------------------------------------------------------------------------------------- #

# 7. compare with the dredge survey value ----

station_df <- read.csv("../data/22565_UNION_FSCS_SVSTA.csv") |>
  filter(DECDEG_BEGLON <= -69, DECDEG_BEGLAT <= 41.7) |>
  mutate(YEAR = substr(CRUISE6, 1, 4)) |>
  select(YEAR, LON = DECDEG_BEGLON, LAT = DECDEG_BEGLAT, DEPTH = AVGDEPTH) |>
  filter(!is.na(DEPTH)) # drop rows with missing coordinates or depth

plot(LAT ~ LON, station_df)

## 7.1 survey stations as spatial points ----
station_sf <- station_df %>%
  st_as_sf(coords = c("LON", "LAT"), crs = 4326, remove = FALSE)

# 7.2 project stations to the same CRS as the final depth raster ----
station_sf <- st_transform(station_sf, crs(depth_mask))

# 7.3. extract raster depth at station locations ----
station_vect <- vect(station_sf)

rast_depth <- extract(depth_mask, station_vect)

# 7.4 attach extracted raster depth back to dataframe ----
station_compare <- station_sf %>%
  mutate(DEPTH_RASTER = rast_depth[, 2]) %>%
  st_drop_geometry()

# 7.5 quick checks ----
summary(station_compare$DEPTH)
summary(station_compare$DEPTH_RASTER)

# 7.6. scatterplot comparison ----
plot(
  station_compare$DEPTH,
  station_compare$DEPTH_RASTER,
  xlab = "Survey depth",
  ylab = "Raster depth"
)
abline(0, 1, col = "red", lwd = 2)

# 7.7 correlation ----
cor(
  station_compare$DEPTH,
  station_compare$DEPTH_RASTER,
  use = "complete.obs"
)

# 7.8 difference ----
station_compare <- station_compare %>%
  mutate(DEPTH_DIFF = DEPTH_RASTER - DEPTH)

summary(station_compare$DEPTH_DIFF)

# ---------------------------------------------------------------------------------------------------------- #
