library(tidyverse)
library(sf)
library(concaveman)
library(gstat)
library(terra)
library(units)

# use dredge survey full station to define a spatial domain for interpolation of sediment data

# 1. Read tow data ----
station_df <- read.csv("../data/22565_UNION_FSCS_SVSTA.csv") |>
  filter(DECDEG_BEGLON <= -69, DECDEG_BEGLAT <= 41.7) |>
  mutate(YEAR = substr(CRUISE6, 1, 4)) |>
  select(YEAR, LON = DECDEG_BEGLON, LAT = DECDEG_BEGLAT) # drop rows with missing coordinates

plot(LAT ~ LON, station_df)

# Convert to sf points (WGS84)
station_sf <- st_as_sf(
  station_df,
  coords = c("LON", "LAT"),
  crs = 4326,
  remove = FALSE
)

# ---------------------------------------------------------------------------------------------------------- #

# 2. Choose a projected CRS (for 1 km grid + variogram/kriging) ----

#    Use UTM zone based on tow centroid

cent <- st_coordinates(st_centroid(st_union(station_sf)))
lon0 <- cent[1]
lat0 <- cent[2]
utm_zone <- floor((lon0 + 180) / 6) + 1
epsg_utm <- if (lat0 >= 0) 32600 + utm_zone else 32700 + utm_zone

message(sprintf(
  "Using projected CRS EPSG:%s (UTM zone %s).",
  epsg_utm,
  utm_zone
))

station_utm <- st_transform(station_sf, crs = epsg_utm)

# ---------------------------------------------------------------------------------------------------------- #

# 3. Build an irregular survey-domain polygon from tow points ----

#    - concave hull, then buffer to include nearby sediment points and avoid edge artifacts

# Create a single multipoint geometry
mp <- st_union(station_utm)

pts <- st_cast(mp, "POINT", warn = FALSE)
dom <- concaveman::concaveman(pts)

# Buffer domain (meters). Adjust if you want tighter/looser.
# This helps kriging near edges and includes sediment points just outside the hull.
dom_buf <- st_buffer(dom, dist = 2 * 1000) # n km buffer

# Visual check of domain and tow points
plot(st_geometry(dom_buf), col = "lightblue")
plot(st_geometry(station_utm), add = TRUE, pch = 16, cex = 0.3)

# Save domain polygon for inspection
st_write(
  st_as_sf(dom_buf),
  "../data/sediment data/NEFSC_survey_domain_polygon.gpkg",
  delete_dsn = TRUE,
  quiet = TRUE
)

st_crs(dom_buf)$epsg

# ---------------------------------------------------------------------------------------------------------- #

# 4. Read sediment points (shapefile) ----

sed_sf <- st_read(
  "../data/sediment data/ecstdb2014/ecstdb2014.shp",
  quiet = TRUE
)

# Ensure it's in WGS84 then project to UTM
# (If already has CRS, st_transform will handle it.)
sed_utm <- st_transform(sed_sf, crs = epsg_utm)

# filter sediment data and remove duplicated locations
sed_utm <- sed_utm %>%
  mutate(MEAN = as.numeric(MEAN)) %>%
  filter(is.finite(MEAN) & MEAN != -9999) |> # remove the -9999 invalid values
  group_by(LONGITUDE, LATITUDE) %>% # group by duplicated location
  arrange(desc(YEAR_COLL), STDEV, .by_group = TRUE) %>% # newest year first; within year keep smaller STDEV
  slice(1) %>% # keep the best record per location
  ungroup()

# Spatial crop to buffered domain (faster variogram/kriging)
sed_utm_crop <- sed_utm[st_intersects(sed_utm, dom_buf, sparse = FALSE), ]

plot(st_geometry(dom_buf), col = "lightblue")
plot(st_geometry(sed_utm_crop), add = TRUE, pch = 16, cex = 0.3) # check cropped points

nrow(sed_utm_crop) # need a minimal of 50 for stable variogram/kriging

# ---------------------------------------------------------------------------------------------------------- #

# 5. Ordinary kriging with spherical variogram (choose best among a few) ----

# gstat prefers sp objects for variogram/kriging
sed_sp <- as(sed_utm_crop, "Spatial")

# Empirical variogram
v_emp <- variogram(MEAN ~ 1, data = sed_sp)

# Fit candidate models; choose lowest SSErr
models <- c("Sph", "Exp", "Gau")
fits <- lapply(models, function(m) {
  fit.variogram(
    v_emp,
    vgm(
      psill = var(sed_utm_crop$MEAN, na.rm = TRUE),
      model = m,
      range = 100000, # 100 km starting guess
      nugget = 0
    )
  )
})

ss <- sapply(fits, function(f) attr(f, "SSErr"))
best_i <- which.min(ss)
best_model <- fits[[best_i]]
message(sprintf(
  "Best variogram model by SSErr: %s (SSErr=%.3f)",
  models[best_i],
  ss[best_i]
))

# ------------------------------------------------------------------------------------------ #

# 6. Build prediction grid INSIDE the irregular domain polygon ----

# Make grid of points (centers) within buffered domain
grid_pts <- st_make_grid(dom_buf, cellsize = 1 * 1000, what = "centers") %>% # n km grid
  st_as_sf() %>%
  st_filter(dom_buf, .predicate = st_within)

# Convert to Spatial for gstat predict
grid_sp <- as(grid_pts, "Spatial")
plot(grid_sp, pch = 16, cex = 0.3)

# Ordinary kriging predictions
kr <- krige(
  formula = MEAN ~ 1,
  locations = sed_sp,
  newdata = grid_sp,
  model = best_model,
  nmax = 50 # limit neighbors for stability/speed; adjust as needed
)

summary(kr$var1.pred) # quick summary
sum(is.na(kr$var1.pred))

# ------------------------------------------------------------------------------------------ #

# 7. Convert kriged points to raster + extract to tow points ----
kr_sf <- st_as_sf(kr)

# Visual check of kriged predictions
plot(kr_sf["var1.pred"])
plot(st_geometry(sed_utm_crop), add = TRUE, pch = 16, cex = 0.3)

# Rasterize using terra
kr_vect <- terra::vect(kr_sf)
# Create an empty raster at same resolution over domain extent as in step 6
r_template <- terra::rast(
  ext(kr_vect),
  resolution = 1000,
  crs = terra::crs(kr_vect)
)
r_phi <- terra::rasterize(kr_vect, r_template, field = "var1.pred", fun = mean)

# Mask raster to the buffered domain polygon (keeps irregular shape)
dom_vect <- terra::vect(st_as_sf(dom_buf))
r_phi_mask <- terra::mask(r_phi, dom_vect)

# Visual check of masked raster and tow points
plot(r_phi_mask)
plot(st_geometry(sed_utm_crop), add = TRUE, pch = 16, cex = 0.2)

# Save raster
terra::writeRaster(
  r_phi_mask,
  "../data/sediment data/kriged.sediment/sediment_phi_kriged_1km.tif",
  overwrite = TRUE
)
