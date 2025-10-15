# 1. Define a global equal-area CRS (Mollweide projection)
library(sf)
library(rnaturalearth)

# 1. Get real world polygon
world <- ne_countries(scale = "medium", returnclass = "sf")

# 2. Define Mollweide CRS (ESRI:54009 is standard)
crs_equal_area <- "ESRI:54009"

# 3. Reproject world to equal-area
world_moll <- st_transform(world, crs_equal_area)
world_moll <- st_simplify(lirrr::world_moll, dTolerance = 1000)  # meters

# 4. Check bbox
st_bbox(world_moll)
plot(st_geometry(world_moll))


usethis::use_data(world_moll, overwrite = TRUE)
