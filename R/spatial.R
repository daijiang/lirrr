#' Make grid cells across the lands of world 
#' 
#' @param cell_size The size of the cell, default is 100 km by 100 km.
#' @param ... Additional arguments to be passed into `sf::st_make_grid()`.
#' @return A list of two sf objectives: `grid_sf` is the grid cells of the world; 
#' `grid_land` only contains the grid cells that intersected with terrestrial lands,
#' which can be further filtered by continent, or country, etc.
#' 
make_world_grids = function(cell_size = 100000, ...){
  grid <- sf::st_make_grid(world_moll, cellsize = cell_size, ...) # 100 km = 100,000 m
  grid_sf <- sf::st_sf(
    cell_id = paste0("cell_", 1:length(grid)),
    geometry = grid
  )
  
  # keep only grid cells that intersect land
  grid_land <- sf::st_intersection(grid_sf, world_moll)
  # grid_land <- sf::st_intersection(grid_sf, dplyr::filter(world_moll, ! continent %in% c("Antarctica")))
  # dplyr::count(st_drop_geometry(grid_land), continent)
  
  list(grid_sf, grid_land)
}


