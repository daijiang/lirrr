#' Make grid cells across the lands of world 
#' 
#' @param cell_size The size of the cell, default is 100 km by 100 km.
#' @param exclude_Antarctica Exclude Antarctica in the land grid cells? Default is TRUE.
#' @param ... Additional arguments to be passed into `sf::st_make_grid()`.
#' @return A list of two sf objectives: `grid_sf` is the grid cells of the world; 
#' `grid_land` only contains the grid cells that intersected with terrestrial lands.
#' @export
#' 
make_world_grids = function(cell_size = 100000, exclude_Antarctica = TRUE, ...){
  world_moll = lirrr::world_moll
  grid <- sf::st_make_grid(world_moll, cellsize = cell_size) # 100 km = 100,000 m
  grid_sf <- sf::st_sf(
    cell_id = paste0("cell_", 1:length(grid)),
    geometry = grid
  )
  
  # keep only grid cells that intersect land
  if(exclude_Antarctica){
    world_moll = dplyr::filter(world_moll, continent != "Antarctica")
  }
  idx <- sf::st_intersects(grid_sf, world_moll, sparse = TRUE)
  has_land <- lengths(idx) > 0
  grid_land <- grid_sf[has_land, ]
  grid_land = unique(grid_land)
  
  # grid_land <- sf::st_intersection(grid_sf, dplyr::filter(world_moll, ! continent %in% c("Antarctica")))
  # dplyr::count(st_drop_geometry(grid_land), continent)
  
  list(grids_world = grid_sf, grids_land = grid_land)
}


