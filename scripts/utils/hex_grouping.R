#%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1 Preliminaries ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%

# Requires `lattice_coords`: a data frame with columns xid and yid,
# one row per lattice location. Set this before sourcing.

library(sf)
library(dplyr)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2 Build hexagonal grouping grid ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 2.1 Convert lattice coordinates to sf points ----
pts <- st_as_sf(distinct(lattice_coords, xid, yid), coords = c("xid", "yid"), crs = NA)

## 2.2 Create hexagonal grid over lattice extent ----
cell_size <- 7.5

hex_grid <- st_make_grid(
  pts,
  cellsize = cell_size,
  square = FALSE
)

hex_grid <- st_sf(
  grid_id = 1:length(hex_grid),
  geometry = hex_grid
)

# Keep only hexes that intersect at least one lattice point
hex_grid <- hex_grid[st_intersects(hex_grid, pts, sparse = FALSE) %>% apply(1, any), ]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3 Assign points to hexes ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 3.1 Spatial join ----
pts_joined <- st_join(pts, hex_grid, join = st_within)

## 3.2 Count points per hex ----
hex_counts <- pts_joined %>%
  st_drop_geometry() %>%
  group_by(grid_id) %>%
  summarise(n = n(), .groups = "drop")

hex_grid <- left_join(hex_grid, hex_counts, by = "grid_id")
hex_grid$n[is.na(hex_grid$n)] <- 0

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4 Merge small hexes ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 4.1 Compute minimum allowed count ----
avg_n <- mean(hex_grid$n)
min_allowed <- avg_n / 2

## 4.2 Iteratively merge undersize cells with largest neighbour ----
neighbors <- st_touches(hex_grid)

repeat {
  small_cells <- which(hex_grid$n < min_allowed)

  if (length(small_cells) == 0) break

  for (i in small_cells) {
    neigh <- neighbors[[i]]

    if (length(neigh) == 0) next

    best_neighbor <- neigh[which.max(hex_grid$n[neigh])]

    hex_grid$geometry[best_neighbor] <- st_union(
      hex_grid$geometry[best_neighbor],
      hex_grid$geometry[i]
    )

    hex_grid$n[best_neighbor] <- hex_grid$n[best_neighbor] + hex_grid$n[i]
    hex_grid$n[i] <- NA
  }

  hex_grid <- hex_grid[!is.na(hex_grid$n), ]
  hex_grid$grid_id <- seq_len(nrow(hex_grid))
  neighbors <- st_touches(hex_grid)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5 Join grid IDs back to lattice ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pts_joined_final <- st_join(pts, hex_grid, join = st_within)

pts_with_grid <- lattice_coords %>%
  left_join(
    pts_joined_final %>% st_drop_geometry() %>% distinct(xid, yid, grid_id),
    by = c("xid", "yid")
  )

#%%%%%%%%%%%%%%%%%%%%%%%%%%
# 6 Save result ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
saveRDS(pts_with_grid, "data/processed/points_grid.rds")
message("Saved data/processed/points_grid.rds")
