library(sf)
library(dplyr)
library(ggplot2)

# ---------------------------
# 1. Convert to sf
# ---------------------------
pts <- st_as_sf(distinct(geomix_setup$lattice_coords,xid,yid), coords = c("xid", "yid"), crs = NA)

# ---------------------------
# 2. Create hex grid
# ---------------------------
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

# Keep only hexes that intersect points
hex_grid <- hex_grid[st_intersects(hex_grid, pts, sparse = FALSE) %>% apply(1, any), ]

# ---------------------------
# 3. Assign points to hexes
# ---------------------------
pts_joined <- st_join(pts, hex_grid, join = st_within)

# Count points per hex
hex_counts <- pts_joined %>%
  st_drop_geometry() %>%
  group_by(grid_id) %>%
  summarise(n = n())

hex_grid <- left_join(hex_grid, hex_counts, by = "grid_id")
hex_grid$n[is.na(hex_grid$n)] <- 0

# ---------------------------
# 4. Compute threshold
# ---------------------------
avg_n <- mean(hex_grid$n)
min_allowed <- avg_n / 2

# ---------------------------
# 5. Merge small hexes with neighbors
# ---------------------------
# Find neighbors
neighbors <- st_touches(hex_grid)

# Iteratively merge small cells
repeat {
  small_cells <- which(hex_grid$n < min_allowed)

  if (length(small_cells) == 0) break

  for (i in small_cells) {
    neigh <- neighbors[[i]]

    if (length(neigh) == 0) next

    # Pick neighbor with largest count
    best_neighbor <- neigh[which.max(hex_grid$n[neigh])]

    # Merge geometries
    hex_grid$geometry[best_neighbor] <- st_union(
      hex_grid$geometry[best_neighbor],
      hex_grid$geometry[i]
    )

    # Update count
    hex_grid$n[best_neighbor] <- hex_grid$n[best_neighbor] + hex_grid$n[i]

    # Mark current cell as removed
    hex_grid$n[i] <- NA
  }

  # Remove merged cells
  hex_grid <- hex_grid[!is.na(hex_grid$n), ]
  hex_grid$grid_id <- seq_len(nrow(hex_grid))
  neighbors <- st_touches(hex_grid)
}

# ---------------------------
# 6. Join back to lattice points
# ---------------------------
pts_joined_final <- st_join(pts, hex_grid, join = st_within)

pts_with_grid <- geomix_setup$lattice_coords %>%
  left_join(
    pts_joined_final %>% st_drop_geometry() %>% distinct(xid, yid, grid_id),
    by = c("xid", "yid")
  )

# ---------------------------
# 7. Save result
# ---------------------------
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
saveRDS(pts_with_grid, "data/processed/points_grid.rds")
