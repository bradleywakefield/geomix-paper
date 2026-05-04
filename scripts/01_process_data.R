#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1 Preliminaries ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 1.1 Libraries ----
library(matrixStats)
library(tidyverse)
library(reshape2)
library(nimble)
library(ggh4x)
library(tictoc)
library(scales)
library(sf)
select <- dplyr::select

## 1.2 Helper functions ----
box_cox <- function(y,lambda = 0.6) (y ^ lambda - 1) / lambda
inv_box_cox <- function(y, lambda=0.6) (y * lambda + 1)^(1 / lambda)
su_colours <-c("yellow", "forestgreen", "maroon", "darkorange3", "darkblue", "red4", "gray30",
               "dodgerblue4", "goldenrod3", "darkolivegreen", "firebrick", "mediumvioletred",
               "chocolate4", "darkslateblue", "slategray4", "seagreen4")

soil_units <- c("GT1","GT2","GT2c","GT3","GT4","GT5","GT5*","GT6")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2 Load raw data ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 2.1 CPT profiles ----
ijv_10cm <- readRDS("data/application/cpt_profiles.rds") %>%
  rename(x = Easting_m, y = Northing_m,
         bathymetry5m = Bathymetry_UHR_5m_LAT,
         bathymetryp5m = Bathymetry_MBES_0p5m_LAT,
         depthBSF = Depth_m_bsf,
         qc = qn_MPa) %>%
  mutate(x = x*100, y = y*100) %>%
  mutate(d = round(depthBSF - bathymetryp5m,1))

## 2.2 CPT stratigraphy annotations ----
sample_wsynth <- readRDS("data/application/cpt_stratigraphy.rds")

## 2.3 Seismic CDP data ----
synth_df <- readRDS("data/application/seismic_cdp.rds")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3 Spatial grid and locations ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 3.1 CDP spatial locations ----
cdp_loctions <- distinct(synth_df,easting,northing) %>%
  mutate(id = 1:n(), Z2_flag = 0,
         easting = easting*100,
         northing = northing*100)

locs <- cdp_loctions %>%
  rename(x=easting,y=northing) %>%
  mutate() %>%
  bind_rows(mutate(distinct(ijv_10cm,x,y),Z2_flag=1,id = 1:n()))
ncpt <- nrow(distinct(ijv_10cm,x,y))

saveRDS(locs,'data/processed/locs.rds')

## 3.2 Build 35 km spatial grid (local CRS) ----
cell_size <- 35000

loc_sf <- st_as_sf(locs,coords=c("x","y"))
cell_grid <- st_make_grid(loc_sf,
                          cellsize = cell_size,
                          square = T) %>%
  st_sf(loc_id = 1:length(.)) %>%
  mutate(centre = st_centroid(.))

coords <- st_coordinates(cell_grid$centre)
dimx <- length(unique(coords[,1]))
dimy <- length(unique(coords[,2]))

cell_grid$xid <- rep(1:dimx,dimy)
cell_grid$yid <- rep(1:dimy,each=dimx)

saveRDS(cell_grid,'data/processed/cell_grid.rds')

## 3.3 Assign locations to grid cells ----
selected_points <- st_join(loc_sf, cell_grid, join = st_within) %>%
  mutate(dist_to_centre = st_distance(geometry, centre, by_element = TRUE)) %>%
  group_by(loc_id) %>%
  arrange(Z2_flag != 1, dist_to_centre) %>%
  slice(1) %>%
  ungroup()

which(!(1:ncpt %in% filter(selected_points,Z2_flag==1)$id))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4 Process geophysics (CDP) data ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

geophys <- synth_df %>%
  filter(depth <= 50) %>%
  mutate(easting = easting*100, northing = northing*100) %>%
  left_join(rename(cdp_loctions,cdp_id = id)) %>%
  filter(cdp_id %in% filter(selected_points,Z2_flag==0)$id) %>%
  select(name = id,cdp_id,x=easting,y=northing,d=depth,soilUnitID,soilUnit) %>%
  left_join(st_set_geometry(select(filter(selected_points,Z2_flag==0),cdp_id=id,loc_id,xid,yid),value=NULL)) %>%
  mutate(cdp_id = cdp_id+100000) %>% relocate(loc_id,xid,yid,.after = cdp_id) %>%
  drop_na(loc_id)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5 Process CPT data ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 5.1 Join stratigraphy to CPT measurements ----
cpt_data <- ijv_10cm %>% left_join(
  sample_wsynth %>%
    select(CPT_name=name,depthBSF=depth,soilUnitID,soilUnit)
) %>% filter(d <= 50) %>%
  left_join(mutate(distinct(ijv_10cm,x,y),cdp_id = 1:n())) %>%
  filter(cdp_id %in% filter(selected_points,Z2_flag==1)$id) %>%
  select(name = CPT_name,cdp_id,x,y,d,qc,soilUnitID,soilUnit) %>%
  mutate(bc_qc = box_cox(qc),.after = qc) %>%
  left_join(st_set_geometry(select(filter(selected_points,Z2_flag==1),cdp_id=id,loc_id,xid,yid),value=NULL)) %>%
  relocate(loc_id,xid,yid,.after = name)

## 5.2 Extend CPT profiles to full depth ----
full_depths <- sample_wsynth %>%
  left_join(select(ijv_10cm,name=CPT_name,depth = depthBSF,d)) %>%
  select(name,depth = d) %>%
  group_by(name) %>%
  complete(depth = seq(min(depth),50,0.1)) %>%
  ungroup() %>% rename(d=depth)

missing <- anti_join(rename(full_depths,CPT_name=name),
          select(ijv_10cm,CPT_name,d)) %>%
  left_join(mutate(distinct(ijv_10cm,CPT_name,x,y,bathymetry5m,bathymetryp5m),cdp_id = 1:n())) %>%
  left_join(distinct(sample_wsynth,CPT_name =name,easting=eastingS,northing=northingS)) %>%
  mutate(d = round(d,1)) %>%
  left_join(select(synth_df,d=depth,easting,northing,soilUnitID,soilUnit)) %>%
  filter(cdp_id %in% filter(selected_points,Z2_flag==1)$id) %>%
  drop_na() %>% select(!c(easting,northing)) %>% rename(name = CPT_name) %>%
  left_join(st_set_geometry(select(filter(selected_points,Z2_flag==1),cdp_id=id,loc_id,xid,yid),value=NULL))  %>%
  relocate(loc_id,xid,yid,.after = name) %>%
  select(any_of(colnames(cpt_data)))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 6 Build combined 3D dataset ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 6.1 Combine geophysics and CPT ----
data <- bind_rows(geophys,cpt_data,missing) %>% arrange(loc_id,d) %>%
  mutate(Z1 = as.integer(soilUnitID)) %>% rename(Z2 = bc_qc) %>%
  filter(d %in% seq(20,100,0.5))

full_df <- data

## 6.2 Train/test split (20% hold-out) ----
set.seed(16)
test_locs <- data %>% drop_na(Z2) %>% distinct(loc_id) %>%
  sample_n(size = ceiling(n()*0.2)) %>% pull(loc_id)
cpt_locs <- data %>% drop_na(Z2) %>% distinct(loc_id) %>%
  filter(!(loc_id %in% test_locs)) %>% pull(loc_id)
data$Z2[which(data$loc_id %in% test_locs)] <- NA

## 6.3 Depth grid and dimensions ----
depth_vec <- sort(distinct(data,d)$d)
dimd <- length(depth_vec)
dims <- c(dimd,dimx,dimy)

## 6.4 Depth groupings for parallel sampling ----
lattice <- expand.grid(d=depth_vec,loc_id = 1:prod(dims[2],dims[3])) %>%
  left_join(select(data,loc_id,d,Z1 = soilUnitID)) %>%
  mutate(Z1 = replace_na(Z1,0)) %>% arrange(loc_id,d)

diagonals <- 0
K <- 8

data$x <- (data$x- min(data$x))/100000
data$y <- (data$y- min(data$y))/100000

groups <- data %>%
  arrange(d) %>%
  group_by(d) %>%
  summarise(n = sum(!is.na(Z2))) %>%
  mutate(n = cumsum(n)) %>%
  ungroup() %>%
  mutate(n = cut_interval(n,length = 250),
         groups = as.numeric(n)) %>%
  select(!n)

data <- data %>%
  left_join(groups)

drop_na(data)$groups %>% table()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 7 Build cross-section line definitions ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 7.1 Mark valid grid cells ----
cell_grid$valid <- cell_grid$loc_id %in% selected_points$loc_id

## 7.2 Helper: bottom y-index for a given x ----
compute_bottom_y <- function(x0){
  cell_grid %>%
    filter(valid) %>%
    filter(xid==30) %>%
    pull(yid) %>% min()
}

## 7.3 Build vertical and diagonal line grids ----
vert_endpoints <- c(30,35,40,45)
diag_endpoints <- c(-8,0,8,16,24,32,48,56,64,72,80,88)
line_grid <- lapply(1:4,function(i) mutate(filter(cell_grid,valid,xid==vert_endpoints[i]),line=i)) %>%
  bind_rows()
for (i in 1:length(diag_endpoints)){
  x0 <- diag_endpoints[i]
  y0 <- compute_bottom_y(x0)
  if(x0 > 40) x1 <- 0 else x1 <- (x0+87)
  npts <- abs(x0-x1)
  new_line_grid <- data.frame(xid=x0:x1, yid=y0:(y0+npts)) %>%
    left_join(cell_grid) %>% filter(valid) %>%
    mutate(line = 4+i)
  line_grid <- bind_rows(line_grid,new_line_grid)
}
line_df <- st_set_geometry(select(line_grid,line,loc_id,xid,yid),value=NULL)

## 7.4 Save line definitions ----
saveRDS(line_df,"data/processed/line_df.rds")
message("Saved data/processed/line_df.rds")

#%%%%%%%%%%%%%%%%%%%%%%%%%%
# 8 Save processed data ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = setdiff(ls(), c("data","lattice","dims","diagonals","K","full_df","test_locs","cpt_locs")))

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
save(data, full_df, test_locs, cpt_locs, dims, K,
     file = "data/processed/data3D.RData")
message("Saved data/processed/data3D.RData")

# Build hexagonal grouping lattice (saves data/processed/points_grid.rds)
# Requires geomix_setup from 02_fit_models.R - run after model setup if needed
# source("scripts/utils/hex_grouping.R")
