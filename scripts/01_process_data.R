library(matrixStats)
library(tidyverse)
library(reshape2)
library(nimble)
library(ggh4x)
library(tictoc)
library(scales)
library(sf)
select <- dplyr::select

box_cox <- function(y,lambda = 0.6) (y ^ lambda - 1) / lambda
inv_box_cox <- function(y, lambda=0.6) (y * lambda + 1)^(1 / lambda)
su_colours <-c("yellow", "forestgreen", "maroon", "darkorange3", "darkblue", "red4", "gray30",
               "dodgerblue4", "goldenrod3", "darkolivegreen", "firebrick", "mediumvioletred",
               "chocolate4", "darkslateblue", "slategray4", "seagreen4")

soil_units <- c("GT1","GT2","GT2c","GT3","GT4","GT5","GT5*","GT6")

ijv_10cm <- readRDS("data/application/cpt_profiles.rds") %>%
  rename(x = Easting_m, y = Northing_m,
         bathymetry5m = Bathymetry_UHR_5m_LAT,
         bathymetryp5m = Bathymetry_MBES_0p5m_LAT,
         depthBSF = Depth_m_bsf,
         qc = qn_MPa) %>%
  mutate(x = x*100, y = y*100) %>%
  mutate(d = round(depthBSF - bathymetryp5m,1))
sample_wsynth <- readRDS("data/application/cpt_stratigraphy.rds")
synth_df <- readRDS("data/application/seismic_cdp.rds")

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

# Now group by loc_id
selected_points <- st_join(loc_sf, cell_grid, join = st_within) %>%
  mutate(dist_to_centre = st_distance(geometry, centre, by_element = TRUE)) %>%
  group_by(loc_id) %>%
  arrange(Z2_flag != 1, dist_to_centre) %>%
  slice(1) %>%
  ungroup()

which(!(1:ncpt %in% filter(selected_points,Z2_flag==1)$id))

geophys <- synth_df %>%
  filter(depth <= 50) %>%
  mutate(easting = easting*100, northing = northing*100) %>%
  left_join(rename(cdp_loctions,cdp_id = id)) %>%
  filter(cdp_id %in% filter(selected_points,Z2_flag==0)$id) %>%
  select(name = id,cdp_id,x=easting,y=northing,d=depth,soilUnitID,soilUnit) %>%
  left_join(st_set_geometry(select(filter(selected_points,Z2_flag==0),cdp_id=id,loc_id,xid,yid),value=NULL)) %>%
  mutate(cdp_id = cdp_id+100000) %>% relocate(loc_id,xid,yid,.after = cdp_id) %>%
  drop_na(loc_id)


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

##############################################################
data <- bind_rows(geophys,cpt_data,missing) %>% arrange(loc_id,d) %>%
  mutate(Z1 = as.integer(soilUnitID)) %>% rename(Z2 = bc_qc) %>%
  filter(d %in% seq(20,100,0.5))

full_df <- data
set.seed(16)
test_locs <- data %>% drop_na(Z2) %>% distinct(loc_id) %>%
  sample_n(size = ceiling(n()*0.2)) %>% pull(loc_id)
cpt_locs <- data %>% drop_na(Z2) %>% distinct(loc_id) %>%
  filter(!(loc_id %in% test_locs)) %>% pull(loc_id)
data$Z2[which(data$loc_id %in% test_locs)] <- NA

depth_vec <- sort(distinct(data,d)$d)
dimd <- length(depth_vec)
dims <- c(dimd,dimx,dimy)

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

rm(list = setdiff(ls(), c("data","lattice","dims","diagonals","K","full_df","test_locs","cpt_locs")))

# Save processed data for use by downstream scripts
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
save(data, full_df, test_locs, cpt_locs, dims, K,
     file = "data/processed/data3D.RData")
message("Saved data/processed/data3D.RData")

# Build hexagonal grouping lattice (saves data/processed/points_grid.rds)
# Requires geomix_setup from 02_fit_models.R - run after model setup if needed
# source("scripts/utils/hex_grouping.R")
