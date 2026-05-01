#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1 Preliminaries ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 1.1 Libraries ----
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyverse)
library(patchwork)
library(concaveman)
library(grid)

# Requires in memory: full_df, geomix_setup
# Run 01_process_data.R and 02_fit_models.R first

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2 Build IJV study area outline ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 2.1 Spatial point data ----
pts <- full_df %>%
  group_by(x,y) %>%
  mutate(Z2flag = if_else(sum(is.na(Z2))==n(),0,1),
         x=x/100,y=y/100,
         Z2flag = factor(Z2flag,levels = c(0,1), labels = c("Common Depth Point","Cone Penetrometer Test Site"))) %>%
  distinct(x,y,loc_id,Z2flag) %>%
  left_join(
    distinct(drop_na(full_df),loc_id) %>%
      mutate(testFlag = !(loc_id %in% distinct(drop_na(geomix_setup$df),loc_id)$loc_id),
             testFlag = factor(testFlag,levels = c(F,T), labels =c("Train CPT","Test CPT")))

  )

## 2.2 Project to WGS84 and build concave hull ----
ijv_sf  <- st_as_sf(pts, coords = c("x","y"), crs = 25831)
ijv_ll  <- st_transform(ijv_sf, 4326)
coords <- st_coordinates(ijv_ll)[, 1:2, drop = FALSE]
hull_xy <- concaveman(coords, concavity = 2)
ijv_outline <- st_sfc(
  st_polygon(list(as.matrix(hull_xy))),
  crs = st_crs(ijv_ll)
)
ijv_outline <- st_buffer(ijv_outline,50)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3 Country basemaps ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nl <- ne_countries(country = "Netherlands", scale = 10L, returnclass = "sf") |>
  st_transform(4326)
bel <- ne_countries(country = "Belgium", scale = 10L, returnclass = "sf") |>
  st_transform(4326)
ger <- ne_countries(country = "Germany", scale = 10L, returnclass = "sf") |>
  st_transform(4326)

## 3.1 Label anchor points ----
nl_cent <- st_coordinates(st_centroid(st_union(nl)))
lab_nl  <- data.frame(lon = nl_cent[1], lat = nl_cent[2], txt = "Netherlands")

ijv_bb   <- st_bbox(ijv_outline)
ijv_cx   <- (ijv_bb["xmin"] + ijv_bb["xmax"]) / 2
ijv_off  <- 0.08 * (ijv_bb["ymax"] - ijv_bb["ymin"])
lab_ijv  <- data.frame(lon = ijv_cx, lat = ijv_bb["ymin"] - ijv_off,
                       txt = "IJmuiden Ver\nWind Farm Zone")

sea_col   <- "#D9ECF5"; land_col <- "#3C6287"; ijv_col <- "#E67E22"

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4 Seismic line layout ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

synth_df <- readRDS("data/application/seismic_cdp.rds")

lines_map <- full_df %>%
  dplyr::filter(d <= 50,xid > 29, xid < 46) %>%
  mutate(easting = x/100, northing = y/100) %>%
  distinct(xid,yid,easting,northing) %>%
  arrange(northing) %>%
  st_as_sf(coords = c("easting","northing"), crs = 25831) %>%
  st_transform(., 4326) %>%
  group_by(xid) %>%
  summarise(geometry = st_union(geometry)) %>%
  st_cast("LINESTRING")

all_lines_map <- synth_df %>%
  mutate(easting = easting, northing = northing) %>%
  distinct(id,easting,northing) %>%
  arrange(northing) %>%
  st_as_sf(coords = c("easting","northing"), crs = 25831) %>%
  st_transform(., 4326) %>%
  group_by(id) %>%
  summarise(geometry = st_union(geometry)) %>%
  st_cast("LINESTRING")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5 Build map-version spatial grid (WGS84) ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cell_size <- 350
locs <- readRDS("data/processed/locs.rds")
locs$x <- locs$x/100; locs$y <- locs$y/100
loc_sf <- st_as_sf(locs,coords=c("x","y"),crs = 25831)
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
cell_grid <- st_transform(cell_grid,crs = 4326)
cell_grid$valid <- cell_grid$loc_id %in% geomix_setup$df$loc_id
cell_grid$z2 <- cell_grid$loc_id %in% drop_na(geomix_setup$df)$loc_id +
  cell_grid$loc_id %in% drop_na(full_df)$loc_id
cell_grid$z2 <- factor(cell_grid$z2,labels = c("No CPT data","Test CPT data","Train CPT data"))

saveRDS(cell_grid,'data/processed/cell_grid_map.rds')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 6 Define cross-section lines ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 6.1 Helper: bottom y-index for a given x ----
compute_bottom_y <- function(x0){
  cell_grid %>%
    filter(valid) %>%
    filter(xid==30) %>%
    pull(yid) %>% min()
}

## 6.2 Build vertical and diagonal line grids ----
vert_endpoints <- c(30,35,40,45)
diag_endpoints <- c(-8,0,8,16,24,32,48,56,64,72,80,88)
line_grid <- lapply(1:4,function(i) mutate(filter(cell_grid,valid,xid==vert_endpoints[i]),line=i)) %>%
  bind_rows()
for (i in 1:length(diag_endpoints)){
  x0 <-  diag_endpoints[i]
  y0 <- compute_bottom_y(x0)
  if(x0 > 40) x1 <- 0 else x1 <- (x0+87)
  npts <- abs(x0-x1)

  new_line_grid <-data.frame(xid=x0:x1,
                             yid=y0:(y0+npts)) %>%
    left_join(cell_grid) %>% filter(valid) %>%
    mutate(line = 4+i)
  line_grid <- bind_rows(line_grid,new_line_grid)
}
line_df <- st_set_geometry(select(line_grid,line,loc_id,xid,yid),value=NULL)
line_grid <- st_as_sf(line_grid)

## 6.3 Save line definitions ----
saveRDS(line_df,"data/processed/line_df.rds")

## 6.4 Line labels and sf objects for plotting ----
line_label <- line_df %>%
  group_by(line) %>%
  summarise(mxid = max(xid), myid = max(yid[which(xid == mxid)])) %>%
  left_join(distinct(full_df,mxid=xid,myid=yid,x,y)) %>%
  mutate(x=x/100 + 500, y = y/100 + 800*(-1)^(myid < 45)) %>%
  st_as_sf(., coords = c("x","y"), crs = 25831) %>%
  st_transform(crs = 4326)

line_sf <- left_join(line_df,distinct(full_df,loc_id,x,y)) %>%
  mutate(x = x/100,y=y/100) %>%
  st_as_sf(., coords = c("x","y"), crs = 25831) %>%
  st_transform(crs = 4326) %>%
  group_by(line) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("LINESTRING")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 7 Map figures ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## 7.1 Overview map: Netherlands + IJV location (mp1) ----
mp1 <- ggplot() +
  theme_void(base_size = 8) +
  theme(panel.background = element_rect(fill = sea_col, color = NA),
        plot.margin = margin(0,0,0,0)) +
  geom_sf(data = bel, fill = 'gray', color = NA) +
  geom_sf(data = ger, fill = 'gray', color = NA) +
  geom_sf(data = nl,  fill = land_col, color = NA) +
  geom_sf(data = ijv_outline, fill = ijv_col, color = NA) +
  annotate("text", x = lab_nl$lon-0.3,  y = lab_nl$lat+0.1,
           label = lab_nl$txt, family = "", fontface = "bold",
           size = 3.3, color = "white", hjust = 0) +
  annotate("text", x = lab_ijv$lon, y = lab_ijv$lat,
           label = lab_ijv$txt, family = "", fontface = "plain",
           size = 2, color = "#1f2e40",vjust = 1.1) +
  annotate("text", x = lab_ijv$lon, y = lab_ijv$lat-0.57,
           label = "Alpha & Beta sites", family = "", fontface = "plain",
           size = 1.5, color = "#1f2e40",vjust = 0) +
  coord_sf(xlim = c(2.6, 7.6), ylim = c(50.6, 54.0), expand = FALSE)

## 7.2 CPT locations map (mp2) ----
bb_sub <- st_bbox(c(
  xmin = 3.418,
  xmax = 3.437,
  ymin = 52.8465,
  ymax = 52.8596
), crs = st_crs(all_lines_map))

mp2 <- ggplot() +
  theme_void() +
  theme(
    plot.background  = element_rect(fill = sea_col, color = NA),
    panel.background = element_rect(fill = sea_col, color = NA),
    legend.position = "bottom",
    legend.title.position = "top",
    legend.direction = "horizontal",
    legend.title.align = 0.5,
    legend.title = element_text(size=8),
    legend.text   = element_text(size = 6, margin = margin(0,0,0,1))
  ) +
  geom_sf(data = filter(cell_grid,valid),fill=NA,show.legend = F,size=0.2, col='#E67E22') +
  geom_sf(data = filter(ijv_ll,Z2flag =="Cone Penetrometer Test Site"),
          aes(col = testFlag), size = 1, fill = 'black',show.legend = TRUE) +
  annotate("rect",xmin = bb_sub['xmin'], xmax = bb_sub['xmax'], ymin = bb_sub['ymin'],
           ymax = bb_sub['ymax'],fill =NA,col='yellow',linewidth=0.5)+
  scale_colour_manual(
    values = c('black','blue'),
    guide  = guide_legend(override.aes = list(size = 3))
  ) +
  scale_shape_manual(values = c(16,25))+
  labs(col = "CPT Locations", shape = "")

## 7.3 Seismic lines map (mp3) ----
mp3 <- ggplot() +
  theme_void() +
  theme(
    plot.background  = element_rect(fill = sea_col, color = NA),
    panel.background = element_rect(fill = sea_col, color = NA),
    legend.position = "bottom",
    legend.title.position = "top",
    legend.direction = "horizontal",
    legend.title.align = 0.5,
    legend.title = element_text(size=8),
    legend.text      = element_text(size = 6),
    legend.key.width = unit(1.4, "lines"),
    legend.key.height= unit(1.2, "lines")
  ) +
  geom_sf(data = all_lines_map , col='#E67E22',size=0.05)+
  geom_sf(data = filter(line_sf,line %in% seq(2,16,2)),aes(col='Example Cross Section Lines'),show.legend = T) +
  geom_sf_text(data = filter(line_label,line %in% seq(2,16,2)),aes(label = line/2),hjust = 0,show.legend = F,size=3) +
  annotate("rect",xmin = bb_sub['xmin'], xmax = bb_sub['xmax'], ymin = bb_sub['ymin'],
           ymax = bb_sub['ymax'],fill =NA,col='yellow',linewidth=0.5)+
  scale_shape_manual(values = c(16,25))+
  scale_color_manual(values = c("Example Cross Section Lines" = "black")) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.2)))+
  labs(col = "Seismic Traces", shape = "")

## 7.4 Inset plots ----
lines_sub <- st_crop(all_lines_map, bb_sub)
grid_sub <- st_crop(filter(cell_grid,valid),bb_sub)
points_sub <- st_crop(filter(ijv_ll,Z2flag =="Cone Penetrometer Test Site"),bb_sub)

mp2p1 <- ggplot() +
  theme_void() +
  theme(
    plot.background  = element_rect(fill = NA, color = NA),
    plot.margin = margin(-1,-1,-1,-1),
    panel.background = element_rect(fill = NA, color = NA)
  ) +
  annotate("rect",xmin = bb_sub['xmin'], xmax = bb_sub['xmax'], ymin = bb_sub['ymin'],
           ymax = bb_sub['ymax'],fill =NA,col='yellow',linewidth=0.5)+
  geom_sf(data = grid_sub,fill=sea_col,show.legend = F,size=0.1, col='#E67E22') +
  geom_sf(data = points_sub,aes(col=testFlag) ,size = 0.5,show.legend = F) +
  scale_colour_manual(
    values = c('black','blue')
  )

mp3p1 <- ggplot() +
  theme_void() +
  theme(
    plot.background  = element_rect(fill = NA, color = NA),
    plot.margin = margin(-1,-1,-1,-1),
    panel.background = element_rect(fill = NA, color = NA)
  ) +
  annotate("rect",xmin = bb_sub['xmin'], xmax = bb_sub['xmax'], ymin = bb_sub['ymin'],
           ymax = bb_sub['ymax'],fill =sea_col,col='yellow',linewidth=0.5)+
  geom_sf(data = lines_sub, col='#E67E22',size=0.1,show.legend = F)

#%%%%%%%%%%%%%%%%%%%%%%%%%%
# 8 Save figures ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("results/figures/map1.pdf",mp1,  width = 4.78/3, height = 2,device = cairo_pdf)
ggsave("results/figures/map2.pdf",mp2,  width = 4.78/3, height = 2,device = cairo_pdf)
ggsave("results/figures/map3.pdf",mp3,  width = 4.78/3, height = 2,device = cairo_pdf)
ggsave("results/figures/map2sub.pdf",mp2p1,  width =1.2, height = 1.2, units = "cm",device = cairo_pdf)
ggsave("results/figures/map3sub.pdf",mp3p1,  width =1.2, height = 1.2, units = "cm",device = cairo_pdf)
