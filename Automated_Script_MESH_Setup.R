#=========================================================================================
# SCRIPT:
# Mountain MESH model setup script
#
# DESCRIPTION:
# This script produces the MESH_drainage_database and MESH_Parameters file (r2c and netcdf) 
# that are required to run both Mountain and Original MESH versions
#
# USER INPUTS:
#   > A Raster file of the modeling domain a value representing each pixcel unique ID
#   > An elevation raster of the NWP climate forcing dataset
#   > An elevation raster of the modeling dataset
#   > A land cover raster
#   > A location of the streamflow gauging station
#
# REQUIREMENTS:
#   > TauDEM
#   > R and its packages
#
# DISCLAIMER:
# This script is made available for other's use on an "as is" basis, with no warranty,
# either expressed or implied, as to its fitness for any particular purpose.
#
# AUTHOR: Zelalem Tesemma, GWSI, University of Saskatchewan, zelalem.tesemma@usask.ca
#
# HISTORY: 12-Jan-2020  Develop and test the preliminary workflow
#         15 - Feb - 2020 furhter development and tested on Bow River at Banff
#=========================================================================================
rm(list = ls())
## loading libraries
# install.packages('ncdf4', repos = "http://cran.us.r-project.org")
library(raster)
#library(RSAGA)
library(sp)
library(maptools)
library(rgeos)
library(rgdal)
library(numbers)
library(shapefiles)
library(ncdf4)
######################################## User Inputs #######################################################################
# Change the following directory path to your working directory path
setwd("C:/PD/Automated_MESH_Setup_Sample")
ncpath <- "C:/PD/Automated_MESH_Setup_Sample"
#
########### Comment either Mountain or Original MESH model version #########################################################
MESHVersion <- "Mountain"
# MESHVersion <- "Original"
######### Set the minimum GRUs fraction in a give cell and the model grid fraction in the watershed boundary ###############
Minimum_Basin_Fraction <- 0.001
Minimum_GRU_Fraction <- 0.008
Minimum_Glacier_Fraction <- 0.0016 # Minimum Glacier area of 0.05 km2
Number_Soil_Layers <- 3
#
# The maximum number of river class you want to have (IAK - River Classification)
IAKmax <- 5
# The maximum channel slope
Min_Chanel_Slope <- 0.001
#
# The minimum threshold flow accumulation area to initiate stream.
MinThresh <- 100 # based on 90 meter resolution dem (0.1km2 = 12, 0.5km2 = 61, 1km2 = 123)
#
###### Name of the drainage basin of study ############################
BasinName <- paste0("BowRiverBasinBanff_GEM_0p125_MinThresh_",MinThresh,"_") 
#
################### Slope and Aspect Classification for GRUs Creations ################################################################################################################
##### Two slope Classes: (0 - 10) and (> 10) and Tow aspects: North and South Facing can be modified to include other orientations ##############
slopeclass <- matrix(c(0, 10, 10, 90, 1, 2), nr = 2, nc = 3)
aspectclass <- matrix(c(0, 90, 270, 90, 270, 360, 1, 2, 1), nr = 3, nc = 3)
# For 4 quadrants (North, East, South and West) (315 < N <= 45), (45 < E <= 135), (135 < S <= 225), (225 < W <= 315) # From GEM
# For 4 quadrants (North, East, South and West) (292.5 < N <= 67.5), (67.5 < E <= 112.5), (112.5 < S <= 247.5), (247.5 < W <= 292.5) # ESRI ArcGIS
########## Set the climate forcing grid domain #################################################################################################
LLXcorner <- -116.75
LLYcorner <- 50.00
NumRow <- 15
NumCol <- 29
XRes <- 0.125
YRes <- 0.125
URXcorner <- (LLXcorner + (NumCol*XRes)) 
URYcorner <- (LLYcorner + (NumRow*YRes)) 
#
### Import the Hydroshed / MERIT DEM, NWP climate forcing elevation and Land cover data for your model ####################
domain_dem <- raster("domain_HydroShed_dem.tif")
#
domain_landcover <- raster("res90_reclass_domain_landcover_30m_buffered_glaciers.tif")
#
#### Creating the basin outlet using the lat and long data ############
streamgauge <- data.frame(lat = 51.173, long = -115.567)
coordinates(streamgauge) <- ~long+lat
crs(streamgauge) <- crs(domain_dem)
raster::shapefile(streamgauge, "approxoutlets.shp", overwrite=TRUE)
#
# nwp_dem <- raster("nwp_dem.tif")
### Import the percent of clay and sand contents for the three soil layers #################
domain_Clay1 <- raster("domain_Clay1.tif")
domain_Clay2 <- raster("domain_Clay2.tif")
domain_Clay3 <- raster("domain_Clay3.tif")
#
domain_Sand1 <- raster("domain_Sand1.tif")
domain_Sand2 <- raster("domain_Sand2.tif")
domain_Sand3 <- raster("domain_Sand3.tif")
#
############################## No Change after this except for GRUs Discritization Section ################################
###########################################################################################################################
NumGrids <- matrix(seq(1,(NumRow*NumCol),1), NumRow, NumCol, byrow = T) 
nwp_grid <- raster(NumGrids)
extent(nwp_grid) <- extent(LLXcorner,(LLXcorner+(NumCol*XRes)),LLYcorner,(LLYcorner+(NumRow*YRes)))
dim(nwp_grid) <- c(NumRow, NumCol)
res(nwp_grid) <- c(XRes, YRes)
crs(nwp_grid) <- crs(domain_dem) # crs(nwp_grid) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
#### Resample the nwp_elevation and the Land cover map into the model_elevation spatial resolution
nwp_zone1 <- resample(nwp_grid, domain_dem, method="ngb")
nwp_zone <- crop(nwp_zone1, nwp_grid)
# nwp_elevn <- resample(nwp_dem, nwp_zone, method="bilinear")
land_cover <- resample(domain_landcover, nwp_zone, method="ngb")
############## Common working variable ####################################################################################
ResFactor <- res(nwp_grid)/res(domain_dem)
ResFactor  <- ResFactor[1]
FdirNumber <- matrix((as.integer(c(0,-1,-1,-1,0,1,1,1,1,1,0,-1,-1,-1,0,1))),8,2)
##################### TauDEM approach (DEM post-processing and Watershed Delineation #######################################
# Pitremove
system("mpiexec -n 8 pitremove -z domain_HydroShed_dem.tif -fel domain_demfel.tif")

# D8 flow directions
system("mpiexec -n 8 D8Flowdir -p domain_demp.tif -sd8 domain_demsd8.tif -fel domain_demfel.tif",show.output.on.console=F,invisible=F)

# Contributing area
system("mpiexec -n 8 AreaD8 -p domain_demp.tif -ad8 domain_demad8.tif -nc")

# Grid Network 
system("mpiexec -n 8 Gridnet -p domain_demp.tif -gord domain_demgord.tif -plen domain_demplen.tif -tlen domain_demtlen.tif")

# Threshold
system("mpiexec -n 8 Threshold -ssa domain_demad8.tif -src domain_demsrc.tif -thresh 100")

# Move Outlets to fall in the flow acculation pixcel
system("mpiexec -n 8 moveoutletstostreams -p domain_demp.tif -src domain_demsrc.tif -o approxoutlets.shp -om outlet.shp")
#
# # Check for the moved outlet
# outpt=read.shp("Outlet.shp")
# approxpt=read.shp("approxoutlets.shp")
# ssa <- raster("domain_demsrc.tif")
# plot(ssa)
# zoom(ssa)
# points(outpt$shp[2],outpt$shp[3],pch=19,col=2)
# points(approxpt$shp[2],approxpt$shp[3],pch=19,col=4)
#
# Contributing area upstream of outlet
system("mpiexec -n 8 Aread8 -p domain_demp.tif -o outlet.shp -ad8 domain_demssa.tif")

# Threshold
system("mpiexec -n 8 Threshold -ssa domain_demssa.tif -src domain_demsrc.tif -thresh 100")

# # Drop Analysis
# system("mpiexec -n 8 Dropanalysis -p domain_demp.tif -fel domain_demfel.tif -ad8 domain_demad8.tif -ssa domain_demssa.tif -drp logandrp.txt -o outlet.shp -par 5 500 10 0")

# Stream Reach and Watershed
system("mpiexec -n 8 Streamnet -fel domain_demfel.tif -p domain_demp.tif -ad8 domain_demad8.tif -src domain_demsrc.tif -o outlet.shp -ord domain_demord.tif -tree domain_demtree.txt -coord domain_demcoord.txt -net domain_demnet.shp -w domain_demw.tif")
# 
## Collecting the drainage files generated from previous steps ###
##########################################################################################################################################
filldem0 = raster("domain_demfel.tif")
filldem1 <- crop(filldem0, nwp_zone)
#
fdir0 = raster("domain_demp.tif")
fdir1 <- crop(fdir0, nwp_zone)
#
facc_at_outlet0 = raster("domain_demssa.tif")
facc_at_outlet1 <- crop(facc_at_outlet0, nwp_zone)
#
slp_fdir0 = raster("domain_demsd8.tif")
slp_fdir1 <- crop(slp_fdir0, nwp_zone)
#
stream_order0 = raster("domain_demord.tif")
stream_order1 <- crop(stream_order0, nwp_zone)

flow_path0 = raster("domain_demplen.tif")
flow_path1 <- crop(flow_path0, nwp_zone)
#
drain_net = raster("domain_demsrc.tif")
#### More work is needed ######################################################################################
## Compute nominal grid area [km2]. Also omitte NA's lie outside of the rastered region
BasinAreaInNWPGrid = aggregate(area(facc_at_outlet1, na.rm=TRUE, weights=FALSE), fact = ResFactor, fun = sum, na.rm=TRUE)
NWPGridFullArea = aggregate(area(filldem1, na.rm=TRUE, weights=FALSE), fact = ResFactor, fun = sum, na.rm=TRUE)
FracBasinAreaInNWPGrid <- (BasinAreaInNWPGrid / NWPGridFullArea)
FracBasinAreaInNWPGrid[FracBasinAreaInNWPGrid <= Minimum_Basin_Fraction] <- NA
FracBasinAreaInNWPGridRes <- resample(FracBasinAreaInNWPGrid, facc_at_outlet1, method="ngb")
facc_at_outlet2 <- mask(facc_at_outlet1, FracBasinAreaInNWPGridRes)
## Masking the maximum flow accumulation generated from the high resolution elevation inside MESH modeling grid
facc1 = aggregate(facc_at_outlet2, fact = ResFactor, fun = max, na.rm=TRUE)
facc2 <- resample(facc1, facc_at_outlet2, method="ngb")
facc3 <- facc_at_outlet2==facc2 & facc2==facc_at_outlet2
facc3[facc3 != 1] <- NA
## Masking the maximum flow accumulation generated from the high resolution elevation inside MESH modeling grid by the minimum elevation of the high resolution elevation
filldem1_facc1 <- mask(filldem1,facc3)
filldem2_facc2 = aggregate(filldem1_facc1, fact = ResFactor, fun = min, na.rm=TRUE)
filldem3_facc3 <- resample(filldem2_facc2, facc_at_outlet1, method="ngb")
filldem4_facc4 <- filldem1_facc1 == filldem3_facc3 & filldem3_facc3 == filldem1_facc1
filldem4_facc4[filldem4_facc4 != 1] <- NA
#
## Masking elevation by the minimum elevation of the maximum flow accumulation cell inside the modeling grid
elevn0 <- mask(filldem1, filldem4_facc4)
#
## River bed elevation at outlet of the grid cell 
Elev <- as.matrix(aggregate(elevn0, fact = ResFactor, fun = min, na.rm=TRUE))
#
## Masking flow direction by the minimum elevation of the maximum flow accumulation inside MESH modeling grid
fdir2 <- as.matrix(mask(fdir1, filldem4_facc4))
## Rank the MESH modeling grid based on outlet invert elevation produced from above
elevn1 <- as.vector(Elev)
rank0 <- rank(-elevn1, ties.method= "first", na.last = "keep")
#
## Grid number - gives order of computation
grid_rank  <- matrix(rank0,NumRow,NumCol)
zz <- which(grid_rank == max(grid_rank,na.rm=TRUE), arr.ind=TRUE)
zz[,2] <- zz[,2] + 1
grid_rank[zz] <- (max(grid_rank,na.rm=TRUE) + 1)
#
rank1 <- raster(grid_rank)
extent(rank1) <- extent(nwp_grid)
dim(rank1) <- dim(nwp_grid)
res(rank1) <- res(nwp_grid)
crs(rank1) <- crs(nwp_grid)
rank2 <- resample(rank1, facc_at_outlet1, method="ngb")
#
## Producing Next of the modeling grid cell from flow direction and rank (Receiving cell number (must be more than N)
#
rank3 <- as.matrix(rank2)
next0 <- matrix(nrow = nrow(rank2), ncol = ncol(rank3))
for (i in 1:8) {
  xy <- which(fdir2 == i, arr.ind=TRUE)
  xz <- xy
  xz[,1] <- (xy[,1] + as.integer(FdirNumber[i,1]))
  xz[,2] <- (xy[,2] + as.integer(FdirNumber[i,2]))
  next0[xy] <- rank3[xz] 
}
#
next1 <- raster(next0)
extent(next1) <- extent(rank2)
dim(next1) <- dim(rank2)
res(next1) <- res(rank2)
crs(next1) <- crs(rank2)
next2 = aggregate(next1, fact = ResFactor, fun = min, na.rm=TRUE)
grid_next <- as.matrix(next2)
#
zz <- which(grid_rank == (max(grid_rank,na.rm=TRUE) - 1), arr.ind=TRUE)
grid_next[zz] <- max(grid_rank,na.rm=TRUE)
#
## Generating flow direction for the modeling grid
grid_fdir <- matrix(NA, NumRow, NumCol)
## Potential to change the line below ##################################
# for (i in 1:8) {
#   xy <- which(!is.na(grid_rank), arr.ind=TRUE)
#   xz <- xy
#   xz[,1] <- (xy[,1] + as.integer(FdirNumber[i,1]))
#   xz[,2] <- (xy[,2] + as.integer(FdirNumber[i,2]))
#   grid_fdir[which(grid_next[xz] == grid_rank[xy], arr.ind=TRUE)] <- i
# }
for (i in 1:nrow(grid_fdir)) {
  for (j in 1:ncol(grid_fdir)) {
    if (is.na(grid_next[i,j])) {
      grid_fdir[i,j] <- NA}
    else if ((!is.na(grid_rank[i,j])) & grid_next[i,j] == 0) {
      grid_fdir[i,j] <- 1}
    else if ((!is.na(grid_rank[i,j+1])) & (grid_next[i,j] == grid_rank[i,j+1])) {
      grid_fdir[i,j] <- 1}
    else if ((!is.na(grid_rank[i-1,j+1])) & (grid_next[i,j] == grid_rank[i-1,j+1])) {
      grid_fdir[i,j] <- 2}
    else if ((!is.na(grid_rank[i-1,j])) & (grid_next[i,j] == grid_rank[i-1,j])) {
      grid_fdir[i,j] <- 3}
    else if ((!is.na(grid_rank[i-1,j-1])) & (grid_next[i,j] == grid_rank[i-1,j-1])) {
      grid_fdir[i,j] <- 4}
    else if ((!is.na(grid_rank[i,j-1])) & (grid_next[i,j] == grid_rank[i,j-1])) {
      grid_fdir[i,j] <- 5}
    else if ((!is.na(grid_rank[i+1,j-1])) & (grid_next[i,j] == grid_rank[i+1,j-1])) {
      grid_fdir[i,j] <- 6}
    else if ((!is.na(grid_rank[i+1,j])) & (grid_next[i,j] == grid_rank[i+1,j])) {
      grid_fdir[i,j] <- 7}
    else if ((!is.na(grid_rank[i+1,j+1])) & (grid_next[i,j] == grid_rank[i+1,j+1])) {
      grid_fdir[i,j] <- 8}
    else {
      grid_fdir[i,j] <- NA}
  }
}
grid_fdir1 <- raster(grid_fdir)
extent(grid_fdir1) <- extent(nwp_grid)
dim(grid_fdir1) <- dim(nwp_grid)
res(grid_fdir1) <- res(nwp_grid)
crs(grid_fdir1) <- crs(nwp_grid)
writeRaster(grid_fdir1, "grid_fdirp.tif", datatype="INT2S", overwrite=TRUE)

#### Generating flow accumulation for all modeling grid cell independently #######
nwp_zone2 <- matrix(NA, nrow(nwp_zone), ncol(nwp_zone))
nwp_zone3 <- matrix(NA, nrow(nwp_zone), ncol(nwp_zone))
nwp_zone4 <- matrix(NA, nrow(nwp_zone), ncol(nwp_zone))
nwp_zone5 <- matrix(NA, nrow(nwp_zone), ncol(nwp_zone))
#
for (i in 1:NumRow) {
  for (j in 1:NumCol) {
    if ((i %% 2 != 0) & (j %% 2 != 0)){
      nwp_zone2[(ResFactor*i-(ResFactor-1)):(ResFactor*i), (ResFactor*j-(ResFactor-1)):(ResFactor*j)] <- 1 }
    if ((i %% 2 != 0) & (j %% 2 == 0)){    
      nwp_zone3[(ResFactor*i-(ResFactor-1)):(ResFactor*i), (ResFactor*j-(ResFactor-1)):(ResFactor*j)] <- 1 }
    if ((i %% 2 == 0) & (j %% 2 != 0)){    
      nwp_zone4[(ResFactor*i-(ResFactor-1)):(ResFactor*i), (ResFactor*j-(ResFactor-1)):(ResFactor*j)] <- 1 }
    if ((i %% 2 == 0) & (j %% 2 == 0)){   
      nwp_zone5[(ResFactor*i-(ResFactor-1)):(ResFactor*i), (ResFactor*j-(ResFactor-1)):(ResFactor*j)] <- 1 }
  }
}
#
fdir3 <- as.matrix(fdir1)*nwp_zone2
fdir4 <- raster(fdir3)
extent(fdir4) <- extent(nwp_zone)
dim(fdir4) <- dim(nwp_zone)
res(fdir4) <- res(nwp_zone)
crs(fdir4) <- crs(nwp_zone)
#
fdir5 <- as.matrix(fdir1)*nwp_zone3
fdir6 <- raster(fdir5)
extent(fdir6) <- extent(nwp_zone)
dim(fdir6) <- dim(nwp_zone)
res(fdir6) <- res(nwp_zone)
crs(fdir6) <- crs(nwp_zone)
#
fdir7 <- as.matrix(fdir1)*nwp_zone4
fdir8 <- raster(fdir7)
extent(fdir8) <- extent(nwp_zone)
dim(fdir8) <- dim(nwp_zone)
res(fdir8) <- res(nwp_zone)
crs(fdir8) <- crs(nwp_zone)
#
fdir9 <- as.matrix(fdir1)*nwp_zone5
fdir10 <- raster(fdir9)
extent(fdir10) <- extent(nwp_zone)
dim(fdir10) <- dim(nwp_zone)
res(fdir10) <- res(nwp_zone)
crs(fdir10) <- crs(nwp_zone)
#
writeRaster(fdir4, "domain_demp_1.tif", datatype="INT2S", overwrite=TRUE)
writeRaster(fdir6, "domain_demp_2.tif", datatype="INT2S", overwrite=TRUE)
writeRaster(fdir8, "domain_demp_3.tif", datatype="INT2S", overwrite=TRUE)
writeRaster(fdir10, "domain_demp_4.tif", datatype="INT2S", overwrite=TRUE)
#
## Calculate Contributing area (km2)
cell_area_weight <- area(fdir1, na.rm=TRUE, weights=FALSE)
#
writeRaster(cell_area_weight, "cell_area_weight.tif", datatype="FLT4S", overwrite=TRUE)
#                                
system("mpiexec -n 8 AreaD8 -p domain_demp_1.tif -ad8 domain_demad8_1.tif -wg cell_area_weight.tif -nc")
system("mpiexec -n 8 AreaD8 -p domain_demp_2.tif -ad8 domain_demad8_2.tif -wg cell_area_weight.tif -nc")
system("mpiexec -n 8 AreaD8 -p domain_demp_3.tif -ad8 domain_demad8_3.tif -wg cell_area_weight.tif -nc")
system("mpiexec -n 8 AreaD8 -p domain_demp_4.tif -ad8 domain_demad8_4.tif -wg cell_area_weight.tif -nc")
#
# Read back the flow accumulation files and merge them together 
facc4 = raster("domain_demad8_1.tif")
facc5 = raster("domain_demad8_2.tif")
facc6 = raster("domain_demad8_3.tif")
facc7 = raster("domain_demad8_4.tif")
facc8 <- merge(facc4, facc5, facc6, facc7)
# writeRaster(facc8, "nwp_grid_facc_km2.tif", datatype="FLT4S", overwrite=TRUE)
#
# Craeting basin boundary mask file
basin_boundary_mask <- as.matrix(facc_at_outlet1/facc_at_outlet1)
basin_boundary_mask1 <- matrix(0, nrow(basin_boundary_mask), ncol(basin_boundary_mask), byrow = T)
for (i in 1:nrow(basin_boundary_mask)) {
  for (j in 1:ncol(basin_boundary_mask)) {
    if (!is.na(basin_boundary_mask[i,j]) & (i-1 >= 1) & (j-1 >= 1) & (i+1 <= nrow(basin_boundary_mask)) & (j+1 <= ncol(basin_boundary_mask))) { 
      if (!is.na(basin_boundary_mask[i-1,j-1]) & !is.na(basin_boundary_mask[i-1,j]) &
          !is.na(basin_boundary_mask[i-1,j+1]) & !is.na(basin_boundary_mask[i,j-1]) &
          !is.na(basin_boundary_mask[i,j+1]) & !is.na(basin_boundary_mask[i+1,j-1]) &
          !is.na(basin_boundary_mask[i+1,j]) & !is.na(basin_boundary_mask[i+1,j+1])) {
        basin_boundary_mask1[i,j] <- 0 }
      else {
        basin_boundary_mask1[i,j] <- 1 } 
    }
  }
}
# Creat a mask file to mask only the edge of the nwp_zone
nwp_zone_edge_mask1 <- matrix(1, nrow(nwp_zone), ncol(nwp_zone))
#
for (i in 1:NumRow) {
  for (j in 1:NumCol) {
    nwp_zone_edge_mask1[(ResFactor*i-(ResFactor-2)):(ResFactor*i-1), (ResFactor*j-(ResFactor-2)):(ResFactor*j-1)] <- 0 
  }
}
#
nwp_zone_edge_mask <- ((nwp_zone_edge_mask1*basin_boundary_mask) + basin_boundary_mask1)
nwp_zone_edge_mask[nwp_zone_edge_mask == 0] <- NA
nwp_zone_edge_mask[nwp_zone_edge_mask > 1] <- 1
#
# Masking the grid flow accumulation by nwp_zone_edge_mask
facc9 <- raster(as.matrix(facc8)*nwp_zone_edge_mask)
extent(facc9) <- extent(nwp_zone)
dim(facc9) <- dim(nwp_zone)
res(facc9) <- res(nwp_zone)
crs(facc9) <- crs(nwp_zone)
#
fdir10 <- (as.matrix(fdir1)*nwp_zone_edge_mask)
#
nwp_zone6 <- (as.matrix(nwp_zone)*basin_boundary_mask)
nwp_zone6[is.na(nwp_zone6)] <- 0
#
for (i in 1:nrow(nwp_zone6)) {
  for (j in 1:ncol(nwp_zone6)) {
    if (!is.na(fdir10[i,j]) & (i-1 >= 1) & (j-1 >= 1) & (i+1 <= nrow(nwp_zone6)) & (j+1 <= ncol(nwp_zone6))) {
      if ((fdir10[i,j] == 1) & (nwp_zone6[i,j] == nwp_zone6[i,j+1])){
        fdir10[i,j] <- NA}
      else if ((fdir10[i,j] == 2) & (nwp_zone6[i,j] == nwp_zone6[i-1,j+1])){
        fdir10[i,j] <- NA}
      else if ((fdir10[i,j] == 3) & (nwp_zone6[i,j] == nwp_zone6[i-1,j])){
        fdir10[i,j] <- NA}
      else if ((fdir10[i,j] == 4) & (nwp_zone6[i,j] == nwp_zone6[i-1,j-1])){
        fdir10[i,j] <- NA}
      else if ((fdir10[i,j] == 5) & (nwp_zone6[i,j] == nwp_zone6[i,j-1])){
        fdir10[i,j] <- NA}
      else if ((fdir10[i,j] == 6) & (nwp_zone6[i,j] == nwp_zone6[i+1,j-1])){
        fdir10[i,j] <- NA}
      else if ((fdir10[i,j] == 7) & (nwp_zone6[i,j] == nwp_zone6[i+1,j])){
        fdir10[i,j] <- NA}
      else if ((fdir10[i,j] == 8) & (nwp_zone6[i,j] == nwp_zone6[i+1,j+1])){ 
        fdir10[i,j] <- NA}
      else { 
        fdir10[i,j] <- fdir10[i,j]}
    }
    else {
      fdir10[i,j] <- NA }
  }
}
#
## Correcting drainage direction based on the neighboring cells
for (i in 1:nrow(nwp_zone6)) {
  for (j in 1:ncol(nwp_zone6)) {
    if (!is.na(fdir10[i,j]) & (i-1 >= 1) & (j-1 >= 1) & (i+1 <= nrow(nwp_zone6)) & (j+1 <= ncol(nwp_zone6))) {
      if ((fdir10[i,j] == 2) & ((nwp_zone6[i,j] + 1) == nwp_zone6[i-1,j+1])) {
        fdir10[i,j] <- 1}
      else if ((fdir10[i,j] == 2) & ((nwp_zone6[i,j] - NumCol) == nwp_zone6[i-1,j+1]))    {
        fdir10[i,j] <- 3}
      else if ((fdir10[i,j] == 4) & ((nwp_zone6[i,j] - 1) == nwp_zone6[i-1,j-1])) {
        fdir10[i,j] <- 5}
      else if ((fdir10[i,j] == 4) & ((nwp_zone6[i,j] - NumCol)  == nwp_zone6[i-1,j-1]))    {
        fdir10[i,j] <- 3}
      else if ((fdir10[i,j] == 6) & ((nwp_zone6[i,j] - 1) == nwp_zone6[i+1,j-1])) {
        fdir10[i,j] <- 5}
      else if ((fdir10[i,j] == 6) & ((nwp_zone6[i,j] + NumCol) == (nwp_zone6[i+1,j-1])))    {
        fdir10[i,j] <- 7}
      else if ((fdir10[i,j] == 8) & ((nwp_zone6[i,j] + 1) == nwp_zone6[i+1,j+1])) {
        fdir10[i,j] <- 1}
      else if ((fdir10[i,j] == 8) & ((nwp_zone6[i,j] + NumCol) == (nwp_zone6[i+1,j+1])))    {
        fdir10[i,j] <- 7}
      else { 
        fdir10[i,j] <- fdir10[i,j]}
    }
    else {
      fdir10[i,j] <- NA }
  }
}
#
fdir11 <- raster(fdir10)
extent(fdir11) <- extent(nwp_zone)
dim(fdir11) <- dim(nwp_zone)
res(fdir11) <- res(nwp_zone)
crs(fdir11) <- crs(nwp_zone)
#
nwp_grid_facc1 <- matrix(0, NumRow*NumCol)
nwp_grid_facc <- matrix(0, NumRow*NumCol,8)
#
for (i in 1:8) {
  fdir12 <- fdir11
  fdir12[fdir12 != i] <- NA
  fdir12[!is.na(fdir12)] <- 1
  facc13 <- mask(facc9,fdir12)
  nwp_grid_facc[,i] = as.vector(aggregate(facc13, fact = ResFactor, fun = sum, na.rm=TRUE))
}
#
nwp_grid_facc[is.na(nwp_grid_facc)] <- 0
# Calculation of Octafurcation area for the eight flow directions
GridFaction1 <- matrix(nwp_grid_facc[,1], NumRow, NumCol, byrow = T)
GridFaction2 <- matrix(nwp_grid_facc[,2], NumRow, NumCol, byrow = T)
GridFaction3 <- matrix(nwp_grid_facc[,3], NumRow, NumCol, byrow = T)
GridFaction4 <- matrix(nwp_grid_facc[,4], NumRow, NumCol, byrow = T)
GridFaction5 <- matrix(nwp_grid_facc[,5], NumRow, NumCol, byrow = T)
GridFaction6 <- matrix(nwp_grid_facc[,6], NumRow, NumCol, byrow = T)
GridFaction7 <- matrix(nwp_grid_facc[,7], NumRow, NumCol, byrow = T)
GridFaction8 <- matrix(nwp_grid_facc[,8], NumRow, NumCol, byrow = T)
#
## Fraction of drainage area 
grid_drain_area <- matrix(0, nrow = NumRow, ncol = NumCol, byrow = T)
#
for (i in 1:nrow(grid_drain_area)) {
  for (j in 1:ncol(grid_drain_area)) {
    
    if ((i-1 >= 1) & (j-1 >= 1) & (i+1 <= NumRow) & (j+1 <= NumCol)) {
      if (is.na(grid_fdir[i,j])) {
        grid_drain_area[i,j] <- 0 }
      else if (grid_fdir[i,j] == 1) {
        grid_drain_area[i,j] <- (GridFaction1[i,j])} 
      else if (grid_fdir[i,j] == 2) {
        grid_drain_area[i,j] <- (GridFaction2[i,j])} 
      else if (grid_fdir[i,j] == 3) {
        grid_drain_area[i,j] <- (GridFaction3[i,j])} 
      else if (grid_fdir[i,j] == 4) {
        grid_drain_area[i,j] <- (GridFaction4[i,j])} 
      else if (grid_fdir[i,j] == 5) {
        grid_drain_area[i,j] <- (GridFaction5[i,j])} 
      else if (grid_fdir[i,j] == 6) {
        grid_drain_area[i,j] <- (GridFaction6[i,j])} 
      else if (grid_fdir[i,j] == 7) {
        grid_drain_area[i,j] <- (GridFaction7[i,j])} 
      else if (grid_fdir[i,j] == 8) {
        grid_drain_area[i,j] <- (GridFaction8[i,j])} 
    }
  }
}
###
for (i in 1:nrow(grid_drain_area)) {
  for (j in 1:ncol(grid_drain_area)) {
    if ((i-1 >= 1) & (j-1 >= 1) & (i+1 <= NumRow) & (j+1 <= NumCol)) {
      if (is.na(grid_fdir[i,j])) {
        grid_drain_area[i,j] <- 0 }
      if (!is.na(grid_fdir[i,j-1]) & grid_fdir[i,j-1] != 1) {
        grid_drain_area[i,j] <- (grid_drain_area[i,j] + GridFaction1[i,j-1])}
      if (!is.na(grid_fdir[i+1,j-1]) & grid_fdir[i+1,j-1] != 2) {
        grid_drain_area[i,j] <- (grid_drain_area[i,j] + GridFaction2[i+1,j-1])}
      if (!is.na(grid_fdir[i+1,j]) & grid_fdir[i+1,j] != 3) {
        grid_drain_area[i,j] <- (grid_drain_area[i,j] + GridFaction3[i+1,j])}
      if (!is.na(grid_fdir[i+1,j+1]) & grid_fdir[i+1,j+1] != 4) {
        grid_drain_area[i,j] <- (grid_drain_area[i,j] + GridFaction4[i+1,j+1])}
      if (!is.na(grid_fdir[i,j+1]) & grid_fdir[i,j+1] != 5) {
        grid_drain_area[i,j] <- (grid_drain_area[i,j] + GridFaction5[i,j+1])}
      if (!is.na(grid_fdir[i-1,j+1]) & grid_fdir[i-1,j+1] != 6) {
        grid_drain_area[i,j] <- (grid_drain_area[i,j] + GridFaction6[i-1,j+1])}
      if (!is.na(grid_fdir[i-1,j]) & grid_fdir[i-1,j] != 7) {
        grid_drain_area[i,j] <- (grid_drain_area[i,j] + GridFaction7[i-1,j])}
      if (!is.na(grid_fdir[i-1,j-1]) & grid_fdir[i-1,j-1] != 8) {
        grid_drain_area[i,j] <- (grid_drain_area[i,j] + GridFaction8[i-1,j-1])}
    }
  }
}

## Grid Drainage Area (FRAC) in m2 (each cell FRAC within the basin) # 1000000 is unit converter from m2 to km2
grid_cell_size <- area(NWPGridFullArea, na.rm=TRUE, weights=FALSE)
grid_cell_size <- grid_cell_size[!is.na(grid_cell_size)]
# compute area [km2]
NominalGrid_Size <- median(grid_cell_size)
#
GridArea <- 1000000 * grid_drain_area * (NominalGrid_Size/as.matrix(NWPGridFullArea))
GridArea1 <- raster(GridArea/1000000)
extent(GridArea1) <- extent(nwp_grid)
dim(GridArea1) <- dim(nwp_grid)
res(GridArea1) <- res(nwp_grid)
crs(GridArea1) <- crs(nwp_grid)
#
writeRaster(GridArea1, "grid_area_weight.tif", datatype="FLT4S", overwrite=TRUE)
#
## Drainage area in km2 
system("mpiexec -n 8 AreaD8 -p grid_fdirp.tif -ad8 grid_da_ad8.tif -wg grid_area_weight.tif -nc")
DA <- as.matrix(raster("grid_da_ad8.tif"))
DA[is.na(DA)] <- 0
#
## Bankfull cross-section area of river channel in m2
Bankfull  <- 0.1667 * DA + 0.1     # BA <- (11.0 + 0.43*(DA))
#
## define river classes (IAK) using drainage area (DA)
IAK <- matrix(0, NumRow, NumCol)
LogDA <- log(as.numeric(DA))
LogDA[LogDA == -Inf] <- NaN
LogDA1 <- as.vector(LogDA)
RR <- integer(IAKmax - 1)
for (i in 1:(IAKmax - 1)) {
RR[i] = quantile(LogDA1, probs = i/IAKmax, na.rm = TRUE)
}
# Apply the background field (filtered by NaN to target active grid cells in the domain).
IAK[!is.nan(LogDA)] <- IAKmax
# Apply IAK values by comparing the transformed values in 'LogDA' to the 'RR' thresholds the (to (IAKmax - 1), IAKmax is applied as the background field above).
for (i in 1:(IAKmax - 1)) {
  IAK[LogDA > RR[i]] = IAKmax - i
}
#
## The internal slope in each grid (Land slope in m/m) 
slp_fdir2 <- mask(slp_fdir1, facc_at_outlet1)
IntSlope <- as.matrix(aggregate(slp_fdir2, fact = ResFactor, fun = mean, na.rm=TRUE))
# IntSlope1 <- zonal(slp_fdir2, nwp_zone, fun='mean', na.rm=TRUE)
# IntSlope <- matrix(IntSlope1[,"mean"], NumRow, NumCol, byrow = T)
#
## Number of channels: Number of equally sized channels traversing an element (No. of channels draining through the cell)
# Can be adjusted based on the stream order file
MaxSord <- maxValue(mask(stream_order1, facc_at_outlet1))
channelclass1 <- matrix(c(MaxSord-100, MaxSord-4, MaxSord-3, MaxSord-2, MaxSord-1, MaxSord-4, MaxSord-3, MaxSord-2, MaxSord-1, MaxSord, 1, 2, 3, 4, 5), nr = 5, nc = 3)
channelclass2 <- matrix(c(0, 1, 2, 3, 4, 1, 2, 3, 4, 5, 5, 4, 3, 2, 1), nr = 5, nc = 3)
#
stord1 <- reclassify(mask(stream_order1, facc_at_outlet1), channelclass1)
stord <- aggregate(stord1, fact = ResFactor, fun = max, na.rm=TRUE) 
#
Chnl <- as.matrix(reclassify(stord, channelclass2))
# # Other option may be to use IAK for channel density as well.
# Chnl <- IAK
#
## Longest flow path for channel length / slope calculation
flow_path_mask0 <- matrix(NA,nrow(fdir2),ncol(fdir2))
for (i in 1:8) {
  xy <- which(fdir2 == i, arr.ind=TRUE)
  xz <- xy
  xz[,1] <- (xy[,1] + as.integer(FdirNumber[i,1]))
  xz[,2] <- (xy[,2] + as.integer(FdirNumber[i,2]))
  flow_path_mask0[xz] <- 1
}
#
flow_path_mask <- raster(flow_path_mask0)
extent(flow_path_mask) <- extent(nwp_zone)
dim(flow_path_mask) <- dim(nwp_zone)
res(flow_path_mask) <- res(nwp_zone)
crs(flow_path_mask) <- crs(nwp_zone)
#
flow_path2 <- mask(flow_path1,filldem4_facc4)
flow_path3 <- mask(flow_path1,flow_path_mask)
ChnlLength1 <- as.matrix(aggregate(flow_path2, fact = ResFactor, fun = mean, na.rm=TRUE))
ChnlLength1[ChnlLength1 == 0] <- 100 # 100 assumed to be the minimum length required to initiate overland flow
ChnlLength2 <- as.matrix(aggregate(flow_path3, fact = ResFactor, fun = mean, na.rm=TRUE))
ChnlLength1[is.na(ChnlLength1)] <- 0
ChnlLength2[is.na(ChnlLength2)] <- 0
ChnlLength <- (ChnlLength1 - ChnlLength2) 
#
# ## Similar to Green Kenue methods but not tested yet ####
# total_length <- 1000 * as.matrix(NWPGridFullArea) * sqrt(NominalGrid_Size) / NominalGrid_Size
# channel_length <- matrix(0, nrow = NumRow, ncol = NumCol, byrow = T)
# channel_elevation = as.matrix(aggregate(filldem1, fact = ResFactor, fun = mean, na.rm=TRUE))
# channel_slope <- matrix(0, nrow = NumRow, ncol = NumCol, byrow = T)
# #
# for (i in 1:nrow(channel_length)) {
#   for (j in 1:ncol(channel_length)) {
#     
#     if ((i-1 >= 1) & (j-1 >= 1) & (i+1 <= NumRow) & (j+1 <= NumCol)) {
#       if (is.na(grid_fdir[i,j])) {
#         channel_length[i,j] <- 0 }
#       else if (grid_fdir[i,j] == 1) {
#         channel_length[i,j] <- (total_length[i,j])
#         channel_slope[i,j] <- ((channel_elevation[i,j] - channel_elevation[i+1,j])/channel_length[i,j])}
#       else if (grid_fdir[i,j] == 2) {
#         channel_length[i,j] <- (sqrt(2)*total_length[i,j])
#         channel_slope[i,j] <- ((channel_elevation[i,j] - channel_elevation[i+1,j+1])/channel_length[i,j])}
#       else if (grid_fdir[i,j] == 3) {
#         channel_length[i,j] <- (total_length[i,j])
#         channel_slope[i,j] <- ((channel_elevation[i,j] - channel_elevation[i,j+1])/channel_length[i,j])}
#       else if (grid_fdir[i,j] == 4) {
#         channel_length[i,j] <- (sqrt(2)*total_length[i,j])
#         channel_slope[i,j] <- ((channel_elevation[i,j] - channel_elevation[i+1,j-1])/channel_length[i,j])}
#       else if (grid_fdir[i,j] == 5) {
#         channel_length[i,j] <- (total_length[i,j])
#         channel_slope[i,j] <- ((channel_elevation[i,j] - channel_elevation[i-1,j])/channel_length[i,j])}
#       else if (grid_fdir[i,j] == 6) {
#         channel_length[i,j] <- (sqrt(2)*total_length[i,j])
#         channel_slope[i,j] <- ((channel_elevation[i,j] - channel_elevation[i-1,j-1])/channel_length[i,j])}
#       else if (grid_fdir[i,j] == 7) {
#         channel_length[i,j] <- (total_length[i,j])
#         channel_slope[i,j] <- ((channel_elevation[i,j] - channel_elevation[i,j-1])/channel_length[i,j])}
#       else if (grid_fdir[i,j] == 8) {
#         channel_length[i,j] <- (sqrt(2)*total_length[i,j])
#         channel_slope[i,j] <- ((channel_elevation[i,j] - channel_elevation[i+1,j-1])/channel_length[i,j])}
#     }
#   }
# }
## River slope in m/m
filldem5 <- mask(filldem1, facc_at_outlet1)
Elevn_flow_path0 <- as.matrix(aggregate(filldem5, fact = ResFactor, fun = max, na.rm=TRUE))
Elevn_flow_path1 <- mask(filldem5,filldem4_facc4)
Elevn_flow_path2 <- mask(filldem5,flow_path_mask)
Elevn_flow_path3 <- as.matrix(aggregate(Elevn_flow_path1, fact = ResFactor, fun = mean, na.rm=TRUE))
Elevn_flow_path4 <- as.matrix(aggregate(Elevn_flow_path2, fact = ResFactor, fun = mean, na.rm=TRUE))
xz <- which(is.na(Elevn_flow_path4), arr.ind=TRUE)
Elevn_flow_path4[xz] <- Elevn_flow_path0[xz]
ChnlSlope <- (Elevn_flow_path4 - Elevn_flow_path3)/ChnlLength
ChnlSlope[is.na(ChnlSlope)] <- 0.001  # Some grid has zero flow length as the grid area is so small and need to be corrected by setting the minimum slope or remove the small basin area in grid area
ChnlSlope[is.infinite(ChnlSlope)] <- 0.001
ChnlSlope[ChnlSlope < Min_Chanel_Slope] <- Min_Chanel_Slope
#
## Reach number for lake, reservoir and/or external routing: user can assign number from 1,2,3...with respect to the number of reservoir / lakes 
Reach  <- matrix(0, nrow = NumRow, ncol = NumCol)
# ###########################################################################
# Collecting the drainage database parameters from the above steps
grid_rank[is.na(grid_rank)] <- 0
#
grid_next[is.na(grid_next)] <- 0
#
DA[is.na(DA)] <- 0
#
ChnlLength[is.na(ChnlLength)]  <- 0
#
Bankfull[is.na(Bankfull)] <- 0
#
ChnlSlope[is.na(ChnlSlope)] <- 0
#
Elev[is.na(Elev)] <- 0
#
IAK[is.na(IAK)] <- 0
#
IntSlope[is.na(IntSlope)] <- 0
#
Chnl[is.na(Chnl)] <- 0
#
Reach[is.na(Reach)] <- 0
#
GridArea[is.na(GridArea)] <- 0
##
drainagedata <- rbind((apply(grid_rank, 2, rev)), (apply(grid_next, 2, rev)), (apply(DA, 2, rev)), (apply(Bankfull, 2, rev)), 
                      (apply(ChnlSlope, 2, rev)), (apply(Elev, 2, rev)), (apply(ChnlLength, 2, rev)), (apply(IAK, 2, rev)),
                      (apply(IntSlope, 2, rev)), (apply(Chnl, 2, rev)), (apply(Reach, 2, rev)), (apply(GridArea, 2, rev)))
####################################################################################################################################
### Model discretization and GRUs creation: Combine the land cover, slope and aspect for mountain MESH version and land cover only
# for Original MESH version ########################################################################################################
####################################################################################################################################
xorigin_value <- paste(":xOrigin                 ", LLXcorner, sep = "")
yorigin_value <- paste(":yOrigin                 ", LLYcorner, sep = "")
xcount_value <- paste(":xCount                   ", NumCol, sep = "")
ycount_value <- paste(":yCount                   ", NumRow, sep = "")
xdelta_value <- paste(":xDelta                   ", XRes, sep = "")
ydelta_value <- paste(":yDelta                   ", YRes, sep = "")
###################################################################################################################################
basin_dem <- mask(filldem1,facc_at_outlet1)
basin_landcover <- mask(land_cover,facc_at_outlet1)
## Calculate Slope and aspect
slope <- terrain(crop(domain_dem, nwp_zone), opt = "slope", unit='degrees', neighbors=8)
aspect <- terrain(crop(domain_dem, nwp_zone), opt = "aspect", unit='degrees', neighbors=8, flatAspect = NA)
#
if (MESHVersion == "Original") {
   grus1 <- basin_landcover }
##################################################################################################################################
if (MESHVersion == "Mountain") {
## Mask the slope and aspect by basin
basin_slope <- mask(slope, facc_at_outlet1)
basin_aspect <- mask(aspect, facc_at_outlet1)
# writeRaster(basin_slope, "R_Slope.tif", datatype="FLT4S", overwrite=TRUE)
# writeRaster(basin_aspect, "R_Aspect.tif", datatype="FLT4S", overwrite=TRUE)
#
## Calculate the Sine and cosine of aspect expressed in radian for the calculation of the zonal GRUs mean aspect
basin_sinaspect <- sin(basin_aspect*pi/180)
basin_cosaspect <- cos(basin_aspect*pi/180)
######################### Reclassifying slope and aspect ####################################################################
#
basin_slopeclass <- reclassify(basin_slope, slopeclass)
basin_aspectclass <- reclassify(basin_aspect, aspectclass)
#
### Model discretization and GRUs creation: Combine the land cover, slope and aspect #############################################
grus1 <- basin_landcover
#
grus1[basin_landcover == 1] <- 1
grus1[basin_landcover == 2] <- 2
#
grus1[basin_landcover == 3 & basin_slopeclass == 2 & basin_aspectclass == 2] <- 3
grus1[basin_landcover == 3 & basin_slopeclass == 2 & basin_aspectclass == 1] <- 4
grus1[basin_landcover == 3 & basin_slopeclass == 1 & basin_aspectclass <= 2] <- 5
grus1[basin_landcover == 3 & basin_slopeclass == 1 & basin_aspectclass <= 1] <- 5
#
grus1[basin_landcover == 4] <- 6
grus1[basin_landcover == 5] <- 7
#
grus1[basin_landcover == 6 & basin_slopeclass == 2 & basin_aspectclass == 2] <- 8
grus1[basin_landcover == 6 & basin_slopeclass == 2 & basin_aspectclass == 1] <- 9
grus1[basin_landcover == 6 & basin_slopeclass == 1 & basin_aspectclass <= 2] <- 10
grus1[basin_landcover == 6 & basin_slopeclass == 1 & basin_aspectclass <= 1] <- 10
#
grus1[basin_landcover == 7] <- 11
#
}
#
# writeRaster(grus1, "GRUs_Produced.tif", datatype="INT2S", overwrite=TRUE)
#
## Polishing GRUs that have less or equal to the minimum land cover fraction in the modeling grid ################################
TotalGrid1 <- zonal(grus1, nwp_zone, 'count', na.rm=TRUE)
TotalGrid <- as.vector(TotalGrid1[,"count"])
GRUGrid_Frac <- matrix(nrow = NumRow*NumCol, ncol = (1 + maxValue(grus1)))
#
for (i in 1:maxValue(grus1)) {
  grus_mask <- grus1
  grus_mask[grus_mask != i] <- NA
  GRU_Frac <- zonal(grus_mask, nwp_zone, 'count', na.rm=TRUE)
  GRUGrid_Frac[,i] <- (as.vector(GRU_Frac[,"count"])/TotalGrid)
}
GRUGrid_Frac[is.na(GRUGrid_Frac)] <- 0
#
PolishedGRUs <- raster(matrix(0, nrow = nrow(grus1), ncol = ncol(grus1), byrow = T))
extent(PolishedGRUs) <- extent(grus1)  
res(PolishedGRUs) <- res(grus1)
crs(PolishedGRUs) <- crs(grus1)
#
for (i in 1:maxValue(grus1)) {
grusmask <- grus1
grusmask[grusmask != i] <- NA
xxGRU_Frac <- GRUGrid_Frac[,i]
  if (i == 2) {                                          # GRU 2 must be glacier otherwise (i == 2) needs changes
     xxGRU_Frac[xxGRU_Frac < Minimum_Glacier_Fraction] <- 0 
     } else {
     xxGRU_Frac[xxGRU_Frac < Minimum_GRU_Fraction] <- 0     
     }
   xxGRU_Frac[xxGRU_Frac != 0] <- i
   yy <- disaggregate(raster(matrix(xxGRU_Frac, nrow = NumRow, ncol = NumCol, byrow = T)), fact = ResFactor, fun = mean, na.rm=TRUE)
   extent(yy) <- extent(grus1)  
   res(yy) <- res(grus1)
   crs(yy) <- crs(grus1)
   PolishedGRUs <- PolishedGRUs + mask(yy, grusmask, inverse=FALSE, updatevalue=0) 
}
#
PolishedGRUs[PolishedGRUs == 0] <- NA
################# Calculating fraction of GRUs in modeling grid ############################################
grus2 <- zonal(PolishedGRUs, nwp_zone, 'count', na.rm=TRUE)
grus3 <- as.vector(grus2[,"count"])
grus4 <- matrix(nrow = NumRow*NumCol, ncol = (1 + maxValue(grus1)))
#
for (i in 1:maxValue(grus1)) {
  grus_mask <- PolishedGRUs
  grus_mask[grus_mask != i] <- NA
  grus5 <- zonal(grus_mask, nwp_zone, 'count', na.rm=TRUE)
  grus6 <- as.vector(grus5[,"count"])
  grus4[,i] <- (grus6/grus3)
}
grus_frac_r2c <- matrix(nrow = NumRow*(1 + maxValue(grus1)), ncol = NumCol)
#
for (i in 1:(maxValue(grus1)+1)) {
  grus7 <- matrix(grus4[,i], nrow(nwp_grid), ncol(nwp_grid), byrow = T)
  grus_frac_r2c[(NumRow*i-(NumRow-1)):(NumRow*i),] <- apply(grus7, 2, rev)
}
grus_frac_r2c[is.na(grus_frac_r2c)] <- 0
##
grus_frac_r2c <- matrix(nrow = nrow(nwp_grid)*(1 + maxValue(grus1)), ncol = ncol(nwp_grid))
for (i in 1:(maxValue(grus1)+1)) {
  grus5 <- matrix(grus4[,i], nrow(nwp_grid), ncol(nwp_grid), byrow = T)
  grus_frac_r2c[(nrow(nwp_grid)*i-(nrow(nwp_grid)-1)):(nrow(nwp_grid)*i),] <- apply(grus5, 2, rev)
}
grus_frac_r2c[is.na(grus_frac_r2c)] <- 0
### Combine the drainage data and fraction of GRUs for MESH_drainage_database file
dranage_database <- rbind(drainagedata,grus_frac_r2c)
########################################################################################################
## Creating MESH drainage database file and save it as MESH_drainage_database.r2c
########################################################################################################
## Add source file name to header and required information for the drainage database file header
##
NominalGrid_Size <- sqrt(NominalGrid_Size)*1000
NominalGrid_Size_Value <- paste(":NominalGridSize_AL   ",NominalGrid_Size, sep = "")
contour_interval <- 1
contour_interval_value <- paste(":ContourInterval   ",contour_interval, sep = "")
Impervious_Area <- 0
Impervious_Area_value <- paste(":ImperviousArea   ",Impervious_Area, sep = "")
Class_Count <- maxValue(grus1) + 1
class_count_value <- paste(":ClassCount   ",Class_Count, sep = "")
NumRiver_Classes <- 1
NumRiver_Classes_value <- paste(":NumRiverClasses   ",NumRiver_Classes, sep = "")
Elev_Conversion <- 1
Elev_Conversion_value <- paste(":ElevConversion   ",Elev_Conversion, sep = "")
Total_NumOfGrids_value <- paste(":TotalNumOfGrids   ",(length(which(grid_rank != 0))), sep = "")  
NumGridsIn_Basin_value <- paste(":NumGridsInBasin   ",(length(which(grid_next != 0))), sep = "") 
Minimum_Slope <- Min_Chanel_Slope
Minimum_Slope_value <- paste(":MinimumSlope   ",Minimum_Slope, sep = "")
ProjectionName <- "LATLONG"
DatumName <- "WGS84"
#####
projection_value <- ":Projection              LATLONG   "
Ellipsoid_value <- ":Ellipsoid               WGS84 "
#
## Write header #####################################################################
header1 <- "#########################################################################"
filetype <- ":FileType r2c  ASCII  EnSim 1.0"
header1 <- c(header1, filetype)
owner <- "# National Research Council Canada (c) 1998-2014"
header1 <- c(header1, owner)
datatype <- "# DataType                 2D Rect Cell"
header1 <- c(header1, datatype, "#")
application <- ":Application    MMESH  "
header1 <- c(header1, application)
version <- ":Version    1.0.0"
header1 <- c(header1, version)
written_by <- ":WrittenBy    MMESHr"
header1 <- c(header1, written_by)
creation_date <- paste(":CreationDate ", date(), sep = "")
header1 <- c(header1, creation_date, "#", "#------------------------------------------------------------------------")
# add source file name to header
sourcefile <- ":SourceFile              DEM from Hydroshed and Land Cover from COE"
header1 <- c(header1, sourcefile)
header1 <- c(header1, NominalGrid_Size_Value)
header1 <- c(header1, contour_interval_value)
header1 <- c(header1, Impervious_Area_value)
header1 <- c(header1, class_count_value)
header1 <- c(header1, NumRiver_Classes_value)
header1 <- c(header1, Elev_Conversion_value)        #  1 if elevations are in S.I. Units (meters), or 0.305 for Imperial Units (feet)
header1 <- c(header1, Total_NumOfGrids_value)
header1 <- c(header1, NumGridsIn_Basin_value)
header1 <- c(header1, Minimum_Slope_value, "#", "#")
#####
header1 <- c(header1, projection_value)
header1 <- c(header1, Ellipsoid_value, "#")
header1 <- c(header1, xorigin_value, yorigin_value, "#")
#####
attributename <- ":AttributeName 1 Rank"
attributetype <- ":AttributeType 1 integer"
header1 <- c(header1, attributename, attributetype)
attributename <- ":AttributeName 2 Next"
attributetype <- ":AttributeType 2 integer"
header1 <- c(header1, attributename, attributetype)
attributename <- ":AttributeName 3 DA"
attributeunit <- ":AttributeUnits 3 km^2"
header1 <- c(header1, attributename, attributeunit)
attributename <- ":AttributeName 4 Bankfull"
header1 <- c(header1, attributename)
attributename <- ":AttributeName 5 ChnlSlope"
header1 <- c(header1, attributename)
attributename <- ":AttributeName 6 Elev"
attributeunit <- ":AttributeUnits 6 m"
header1 <- c(header1, attributename, attributeunit)
attributename <- ":AttributeName 7 ChnlLength"
attributeunit <- ":AttributeUnits 7 m"
header1 <- c(header1, attributename, attributeunit)
attributename <- ":AttributeName 8 IAK"
header1 <- c(header1, attributename)
attributename <- ":AttributeName 9 IntSlope"
header1 <- c(header1, attributename)
attributename <- ":AttributeName 10 Chnl"
header1 <- c(header1, attributename)
attributename <- ":AttributeName 11 Reach"
header1 <- c(header1, attributename)
attributename <- ":AttributeName 12 GridArea"
attributeunit <- ":AttributeUnits 12 m^2"
header1 <- c(header1, attributename, attributeunit)
##
for (i in 1:(1+maxValue(grus1))) {
  j = i + 12
  attributename <- paste(":AttributeName",j,"Class",i,"(",i,")")
  attributeunit <- paste(":AttributeUnits ",j," 0-1")
  header1 <- c(header1, attributename, attributeunit)
}
###
header1 <- c(header1, "#", xcount_value, ycount_value, xdelta_value, ydelta_value, "#")
header1 <- c(header1, ":EndHeader")
##
eol <- "\n"
con <- file(paste0(BasinName,MESHVersion,"_MESH_drainage_database.r2c"), open = "w")
writeLines(header1, con = con, sep = eol)
close(con)

for (row in 1:nrow(dranage_database)) {
  dranage_database1 <- formatC(dranage_database[row, ],
                               digits = 6, width = 1,
                               format = "f")
  cat(dranage_database1, "\n", sep = "  ", file = paste0(BasinName,MESHVersion,"_MESH_drainage_database.r2c"), append = TRUE)
 }
#############################################################################################################################
## Calculate Slope, aspect, delta (difference in elevation between the nwp_elevation and model_elevation)
domain_dem1 <- mask(crop(domain_dem, nwp_zone), PolishedGRUs)
#
elev_max1 <- aggregate(domain_dem1, fact = ResFactor, fun = max, na.rm=TRUE)
elev_max <- resample(elev_max1, domain_dem1, method="ngb")
#
elev_mean1 <- aggregate(domain_dem1, fact = ResFactor, fun = mean, na.rm=TRUE)
elev_mean <- resample(elev_mean1, domain_dem1, method="ngb")
#
# elev_max1 <- zonal(domain_dem1, nwp_zone, 'max', na.rm=TRUE)
# elev_max2 <- matrix(elev_max1[,"max"] ,NumRow, NumCol, byrow = T)
# elev_max3 <- raster(elev_max2)
# extent(elev_max3) <- extent(LLXcorner,(LLXcorner+(NumCol*XRes)),LLYcorner,(LLYcorner+(NumRow*YRes)))
# dim(elev_max3) <- c(NumRow, NumCol)
# res(elev_max3) <- c(XRes, YRes)
# crs(elev_max3) <- crs(domain_dem1) # crs(nwp_grid) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
# elev_max <- resample(elev_max3, domain_dem1, method="ngb")
# #
# elev_mean1 <- zonal(domain_dem1, nwp_zone, 'mean', na.rm=TRUE)
# elev_mean2 <- matrix(elev_mean1[,"mean"] ,NumRow, NumCol, byrow = T)
# elev_mean3 <- raster(elev_mean2)
# extent(elev_mean3) <- extent(LLXcorner,(LLXcorner+(NumCol*XRes)),LLYcorner,(LLYcorner+(NumRow*YRes)))
# dim(elev_mean3) <- c(NumRow, NumCol)
# res(elev_mean3) <- c(XRes, YRes)
# crs(elev_mean3) <- crs(domain_dem1) # crs(nwp_grid) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
# elev_mean <- resample(elev_mean3, domain_dem1, method="ngb")
#
# delta <- (domain_dem - model_nwp_dem)                            # check and fix the extent issues
delta <- (domain_dem1 - elev_mean)
delta_elev_max <- (delta / elev_max)
#
## Weight for grid curvature calculation for wind speed corrections
weight1 = matrix(c(0,1,0,1,0,1,0,1,0), nrow=3)
weight2 = matrix(c(1,0,1,0,0,0,1,0,1), nrow=3)
curve1 <- 0.5*focal(crop(domain_dem, nwp_zone), w=weight1,fun=sum, na.rm=FALSE)
curve2 <- 0.5*focal(crop(domain_dem, nwp_zone), w=weight1,fun=sum, na.rm=FALSE) 
dem_cell_area <- area(crop(domain_dem, nwp_zone), na.rm=TRUE, weights=FALSE)
dem_cell_area <- dem_cell_area[!is.na(dem_cell_area)]
dem_cell_length <- 1000*sqrt(median(dem_cell_area))
dem_cell_length <- min(dem_cell_length[!is.na(dem_cell_length)],300)
curve <- 0.25 * (((2.0 * sqrt(2.0) + 2) * crop(domain_dem, nwp_zone)) - (sqrt(2.0)*curve1) - (curve2)) / (2.0 * sqrt(2.0) * dem_cell_length)
#
basin_delta <- mask(delta,PolishedGRUs)
basin_delta_elev_max <- mask(delta_elev_max,PolishedGRUs)
basin_slope0 <- mask(slope,PolishedGRUs)
basin_aspect0 <- mask(aspect,PolishedGRUs)
basin_curve <- mask(curve,PolishedGRUs)
#### Calculate the Sine and cosine of aspect expressed in radian for the calculation of the zonal GRUs mean aspect
basin_sinaspect <- sin(basin_aspect0*pi/180)
basin_cosaspect <- cos(basin_aspect0*pi/180)
#
### Calculate the weighted average values of elevation, slope, aspect, delta and curvature for MESH_Parameters.r2c file
basin_elevn1  <- matrix(nrow = NumRow*maxValue(grus1), ncol = NumCol)
basin_delta1 <- matrix(nrow = NumRow*maxValue(grus1), ncol = NumCol)
basin_delta_elev_max1 <- matrix(nrow = NumRow*maxValue(grus1), ncol = NumCol)
basin_slope1 <- matrix(nrow = NumRow*maxValue(grus1), ncol = NumCol)
basin_aspect1 <- matrix(nrow = NumRow*maxValue(grus1), ncol = NumCol)
basin_curve1 <- matrix(nrow = NumRow*maxValue(grus1), ncol = NumCol)

#
for (i in 1 : maxValue(grus1)) {
  grus_mask <- PolishedGRUs
  grus_mask[grus_mask != i] <- NA
  #
  basin_elevn2 <- mask(domain_dem1, grus_mask)
  basin_delta2 <- mask(basin_delta, grus_mask)
  basin_delta_elev_max2 <- mask(basin_delta_elev_max, grus_mask)
  basin_slope2 <- mask(basin_slope, grus_mask)
  basin_sinaspect2 <- mask(basin_sinaspect, grus_mask)
  basin_cosaspect2 <- mask(basin_cosaspect, grus_mask)
  basin_curve2 <- mask(basin_curve, grus_mask)
#
  basin_elevn3 <- as.matrix(aggregate(basin_elevn2, fact = ResFactor, fun = mean, na.rm=TRUE))
  basin_slope3 <- as.matrix(aggregate(basin_slope2, fact = ResFactor, fun = mean, na.rm=TRUE))
  basin_delta3 <- as.matrix(aggregate(basin_delta2, fact = ResFactor, fun = mean, na.rm=TRUE))
  basin_delta3 <- as.matrix(aggregate(basin_delta2, fact = ResFactor, fun = mean, na.rm=TRUE))
  basin_delta_elev_max3 <- as.matrix(aggregate(basin_delta_elev_max2, fact = ResFactor, fun = mean, na.rm=TRUE))
  basin_sinaspect3 <- as.matrix(aggregate(basin_sinaspect2, fact = ResFactor, fun = mean, na.rm=TRUE))
  basin_cosaspect3 <- as.matrix(aggregate(basin_cosaspect2, fact = ResFactor, fun = mean, na.rm=TRUE))
  basin_aspect2 <- ((360 + (atan2(basin_sinaspect3,basin_cosaspect3))*(180/pi)) %% 360)
  basin_curve3 <- as.matrix(aggregate(basin_curve2, fact = ResFactor, fun = mean, na.rm=TRUE))
# #  
#   basin_elevn3 <- zonal(basin_elevn2, nwp_zone, 'mean', na.rm=TRUE)
#   basin_elevn3 <- matrix(basin_elevn3[,"mean"] ,NumRow, NumCol, byrow = T)
#   basin_delta3 <- zonal(basin_delta2, nwp_zone, 'mean', na.rm=TRUE)
#   basin_delta3 <- matrix(basin_delta3[,"mean"], NumRow, NumCol, byrow = T)
#   basin_delta_elev_max3 <- zonal(basin_delta_elev_max2, nwp_zone, 'mean', na.rm=TRUE)
#   basin_delta_elev_max3 <- matrix(basin_delta_elev_max3[,"mean"], NumRow, NumCol, byrow = T)    
#   basin_slope3 <- zonal(basin_slope2, nwp_zone, 'mean', na.rm=TRUE)
#   basin_slope3 <- matrix(basin_slope3[,"mean"], NumRow, NumCol, byrow = T)
#   basin_sinaspect3 <- zonal(basin_sinaspect2, nwp_zone, 'mean', na.rm=TRUE)
#   basin_sinaspect3 <- matrix(basin_sinaspect3[,"mean"], NumRow, NumCol, byrow = T)
#   basin_cosaspect3 <- zonal(basin_cosaspect2, nwp_zone, 'mean', na.rm=TRUE)
#   basin_cosaspect3 <- matrix(basin_cosaspect3[,"mean"], NumRow, NumCol, byrow = T)
#   basin_aspect2 <- ((360 + (atan2(basin_sinaspect3,basin_cosaspect3))*(180/pi)) %% 360)
#   #basin_aspect2 <- ((360 + (atan(basin_sinaspect3/basin_cosaspect3))*(180/pi)) %% 360)
#   basin_curve3 <- zonal(basin_curve2, nwp_zone, 'mean', na.rm=TRUE)
#   basin_curve3 <- matrix(basin_curve3[,"mean"], NumRow, NumCol, byrow = T)
#
  basin_elevn1[(NumRow*i-(NumRow-1)):(NumRow*i),] <- apply(basin_elevn3, 2, rev)
  basin_delta1[(NumRow*i-(NumRow-1)):(NumRow*i),] <- apply(basin_delta3, 2, rev)
  basin_delta_elev_max1[(NumRow*i-(NumRow-1)):(NumRow*i),] <- apply(basin_delta_elev_max3, 2, rev)
  basin_slope1[(NumRow*i-(NumRow-1)):(NumRow*i),] <- apply(basin_slope3, 2, rev)
  basin_aspect1[(NumRow*i-(NumRow-1)):(NumRow*i),] <- apply(basin_aspect2, 2, rev)
  basin_curve1[(NumRow*i-(NumRow-1)):(NumRow*i),] <- apply(basin_curve3, 2, rev)
}  
#
elevnslopeaspectdeltacurve <- rbind(basin_elevn1,basin_slope1,basin_aspect1,basin_delta1,basin_delta_elev_max1,basin_curve1)
elevnslopeaspectdeltacurve[is.na(elevnslopeaspectdeltacurve)] <- 0
#
## Soil data for MESH_parameters file
Clay1 <- as.matrix(aggregate(mask(domain_Clay1, facc_at_outlet1), fact = ResFactor, fun = mean, na.rm=TRUE))
Clay2 <- as.matrix(aggregate(mask(domain_Clay2, facc_at_outlet1), fact = ResFactor, fun = mean, na.rm=TRUE))
Clay3 <- as.matrix(aggregate(mask(domain_Clay3, facc_at_outlet1), fact = ResFactor, fun = mean, na.rm=TRUE))
#
Sand1 <- as.matrix(aggregate(mask(domain_Sand1, facc_at_outlet1), fact = ResFactor, fun = mean, na.rm=TRUE))
Sand2 <- as.matrix(aggregate(mask(domain_Sand2, facc_at_outlet1), fact = ResFactor, fun = mean, na.rm=TRUE))
Sand3 <- as.matrix(aggregate(mask(domain_Sand3, facc_at_outlet1), fact = ResFactor, fun = mean, na.rm=TRUE))
gridslope0 <- as.matrix(aggregate(basin_slope, fact = ResFactor, fun = mean, na.rm=TRUE))
#
# Clay1 <- zonal(domain_Clay1, nwp_zone, 'mean', na.rm=TRUE)
# Clay1 <- matrix(Clay1[,"mean"], NumRow, NumCol, byrow = T)
# Clay2 <- zonal(domain_Clay2, nwp_zone, 'mean', na.rm=TRUE)
# Clay2 <- matrix(Clay2[,"mean"], NumRow, NumCol, byrow = T)
# Clay3 <- zonal(domain_Clay3, nwp_zone, 'mean', na.rm=TRUE)
# Clay3 <- matrix(Clay3[,"mean"], NumRow, NumCol, byrow = T)
# #
# Sand1 <- zonal(domain_Sand1, nwp_zone, 'mean', na.rm=TRUE)
# Sand1 <- matrix(Sand1[,"mean"], NumRow, NumCol, byrow = T)
# Sand2 <- zonal(domain_Sand2, nwp_zone, 'mean', na.rm=TRUE)
# Sand2 <- matrix(Sand2[,"mean"], NumRow, NumCol, byrow = T)
# Sand3 <- zonal(domain_Sand3, nwp_zone, 'mean', na.rm=TRUE)
# Sand3 <- matrix(Sand3[,"mean"], NumRow, NumCol, byrow = T)
#
# Prepare grided drainage density
drain_net[drain_net == 0] <- NA
basin_DN <- mask(crop(drain_net, nwp_zone), facc_at_outlet1)
basin_nwp_zone <- mask(nwp_zone, facc_at_outlet1)
#
dn_cell_length <- sqrt(2) * sqrt(area(basin_DN, na.rm=TRUE, weights=FALSE))
grid_area_km2 <- area(basin_nwp_zone, na.rm=TRUE, weights=FALSE)
# zonal_drnr_km <- zonal(dn_cell_length, nwp_zone, 'sum', na.rm=TRUE)
# zonal_area_km2 <- zonal(grid_area_km2, nwp_zone, 'sum', na.rm=TRUE)
zonal_drnr_km <- aggregate(dn_cell_length, fact = ResFactor, fun = sum, na.rm=TRUE)
zonal_drnr_km[is.na(zonal_drnr_km)] <- 0.1  # when the grid has no stream use the minimum length of 0.1 km (100 meter) for overland flow.
zonal_area_km2 <- aggregate(grid_area_km2, fact = ResFactor, fun = sum, na.rm=TRUE)
drainge_dens1 <- zonal_drnr_km / zonal_area_km2
drainge_dens1[is.na(drainge_dens1)] <- 0
drainge_dens <- apply(as.matrix(drainge_dens1), 2, rev)
#
Clay11 <- apply(Clay1, 2, rev)
Clay22 <- apply(Clay2, 2, rev)
Clay33 <- apply(Clay3, 2, rev)
#
Sand11 <- apply(Sand1, 2, rev)
Sand22 <- apply(Sand2, 2, rev)
Sand33 <- apply(Sand3, 2, rev)
#
soildatabase <- rbind(Clay11,Clay22,Clay33,Sand11,Sand22,Sand33)
soildatabase[is.na(soildatabase)] <- 0
#
gridslope <- apply(gridslope0, 2, rev)
#
#######################################################################################################
## Creating the Mountain MESH Parameter file and save it as MESH_Parameter.r2c
######################################################################################################
if (MESHVersion == "Mountain") {
### Combine soil, elevation, slope, aspect, delta, delta_elev_max, curve and drainage density for MESH_parameters file
soilelevnslopeaspectdeltacurvedraindens <- rbind(soildatabase,elevnslopeaspectdeltacurve,drainge_dens)
} else {
### Combine Soil, slope and drainage density for MESH_parameters file
soilelevnslopeaspectdeltacurvedraindens <- rbind(soildatabase,gridslope,drainge_dens)
}
### write header ###
header1 <- "#########################################################################"
filetype <- ":FileType r2c  ASCII  EnSim 1.0"
header1 <- c(header1, filetype, "#")
owner <- "# National Research Council Canada (c) 1998-2014"
header1 <- c(header1, owner)
datatype <- "# DataType                 2D Rect Cell"
header1 <- c(header1, datatype, "#")
application <- ":Application    MMESH  "
header1 <- c(header1, application)
version <- ":Version    1.0.0"
header1 <- c(header1, version)
written_by <- ":WrittenBy    MMESHr"
header1 <- c(header1, written_by)
creation_date <- paste(":CreationDate ", date(), sep = "")
header1 <- c(header1, creation_date, "#", "#------------------------------------------------------------------------", "#", "#")
header1 <- c(header1, projection_value)
header1 <- c(header1, Ellipsoid_value, "#")
header1 <- c(header1, xorigin_value, yorigin_value, "#")
##
if (MESHVersion == "Mountain") {
#  
for (i in 1:Number_Soil_Layers) {
  j <- i 
  attributename <- paste(":AttributeName ", j, " Clay	", i, "")
  header1 <- c(header1, attributename)
}
#
for (i in 1:Number_Soil_Layers) {
  j <- j + 1
  attributename <- paste(":AttributeName ", j, " Sand	", i, "")
  header1 <- c(header1, attributename)
}
#
for (i in 1:maxValue(grus1)) {
  j <- j + 1
  attributename <- paste(":AttributeName ", j, " elevation	", i, "")
  header1 <- c(header1, attributename)
}
#
for (i in 1:maxValue(grus1)) {
  j <- j + 1
  attributename <- paste(":AttributeName ", j, " slope	", i, "") 
  header1 <- c(header1, attributename)
}
#
for (i in 1:maxValue(grus1)) {
  j <- j + 1
  attributename <- paste(":AttributeName ", j, " aspect	", i, "")
  header1 <- c(header1, attributename)
}
#
for (i in 1:maxValue(grus1)) {
  j <- j + 1
  attributename <- paste(":AttributeName ", j, " delta	", i, "") 
  header1 <- c(header1, attributename)
}
#
for (i in 1:maxValue(grus1)) {
  j <- j + 1
  attributename <- paste(":AttributeName ", j, " delta_elevmax	", i, "") 
  header1 <- c(header1, attributename)
}
#
for (i in 1:maxValue(grus1)) {
  j <- j + 1
  attributename <- paste(":AttributeName ", j, " curvature	", i, "") 
  header1 <- c(header1, attributename)
}
#
j <- j + 1
attributename <- paste(":AttributeName ", j, " dd	")
header1 <- c(header1, attributename)
attributeunits <- paste(":AttributeUnits ", j, " km km-2 ")
header1 <- c(header1, attributeunits)
#
header1 <- c(header1, "#", xcount_value, ycount_value, xdelta_value, ydelta_value, "#")
header1 <- c(header1, ":EndHeader")

eol <- "\n"
con <- file(paste0(BasinName,MESHVersion,"_MESH_parameters.r2c"), open = "w")
writeLines(header1, con = con, sep = eol)
close(con)

for (row in 1:nrow(soilelevnslopeaspectdeltacurvedraindens)) {
  mountainparams <- formatC(soilelevnslopeaspectdeltacurvedraindens[row, ],
                            digits = 6, width = 1,
                            format = "f")
  cat(mountainparams, "\n", sep = "  ", file = paste0(BasinName,MESHVersion,"_MESH_parameters.r2c"), append = TRUE)
}
} else {
#
for (i in 1:Number_Soil_Layers) {
  j <- i 
  attributename <- paste(":AttributeName ", j, " Clay	", i, "")
  header1 <- c(header1, attributename)
}
#
for (i in 1:Number_Soil_Layers) {
  j <- j + 1
  attributename <- paste(":AttributeName ", j, " Sand	", i, "")
  header1 <- c(header1, attributename)
}
#
j <- j + 1
attributename <- paste(":AttributeName ", j, " xslp	")
header1 <- c(header1, attributename)
attributeunits <- paste(":AttributeUnits ", j, " m/m ")
header1 <- c(header1, attributeunits)
#
j <- j + 1
attributename <- paste(":AttributeName ", j, " dd	")
header1 <- c(header1, attributename)
attributeunits <- paste(":AttributeUnits ", j, " km km-2 ")
header1 <- c(header1, attributeunits)
#
header1 <- c(header1, "#", xcount_value, ycount_value, xdelta_value, ydelta_value, "#")
header1 <- c(header1, ":EndHeader")

eol <- "\n"
con <- file(paste0(BasinName,MESHVersion,"_MESH_parameters.r2c"), open = "w")
writeLines(header1, con = con, sep = eol)
close(con)

for (row in 1:nrow(soilgridslopedraingedens)) {
  mountainparams <- formatC(soilgridslopedraingedens[row, ],
                            digits = 6, width = 1,
                            format = "f")
  cat(mountainparams, "\n", sep = "  ", file = paste0(BasinName,MESHVersion,"_MESH_parameters.r2c"), append = TRUE)
 }
}
###########################################################################################################################
###########################################################################################################################































































# ####################################################################################################################################
# #################### Creating netcdf file format of the MESH drainage database file (MESH_drainage_database.nc) ####################
# #
# grus_frac_netcdf <- array(NA, dim=c(ncol(nwp_grid), nrow(nwp_grid), (1 + maxValue(grus1))))
# for (i in 1:(maxValue(grus1)+1)) {
#   grus6 <- matrix(grus4[,i], nrow(nwp_grid), ncol(nwp_grid), byrow = T)
#   grus_frac_netcdf[,,i] <- array(t(apply(grus6, 2, rev)))
# }
# grus_frac_netcdf[is.na(grus_frac_netcdf)] <- 0
# ##
# # set path and file name
# ncname <- paste0(MESHVersion,"_MESH_drainage_database")   #"MESH_drainage_database"  
# ncfname <- paste(ncpath, ncname, ".nc", sep="")
# dname <- "mddf"  # note: ddf means MESH drainage database file (not temporary)
# # generate lons, lats and set time
# lon <- as.array(seq((LLXcorner + (0.5*XRes)),(URXcorner - (0.5*XRes)),XRes))
# lat <- as.array(seq((LLYcorner + (0.5*YRes)),(URYcorner - (0.5*YRes)),YRes))
# gru <- as.array(seq(1,Class_Count,1))
# nval <- as.array(1)
# ##
# nlon <- length(lon)
# nlat <- length(lat)
# ngru <- length(gru)
# # create and write the netCDF file -- ncdf4 version
# # define dimensions
# londim <- ncdim_def("lon","degrees_east",as.double(lon)) 
# latdim <- ncdim_def("lat","degrees_north",as.double(lat))
# grudim <- ncdim_def("gru","number_of_grus",as.double(gru))
# nvaldim <- ncdim_def("nval","parameter",as.double(nval))
# 
# # define variables
# FillValue <- -1.0
# #
# dlname <- "Rank" 
# rank_def <- ncvar_def("Rank","integer",list(londim,latdim),fillvalue,dlname,prec="single")
# dlname <- "Next" 
# next_def <- ncvar_def("Next","integer",list(londim,latdim),fillvalue,dlname,prec="single")
# dlname <- "DA" 
# da_def <- ncvar_def("DA","km^2",list(londim,latdim),fillvalue,dlname,prec="single")
# dlname <- "Bankfull"
# bankful_def <- ncvar_def("Bankfull","m^2",list(londim,latdim),fillvalue,dlname,prec="single")
# dlname <- "ChnlSlope"
# chnlslope_def <- ncvar_def("ChnlSlope","m/m",list(londim,latdim),fillvalue,dlname,prec="single")
# dlname <- "Elev"  
# elev_def <- ncvar_def("Elev","m",list(londim,latdim),fillvalue,dlname,prec="single")
# dlname <- "ChnlLength"
# chnlength_def <- ncvar_def("ChnlLength","m",list(londim,latdim),fillvalue,dlname,prec="single")
# dlname <- "IAK"
# iak_def <- ncvar_def("IAK","integer",list(londim,latdim),fillvalue,dlname,prec="single")
# dlname <- "IntSlope"
# intslope_def <- ncvar_def("IntSlope","m/m",list(londim,latdim),fillvalue,dlname,prec="single")
# dlname <- "Chnl"
# chnl_def <- ncvar_def("Chnl","integer",list(londim,latdim),fillvalue,dlname,prec="single")
# dlname <- "Reach"
# reach_def <- ncvar_def("Reach","integer",list(londim,latdim),fillvalue,dlname,prec="single")
# dlname <- "GridArea" 
# gridarea_def <- ncvar_def("GridArea","m^2",list(londim,latdim),fillvalue,dlname,prec="single")
# dlname <- "GRUs"
# grus_def <- ncvar_def("GRUs","fraction",list(londim,latdim,grudim),fillvalue,dlname,prec="single")
# dlname <- "Basin_Parameters"
# basin_def <- ncvar_def("Basin_Parameters","1",NULL,NULL,longname=dlname,prec="char")
# 
# # create netCDF file and put arrays
# ncout <- nc_create(ncfname,list(rank_def,next_def,da_def,bankful_def,chnlslope_def,elev_def,chnlength_def,iak_def,intslope_def,chnl_def,reach_def,gridarea_def,grus_def,basin_def),force_v4=TRUE)
# # put variables
# ncvar_put(ncout,rank_def,t(apply(grid_rank, 2, rev)))
# ncvar_put(ncout,next_def,t(apply(grid_next, 2, rev)))
# ncvar_put(ncout,da_def,t(apply(DA, 2, rev)))
# ncvar_put(ncout,bankful_def,t(apply(Bankfull, 2, rev)))
# ncvar_put(ncout,chnlslope_def,t(apply(ChnlSlope, 2, rev)))
# ncvar_put(ncout,elev_def,t(apply(Elev, 2, rev)))
# ncvar_put(ncout,chnlength_def,t(apply(ChnlLength, 2, rev)))
# ncvar_put(ncout,iak_def,t(apply(IAK, 2, rev)))
# ncvar_put(ncout,intslope_def,t(apply(IntSlope, 2, rev)))
# ncvar_put(ncout,chnl_def,t(apply(Chnl, 2, rev)))
# ncvar_put(ncout,reach_def,t(apply(Reach, 2, rev)))
# ncvar_put(ncout,gridarea_def,t(apply(GridArea, 2, rev)))
# ncvar_put(ncout,grus_def,grus_frac_netcdf)
# 
# # put the MESH_drainage_database attributes
# ncatt_put(ncout,"Basin_Parameters","NominalGridSize_AL",NominalGrid_Size)
# ncatt_put(ncout,"Basin_Parameters","ContourInterval",contour_interval)
# ncatt_put(ncout,"Basin_Parameters","ImperviousArea",Impervious_Area)
# ncatt_put(ncout,"Basin_Parameters","ClassCount",Class_Count)
# ncatt_put(ncout,"Basin_Parameters","NumRiverClasses",NumRiver_Classes)
# ncatt_put(ncout,"Basin_Parameters","ElevConversion",Elev_Conversion)
# ncatt_put(ncout,"Basin_Parameters","TotalNumOfGrids",(length(which(!is.na(grid_rank)))))
# ncatt_put(ncout,"Basin_Parameters","NumGridsInBasin",(length(which(!is.na(grid_next)))))
# ncatt_put(ncout,"Basin_Parameters","MinimumSlope",Minimum_Slope)
# ncatt_put(ncout,"Basin_Parameters","Projection",ProjectionName)
# ncatt_put(ncout,"Basin_Parameters","Ellipsoid",DatumName)
# ncatt_put(ncout,"Basin_Parameters","xOrigin",LLXcorner)
# ncatt_put(ncout,"Basin_Parameters","yOrigin",LLYcorner)
# ncatt_put(ncout,"Basin_Parameters","xCount",NumCol)
# ncatt_put(ncout,"Basin_Parameters","yCount",NumRow)
# ncatt_put(ncout,"Basin_Parameters","xDelta",XRes)
# ncatt_put(ncout,"Basin_Parameters","yDelta",YRes)
# # add global attributes
# titlename <- paste("MESH_drainage_database_file")
# ncatt_put(ncout,0,"title",titlename)
# institutionweb <- paste("https://gwf.usask.ca/")
# ncatt_put(ncout,0,"institution",institutionweb)
# datasource <- paste("HydroShed for DEM and COE for land cover")
# ncatt_put(ncout,0,"source",datasource)
# referencepaper <- paste("xx et al. (xxxx) journal xx:xx-xx")
# ncatt_put(ncout,0,"references",referencepaper)
# history <- paste("Z.K. Tesemma", date(), sep=", ")
# ncatt_put(ncout,0,"history",history)
# ConventionsCode <- paste("CF-1.0")
# ncatt_put(ncout,0,"Conventions",ConventionsCode)
# # Get a summary of the created file:
# ncout
# # close the file, writing data to disk
# nc_close(ncout)
# #################################################################################################################################
