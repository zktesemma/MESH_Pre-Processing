#=========================================================================================
# SCRIPT:
# Automated MESH Setup Script
#
# DESCRIPTION:
# This script produces the MESH_drainage_database and MESH_Parameters file (r2c and netcdf) 
# that are required to run both Mountain and Original MESH versions
#
# USER INPUTS:
#   > The modeling domain box that is consistent with the climate forcing data to be used to drive the model
#   > High resolution elevation (dem).
#   > Elevation raster of the NWP climate forcing data if available.
#   > Land cover raster map that covers the study domain not the basin to avoid missing small fractions
#   > streamflow gauging stations locations [Longitude and Latitude].
#   > Reservoirs outlet locations [Longitude and Latitude].
#
# SOFTWARE REQUIREMENTS:
#   > TauDEM
#   > R and the required packages
#
# DISCLAIMER:
#   > This script is made available for others use on an "as is" basis, with no warranty, either expressed or implied, 
#     as to its fitness for any particular purpose.
# AUTHOR: Zelalem Tesemma, GWSI, University of Saskatchewan, zelalem.tesemma@usask.ca
#
# HISTORY: 12-Jan-2020: Developed and tested on Bow River Basin at Banff to check the preliminary work flow
#          10-Feb-2021: Improved and included various options. Tested on Bow River Basin at Banff
#=============================================================================================================
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
library(geosphere)
######################################## User Inputs #######################################################################
# Change the following directory path to your working directory path
setwd("C:/PD/Automated_MESH_Setup_Sample")
ncpath <- "C:/PD/Automated_MESH_Setup_Sample"
#
########### Comment either Mountain or Original MESH model version #########################################################
MESHVersion <- "Mountain"
# MESHVersion <- "Original"
#
InsertedReservoir <- "NoReservoir"
# InsertedReservoir <- "Reservoir"
#
UpdateIrrigation <- "NoIrrigation"
#InsertedIrrigation <- "Irrigation"
#
########### Prepare the location file for Reservoir insertion ############################################################
if (InsertedReservoir == "Reservoir") {
#### Reading Reservoir location shapefile or csv file (reservoirs.shp / reservoirs.csv ) of lat-long of resevoir outlet ##
  reservoirlocation <- read.csv("reservoirs.csv")
  coordinates(reservoirlocation) = ~ LONGITUDE+LATITUDE
  # reservoirlocation <- readOGR(dsn = ".", layer = "approxoutlets")
  # reservoirlocation = SpatialPoints(reservoirlocation)
}
#
######### Set the minimum GRUs fraction in a give cell and the model grid fraction in the watershed boundary ###############
Minimum_Basin_Fraction <- 0.001 # Minimum basin area in modeling grid to be polished (0.1%)
Minimum_GRU_Fraction <- 0.05   # Minimum GRUs fraction in grid area to be polished (1%) (5%)
Minimum_Glacier_Fraction <- 0.0016 # Minimum Glacier GRUs area in the modeling grid to be polished (0.05 km2)
#
# The maximum number of river class you want to have (IAK - River Classification)
IAKmax <- 5
# The maximum channel slope
Min_Chanel_Slope <- 0.001
#
# The minimum threshold flow accumulation area to initiate stream.
MinThresh <- 123 # based on 90 meter resolution dem (0.1km2 = 12, 0.5km2 = 61, 1km2 = 123) # 180 from stream drop analysis
#
###### Name of the drainage basin of study ############################
BasinName <- paste0("BowRiverBasinBanff_GEM_0p125_MinThresh_",MinThresh,"_")
# BasinName <- paste0("BowElbowRiverBasins_GEM_0p125_MinThresh_",MinThresh,"_") 
#
####### Flow direction numbers
FdirNumber <- matrix((as.integer(c(0,-1,-1,-1,0,1,1,1,1,1,0,-1,-1,-1,0,1))),8,2)
#
##### Spatial filtering for aspect to reduce noise spike.
SpatialFilter = matrix(1/9, nrow = 3, ncol = 3)
# SpatialFilter = matrix(1/121, nrow = 11, ncol = 11)
################### Land cover, Slope and Aspect Classification for GRUs Creations ################################################################################################################
##### Land cover and two slope Classes: (0 - 10) and (> 10) and two aspects: North and South Facing can be modified to include other orientations ##############
# # For 4 quadrants (North, East, South and West) (315 < N <= 45), (45 < E <= 135), (135 < S <= 225), (225 < W <= 315) # From GEM
# # For 4 quadrants (North, East, South and West) (292.5 < N <= 67.5), (67.5 < E <= 112.5), (112.5 < S <= 247.5), (247.5 < W <= 292.5) # ESRI ArcGIS
LandCoverClass <- matrix(c(seq(0,6),
                           seq(1,7),
                           1,2,3,4,5,6,7), nr = 7, nc = 3)
###### The value of glacier in the reclassified land cover file 
GlacierInLandCover <- c(2,999,999)
#
SlopeClass <- matrix(c(-1, 10, 
                       10, 90, 
                       1, 2), nr = 2, nc = 3)
#
# ElevnClass <- matrix(c(-99, 0, 
#                        0, 1000, 
#                        1000, 2000,
#                        2000, 3000), nr = 4, nc = 3)
# 
AspectClass <- matrix(c(-1, 90, 270,
                        90, 270, 360,
                        1, 2, 1), nr = 3, nc = 3)
#
# AspectClass <- matrix(c(-1, 315, 135, 45, 225,
#                         45, 360, 225, 135, 315,
#                         1, 1, 2, 3, 4), nr = 5, nc = 3)
# #
GRUsClass <- matrix(c(seq(0,20),
                      seq(1,21),
                      1,1,1,2,2,2,3,4,5,6,6,6,7,7,7,8,9,10,11,11,11), nr = 21, nc = 3)
#                      1,1,1,1,1,2,2,2,2,2,3,4,5,6,7,8,8,8,8,8,9,9,9,9,9,10,11,12,13,14,15,15,15,15,15), nr = 35, nc = 3)
########## Set the climate forcing grid domain ###################################################################################################
LLXcorner <- -116.75
LLYcorner <- 50.00
NumRow <- 15
NumCol <- 29
XRes <- 0.125
YRes <- 0.125
URXcorner <- (LLXcorner + (NumCol*XRes)) 
URYcorner <- (LLYcorner + (NumRow*YRes))
#
### Import the Hydroshed / MERIT DEM, NWP climate forcing elevation data for model domain ####################
domain_dem <- raster("domain_HydroShed_dem.tif")
# domain_dem <- raster("domain_MERIT_dem.tif")
#### Produce a raster for the NWP grid which is equal to the modeling grid #########################################
NumGrids <- matrix(seq(1,(NumRow*NumCol),1), NumRow, NumCol, byrow = T) 
nwp_grid <- raster(NumGrids)
extent(nwp_grid) <- extent(LLXcorner,(LLXcorner+(NumCol*XRes)),LLYcorner,(LLYcorner+(NumRow*YRes)))
dim(nwp_grid) <- c(NumRow, NumCol)
res(nwp_grid) <- c(XRes, YRes)
crs(nwp_grid) <- crs(domain_dem) # crs(nwp_grid) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
#
#### Resolving spatial resolution issues ##################################################################################
ResFactor <- round(res(nwp_grid)/res(domain_dem))
ResFactor  <- ResFactor[1]
#
nwp_zone  <- disaggregate(nwp_grid, fact = ResFactor)            # disaggregate nwp_grid
domain_dem1 <- resample(domain_dem, nwp_zone, method="bilinear")
writeRaster(domain_dem1, "resample_domain_dem.tif", datatype="FLT4S", overwrite=TRUE)
#
### Import the Land cover data for the model domain and preprocess it for future use ####################
# domainlandcover <- raster("domain_landcover.tif")
# domain_landcover <- reclassify(domainlandcover, LandCoverClass)
# domainlandcover[domainlandcover == 2] <- 99
# #
# Dissolved_RGI_Glacier_Map  <- readOGR('Dissolved_RGI_Glacier_Map.shp')
# domainglaciercover <- rasterize(Dissolved_RGI_Glacier_Map, domainlandcover)
# domainglaciercover[domainglaciercover != 2] <- 2
# writeRaster(domainglaciercover, "rgi60_glaciers.tif", datatype="INT2S", overwrite=TRUE)
# #
# domainlandcover[domainglaciercover == 2] <- 2
# domainlandcover[domainlandcover == 99] <- 3
# 
# land_cover <- resample(domain_landcover, nwp_zone, method="ngb")
# writeRaster(land_cover, "res30_reclass_domain_landcover_glacier_corrected.tif", datatype="INT2S", overwrite=TRUE)
land_cover <- raster("res30_reclass_domain_landcover_glacier_corrected.tif")
#
#### Read irrigation fraction and segregate the cropland into Irrigated and Non-Irrigated land use ################
if (UpdateIrrigation == "Irrigation") {
domain_irr <- raster("domain_irr_tot.tif")                               # read irrigation raster
domain_irr[domain_irr <= 0] <- NA                     
domain_irrRes <- resample(domain_irr, nwp_zone, method="ngb")            # make it similar scale to land cover
land_cover[land_cover == 4 & domain_irrRes > 0] <- 5                    # assign the irrigated landcover as 20
## The Order of crop should be first Nor irrigated then irrigated
CropInLandCover <- c(4,999,999)           # Non Irrigated
IrrigInLandCover <- c(5,999,999)          # Irrigated
}
#
#### Creating the basin outlet using the lat and long data ##############################################################
# streamgauge <- data.frame(LONGITUDE = c(-115.5717), LATITUDE = c(51.17222))
streamgauge <- data.frame(LONGITUDE = c(-115.5717,-114.501), LATITUDE = c(51.17222,50.996))
coordinates(streamgauge) = ~ LONGITUDE+LATITUDE
crs(streamgauge) <- crs(domain_dem)
raster::shapefile(streamgauge, "approxoutlets.shp", overwrite=TRUE)
#
if (MESHVersion == "Mountain") {
###### Read and preprocess NWP climate forcing elevation data for model domain
# nwp_dem <- raster("domain_GEM_dem.tif")
# nwp_dem <- raster("bow_domain_10km_GEM_Elevation_20031218.tif")
# model_nwp_dem <- resample(nwp_dem, nwp_zone, method="bilinear")
#### Calculate the maximum and mean elevation and delta (difference in elevation between the nwp_elevation (if available) and the high resolution dem)
 elev_max1 <- aggregate(domain_dem1, fact = ResFactor, fun = max, na.rm=TRUE)
 elev_max <- disaggregate(elev_max1, fact = ResFactor, fun = max, na.rm=TRUE)
#  
 elev_mean1 <- aggregate(domain_dem1, fact = ResFactor, fun = mean, na.rm=TRUE)
 elev_mean <- disaggregate(elev_mean1, fact = ResFactor, fun = mean, na.rm=TRUE)
 delta <- (domain_dem1 - resample(elev_mean1, nwp_zone, method="bilinear"))
 delta_elev_max <- ((domain_dem1 - elev_mean) / elev_max)
#  
#### delta if the elevation of the nwp model to be used
# elev_mean <- disaggregate(nwp_dem, fact = ResFactor, fun = mean, na.rm=TRUE)
# delta <- (domain_dem1 - model_nwp_dem)
# delta_elev_max <- ((domain_dem1 - elev_mean) / elev_max)
##############################################################################################################
# elev_band <- domain_dem1
# elev_band[domain_dem1 <= elev_mean] <- 1
# elev_band[domain_dem1 > elev_mean] <- 2
}
#
### Import the percent of clay, sand and Organic Matter contents for the required soil layers #################
NumbSoilData <- 2 # (2 if only %Clay and %Sand) and (3 if %Clay, %Sand and %OM) 
#
##### Read soil data (%Clay, %Sand and %OM) from geotiff files in the following order where Clay1,Clay2,...then Sand1,Sand2,....) 
# rlist <- list.files(path = getwd(), pattern = "domain_Clay.*.tif$")
# for(i in rlist) { assign(unlist(strsplit(i, "[.]"))[1], raster(i)) }
# rlist <- list.files(path = getwd(), pattern = "domain_Sand.*.tif$")
# for(i in rlist) { assign(unlist(strsplit(i, "[.]"))[1], raster(i)) }
# CLAYSANDSOM <- resample(stack(domain_Clay1, domain_Clay2, domain_Clay3, domain_Sand1, domain_Sand2, domain_Sand3), nwp_zone, method="ngb")
#
####### Read soil data (%Clay, %Sand and %OM) from NetCDF file in the following order where Clay1,Clay2,...then Sand1,Sand2,....)
####### DateSet from "Wei Shangguan1 et al., 2014. A global soil data set for earth system modeling"
CLAY <- brick('CLAY.nc', varname="CLAY")
SAND <- brick('SAND.nc', varname="SAND")
#
CLAY1 <- ((0.045/0.1)*CLAY[[1]] + (0.046/0.1)*CLAY[[2]] + (0.009/0.1)*CLAY[[3]])
CLAY2 <- ((0.066/0.25)*CLAY[[3]] + (0.123/0.25)*CLAY[[4]] + (0.061/0.25)*CLAY[[5]])
CLAY3 <- ((0.143/1.946)*CLAY[[5]] + (0.336/1.946)*CLAY[[6]] + (0.554/1.946)*CLAY[[7]] + (0.913/1.946)*CLAY[[8]]) 
#
SAND1 <- ((0.045/0.1)*SAND[[1]] + (0.046/0.1)*SAND[[2]] + (0.009/0.1)*SAND[[3]])
SAND2 <- ((0.066/0.25)*SAND[[3]] + (0.123/0.25)*SAND[[4]] + (0.061/0.25)*SAND[[5]])
SAND3 <- ((0.143/1.946)*SAND[[5]] + (0.336/1.946)*SAND[[6]] + (0.554/1.946)*SAND[[7]] + (0.913/1.946)*SAND[[8]])
#
if (NumbSoilData > 2) {
SOM <- (1/0.58)*0.01*brick('OC.nc', varname="OC")  # To convert Soil organic carbon into Soil organic matter multiply it by a factor of (1/0.58) [Pribyl, 2010].
SOM1 <- ((0.045/0.1)*SOM[[1]] + (0.046/0.1)*SOM[[2]] + (0.009/0.1)*SOM[[3]])
SOM2 <- ((0.066/0.25)*SOM[[3]] + (0.123/0.25)*SOM[[4]] + (0.061/0.25)*SOM[[5]])
SOM3 <- ((0.143/1.946)*SOM[[5]] + (0.336/1.946)*SOM[[6]] + (0.554/1.946)*SOM[[7]] + (0.913/1.946)*SOM[[8]]) 
CLAYSANDSOM <- resample(brick(CLAY1,CLAY2,CLAY3,SAND1,SAND2,SAND3,SOM1,SOM2,SOM3), nwp_zone, method="ngb") 
} else {
CLAYSANDSOM <- resample(brick(CLAY1,CLAY2,CLAY3,SAND1,SAND2,SAND3), nwp_zone, method="ngb")  
}
#
####### Produce the lat long of the center of the forcing grid ####################
YLat <- matrix(0, NumRow, NumCol, byrow = T)
XLon <- matrix(0, NumRow, NumCol, byrow = T)
#
for (i in 1:NumRow) {
  for (j in 1:NumCol) {
    XLon[i,j] <- LLXcorner + (XRes/2) + (j-1)*XRes
    YLat[i,j] <- LLYcorner + (YRes/2) + (i-1)*XRes
  }
}
#
## Average side-length of the routing element (using the WATFLOOD square-grid approach).	
Xlength <- distm (c(XLon[(NumRow/2),1], YLat[(NumRow/2),1]), c(XLon[(NumRow/2),NumCol], YLat[(NumRow/2),NumCol]), fun = distHaversine)
Ylength <- distm (c(XLon[1,(NumCol/2)], YLat[1,(NumCol/2)]), c(XLon[NumRow,(NumCol/2)], YLat[NumRow,(NumCol/2)]), fun = distHaversine)
NominalGridSize <- 0.5 * ((Xlength / (NumCol - 1)) + (Ylength / (NumRow - 1)))
#
###################################################################################################################################################
##################### Few Change after this except for GRUs discretization Section ################################
##################### TauDEM approach (DEM post-processing and Watershed Delineation ##############################
###################################################################################################################################################
# Pitremove 
system("mpiexec -n 12 pitremove -z resample_domain_dem.tif -fel domain_demfel.tif")

# D8 flow directions
system("mpiexec -n 12 D8Flowdir -p domain_demp.tif -sd8 domain_demsd8.tif -fel domain_demfel.tif") #,show.output.on.console=F,invisible=F)

# Contributing area
system("mpiexec -n 12 AreaD8 -p domain_demp.tif -ad8 domain_demad8.tif -nc")

# Grid Network 
system("mpiexec -n 12 Gridnet -p domain_demp.tif -gord domain_demgord.tif -plen domain_demplen.tif -tlen domain_demtlen.tif")

# Threshold
cmd=paste("mpiexec -n 12 Threshold -ssa domain_demad8.tif -src domain_demsrc.tif -thresh ",MinThresh)
system(cmd)

# Move Outlets to fall in the flow accumulation pixels
system("mpiexec -n 12 moveoutletstostreams -p domain_demp.tif -src domain_demsrc.tif -o approxoutlets.shp -om outlet.shp")

# maxdist <- ResFactor 
# cmd=paste("mpiexec -n 12 moveoutletstostreams -p domain_demp.tif -src domain_demsrc.tif -o approxoutlets.shp -om outlet.shp -md",maxdist)
# system(cmd)
####### Check for the moved outlet ########################################################################################################
# outpt=read.shp("Outlet.shp")
# approxpt=read.shp("approxoutlets.shp")
# ssa <- raster("domain_demsrc.tif")
# plot(ssa)
# zoom(ssa)
# points(outpt$shp[2],outpt$shp[3],pch=19,col=2)
# points(approxpt$shp[2],approxpt$shp[3],pch=19,col=4)
#
# Contributing area upstream of outlet
system("mpiexec -n 12 AreaD8 -p domain_demp.tif -o outlet.shp -ad8 domain_demssa.tif -nc")

# Threshold
cmd=paste("mpiexec -n 12 Threshold -ssa domain_demssa.tif -src domain_demsrc.tif -thresh ",MinThresh)
system(cmd)

# # Drop Analysis
# system("mpiexec -n 8 Dropanalysis -p domain_demp.tif -fel domain_demfel.tif -ad8 domain_demad8.tif -ssa domain_demssa.tif -drp logandrp.txt -o outlet.shp -par 5 500 10 0")

# Stream Reach and Watershed
system("mpiexec -n 8 Streamnet -fel domain_demfel.tif -p domain_demp.tif -ad8 domain_demad8.tif -src domain_demsrc.tif -o outlet.shp -ord domain_demord.tif -tree domain_demtree.txt -coord domain_demcoord.txt -net domain_demnet.shp -w domain_demw.tif")
# 
# # Calculate Height Above Nearest Drainage (HAND)
# # DInf flow directions
# system("mpiexec -n 8 DinfFlowdir -ang domain_demang.tif -fel domain_demfel.tif",show.output.on.console=F,invisible=F)
# system("mpiexec -n 8 DinfDistDown -ang domain_demang.tif -fel domain_demfel.tif -src domain_demsrc.tif -m ave v -dd domain_demhand.tif -nc",show.output.on.console=F,invisible=F)
# #
######## Collecting the drainage files generated from previous steps ###
##########################################################################################################################################
filldem1 = raster("domain_demfel.tif")
fdir1 = raster("domain_demp.tif")
facc_at_outlet1 = raster("domain_demssa.tif")
stream_order1 = raster("domain_demord.tif")
flow_path1 = raster("domain_demplen.tif")
drain_net = raster("domain_demsrc.tif")
#
#### Generating flow accumulation for all modeling grid cell independently ################################
nwp_zone1 <- matrix(NA, nrow(nwp_zone), ncol(nwp_zone))
nwp_zone2 <- matrix(NA, nrow(nwp_zone), ncol(nwp_zone))
nwp_zone3 <- matrix(NA, nrow(nwp_zone), ncol(nwp_zone))
nwp_zone4 <- matrix(NA, nrow(nwp_zone), ncol(nwp_zone))
#
for (i in 1:NumRow) {
  for (j in 1:NumCol) {
    if ((i %% 2 != 0) & (j %% 2 != 0)){
      nwp_zone1[(ResFactor*i-(ResFactor-1)):(ResFactor*i), (ResFactor*j-(ResFactor-1)):(ResFactor*j)] <- 1 }
    if ((i %% 2 != 0) & (j %% 2 == 0)){    
      nwp_zone2[(ResFactor*i-(ResFactor-1)):(ResFactor*i), (ResFactor*j-(ResFactor-1)):(ResFactor*j)] <- 1 }
    if ((i %% 2 == 0) & (j %% 2 != 0)){    
      nwp_zone3[(ResFactor*i-(ResFactor-1)):(ResFactor*i), (ResFactor*j-(ResFactor-1)):(ResFactor*j)] <- 1 }
    if ((i %% 2 == 0) & (j %% 2 == 0)){   
      nwp_zone4[(ResFactor*i-(ResFactor-1)):(ResFactor*i), (ResFactor*j-(ResFactor-1)):(ResFactor*j)] <- 1 }
  }
}
#
fdir2 <- raster(as.matrix(fdir1)*nwp_zone1)
extent(fdir2) <- extent(nwp_zone)
dim(fdir2) <- dim(nwp_zone)
res(fdir2) <- res(nwp_zone)
crs(fdir2) <- crs(nwp_zone)
#
fdir3 <- raster(as.matrix(fdir1)*nwp_zone2)
extent(fdir3) <- extent(nwp_zone)
dim(fdir3) <- dim(nwp_zone)
res(fdir3) <- res(nwp_zone)
crs(fdir3) <- crs(nwp_zone)
#
fdir4 <- raster(as.matrix(fdir1)*nwp_zone3)
extent(fdir4) <- extent(nwp_zone)
dim(fdir4) <- dim(nwp_zone)
res(fdir4) <- res(nwp_zone)
crs(fdir4) <- crs(nwp_zone)
#
fdir5 <- raster(as.matrix(fdir1)*nwp_zone4)
extent(fdir5) <- extent(nwp_zone)
dim(fdir5) <- dim(nwp_zone)
res(fdir5) <- res(nwp_zone)
crs(fdir5) <- crs(nwp_zone)
#
writeRaster(fdir2, "domain_demp_1.tif", datatype="INT2S", overwrite=TRUE)
writeRaster(fdir3, "domain_demp_2.tif", datatype="INT2S", overwrite=TRUE)
writeRaster(fdir4, "domain_demp_3.tif", datatype="INT2S", overwrite=TRUE)
writeRaster(fdir5, "domain_demp_4.tif", datatype="INT2S", overwrite=TRUE)
#
##### Calculate Local (within the modeling grid) Contributing area calculation using weight and without weight
cell_area_weight <- area(nwp_zone, na.rm=TRUE, weights=FALSE)
writeRaster(cell_area_weight, "cell_area_weight.tif", datatype="FLT4S", overwrite=TRUE)
system("mpiexec -n 8 AreaD8 -p domain_demp_1.tif -ad8 domain_demad8_1.tif -wg cell_area_weight.tif -nc")
system("mpiexec -n 8 AreaD8 -p domain_demp_2.tif -ad8 domain_demad8_2.tif -wg cell_area_weight.tif -nc")
system("mpiexec -n 8 AreaD8 -p domain_demp_3.tif -ad8 domain_demad8_3.tif -wg cell_area_weight.tif -nc")
system("mpiexec -n 8 AreaD8 -p domain_demp_4.tif -ad8 domain_demad8_4.tif -wg cell_area_weight.tif -nc")
facc1 <- merge(raster("domain_demad8_1.tif"), raster("domain_demad8_2.tif"), raster("domain_demad8_3.tif"), raster("domain_demad8_4.tif"))
#
system("mpiexec -n 8 AreaD8 -p domain_demp_1.tif -ad8 domain_demad8_11.tif -nc")
system("mpiexec -n 8 AreaD8 -p domain_demp_2.tif -ad8 domain_demad8_22.tif -nc")
system("mpiexec -n 8 AreaD8 -p domain_demp_3.tif -ad8 domain_demad8_33.tif -nc")
system("mpiexec -n 8 AreaD8 -p domain_demp_4.tif -ad8 domain_demad8_44.tif -nc")
facc2 <- merge(raster("domain_demad8_11.tif"), raster("domain_demad8_22.tif"), raster("domain_demad8_33.tif"), raster("domain_demad8_44.tif"))
#
#######################################################################################################################################
## Compute the nominal grid area [km2] and fraction of the basin in every modeling grid cell.
NWPGridFullArea <- aggregate(area(nwp_zone, na.rm=TRUE, weights=FALSE), fact = ResFactor, fun = sum, na.rm=TRUE)
BasinAreaInNWPGrid <- aggregate(area(facc_at_outlet1, na.rm=TRUE, weights=FALSE), fact = ResFactor, fun = sum, na.rm=TRUE)
FracBasinAreaInNWPGrid <- (BasinAreaInNWPGrid / NWPGridFullArea)
FracBasinAreaInNWPGrid[FracBasinAreaInNWPGrid <= Minimum_Basin_Fraction] <- NA
FracBasinAreaInNWPGridRes <- disaggregate(FracBasinAreaInNWPGrid, fact = ResFactor, fun = max, na.rm=TRUE)
facc_at_outlet2 <- mask(facc_at_outlet1, FracBasinAreaInNWPGridRes)
#
## Masking the maximum flow accumulation which generated from the high resolution elevation within the MESH modeling grid cell
## Major river: Uses the flow accumulation over the whole basin 
MajorRiv1 <- aggregate(facc_at_outlet2, fact = ResFactor, fun = max, na.rm=TRUE)
MajorRiv2 <- disaggregate(MajorRiv1, fact = ResFactor, fun = max, na.rm=TRUE)
MajorRiv3 <- (MajorRiv2 == facc_at_outlet2)
MajorRiv3[MajorRiv3 != 1] <- NA
#
# writeRaster(MajorRiv3, "MESH_MajorRiver.tif", datatype="INT2S", overwrite=TRUE)
#
## Masking the longest flow path length using the Major River outlet points from modeling grid cell
flow_path_MajorRiv1 <- mask(flow_path1, MajorRiv3)
flow_path_MajorRiv2 <- aggregate(flow_path_MajorRiv1, fact = ResFactor, fun = max, na.rm=TRUE)
flow_path_MajorRiv3 <- disaggregate(flow_path_MajorRiv2, fact = ResFactor, fun = max, na.rm=TRUE)
MajorRiv4 <- (flow_path_MajorRiv3 == flow_path_MajorRiv1)
MajorRiv4[MajorRiv4 != 1] <- NA
MajorRiv <- mask(facc2, MajorRiv4)
MajorRiv[MajorRiv < as.integer(ResFactor/2)] <- NA
##
MajorRiv[!is.na(MajorRiv)] <- 2
MajorRiv[is.na(MajorRiv)] <- 0
#
## Minor river: Uses the local (within the modeling grid cell) flow accumulation
MinorRiv1 <- aggregate(mask(facc1, facc_at_outlet2), fact = ResFactor, fun = max, na.rm=TRUE)
MinorRiv2 <- disaggregate(MinorRiv1, fact = ResFactor, fun = max, na.rm=TRUE)
MinorRiv3 <- mask(facc1, facc_at_outlet2)
MinorRiv4 <- (MinorRiv2 == MinorRiv3)
MinorRiv4[MinorRiv4 != 1] <- NA
#
## Masking the longest flow path length using the Minor River outlet points from modeling grid cell
flow_path_MinorRiv1 <- mask(flow_path1, MinorRiv4)
flow_path_MinorRiv2 <- aggregate(flow_path_MinorRiv1, fact = ResFactor, fun = max, na.rm=TRUE)
flow_path_MinorRiv3 <- disaggregate(flow_path_MinorRiv2, fact = ResFactor, fun = max, na.rm=TRUE)
MinorRiv <- (flow_path_MinorRiv3 == flow_path_MinorRiv1)
MinorRiv[is.na(MinorRiv)] <- 0
##
MajorMinorRiv1 <- (MajorRiv + MinorRiv)
MajorMinorRiv1[MajorMinorRiv1 == 0] <- NA
#
#########################################################################################################################
MajorMinorRiv2 <- aggregate(MajorMinorRiv1, fact = ResFactor, fun = max, na.rm=TRUE)
MajorMinorRiv3 <- disaggregate(MajorMinorRiv2, fact = ResFactor, fun = max, na.rm=TRUE)
MajorMinorRiv4 <- (MajorMinorRiv3 == MajorMinorRiv1)
MajorMinorRiv4[MajorMinorRiv4 != 1] <- NA
MajorMinorRiv <- mask(MajorMinorRiv3, MajorMinorRiv4)
##
### Selecting a single outlet per modeling grid cell in case there are multiple outlet of the same flow accumulation and longest flow path length
Riv_nwp_zone <- as.matrix(mask(nwp_zone, MajorMinorRiv))
xx <- which(as.matrix(MajorMinorRiv) > 0, arr.ind=TRUE)
yy <- nwp_zone[xx]
zz <- yy[duplicated(yy)]
#
if (length(zz) > 0) {
for (i in 1 : length(zz)) {
yz <- which(Riv_nwp_zone == zz[i], arr.ind=TRUE)
yz[1,] <- NA
MajorMinorRiv[yz] <- NA
} }
#
writeRaster(MajorMinorRiv, "MajorMinorRiv.tif", datatype="INT2S", overwrite=TRUE)
#
## Masking elevation by the longest flow path length of the maximum flow accumulation cell (either major or minor river) inside the modeling grid
### and add an extra elevation grid cells for final outlet points that have elevation values lower than the basins minimum 
elevn1 <- mask(filldem1, MajorMinorRiv)
#
streamgauge = SpatialPoints(readOGR(dsn = ".", layer = "outlet"))
outletcelnum <- extract(nwp_zone, streamgauge, cellnumbers=T)
xx <- rowColFromCell(nwp_zone, outletcelnum[,1])
yy <- xx
for (i in 1 : length(xx[,1])) {
for (j in 1 : (ResFactor*ResFactor)) {
yy[i,1] <- (yy[i,1] + as.integer(FdirNumber[fdir1[t(yy[i,])],1]))
yy[i,2] <- (yy[i,2] + as.integer(FdirNumber[fdir1[t(yy[i,])],2]))
if (nwp_zone[t(xx[i,])] != nwp_zone[t(yy[i,])]) {
    break
} } }
#
elevn1[yy] <- filldem1[yy]
elevn2 <- elevn1 
elevn2[yy] <- minValue(elevn1) - 1
#
## Keep the location of the basin outlets and their next of outlets for later use
basin_outlet <- xx
basin_outlet_next <- yy
#
## River bed elevation at the outlet of the local river within the modeling grid cell 
Elevn <- as.matrix(aggregate(elevn1, fact = ResFactor, fun = min, na.rm=TRUE))
#
## Rank the MESH modeling grid cell based on invert elevation produced from above step
elevn3 <- as.vector(as.matrix(aggregate(elevn2, fact = ResFactor, fun = min, na.rm=TRUE)))
rank1 <- rank(-elevn3, ties.method= "first", na.last = "keep")
#
## Produce the modeling grid cell Rank number based on the invert elevation - that gives order of computation
grid_rank  <- matrix(rank1, NumRow, NumCol)
#
rank2 <- raster(grid_rank)
extent(rank2) <- extent(nwp_grid)
dim(rank2) <- dim(nwp_grid)
res(rank2) <- res(nwp_grid)
crs(rank2) <- crs(nwp_grid)
rank3 = disaggregate(rank2, fact = ResFactor, fun = max, na.rm=TRUE)
writeRaster(rank2, "MESH_Rank.tif", datatype="INT2S", overwrite=TRUE)
#
###### Creating basin boundary mask file ###########################################################
basin_boundary_mask <- as.matrix(facc_at_outlet1)
basin_boundary_mask[basin_boundary_mask > 0] <- 1
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
#
########### Creating a mask file to mask only the edge of the nwp_zone grid cell
nwp_zone_edge_mask1 <- matrix(1, nrow(nwp_zone), ncol(nwp_zone))
for (i in 1:NumRow) {
  for (j in 1:NumCol) {
    nwp_zone_edge_mask1[(ResFactor*i-(ResFactor-2)):(ResFactor*i-1), (ResFactor*j-(ResFactor-2)):(ResFactor*j-1)] <- 0 
} }
#
nwp_zone_edge_mask <- ((nwp_zone_edge_mask1*basin_boundary_mask) + basin_boundary_mask1)
nwp_zone_edge_mask[nwp_zone_edge_mask == 0] <- NA
nwp_zone_edge_mask[nwp_zone_edge_mask > 1] <- 1
#
########### Masking the local (within the grid) flow direction by edge of the nwp mask file (nwp_zone_edge_mask)
fdir6 <- (as.matrix(fdir1)*nwp_zone_edge_mask)
########### Masking out the flow directions cells that are not leaving the modeling grid cell
xx <- which(fdir6 > 0, arr.ind=TRUE)
yy <- fdir1[xx]
zz <- xx  + FdirNumber[yy,]
for (i in 1 : length(xx[,1])) {
if (nwp_zone[t(xx[i,])] == nwp_zone[t(zz[i,])]) {
fdir6[t(xx[i,])] <- NA }
}
fdir6[basin_outlet] <- fdir1[basin_outlet]
#
## Local Rivers: Produce the Next of the modeling grid cell from the flow direction, Rank and length threshold above which the diagonal flow direction will be assigned
next1 <- fdir6
nwp_zone_rank <- as.matrix(nwp_zone)
nwp_zone_next <- matrix(nrow = nrow(nwp_zone), ncol = ncol(nwp_zone))
#
xx <- which(fdir6 > 0, arr.ind = TRUE)
yy <- xx
#
ListResFactor <- seq(1, as.integer(4*ResFactor), by = as.integer(ResFactor/2))
#
for (i in 1 : 5) {
###
for (j in ListResFactor[i] : ListResFactor[i + 1] - 1) {
zz <- fdir1[yy]
yy <- yy + FdirNumber[zz,]
}
#
next1[xx] <- rank3[yy]
nwp_zone_next[xx] <- nwp_zone_rank[yy]
zz <- which(nwp_zone_next == nwp_zone_rank, arr.ind = TRUE)
if (length(zz[,1]) < 1) {
  break }
zx <- cellFromRowCol(nwp_zone, xx[,1], xx[,2])
zy <- cellFromRowCol(nwp_zone, zz[,1], zz[,2])
yy <- yy[zx %in% zy,]
xx <- zz
###
}
#
## Major River: Produce the Next of the modeling grid cell from the flow direction, Rank and length threshold above which the diagonal flow direction will be assigned
#
xx <- which(as.matrix(MajorMinorRiv) > 0, arr.ind = TRUE)
yy <- xx

for (i in 1 : 2) {
  ###
  for (j in ListResFactor[i] : ListResFactor[i + 1] - 1) {
    zz <- fdir1[yy]
    yy <- yy + FdirNumber[zz,]
  }
next1[xx] <- rank3[yy]
nwp_zone_next[xx] <- nwp_zone_rank[yy]
zz <- which(as.matrix(MajorMinorRiv3[yy]) == 1, arr.ind = TRUE)
if (length(zz[,1]) < 1) {
  break }
zx <- cellFromRowCol(nwp_zone, xx[,1], xx[,2])
zy <- cellFromRowCol(nwp_zone, zz[,1], zz[,2])
yy <- yy[zx %in% zy,]
xx <- zz
###
}
### Correcting flow direction of the outlets based on the grid flow direction values.
next1[basin_outlet] <- rank3[basin_outlet_next]
nwp_zone_next[basin_outlet] <- nwp_zone_rank[basin_outlet_next]
##
##### Generating flow direction for the modeling grid cell ############################################
grid_fdir1 <- nwp_zone   #matrix(nrow = nrow(nwp_zone), ncol = ncol(nwp_zone))
grid_fdir1[!is.na(grid_fdir1)] <- NA
#
zz <- which(nwp_zone_next == (nwp_zone_rank + 1), arr.ind = TRUE)
grid_fdir1[zz] <- 1

zz <- which(nwp_zone_next == (nwp_zone_rank + 1 - NumCol), arr.ind = TRUE)
grid_fdir1[zz] <- 2

zz <- which(nwp_zone_next == (nwp_zone_rank - NumCol), arr.ind = TRUE)
grid_fdir1[zz] <- 3

zz <- which(nwp_zone_next == (nwp_zone_rank - NumCol - 1), arr.ind = TRUE)
grid_fdir1[zz] <- 4

zz <- which(nwp_zone_next == (nwp_zone_rank - 1), arr.ind = TRUE)
grid_fdir1[zz] <- 5

zz <- which(nwp_zone_next == (nwp_zone_rank + NumCol - 1), arr.ind = TRUE)
grid_fdir1[zz] <- 6

zz <- which(nwp_zone_next == (nwp_zone_rank + NumCol), arr.ind = TRUE)
grid_fdir1[zz] <- 7

zz <- which(nwp_zone_next == (nwp_zone_rank + NumCol + 1), arr.ind = TRUE)
grid_fdir1[zz] <- 8

zz <- which(nwp_zone_next == nwp_zone_rank, arr.ind = TRUE)
grid_fdir1[zz] <- fdir1[zz]
#
writeRaster(grid_fdir1, "MESH_FDir4.tif", datatype="INT2S", overwrite=TRUE)
##
### Produce the Next of the Major River modeling grid cell from the flow direction, Rank and length threshold below which the diagonal flow direction will be assigned
next2 <- MajorMinorRiv
grid_fdir2 <- next2
xx <- which(as.matrix(MajorMinorRiv) > 0, arr.ind=TRUE)
next2[xx] <- next1[xx]
grid_fdir2[xx] <- grid_fdir1[xx]
#
grid_next <- aggregate(next2, fact = ResFactor, fun = max, na.rm=TRUE)
writeRaster(grid_next, "MESH_Next.tif", datatype="INT2S", overwrite=TRUE)
grid_next <- as.matrix(grid_next)
#
grid_fdir <- aggregate(grid_fdir2, fact = ResFactor, fun = min, na.rm=TRUE)
writeRaster(grid_fdir, "grid_fdirp.tif", datatype="INT2S", overwrite=TRUE)
grid_fdir <- as.matrix(grid_fdir)
#
####### Mask and sum the flow accumulation of within the modeling grid cell into the eight direction using flow direction as a mask file 
### (Calculation of the local drainage area in all eight flow directions) 
nwp_grid_facc <- matrix(0, NumRow*NumCol, 8)
#
for (i in 1 : 8) {
fdir7 <- grid_fdir1
fdir7[fdir7 != i] <- NA
facc3 <- mask(facc1, fdir7)
nwp_grid_facc[,i] = as.vector(aggregate(facc3, fact = ResFactor, fun = sum, na.rm=TRUE))
}
#
nwp_grid_facc[is.na(nwp_grid_facc)] <- 0
#
GridFaction1 <- matrix(nwp_grid_facc[,1], NumRow, NumCol, byrow = T)
GridFaction2 <- matrix(nwp_grid_facc[,2], NumRow, NumCol, byrow = T)
GridFaction3 <- matrix(nwp_grid_facc[,3], NumRow, NumCol, byrow = T)
GridFaction4 <- matrix(nwp_grid_facc[,4], NumRow, NumCol, byrow = T)
GridFaction5 <- matrix(nwp_grid_facc[,5], NumRow, NumCol, byrow = T)
GridFaction6 <- matrix(nwp_grid_facc[,6], NumRow, NumCol, byrow = T)
GridFaction7 <- matrix(nwp_grid_facc[,7], NumRow, NumCol, byrow = T)
GridFaction8 <- matrix(nwp_grid_facc[,8], NumRow, NumCol, byrow = T)
#
## Calculate the fraction of drainage area in all modeling grid cells
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
} } }
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
} } }
#
######### compute the grid area [each cell FRAC within the basin] which is the product of drainage area fraction and the nominal grid area [m2]
NominalGrid_Area <- mean(NWPGridFullArea[!is.na(NWPGridFullArea)])
GridArea <- 1000000 * NominalGrid_Area * grid_drain_area / as.matrix(NWPGridFullArea)
# GridArea <- 1000000 * grid_drain_area * NominalGrid_Area / (ResFactor*ResFactor)
GridArea1 <- raster(GridArea/1000000)
extent(GridArea1) <- extent(nwp_grid)
dim(GridArea1) <- dim(nwp_grid)
res(GridArea1) <- res(nwp_grid)
crs(GridArea1) <- crs(nwp_grid)
#
writeRaster(GridArea1, "grid_area_weight.tif", datatype="FLT4S", overwrite=TRUE)
#
############ Drainage area in of the modeling grid [km2]
system("mpiexec -n 8 AreaD8 -p grid_fdirp.tif -ad8 grid_da_ad8.tif -wg grid_area_weight.tif -nc")
DA <- as.matrix(raster("grid_da_ad8.tif"))
DA[is.na(DA)] <- 0
#
########### Bankfull cross-section area of river channel in [m2]
Bankfull  <- 0.043 * DA + 1.1         # Default watflood equation
# Bankfull  <- 0.1667 * DA + 0.1     # best fit regression equation from Green Kenue
#
############## Define river classes (IAK) using drainage area (DA)
IAK <- matrix(0, NumRow, NumCol)
LogDA <- log(as.numeric(DA))
LogDA[LogDA == -Inf] <- NaN
LogDA1 <- as.vector(LogDA)
RR <- integer(IAKmax - 1)
for (i in 1:(IAKmax - 1)) {
  RR[i] = quantile(LogDA1, probs = i/IAKmax, na.rm = TRUE)
}
############# Apply the background field (filtered by NaN to target active grid cells in the domain).
IAK[!is.nan(LogDA)] <- IAKmax
########## Apply IAK values by comparing the transformed values in 'LogDA' to the 'RR' thresholds the (to (IAKmax - 1), IAKmax is applied as the background field above).
for (i in 1:(IAKmax - 1)) {
  IAK[LogDA > RR[i]] = IAKmax - i
}
#
############# The internal slope in each grid (Land slope in m/m) #####################
slope <- terrain(domain_dem1, opt = "slope", unit='degrees')
# slope <- terrain(domain_dem1, opt = "slope", unit='degrees', neighbors=8)
IntSlope1 <- mask(tan(slope*pi/180), facc_at_outlet1)
IntSlope <- as.matrix(aggregate(IntSlope1, fact = ResFactor, fun = mean, na.rm=TRUE))
#
############# Number of channels: Number of equally sized channels traversing an element (No. of channels draining through the cell)
########### It can be adjusted based on the stream order file
MaxSord <- maxValue(mask(stream_order1, facc_at_outlet1))
channelclass1 <- matrix(c(MaxSord-100, MaxSord-4, MaxSord-3, MaxSord-2, MaxSord-1, MaxSord-4, MaxSord-3, MaxSord-2, MaxSord-1, MaxSord, 1, 2, 3, 4, 5), nr = 5, nc = 3)
channelclass2 <- matrix(c(0, 1, 2, 3, 4, 1, 2, 3, 4, 5, 5, 4, 3, 2, 1), nr = 5, nc = 3)
#
stord1 <- reclassify(mask(stream_order1, facc_at_outlet1), channelclass1)
stord <- aggregate(stord1, fact = ResFactor, fun = max, na.rm=TRUE) 
#
Chnl <- as.matrix(reclassify(stord, channelclass2))
#
####### Channel length similar to Green Kenue methods to keep the WATFLOOD assumption of channel length from center to center of the modeling grid ####
channel_length <- matrix(0, nrow = NumRow, ncol = NumCol, byrow = T)
channel_slope <- matrix(0, nrow = NumRow, ncol = NumCol, byrow = T)
# #
for (i in 1:nrow(channel_length)) {
  for (j in 1:ncol(channel_length)) {
    
    if ((i-1 >= 1) & (j-1 >= 1) & (i+1 <= NumRow) & (j+1 <= NumCol)) {
      if (is.na(grid_fdir[i,j])) {
        channel_length[i,j] <- 0 }
      else if (grid_fdir[i,j] == 1) {
        channel_length[i,j] <- distm (c(XLon[i,j], YLat[i,j]), c(XLon[i,j+1], YLat[i,j+1]), fun = distHaversine)
        channel_slope[i,j] <- ((Elevn[i,j] - Elevn[i,j+1])/channel_length[i,j])}
      else if (grid_fdir[i,j] == 2) {
        channel_length[i,j] <- distm (c(XLon[i,j], YLat[i,j]), c(XLon[i-1,j+1], YLat[i-1,j+1]), fun = distHaversine)
        channel_slope[i,j] <- ((Elevn[i,j] - Elevn[i-1,j+1])/channel_length[i,j])}
      else if (grid_fdir[i,j] == 3) {
        channel_length[i,j] <- distm (c(XLon[i,j], YLat[i,j]), c(XLon[i-1,j], YLat[i-1,j]), fun = distHaversine)
        channel_slope[i,j] <- ((Elevn[i,j] - Elevn[i-1,j])/channel_length[i,j])}
      else if (grid_fdir[i,j] == 4) {
        channel_length[i,j] <- distm (c(XLon[i,j], YLat[i,j]), c(XLon[i-1,j-1], YLat[i-1,j-1]), fun = distHaversine)
        channel_slope[i,j] <- ((Elevn[i,j] - Elevn[i-1,j-1])/channel_length[i,j])}
      else if (grid_fdir[i,j] == 5) {
        channel_length[i,j] <- distm (c(XLon[i,j], YLat[i,j]), c(XLon[i,j-1], YLat[i,j-1]), fun = distHaversine)
        channel_slope[i,j] <- ((Elevn[i,j] - Elevn[i,j-1])/channel_length[i,j])}
      else if (grid_fdir[i,j] == 6) {
        channel_length[i,j] <- distm (c(XLon[i,j], YLat[i,j]), c(XLon[i+1,j-1], YLat[i+1,j-1]), fun = distHaversine)
        channel_slope[i,j] <- ((Elevn[i,j] - Elevn[i+1,j-1])/channel_length[i,j])}
      else if (grid_fdir[i,j] == 7) {
        channel_length[i,j] <- distm (c(XLon[i,j], YLat[i,j]), c(XLon[i+1,j], YLat[i+1,j]), fun = distHaversine)
        channel_slope[i,j] <- ((Elevn[i,j] - Elevn[i+1,j])/channel_length[i,j])}
      else if (grid_fdir[i,j] == 8) {
        channel_length[i,j] <- distm (c(XLon[i,j], YLat[i,j]), c(XLon[i+1,j+1], YLat[i+1,j+1]), fun = distHaversine)
        channel_slope[i,j] <- ((Elevn[i,j] - Elevn[i+1,j+1])/channel_length[i,j])}
    }
  }
}
channel_slope[is.na(channel_slope)] <- 0
channel_slope[channel_slope < Min_Chanel_Slope] <- Min_Chanel_Slope
#
## Reach number for lake, reservoir and/or external routing: user can assign number from 1,2,3...with respect to the number of reservoir / lakes 
nwp_reach  <- (nwp_grid - nwp_grid)
if (InsertedReservoir == "Reservoir") {
  rescelnum <- extract(nwp_reach, reservoirlocation, cellnumbers=T)
  for (i in 1:length(rescelnum)) { 
    nwp_reach[rescelnum[i]] <- i
  }
}
#
Reach <- as.matrix(nwp_reach)
#
################## Collecting the drainage database parameters from the above steps
grid_rank[is.na(grid_rank)] <- 0
#
grid_next[is.na(grid_next)] <- 0
#
DA[is.na(DA)] <- 0
#
channel_length[is.na(channel_length)]  <- 0
#
Bankfull[is.na(Bankfull)] <- 0
#
channel_slope[is.na(channel_slope)] <- 0
#
Elevn[is.na(Elevn)] <- 0
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
                      (apply(channel_slope, 2, rev)), (apply(Elevn, 2, rev)), (apply(channel_length, 2, rev)), (apply(IAK, 2, rev)),
                      (apply(IntSlope, 2, rev)), (apply(Chnl, 2, rev)), (apply(Reach, 2, rev)), (apply(GridArea, 2, rev)))
#
rm(list=setdiff(ls(), c("nwp_zone", "drainagedata", "facc_at_outlet1", "land_cover", "domain_dem1", "drain_net", "CLAYSANDSOM", "GlacierInLandCover",
                        "MESHVersion", "NominalGrid_Area", "SpatialFilter", "Minimum_GRU_Fraction", "Minimum_Glacier_Fraction", "UpdateIrrigation", "domain_irrRes",
                        "NumRow", "NumCol", "ResFactor" ,"GRUsClass", "AspectClass", "SlopeClass", "Min_Chanel_Slope","BasinName", "NominalGridSize",
                        "LLXcorner", "LLYcorner", "XRes", "YRes", "grid_rank", "grid_next", "IAKmax", "NumbSoilData", "delta", "delta_elev_max", "elev_band")))
###################################################################################################################################################################
###################################################################################################################################################################
#### Masking the land cover by the river basin ######################
basin_landcover <- mask(land_cover, facc_at_outlet1)
##### Calculate Slope 
slope <- terrain(domain_dem1, opt = "slope", unit='degrees')
# slope <- terrain(domain_dem1, opt = "slope", unit='degrees', neighbors=8)
basin_slope <- mask(slope, facc_at_outlet1)
#
if (MESHVersion == "Original") {
  grus1 <- basin_landcover }
#
if (MESHVersion == "Mountain") {
  #
  ## Calculate aspect
  aspect <- focal(terrain(domain_dem1, opt = "aspect", unit='degrees'), SpatialFilter)
  # aspect <- focal(terrain(domain_dem1, opt = "aspect", unit='degrees', neighbors=8, flatAspect = NA), SpatialFilter)
  ## Mask the aspect by basin
  basin_aspect <- mask(aspect, facc_at_outlet1)
  #
  #  writeRaster(slope, "R_Slope.tif", datatype="FLT4S", overwrite=TRUE)
  #  writeRaster(aspect, "R_Aspect.tif", datatype="FLT4S", overwrite=TRUE)
  #
  ###### GRUs Creation and model discretization by Combining land cover, slope and aspect################################################## 
  ##### for mountain MESH version and land cover only for Original MESH version
  ######################### Reclassifying slope and aspect ################################################################################
  
  basin_landcover <- land_cover
  basin_slope <- slope
  basin_aspect <- aspect
  
  basin_slopeclass <- reclassify(basin_slope, SlopeClass)
  basin_aspectclass <- reclassify(basin_aspect, AspectClass)
  # writeRaster(basin_aspectclass, "NR_basin_aspectclass.tif", datatype="FLT4S", overwrite=TRUE)
  #
  ### Model discretization and GRUs creation: Combine the land cover, slope and aspect ####################################################
  grus1 <- basin_landcover
  #
  for (i in 1 : maxValue(basin_landcover)) {
    ### GRUs discretization when there are two slope aspect classes #########
    grus1[basin_landcover == i & basin_slopeclass == 2 & basin_aspectclass == 2] <- 3*i-2
    grus1[basin_landcover == i & basin_slopeclass == 2 & basin_aspectclass == 1] <- 3*i-1
    grus1[basin_landcover == i & basin_slopeclass != 2] <- 3*i
    #
    # ### GRUs discretization when there are four slope aspect classes #########
    # grus1[basin_landcover == i & basin_slopeclass == 2 & basin_aspectclass == 2] <- 4*i-3 #1
    # grus1[basin_landcover == i & basin_slopeclass == 2 & basin_aspectclass == 1] <- 4*i-2 #2
    # grus1[basin_landcover == i & basin_slopeclass == 2 & basin_aspectclass == 3] <- 4*i-1 #3
    # grus1[basin_landcover == i & basin_slopeclass == 2 & basin_aspectclass == 4] <- 4*i   #4
    # grus1[basin_landcover == i & basin_slopeclass != 2] <- 4*i                            #4
    #
    # #### GRUs discretization when there are two slope aspect classes and two elevation bands #########
    # grus1[basin_landcover == i & elev_band == 2 & basin_slopeclass == 2 & basin_aspectclass != 2] <- 5*i-4
    # grus1[basin_landcover == i & elev_band == 2 & basin_slopeclass == 2 & basin_aspectclass != 1] <- 5*i-3
    # grus1[basin_landcover == i & elev_band == 2 & basin_slopeclass != 2] <- 5*i
    # grus1[basin_landcover == i & elev_band == 1 & basin_slopeclass == 2 & basin_aspectclass != 2] <- 5*i-2
    # grus1[basin_landcover == i & elev_band == 1 & basin_slopeclass == 2 & basin_aspectclass != 1] <- 5*i-1
    # grus1[basin_landcover == i & elev_band == 1 & basin_slopeclass != 2] <- 5*i
    #
  }
  #
  grus <- reclassify(grus1, GRUsClass)
  grus1 <- grus
  #
}
# writeRaster(grus1, "M_GRUs_Produced.tif", datatype="INT2S", overwrite=TRUE)
#
##### Polishing GRUs that have less or equal to the minimum land cover fraction in the modeling grid ################################
grus_mask <- grus1
grus_mask[grus_mask != 1] <- 1
TotalGrid <- as.vector(aggregate(grus_mask, fact = ResFactor, fun = sum, na.rm=TRUE))
#
GRUGrid_Frac <- matrix(nrow = NumRow*NumCol, ncol = maxValue(grus1))
#
for (i in 1 : maxValue(grus1)) {
  grus_mask <- grus1
  grus_mask[grus_mask != i] <- NA
  grus_mask[grus_mask == i] <- 1
  GRU_Frac <- as.vector(aggregate(grus_mask, fact = ResFactor, fun = sum, na.rm=TRUE))
  GRUGrid_Frac[,i] <- (GRU_Frac/TotalGrid)
}
GRUGrid_Frac[is.na(GRUGrid_Frac)] <- 0
####################################################################################################
PolishedGRUs <- raster(matrix(0, nrow = nrow(grus1), ncol = ncol(grus1), byrow = T))
extent(PolishedGRUs) <- extent(grus1)
res(PolishedGRUs) <- res(grus1)
crs(PolishedGRUs) <- crs(grus1)
#
for (i in 1:maxValue(grus1)) {
  grusmask <- grus1
  grusmask[grusmask != i] <- NA
  xxGRU_Frac <- GRUGrid_Frac[,i]
  if (is.element(i, GlacierInLandCover)) {
    xxGRU_Frac[xxGRU_Frac < Minimum_Glacier_Fraction] <- 0
    GLGRU_Frac <- xxGRU_Frac
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
# writeRaster(PolishedGRUs, "PolishedGRUs.tif", datatype="INT2S", overwrite=TRUE)
##### Calculating fraction of GRUs in the modeling grid ############################################
grus_mask <- PolishedGRUs
grus_mask[grus_mask != 1] <- 1
grus2 <- as.vector(aggregate(grus_mask, fact = ResFactor, fun = sum, na.rm=TRUE))
#
grus3 <- matrix(nrow = NumRow*NumCol, ncol = (1 + maxValue(grus1)))
#
for (i in 1 : maxValue(grus1)) {
  grus_mask <- PolishedGRUs
  grus_mask[grus_mask != i] <- NA
  grus_mask[grus_mask == i] <- 1
  grus3[,i] <- as.vector(aggregate(grus_mask, fact = ResFactor, fun = sum, na.rm=TRUE)) / grus2
}
#
grus3[is.na(grus3)] <- 0
#
#######################################################################################################
# for (i in 1 : maxValue(grus1)) {
#   if (is.element(i, GlacierInLandCover)) {
#     grus3[,i] <- 0
#   }
# }
# #
# grus4 <- (grus3 * (1 - GLGRU_Frac) / rowSums(grus3, na.rm = TRUE, dims = 1))
# grus3 <- grus4
# #
# for (i in 1 : maxValue(grus1)) {
#   if (is.element(i, GlacierInLandCover)) {
#     grus3[,i] <- GLGRU_Frac
#   }
# }
# #
# grus3[is.na(grus3)] <- 0
##########################################################################################################
#### Updating the irrigation fraction of the modeling grid ###############################################
if (UpdateIrrigation == "Irrigation") {
  Irrig_mask <-  mask(domain_irrRes, facc_at_outlet1)
  IrrigFracSum <- as.vector(aggregate(Irrig_mask, fact = ResFactor, fun = sum, na.rm=TRUE))
  Irrig_mask[Irrig_mask != 1] <- 1
  IrrigTotal <- as.vector(aggregate(Irrig_mask, fact = ResFactor, fun = sum, na.rm=TRUE))
  IrrigFract <- (IrrigFracSum/IrrigTotal)
  IrrigFract[is.na(IrrigFract)] <- 0
  #  
  grus_frac_r2c <- matrix(nrow = NumRow*(1 + maxValue(grus1)), ncol = NumCol)
  #
  for (i in 1 : (1 + maxValue(grus1))) {
    
    if (is.element(i, CropInLandCover)) {
      cropgrus <- ((1 - IrrigFract)*(grus3[,i] + grus3[,i + 1]))
      grus_frac_r2c[(NumRow*i-(NumRow-1)):(NumRow*i),] <- apply((matrix(cropgrus, NumRow, NumCol, byrow = T)), 2, rev)
    } 
    else if (is.element(i, IrrigInLandCover)) {
      Irriggrus <- (IrrigFract*(grus3[,i] + grus3[,i - 1]))
      grus_frac_r2c[(NumRow*i-(NumRow-1)):(NumRow*i),] <- apply((matrix(Irriggrus, NumRow, NumCol, byrow = T)), 2, rev)
      
    } else {    
      grus_frac_r2c[(NumRow*i-(NumRow-1)):(NumRow*i),] <- apply((matrix(grus3[,i], NumRow, NumCol, byrow = T)), 2, rev)
    }
  }
} else {
  ############################################################################################################################
  grus_frac_r2c <- matrix(nrow = NumRow*(1 + maxValue(grus1)), ncol = NumCol)
  #
  for (i in 1 : (1 + maxValue(grus1))) {
    grus_frac_r2c[(NumRow*i-(NumRow-1)):(NumRow*i),] <- apply((matrix(grus3[,i], NumRow, NumCol, byrow = T)), 2, rev)
  }
}
grus_frac_r2c[is.na(grus_frac_r2c)] <- 0
##
### Combine the drainage data and fraction of GRUs for MESH_drainage_database file
dranage_database <- rbind(drainagedata, grus_frac_r2c)
#
######### common data for MESH_drainage_database and MESH_parameters files ########################
xorigin_value <- paste(":xOrigin                 ", LLXcorner, sep = "")
yorigin_value <- paste(":yOrigin                 ", LLYcorner, sep = "")
xcount_value <- paste(":xCount                   ", NumCol, sep = "")
ycount_value <- paste(":yCount                   ", NumRow, sep = "")
xdelta_value <- paste(":xDelta                   ", XRes, sep = "")
ydelta_value <- paste(":yDelta                   ", YRes, sep = "")
########################################################################################################
## Creating MESH drainage database file and save it as MESH_drainage_database.r2c
########################################################################################################
## Add source file name to header and required information for the drainage database file header
NominalGrid_Area_Value <- paste(":NominalGridSize_AL   ",round(NominalGridSize, 2), sep = "")
contour_interval_value <- paste(":ContourInterval   ",1, sep = "")
Impervious_Area_value <- paste(":ImperviousArea   ",0, sep = "")
class_count_value <- paste(":ClassCount   ",maxValue(grus1) + 1, sep = "")
NumRiver_Classes_value <- paste(":NumRiverClasses   ",IAKmax, sep = "")
Elev_Conversion_value <- paste(":ElevConversion   ",1, sep = "")
Total_NumOfGrids_value <- paste(":TotalNumOfGrids   ",(length(which(grid_rank != 0))), sep = "")  
NumGridsIn_Basin_value <- paste(":NumGridsInBasin   ",(length(which(grid_next != 0))), sep = "") 
Minimum_Slope_value <- paste(":MinimumSlope   ",Min_Chanel_Slope, sep = "")
ProjectionName <- "LATLONG"
DatumName <- "WGS84"
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
## add source file name to header
sourcefile <- ":SourceFile              DEM from Hydroshed and Land Cover from CEC"
header1 <- c(header1, sourcefile)
header1 <- c(header1, NominalGrid_Area_Value)
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
for (i in 1:(1 + maxValue(grus1))) {
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
###################################################################################################################################
#######################################################################################################
## Creating the Mountain MESH Parameter file and save it as MESH_Parameter.r2c
######################################################################################################
######## Soil data for MESH_parameters file
ClaySandSomMasked <- as.array(aggregate(mask(CLAYSANDSOM, facc_at_outlet1), fact = ResFactor, fun = mean, na.rm=TRUE))
ClaySandSomData <- apply(ClaySandSomMasked, c(2,3), rev)
soildatabase <- matrix(nrow = NumRow*dim(ClaySandSomData)[3], ncol = NumCol)
for (i in 1 : dim(ClaySandSomData)[3]) {
  soildatabase[(NumRow*i-(NumRow-1)):(NumRow*i),] <- ClaySandSomData[,,i]
}
soildatabase[is.na(soildatabase)] <- 0
Number_Soil_Layers <- as.integer(dim(ClaySandSomData)[3]/NumbSoilData)
#
############### compute drainage density in the modeling grid ###
drain_net[drain_net == 0] <- NA
basin_DN <- mask(drain_net, facc_at_outlet1)
basin_nwp_zone <- mask(nwp_zone, facc_at_outlet1)
#
dn_cell_length <- sqrt(2) * sqrt(area(basin_DN, na.rm=TRUE, weights=FALSE))
grid_area_km2 <- area(basin_nwp_zone, na.rm=TRUE, weights=FALSE)
zonal_drnr_km <- aggregate(dn_cell_length, fact = ResFactor, fun = sum, na.rm=TRUE)
zonal_drnr_km[is.na(zonal_drnr_km)] <- 0.1  # when the grid has no stream use the minimum length of 0.1 km (100 meter) for overland flow.
zonal_area_km2 <- aggregate(grid_area_km2, fact = ResFactor, fun = sum, na.rm=TRUE)
drainge_dens1 <- zonal_drnr_km / zonal_area_km2
drainge_dens1[is.na(drainge_dens1)] <- 0
drainge_dens <- apply(as.matrix(drainge_dens1), 2, rev)
############################################################################################################
#
if (MESHVersion == "Mountain") {
  #### Calculate Mountain MESH model parameters
  ##### Weight for grid curvature that will be used to calculate wind speed correction parameters
  weight1 = matrix(c(0,1,0,1,0,1,0,1,0), nrow=3)
  weight2 = matrix(c(1,0,1,0,0,0,1,0,1), nrow=3)
  curve1 <- 0.5*focal(domain_dem1, w=weight1,fun=sum, na.rm=FALSE)
  curve2 <- 0.5*focal(domain_dem1, w=weight1,fun=sum, na.rm=FALSE) 
  dem_cell_area <- area(domain_dem1, na.rm=TRUE, weights=FALSE)
  dem_cell_area <- dem_cell_area[!is.na(dem_cell_area)]
  dem_cell_length <- 1000*sqrt(median(dem_cell_area))
  dem_cell_length <- min(dem_cell_length[!is.na(dem_cell_length)],300)
  curve <- 0.25 * (((2.0 * sqrt(2.0) + 2) * domain_dem1) - (sqrt(2.0)*curve1) - (curve2)) / (2.0 * sqrt(2.0) * dem_cell_length)
  curve_Max <- 2.0 * disaggregate(aggregate(abs(curve), fact = ResFactor, fun = max, na.rm=TRUE), fact = ResFactor, fun = max, na.rm=TRUE)
  #
  basin_elevn <- mask(domain_dem1, PolishedGRUs)
  basin_delta <- mask(delta, PolishedGRUs)
  basin_delta_elev_max <- mask(delta_elev_max, PolishedGRUs)
  basin_slope0 <- mask(slope, PolishedGRUs)
  basin_aspect0 <- mask(aspect, PolishedGRUs)
  #  basin_curve <- mask(curve, PolishedGRUs)
  basin_curve <- mask((curve/curve_Max), PolishedGRUs)
  #
  ###### Calculate the Sine and cosine of aspect expressed in radian for the calculation of the zonal GRUs mean aspect
  basin_sinaspect <- sin(basin_aspect0*pi/180)
  basin_cosaspect <- cos(basin_aspect0*pi/180)
  #
  ###### Calculate the weighted average values of elevation, slope, aspect, delta and curvature for MESH_Parameters.r2c file
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
    basin_elevn2 <- mask(basin_elevn, grus_mask)
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
    basin_sinaspect3 <- as.matrix(aggregate(basin_sinaspect2, fact = ResFactor, fun = sum, na.rm=TRUE))
    basin_cosaspect3 <- as.matrix(aggregate(basin_cosaspect2, fact = ResFactor, fun = sum, na.rm=TRUE))
    basin_aspect2 <- ((360 + (atan2(basin_sinaspect3,basin_cosaspect3))*(180/pi)) %% 360)
    basin_curve3 <- as.matrix(aggregate(basin_curve2, fact = ResFactor, fun = mean, na.rm=TRUE))
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
  ######### Combine soil, elevation, slope, aspect, delta, delta_elev_max, curve and drainage density for MESH_parameters file
  soilelevnslopeaspectdeltacurvedraindens <- rbind(soildatabase,elevnslopeaspectdeltacurve,drainge_dens)
} else {
  ######## Combine Soil, slope and drainage density for MESH_parameters file
  gridslope0 <- as.matrix(aggregate(tan(basin_slope*pi/180), fact = ResFactor, fun = mean, na.rm=TRUE))
  gridslope <- apply(gridslope0, 2, rev)
  gridslope[is.na(gridslope)] <- 0
  soilelevnslopeaspectdeltacurvedraindens <- rbind(soildatabase,gridslope,drainge_dens)
  #
}
############# Write header for the MESH parameter r2c file ###
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
  for (i in 1 : Number_Soil_Layers) {
    j <- i 
    attributename <- paste(":AttributeName ", j, " Clay	", i, "")
    header1 <- c(header1, attributename)
    attributeunits <- paste(":AttributeUnits ", j, " % ")
    header1 <- c(header1, attributeunits)
    
  }
  #
  for (i in 1 : Number_Soil_Layers) {
    j <- j + 1
    attributename <- paste(":AttributeName ", j, " Sand	", i, "")
    header1 <- c(header1, attributename)
    attributeunits <- paste(":AttributeUnits ", j, " % ")
    header1 <- c(header1, attributeunits)
  }
  #
  if (NumbSoilData > 2) {
    for (i in 1 : Number_Soil_Layers) {
      j <- j + 1
      attributename <- paste(":AttributeName ", j, " orgm	", i, "")
      header1 <- c(header1, attributename)
      attributeunits <- paste(":AttributeUnits ", j, " % ")
      header1 <- c(header1, attributeunits)
    }
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
  #  
  for (row in 1:nrow(soilelevnslopeaspectdeltacurvedraindens)) {
    meshparams <- formatC(soilelevnslopeaspectdeltacurvedraindens[row, ],
                          digits = 6, width = 1,
                          format = "f")
    cat(meshparams, "\n", sep = "  ", file = paste0(BasinName,MESHVersion,"_MESH_parameters.r2c"), append = TRUE)
  }
} else {
  #
  for (i in 1 : Number_Soil_Layers) {
    j <- i
    attributename <- paste(":AttributeName ", j, " Clay	", i, "")
    header1 <- c(header1, attributename)
    attributeunits <- paste(":AttributeUnits ", j, " % ")
    header1 <- c(header1, attributeunits)
  }
  #
  for (i in 1 : Number_Soil_Layers) {
    j <- j + 1
    attributename <- paste(":AttributeName ", j, " Sand	", i, "")
    header1 <- c(header1, attributename)
    attributeunits <- paste(":AttributeUnits ", j, " % ")
    header1 <- c(header1, attributeunits)
  }
  #
  if (NumbSoilData > 2) {
    for (i in 1 : Number_Soil_Layers) {
      j <- j + 1
      attributename <- paste(":AttributeName ", j, " orgm	", i, "")
      header1 <- c(header1, attributename)
      attributeunits <- paste(":AttributeUnits ", j, " % ")
      header1 <- c(header1, attributeunits)
    }
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
  
  for (row in 1:nrow(soilelevnslopeaspectdeltacurvedraindens)) {
    meshparams <- formatC(soilelevnslopeaspectdeltacurvedraindens[row, ],
                          digits = 6, width = 1,
                          format = "f")
    cat(meshparams, "\n", sep = "  ", file = paste0(BasinName,MESHVersion,"_MESH_parameters.r2c"), append = TRUE)
  }
}
###########################################################################################################################
###########################################################################################################################
