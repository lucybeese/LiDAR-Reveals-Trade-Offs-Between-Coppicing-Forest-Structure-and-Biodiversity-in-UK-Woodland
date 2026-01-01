##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### Bontuchel Biodiversity Analysis ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

##%%%%%%%%%%%%%%%%%%%%%##

#### Data Processing ####

##%%%%%%%%%%%%%%%%%%%%%##

# Load libraries
library(terra)
library(MultiscaleDTM)
library(viridisLite)
library(spatstat.geom)
library(sf)
library(dplyr)
library(tidyr)
library(stringr)
library(lidR)
library(sp)
library(ggplot2)
library(weights)
library(corrplot)
library(randomForest) 
library(caret)
library(car)
library(leafR)
library(RColorBrewer)
library(RCSF)
library(plotrix)
library(car)      
library(MuMIn)    
library(leaps)    
library(mgcv)
library(MASS)
library(vegan)
library(lme4)

##%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### STEP 1: PROCESS DATA ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%##

# Set working directory
setwd("~/Library/CloudStorage/OneDrive-UniversityofBristol/Documents/PhD/Project 3_Bontuchel/Data")

#--------------------------------------------------#
# Load in shapefile of site and create a dataframe #
#--------------------------------------------------#

bontuchel<-vect('Bontuchel_Coppice_Panels.shp')
bon_df<- as.data.frame(bontuchel)

# Determine total area of coppicing panels
bon_df$Area<-as.numeric(bon_df$Area)
sum(bon_df$Area,na.rm = TRUE)
summary(bon_df$Area)
ext(bontuchel)

# Remove unnecessary columns from the data set and clean columns
bon_df$Id<-NULL
bon_df$DISS<-NULL
names(bon_df)[names(bon_df) == "box_yr"] <- "Coppicing_year"
names(bon_df)[names(bon_df) == "ID_No"] <- "Coppicing_Plot"

#-----------------------------#
# Load in canopy height model #
#-----------------------------#

bon_CHM_2020<-rast('Bont_CHM.tif')

# Replace negative values with zero
bon_CHM_2020
bon_CHM_2020[bon_CHM_2020 < 0] <- 0

#-------------#
# Load in DTM #
#-------------#

bon_DTM_2020<- rast('Bont_DTM_Merged.tif')

# Replace error values with NA
min(bon_DTM_2020)
bon_DTM_2020[bon_DTM_2020 == -999999] <- NA

#------------------------#
# Load in Dormouse Boxes #
#------------------------#

boxes<-vect('Boxcoordinates.shp')
boxes_df<- as.data.frame(boxes)

# Remove unnecessary columns from the data set and clean columns
boxes_df$Object_ID<- NULL
names(boxes_df)[names(boxes_df) == "Descriptio"] <- "Dor_Box"

#----------------------------#
# Project to the correct CRS #
#----------------------------#

# epsg for British National Grid coordinate system
terra::crs(bon_CHM_2020)
set.crs(bon_CHM_2020, "epsg:27700") 

terra::crs(bon_DTM_2020)
set.crs(bon_DTM_2020, "epsg:27700")

terra::crs(bontuchel)
terra::crs(boxes)

#------#
# Crop #
#------#

bon_CHM_2020 <-crop(bon_CHM_2020, ext(bontuchel)+20)
bon_DTM_2020 <-crop(bon_DTM_2020, ext(bontuchel)+20)

# Save
writeRaster(bon_DTM_2020, "bon_DTM_2020_rep.tif", overwrite=T)
writeRaster(bon_CHM_2020, "bon_CHM_2020_rep.tif", overwrite=T)
write.csv(bon_df,"bon_df.csv",row.names = F)
write.csv(boxes_df,"boxes_df.csv",row.names = F)

#-----------------------#
##### Load all data #####
#-----------------------#

bon_DTM_2020<- rast("bon_DTM_2020_rep.tif")
bon_CHM_2020<- rast("bon_CHM_2020_rep.tif")
bon_df<-read.csv('bon_df.csv')
boxes_df<-read.csv('boxes_df.csv')
bontuchel<-vect('Bontuchel_Coppice_Panels.shp')
boxes<-vect('Boxcoordinates.shp')

#------#
# Plot #
#------#

# Open a PNG device
pdf("~/Library/CloudStorage/OneDrive-UniversityofBristol/GitHub/Thesis/Chapter4/Figures/CHM_2020.pdf")

# Plot
plot(bon_CHM_2020)
plot(bontuchel, add=TRUE, border='#810847', lwd= 2)

# Close the device
dev.off()

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### STEP 2: LOAD FUNCTIONS ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

# Create and load updated ForestGapR Functions to work with terra

#-------------------#
### GetForestGaps ###
#-------------------#

getForestGaps <- function(chm_layer, threshold = 10, size = c(1, 10^4)) {
  chm_layer[chm_layer > threshold] <- NA
  chm_layer[chm_layer <= threshold] <- 1
  gaps <- terra::patches(chm_layer, directions = 8, allowGaps = FALSE)
  rcl <- terra::freq(gaps)
  rcl$layer<-NULL
  rcl[, 2] <- rcl[, 2] * terra::res(chm_layer)[1]^2
  z <- terra::classify(gaps, rcl = rcl, right = FALSE)
  z[is.na(gaps)] <- NA
  gaps[z > size[2]] <- NA
  gaps[z < size[1]] <- NA
  gaps <- terra::patches(gaps, directions = 8, allowGaps = FALSE)
  names(gaps) <- "gaps"
  return(gaps)
}

#-------------#
### GapSPDF ###
#-------------#

GapSPDF <- function(gap_layer){
  gaps_poly <- terra::as.polygons(gap_layer, dissolve=TRUE, na.rm=TRUE, values=TRUE)
  names(gaps_poly) <- "gap_id"
  gaps_sf <- sf::st_as_sf(gaps_poly)
  gaps_df <- terra::as.data.frame(terra::centroids(gaps_poly), geom="XY")
  sf_polys<- sf::as_Spatial(gaps_sf)
  gaps_spdf<- sp::SpatialPolygonsDataFrame(sf_polys, gaps_df)
  return(gaps_spdf)
}

#--------------#
### GapStats ###
#--------------#

GapStats <- function(gap_layer, chm_layer) {
  GiniCoeff <- function(x, finite.sample = TRUE, na.rm = TRUE) {
    if (!na.rm && any(is.na(x))) {
      return(NA_real_)
    }
    x <- as.numeric(stats::na.omit(x))
    n <- length(x)
    x <- sort(x)
    G <- 2 * sum(x * 1L:n) / sum(x) - (n + 1L)
    if (finite.sample) {
      GC <- G / (n - 1L)
    } else {
      GC <- G / n
    }
    return(GC)
  }
  
  Range_Func <- function(x) {
    max_val <- max(x, na.rm = TRUE)
    min_val <- min(x, na.rm = TRUE)
    range_val <- max_val - min_val
    return(range_val)
  }
  
  gap_list <- data.frame(terra::freq(gap_layer))
  gap_list$layer<-NULL
  gap_list$count <-gap_list$count * terra::res(gap_layer)[1]^2
  gap_list <- gap_list[!is.na(gap_list[, 1]), ]
  
  chm_max <- stats::aggregate(chm_layer[], by = list(gap_layer[]), FUN = max)
  chm_min <- stats::aggregate(chm_layer[], by = list(gap_layer[]), FUN = min)
  chm_mean <- round(stats::aggregate(chm_layer[], by = list(gap_layer[]), FUN = mean), 2)
  chm_sd <- round(stats::aggregate(chm_layer[], by = list(gap_layer[]), FUN = sd), 2)
  chm_gini <- round(stats::aggregate(chm_layer[], by = list(gap_layer[]), GiniCoeff), 2)
  chm_range <- round(stats::aggregate(chm_layer[], by = list(gap_layer[]), Range_Func), 2)
  
  # Rename columns
  colnames(gap_list) <- c("gap_id", "gap_area")
  colnames(chm_max) <- c("gap_id", "chm_max")
  colnames(chm_min) <- c("gap_id", "chm_min")
  colnames(chm_mean) <- c("gap_id", "chm_mean")
  colnames(chm_sd) <- c("gap_id", "chm_sd")
  colnames(chm_gini) <- c("gap_id", "chm_gini")
  colnames(chm_range)<- c("gap_id", "chm_range")
  
  # Merge all results into one data frame
  gap_list <- merge(gap_list, chm_max, by = "gap_id")
  gap_list <- merge(gap_list, chm_min, by = "gap_id")
  gap_list <- merge(gap_list, chm_mean, by = "gap_id")
  gap_list <- merge(gap_list, chm_sd, by = "gap_id")
  gap_list <- merge(gap_list, chm_gini, by = "gap_id")
  gap_list <- merge(gap_list, chm_range, by = "gap_id")
  
  return(gap_list)
}

#-------------------#
### Canopy Volume ###
#-------------------#

canopy_volume <- function(chm) {
  volume <- sum(values(chm), na.rm = TRUE) * res(chm)[1] * res(chm)[2]
  return(volume)
}

#------------------------#
### Proportion of gaps ###
#------------------------#

calculate_gap_proportion <- function(gap_layer, plot_area) {
  gap_area <- global(gap_layer, "sum", na.rm = TRUE)[1, 1] * res(gap_layer)[1] * res(gap_layer)[2]
  proportion <- gap_area / plot_area
  return(proportion)
}
# m3

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### STEP 3: CREATE TOPOGRAPHIC RASTERS ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#-------------------------------#
### Calculate terrain metrics ###
#-------------------------------#

# Calculate terrain metrics
bon_Slope<- terrain(bon_DTM_2020, v="slope", unit='degrees')
bon_TPI <- TPI(bon_DTM_2020, w = 10, shape = "circle", unit= "map") # Scale 10m radii
bon_TRI <- terrain(bon_DTM_2020, v="TRI")
bon_Aspect<- terrain(bon_DTM_2020, v="aspect", unit='degrees')
bon_Slope_rad<- terrain(bon_DTM_2020, v="slope", unit='radians')
bon_Aspect_rad<- terrain(bon_DTM_2020, v="aspect", unit='radians')
bon_DTM_2020_hillshade <- shade(bon_Slope_rad, bon_Aspect_rad, 40, 270)
bon_Roughness<- terrain(bon_DTM_2020, v = "roughness")
bon_flowdir<- terrain(bon_DTM_2020, v="flowdir")
bon_planc<- Qfit(bon_DTM_2020, metrics = c("profc")) # how the surface bends along the direction of the slope.
bon_profc<- Qfit(bon_DTM_2020, metrics = c("planc")) # the curvature of the terrain perpendicular to the slope.
bon_meanc<- Qfit(bon_DTM_2020, metrics = c("meanc")) # general measure of the curvature of the surface, incorporating both the profile and plan curvatures into a single value.

#---------------------------------#
# Topographic Wetness Index (TWI) #
#---------------------------------#

flow_accum <- terrain(bon_DTM_2020, v = "flowdir", unit = "radians")
slope_radians <- terrain(bon_DTM_2020, v = "slope", unit = "radians")

# Calculate TWI (TWI = log(flow accumulation / tan(slope)))
bon_TWI <- log(flow_accum / tan(slope_radians))
bon_TWI[!is.finite(bon_TWI)] <- NA

#------------------------#
# Solar Radiation (PISR) #
#------------------------#

# Normalize hillshade values between 0 and 1 to approximate relative solar exposure
bon_solar <- bon_DTM_2020_hillshade / max(values(bon_DTM_2020_hillshade), na.rm = TRUE)

#-----------------------#
### Check the Rasters ###
#-----------------------#

summary(bon_Slope)
summary(bon_TPI)
summary(bon_TRI)
summary(bon_Aspect)
summary(bon_DTM_2020_hillshade)
summary(bon_Roughness)
summary(bon_flowdir)
summary(bon_planc)
summary(bon_profc)
summary(bon_meanc)
summary(bon_TWI)
summary(bon_solar)

#----------------------#
### Plot the Rasters ###
#----------------------#

plot(bon_Slope)
plot(bon_TPI)
plot(bon_TRI)
plot(bon_Aspect)
plot(bon_DTM_2020_hillshade)
plot(bon_Roughness)
plot(bon_flowdir)
plot(bon_planc)
plot(bon_profc)
plot(bon_meanc)
bon_planc_clamped <- clamp(bon_planc, lower = -0.5, upper = 0.5, values = TRUE)
plot(bon_planc_clamped)
bon_profc_clamped <- clamp(bon_profc, lower = -0.5, upper = 0.5, values = TRUE)
plot(bon_profc_clamped)
bon_meanc_clamped <- clamp(bon_meanc, lower = -0.5, upper = 0.5, values = TRUE)
plot(bon_meanc_clamped)
plot(bon_TWI)
plot(bon_solar)

#----------------------#
### Save the Rasters ###
#----------------------#

writeRaster(bon_Slope, "bon_Slope.tif", overwrite=T)
writeRaster(bon_TPI, "bon_TPI.tif", overwrite=T)
writeRaster(bon_TRI, "bon_TRI.tif", overwrite=T)
writeRaster(bon_Aspect, "bon_Aspect.tif", overwrite=T)
writeRaster(bon_DTM_2020_hillshade, "bon_Hillshade.tif",overwrite=T)
writeRaster(bon_flowdir, "bon_Flowdir.tif",overwrite=T)
writeRaster(bon_Roughness, "bon_Roughness.tif",overwrite=T)
writeRaster(bon_planc, "bon_Planc.tif",overwrite=T)
writeRaster(bon_profc, "bon_Profc.tif",overwrite=T)
writeRaster(bon_meanc, "bon_Meanc.tif",overwrite=T)
writeRaster(bon_TWI, "bon_TWI.tif",overwrite=T)
writeRaster(bon_solar, "bon_Solar.tif",overwrite=T)

#--------------------------#
##### Load the Rasters #####
#--------------------------#

bon_Slope <- rast('bon_Slope.tif')
bon_TPI<- rast('bon_TPI.tif')
bon_TRI<- rast('bon_TRI.tif')
bon_Aspect<- rast('bon_Aspect.tif')
bon_Roughness<- rast('bon_Roughness.tif')
bon_flowdir<- rast('bon_Flowdir.tif')
bon_planc<- rast('bon_Planc.tif')
bon_profc<- rast('bon_Profc.tif')
bon_DTM_2020_hillshade<- rast('bon_Hillshade.tif')
bon_meanc<- rast('bon_Meanc.tif')
bon_TWI<- rast('bon_TWI.tif')
bon_solar<- rast('bon_Solar.tif')

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### STEP 4: EXTRACT METRICS FROM POINT CLOUD ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#---------------#
# First Returns #
#---------------#

# Load LAS file 
pt_cloud <- readLAS("Dormouse_Normalised.laz")

# Set the CRS to EPSG:27700 (British National Grid)
projection(pt_cloud) <- "EPSG:27700"

summary(pt_cloud$Z)
plot(pt_cloud)

# Set any first returns below 0 to zero
pt_cloud@data$Z[pt_cloud@data$Z<0]<-0

# Filter the first returns (Class 1 corresponds to first returns)
first_returns <- filter_poi(pt_cloud, ReturnNumber == 1)

min(first_returns$Z)
max(first_returns$Z)
mean(first_returns$Z)
hist(first_returns$Z)

#-----------#
# Intensity #
#-----------#

# Filter the for Intensity
pt_intensity <- pt_cloud$Intensity

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### STEP 5: EXTRACT METRICS FOR EACH PLOT ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#-------------------------------------------------------#
### Extract canopy and topographic data for each plot ###
#-------------------------------------------------------#

# Create data frame to store data
plot_data <- data.frame(matrix(NA, nrow=dim(bontuchel)[1], ncol=27))
names(plot_data)<-c("Coppicing_Plot", "plotvalue_ID", "Coppicing_year","Coppicing_area", "Mean_Canopy_Height_2020", "Perc_under_2m_2020", 
                    "Can_cover_2020", "Height_cv_2020",
                    "Elevation","Slope","TRI","Max_Canopy_Height_2020", "FT_5_10m","FT_1.37m",
                    "FT_10_20m","Intensity","Aspect","TPI", 'Roughness', "LAI", "VCI", "Canopy_Volume_2020",
                    "FlowDir", "Plane_Curve", "Profile_Curve", "Mean_Curve", "Solar_Radiation")

plot_data$Coppicing_Plot<-bontuchel$ID_No
plot_data$Coppicing_year<-bontuchel$box_yr
plot_data$Coppicing_area<-bontuchel$Area

# Run loop across each plot and extract data
for (i in 1:dim(plot_data)[1]){
  
  #//////////#
  # CHM 2020 #
  #//////////#
  
  # Crop
  bon_2020_CHM_crop<-crop(bon_CHM_2020, bontuchel[i,])
  
  #--------------------------------#
  # Mean Top of Canopy Height 2020 #
  #--------------------------------#
  plot_data[i,"Mean_Canopy_Height_2020"]<-global(bon_2020_CHM_crop, fun='mean', na.rm=TRUE) 
  
  #-----------------------------------------#
  # Coefficient of variation in Height 2020 #
  #-----------------------------------------#
  # Coefficient of variation 2020
  plot_data[i,"Height_cv_2020"]<-global(bon_2020_CHM_crop, fun='sd', na.rm=TRUE)/global(bon_2020_CHM_crop, fun='mean', na.rm=TRUE) 
  
  #----------------------------------#
  # % pixels below 2m threshold 2020 #
  #----------------------------------#
  plot_data[i,"Perc_under_2m_2020"]<-100*(length(which(as.vector(na.omit(as.vector(bon_2020_CHM_crop)))<2))/length(as.vector(na.omit(as.vector(bon_2020_CHM_crop))))) 
  
  #----------------------------------#
  # % pixels above 2m threshold 2020 #
  #----------------------------------#
  plot_data[i,"Can_cover_2020"]<-100- plot_data[i,"Perc_under_2m_2020"]
  
  #------------------------#
  # Max Canopy Height 2020 #
  #------------------------#
  values_MCH_crop_2020 <- values(bon_2020_CHM_crop, na.rm = TRUE)
  Max_Canopy_Height_2020<-stats::quantile(values_MCH_crop_2020, probs=0.99, na.rm = TRUE)
  plot_data[i,"Max_Canopy_Height_2020"]<-Max_Canopy_Height_2020
  
  #---------------#
  # Canopy volume #
  #---------------#
  plot_data[i,"Canopy_Volume_2020"]<-canopy_volume(bon_2020_CHM_crop)
  
  #/////////////#
  # Point Cloud #
  #/////////////#
  
  # Convert the shapefile to a spatial object compatible with lidR
  polygon_ <- as(bontuchel[i,], "Spatial")
  target <- st_crs(27700)  # Define the target CRS using sf
  
  # Convert the target CRS to sp's CRS format
  target_sp <- CRS(st_as_text(target))
  
  # Transform the polygon_lidar using spTransform()
  polygon_ <- spTransform(polygon_, target_sp)
  
  # Crop point cloud
  fr_pt_crop<-clip_roi(first_returns, polygon_)
  
  #-------------------------------------------------------------------#
  # Calculate the percentage of first returns between 5 and 10 meters #
  #-------------------------------------------------------------------#
  
  # Filter
  first_returns_5_10m <- filter_poi(fr_pt_crop,(Z >= 5 & Z <= 10))
  
  # Get the total number of first returns and those between 5 and 10 meters
  total_first_returns <- nrow(fr_pt_crop)
  first_returns_5_10m_count <- nrow(first_returns_5_10m)
  
  # Calculate the percentage
  percent_5_10m <- (first_returns_5_10m_count / total_first_returns) * 100
  plot_data[i,"FT_5_10m"]<-percent_5_10m
  
  #-------------------------------------------------------------#
  # Calculate the percentage of first returns above 1.37 meters #
  #-------------------------------------------------------------#
  
  # Filter
  first_returns_1.37m <- filter_poi(fr_pt_crop,(Z >= 1.37))
  
  # Get the total number of first returns and those above 1.37 meters
  first_returns_1.37m_count <- nrow(first_returns_1.37m)
  
  # Calculate the percentage
  percent_1.37m <- (first_returns_1.37m_count / total_first_returns) * 100
  plot_data[i,"FT_1.37m"]<-percent_1.37m
  
  #-------------------------------------------------------------#
  # Calculate the percentage of first returns above 10-20m meters #
  #-------------------------------------------------------------#
  
  # Calculate the percentage of first returns 10-20m meters
  first_returns_10_20m <- filter_poi(fr_pt_crop,(Z >= 10 & Z <= 20))
  
  # Get the total number of first returns and those above 1.37 meters
  first_returns_10_20m_count <- nrow(first_returns_10_20m)
  
  # Calculate the percentage
  percent_10_20m <- (first_returns_10_20m_count / total_first_returns) * 100
  plot_data[i,"FT_10_20m"]<-percent_10_20m
  
  #-----------#
  # Intensity #
  #-----------#
  
  # Crop
  pt_crop<-clip_roi(pt_cloud, polygon_)
  
  # Mean intensity
  plot_data[i,"Intensity"] <- mean(pt_crop$Intensity)
  
  #-----#
  # LAI #
  #-----#
  
  # Save the cropped point cloud to a temporary file
  temp_file <- tempfile(fileext = ".laz")
  writeLAS(pt_crop, temp_file)
  
  # Calculate Leaf Area Density (LAD) from voxelization
  VOXELS_LAD <- lad.voxels(temp_file, grain.size = 2)
  
  # Calculate the LAD profile
  lad_profile <- lad.profile(VOXELS_LAD)
  
  # Calculate LAI derived from LAD profile
  lidar_lai <- lai(lad_profile)
  
  # Add LAI to dataframe
  plot_data[i, "LAI"] <- lidar_lai
  
  #-----#
  # VCI #
  #-----#
  
  ??lai
  
  # Calculate Vertical Complexity Index (VCI)
  vci <- VCI(pt_crop@data$Z, by = 1, zmax = 20)
  
  # Add VCI to dataframe
  plot_data[i, "VCI"] <- vci
  
  #/////////#
  # Terrain #
  #/////////#
  
  #-----------#
  # Elevation #
  #-----------#
  
  bon_DTM_2020_crop<-crop(bon_DTM_2020,bontuchel[i,])
  plot_data[i,"Elevation"]<-global(bon_DTM_2020_crop, fun='mean', na.rm=TRUE)
  
  #--------#
  # Aspect #
  #--------#
  
  Aspect_crop<-crop(bon_Aspect,bontuchel[i,])
  plot_data[i,"Aspect"]<-global(Aspect_crop, fun='mean', na.rm=TRUE)
  
  #-----#
  # TPI #
  #-----#
  
  TPI_crop<-crop(bon_TPI,bontuchel[i,])
  plot_data[i,"TPI"]<-global(TPI_crop, fun='mean', na.rm=TRUE)
  
  #-------#
  # Slope #
  #-------#
  
  Slope_crop<-crop(bon_Slope,bontuchel[i,])
  plot_data[i,"Slope"]<-global(Slope_crop, fun='mean', na.rm=TRUE)
  
  #-----#
  # TRI #
  #-----#
  ??global
  TRI_crop<-crop(bon_TRI,bontuchel[i,])
  plot_data[i,"TRI"]<-global(TRI_crop, fun='mean', na.rm=TRUE)
  
  #-----------#
  # Roughness #
  #-----------#
  
  roughness_crop<-crop(bon_Roughness,bontuchel[i,])
  plot_data[i,"Roughness"]<-global(roughness_crop, fun='mean', na.rm=TRUE)
  
  #---------#
  # FlowDir #
  #---------#
  
  flowdir_crop<-crop(bon_flowdir,bontuchel[i,])
  plot_data[i,"FlowDir"]<-global(flowdir_crop, fun='mean', na.rm=TRUE)
  
  #-----------------#
  # Plane Curvature #
  #-----------------#
  
  planc_crop<-crop(bon_planc,bontuchel[i,])
  plot_data[i,"Plane_Curve"]<-global(planc_crop, fun='mean', na.rm=TRUE)
  
  #-------------------#
  # Profile Curvature #
  #-------------------#
  
  profc_crop<-crop(bon_profc,bontuchel[i,])
  plot_data[i,"Profile_Curve"]<-global(profc_crop, fun='mean', na.rm=TRUE)
  
  #----------------#
  # Mean Curvature #
  #----------------#
  
  meanc_crop<-crop(bon_meanc,bontuchel[i,])
  plot_data[i,"Mean_Curve"]<-global(meanc_crop, fun='mean', na.rm=TRUE)
  
  #-----------------#
  # Solar Radiation #
  #-----------------#
  
  solar_crop<-crop(bon_solar,bontuchel[i,])
  plot_data[i,"Solar_Radiation"]<-global(solar_crop, fun='mean', na.rm=TRUE)
  
  # Give each plotvalue a unique ID
  plot_data[i, "plotvalue_ID"] <- paste0("pv_", i)
  
  # Progress (% of plots completed)
  print((i/dim(plot_data)[1])*100)
  
}

# Minimum DTM Value (for Height above lowest point, HALP)
min_DTM <- min(values(bon_DTM_2020), na.rm = TRUE)
plot_data$HALP<-plot_data$Elevation-min_DTM

# Aspect (cos transformed)
plot_data$Aspect_Cos<-cos(plot_data$Aspect)

#------------#
# Check data #
#------------#

summary(plot_data)

# Canopy Structure (CHM 2020):
#   Mean_Canopy_Height_2020:
#   6.5–17.3 m: Reasonable for young to mid-aged secondary woodland or coppice regrowth.
# Perc_under_2m_2020:
#   ~0–16.5%:
#   Very few areas have canopy below 2 m — suggesting mostly mature regrowth, not super young scrub.
# Can_cover_2020:
#   ~83–99.9%:
#   Mostly very dense canopy — matches a mature or nearly closed forest/coppice structure.
# Height_cv_2020 (coefficient of variation in height):
#   0.2–0.7:
#   Moderate variability — again typical of heterogeneous canopy recovery.
# Max_Canopy_Height_2020:
#   11.2–32.3 m:
#   Looks realistic. 30+ m suggests some taller veteran trees or parts of plots that weren't fully cleared during coppicing.
# Canopy_Volume_2020:
# 10,597–89,148 (units not shown):
# Range seems fine assuming voxel volumes or canopy volume units like m³/ha (depends on your function canopy_volume() scaling).

# LiDAR Point Cloud Metrics:
#   FT_5_10m:
#   2–74%:
#   Some plots dominated by midstory (5–10 m) growth; others less so.
# FT_1.37m:
#   85–99.9%:
#   Very high proportion of returns above breast height (1.37 m), as expected in woody vegetation.
# FT_10_20m:
#   4.5–89.8%:
#   Again a wide variation — reasonable given a gradient of regrowth/maturity.
# Intensity:
#   ~365–770:
#   Values look plausible for laser return intensities (depends on your sensor and calibration, but these are typical numbers).
# LAI (Leaf Area Index):
#   2.1–5.3:
#   Matches dense temperate forests (LAI in mature broadleaf forests often around 3–6).
# VCI (Vertical Complexity Index):
#   0.74–0.99:
#   Quite high, suggesting structurally complex canopies — plausible in mixed aged coppice regrowth.

# Terrain Variables (DTM derived):
#   Elevation:
#   128–171 m:
#   Moderately lowland topography.
# Slope:
#   ~0.14–0.4 radians (~8–23°):
#   Slopes are gentle to moderate.
# Aspect:
#   ~2–5 radians (around 114°–286° degrees) → aspect mainly SE to W.
# Aspect_Cos transformation: ranges from -0.99 to +0.32 — this is expected.
# TPI, TRI, Roughness:
#   All ranges make sense; small, rugged variation typical of wooded lowlands.
# FlowDir, Plane_Curve, Profile_Curve, Mean_Curve:
#   All values small and centered around 0, which matches natural microtopography.
# Solar_Radiation:
#   0.46–0.85 (relative units?) — moderate values, possibly representing normalized solar radiation exposure.
# HALP (Height Above Lowest Point):
#   6.75–50.35 m:
#   A 40 m vertical range is reasonable given the elevation range and terrain.

# List of columns you WANT to be numeric
cols_to_numeric <- c("Mean_Canopy_Height_2020", "Perc_under_2m_2020", "Can_cover_2020", 
                     "Height_cv_2020", "Elevation", "Slope", "TRI", "Max_Canopy_Height_2020",
                     "FT_5_10m", "FT_1.37m", "FT_10_20m", "Intensity", "Aspect",
                     "TPI", "Roughness", "LAI", "VCI", "Canopy_Volume_2020", "FlowDir",
                     "Plane_Curve", "Profile_Curve", "Mean_Curve", "Solar_Radiation",
                     "HALP", "Aspect_Cos")

# Now, convert these columns to numeric
plot_data[cols_to_numeric] <- lapply(plot_data[cols_to_numeric], function(x) as.numeric(as.character(x)))

# Select only numeric columns
numeric_data <- plot_data[sapply(plot_data, is.numeric)]

# Turn data long
numeric_data_long <- numeric_data %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value")

# Plot all
ggplot(numeric_data_long, aes(x = Value)) +
  geom_histogram(fill = "pink", color = "black", bins = 20) +
  facet_wrap(~ Variable, scales = "free", ncol = 4) +
  theme_minimal()

# Plot all
ggplot(numeric_data_long, aes(x = Value)) +
  geom_boxplot(fill = "pink", color = "black") +
  facet_wrap(~ Variable, scales = "free", ncol = 4) +
  theme_minimal()

# Save the dataframe
write.csv(plot_data,"Topo_metrics.csv",row.names = F)

# Load
plot_data<- read.csv("Topo_metrics.csv")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### STEP 7: DORMOUSE DATA ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%##

# Load in data and create a dataframe 
dormice<-vect('Individual Species - Muscardinus avellanarius_points.shp')
dor_df<- as.data.frame(dormice)

dormice2<-vect('Coed Fron Wyllt, Bontuchel - Dormouse Monitoring Results May21 to Oct24.xlsx', layer = "Records")
dor_df_2<- as.data.frame(dormice2)

dormice3<-vect('Bontuchel NDMP.xlsx', layer='NDMP_box')
dor_df_3<- as.data.frame(dormice3)

dormice4<-vect('E130908_new export inc custom fields.xlsx', layer='Records')
dor_df_4<- as.data.frame(dormice4)

# Fix formatting
colnames(dor_df_2) <- dor_df_2[1, ]
dor_df_2 <- dor_df_2[-1, ]
dor_df_2 <- dor_df_2[, !apply(dor_df_2, 2, function(x) all(is.na(x)))]

colnames(dor_df_3) <- dor_df_3[1, ]
dor_df_3 <- dor_df_3[-1, ]
dor_df_3 <- dor_df_3[, !apply(dor_df_3, 2, function(x) all(is.na(x)))]

colnames(dor_df)
colnames(dor_df_2)
colnames(dor_df_3)
colnames(dor_df_4)

#---------------------#
# Clean up site names #
#---------------------#

#////////#
# dor_df #
#////////#

# Remove leading/trailing spaces
dor_df$SITE_NAME <- trimws(dor_df$SITE_NAME)

#Replace hyphens, underscores, or other separators with spaces
dor_df$SITE_NAME  <- gsub("[-_]", " ", dor_df$SITE_NAME )

# Capitalize the first letter of each word (e.g., Site A, Site B)
dor_df$SITE_NAME  <- tools::toTitleCase(dor_df$SITE_NAME )

# Remove particular words (case-sensitive)
dor_df$SITE_NAME <- gsub("Bontuchel:", "", dor_df$SITE_NAME)
dor_df$SITE_NAME <- gsub("Coed Fron Wyllt,", "", dor_df$SITE_NAME)
dor_df$SITE_NAME <- gsub("Bontuchel, Box ", "", dor_df$SITE_NAME)
dor_df$SITE_NAME <- gsub("Box ", "", dor_df$SITE_NAME)
dor_df$SITE_NAME[grepl("Bontuchel", dor_df$SITE_NAME)] <- NA
dor_df$SITE_NAME[grepl("Whole Site", dor_df$SITE_NAME)] <- NA

# Replace entries without 'Canopy' with NA, create a new column for canopy boxes separate to ground boxes
dor_df$SITE_NAME_CANOPY_BOXES<-dor_df$SITE_NAME
dor_df$SITE_NAME_CANOPY_BOXES[!grepl("Canopy", dor_df$SITE_NAME_CANOPY_BOXES)] <- NA

# Replace entries with 'Canopy' as NA
dor_df$SITE_NAME[grepl("Canopy", dor_df$SITE_NAME)] <- NA

dor_df$SITE_NAME <- trimws(dor_df$SITE_NAME)
unique(dor_df$SITE_NAME)

dor_df$SITE_NAME[grepl("Bont Uchel", dor_df$SITE_NAME)] <- NA
dor_df$SITE_NAME[grepl("Coed Fron Wyllt", dor_df$SITE_NAME)] <- NA
dor_df$SITE_NAME[grepl("Coed y Fron Wyllt", dor_df$SITE_NAME)] <- NA

#------------------#
# Clean up dataset #
#------------------#

# Create data frame to store data
dormice_df<-data.frame(matrix(dim(dor_df)[1],6,data=NA))
names(dormice_df)<-c("Longitude", "Latitude", "Dor_Box","Year_DM_Recorded", "No_of_Dormice","Canopy_Box")

dormice_df$Longitude<-dor_df$EASTING
dormice_df$Latitude<-dor_df$NORTHING
dormice_df$Dor_Box<-dor_df$SITE_NAME
dormice_df$Year_DM_Recorded<-dor_df$YEAR
dormice_df$No_of_Dormice<-dor_df$ABUNDANCE
dormice_df$Canopy_Box<-dor_df$SITE_NAME_CANOPY_BOXES

dormice<-vect(dormice_df,geom=c("Longitude", "Latitude"))
plot(dormice)

#--------------#
# DM Abundance #
#--------------#

# Where no abundance is specified it should be assumed at least one was found (from notes)
dormice_df <- dormice_df %>%
  mutate(No_of_Dormice = ifelse(is.na(No_of_Dormice), 1, No_of_Dormice))

unique(dormice_df$No_of_Dormice)

#//////////#
# dor_df_2 #
#//////////#

colnames(dor_df_2)

# Pivot data
df2_long <- pivot_longer(dor_df_2, cols = -`Box No.`, names_to = "Date_DM_Recorded", values_to = "No_of_Dormice")
df2_long$Date_DM_Recorded <- as.Date(df2_long$Date_DM_Recorded, format = "%Y/%m/%d")
unique(df2_long$Date_DM_Recorded )

# Fix typo date
df2_long$Date_DM_Recorded[df2_long$Date_DM_Recorded == "0025-09-20"] <- "2020-09-25"
df2_long$Date_DM_Recorded <- as.Date(df2_long$Date_DM_Recorded, format = "%Y/%m/%d")
df2_long$Year_DM_Recorded <- format(df2_long$Date_DM_Recorded, "%Y")

# Rename
df2_long <- df2_long %>% rename(Dor_Box = `Box No.`)

#//////////#
# dor_df_3 #
#//////////#

colnames(dor_df_3)

# Remove unneccessary columns
dor_df_3$'NA'<-NULL
dor_df_3$Sex<-NULL
dor_df_3$Age<-NULL
dor_df_3$Weight<-NULL
dor_df_3$Torpid<-NULL
dor_df_3$BreedingCon<-NULL
dor_df_3$Notes<-NULL
dor_df_3$NA.1<-NULL

# Rename
dor_df_3 <- dor_df_3 %>% rename(Dor_Box = `BoxID`)
dor_df_3 <- dor_df_3 %>% rename(Date_DM_Recorded = `date`)
dor_df_3 <- dor_df_3 %>% rename(No_of_Dormice = `NoPerBox`)

# Add year
dor_df_3$Date_DM_Recorded <- as.Date(dor_df_3$Date_DM_Recorded, format = "%Y/%m/%d")
dor_df_3$Year_DM_Recorded <- format(dor_df_3$Date_DM_Recorded, "%Y")

#//////////#
# dor_df_4 #
#//////////#

# Remove leading/trailing spaces
dor_df_4$'Site Name' <- trimws(dor_df_4$'Site Name')

# Replace hyphens, underscores, or other separators with spaces
dor_df_4$'Site Name'  <- gsub("[-_]", " ", dor_df_4$'Site Name' )

# Capitalize the first letter of each word (e.g., Site A, Site B)
dor_df_4$'Site Name'  <- tools::toTitleCase(dor_df_4$'Site Name' )

# Remove particular words (case-sensitive)
dor_df_4$'Site Name' <- gsub("Bontuchel:", "", dor_df_4$'Site Name')
dor_df_4$'Site Name' <- gsub("Coed Fron Wyllt,", "", dor_df_4$'Site Name')
dor_df_4$'Site Name' <- gsub("Bontuchel, Box ", "", dor_df_4$'Site Name')
dor_df_4$'Site Name' <- gsub("Box ", "", dor_df_4$'Site Name')
dor_df_4$'Site Name'[grepl("Bontuchel", dor_df_4$'Site Name')] <- NA
dor_df_4$'Site Name'[grepl("Whole Site", dor_df_4$'Site Name')] <- NA
dor_df_4$'Site Name'[grepl("Weir Site 1", dor_df_4$'Site Name')] <- NA
dor_df_4$'Site Name'[grepl("Weir Site 2", dor_df_4$'Site Name')] <- NA
dor_df_4$'Site Name'[grepl("Weir Site 3", dor_df_4$'Site Name')] <- NA
dor_df_4$'Site Name'[grepl("Weir Site 4", dor_df_4$'Site Name')] <- NA
dor_df_4$'Site Name'[grepl("Weir Site 5", dor_df_4$'Site Name')] <- NA
dor_df_4$'Site Name'[grepl("Weir Site 6", dor_df_4$'Site Name')] <- NA
dor_df_4$'Site Name'[grepl("Weir Site 7", dor_df_4$'Site Name')] <- NA
dor_df_4$'Site Name'[grepl("Weir Site 8", dor_df_4$'Site Name')] <- NA
dor_df_4$'Site Name'[grepl("Coed y Fron Wyllt", dor_df_4$'Site Name')] <- NA
dor_df_4$'Site Name'[grepl("Coed Fron Wyllt", dor_df_4$'Site Name')] <- NA
dor_df_4$'Site Name'[grepl("Bont Uchel", dor_df_4$'Site Name')] <- NA

colnames(dor_df_4)
unique(dor_df_4$'Site Name')

# Replace entries without 'Canopy' with NA, create a new column for canopy boxes separate to ground boxes
dor_df_4$SITE_NAME_CANOPY_BOXES<-dor_df_4$'Site Name'
dor_df_4$SITE_NAME_CANOPY_BOXES[!grepl("Canopy", dor_df_4$SITE_NAME_CANOPY_BOXES)] <- NA

# Replace entries with 'Canopy' as NA
dor_df_4$'Site Name'[grepl("Canopy", dor_df_4$'Site Name')] <- NA

dor_df_4$'Site Name' <- trimws(dor_df_4$'Site Name')
unique(dor_df_4$'Site Name')

#------------------#
# Clean up dataset #
#------------------#

colnames(dor_df_4)

# Create data frame to store data
dormice_df_4<-data.frame(matrix(dim(dor_df_4)[1],7,data=NA))
names(dormice_df_4)<-c("Longitude", "Latitude", "Dor_Box","Year_DM_Recorded", "No_of_Dormice", "Negative_Rec", "Canopy_Box")

dormice_df_4$Longitude<-dor_df_4$Easting
dormice_df_4$Latitude<-dor_df_4$Northing
dormice_df_4$Dor_Box<-dor_df_4$'Site Name'
dormice_df_4$Year_DM_Recorded<-dor_df_4$Year
dormice_df_4$No_of_Dormice<-dor_df_4$Abundance
dormice_df_4$Negative_Rec<-dor_df_4$"Negative Record"
dormice_df_4$Canopy_Box<-dor_df_4$SITE_NAME_CANOPY_BOXES

dormice_4<-vect(dormice_df_4,geom=c("Longitude", "Latitude"))
plot(dormice_4, col='pink', add=TRUE)

#--------------#
# DM Abundance #
#--------------#

# Fix abundance based on Negative Record logic
dormice_df_4 <- dormice_df_4 %>%
  # If Negative Record is "Yes" and No_of_Dormice is NA, set it to 0
  mutate(No_of_Dormice = ifelse(Negative_Rec == "Yes" & is.na(No_of_Dormice), 0, No_of_Dormice)) %>%
  # If Negative Record is "No" and No_of_Dormice is still NA, leave as NA (they found dormice but didn't record how many)
  mutate(No_of_Dormice = ifelse(Negative_Rec == "No" & is.na(No_of_Dormice), NA, No_of_Dormice)) %>%
  # Where abundance is still NA (and not a negative record), assume at least 1 dormouse was found
  mutate(No_of_Dormice = ifelse(is.na(No_of_Dormice), 1, No_of_Dormice))

# Check result
unique(dormice_df_4$No_of_Dormice)

#-------#
# Merge #
#-------#

head(boxes_df)
head(dormice_df)
head(df2_long)
head(dor_df_3)
head(dormice_df_4)

nrow(dormice_df)
nrow(df2_long)
nrow(dor_df_3)
nrow(dormice_df_4)

# Convert all columns to characters
dormice_df <- dormice_df %>% mutate_all(as.character)
df2_long <- df2_long %>% mutate_all(as.character)
dor_df_3 <- dor_df_3 %>% mutate_all(as.character)
dormice_df_4 <- dormice_df_4 %>% mutate_all(as.character)

# Add a column indicating the dataset source
dormice_df <- dormice_df %>% mutate(Source = "dormice_df", Num_Cols = ncol(dormice_df))
df2_long <- df2_long %>% mutate(Source = "df2_long", Num_Cols = ncol(df2_long))
dor_df_3 <- dor_df_3 %>% mutate(Source = "dor_df_3", Num_Cols = ncol(dor_df_3))
dormice_df_4 <- dormice_df_4 %>% mutate(Source = "dormice_df_4", Num_Cols = ncol(dormice_df_4))

# Merge datasets
merged_df <- bind_rows(dormice_df, df2_long, dor_df_3, dormice_df_4)

# Identify common columns
common_cols <- Reduce(intersect, list(names(dormice_df), names(df2_long), names(dor_df_3), names(dormice_df_4)))

# Remove duplicates while keeping the record with:
# 1. The most columns
# 2. The fewest NA values (as a tie-breaker)
deduplicated_df <- merged_df %>%
  arrange(desc(Num_Cols), desc(rowSums(!is.na(.)))) %>%  # Sort first by Num_Cols, then by non-NA count
  distinct(across(all_of(common_cols)), .keep_all = TRUE)   # Keep only the best record per duplicate set


# Check if exact duplicates remain and remove them
dormice_df <- deduplicated_df %>% distinct()

dormice_df <-dormice_df %>%dplyr::select(-Source, -Num_Cols)  # Remove helper columns

# Check the final number of rows
nrow(dormice_df) # 4869 obvs

#---------------------------------------------------#
# Determine which Dormouse boxes are in which plots #
#---------------------------------------------------#

# Merge data dormice records with the associated latitude and longitude
dormice_df$Dor_Box<-as.numeric(dormice_df$Dor_Box)
dormice_df$Longitude<-NULL
dormice_df$Latitude<-NULL
boxes_df$Dor_Box<-as.numeric(boxes_df$Dor_Box)

head(boxes_df)
head(dormice_df)
merged_df<- merge(dormice_df, boxes_df, by= 'Dor_Box', all = TRUE)
head(merged_df)
nrow(dormice_df)
nrow(merged_df)

# Create a vector of this and save
merged<- vect(merged_df, geom=c('Longitude', 'Latitude'))
plot(merged)
set.crs(merged, "epsg:27700")
writeVector(merged, 'DM_and_boxes.shp', filetype = "ESRI Shapefile", overwrite=TRUE)

# See which Dormice records intersect which plots
DM_in_plots <- terra::intersect(merged, bontuchel)
DM_in_plots_df<-as.data.frame(DM_in_plots)
head(DM_in_plots_df)
nrow(DM_in_plots_df) # Only 1336 dormice records out of 4869 were in our plots

# Tidy data
DM_in_plots_df$DISS<-NULL
DM_in_plots_df$Id<-NULL

# Identify cases where a single Dor_Box has multiple Lat/Long values
lat_long_check <- merged_df %>%
  group_by(Dor_Box) %>%
  summarise(n_unique_lat = n_distinct(Latitude), n_unique_long = n_distinct(Longitude)) %>%
  filter(n_unique_lat > 1 | n_unique_long > 1)

# Print problematic Dor_Box entries (if any)
print(lat_long_check)

merged_df_unique <- merged_df %>%
  dplyr::select(Dor_Box, Latitude, Longitude) %>%
  distinct(Dor_Box, .keep_all = TRUE)  # Ensures each Dor_Box has one Lat/Long pair

DM_in_plots_df <- merge(DM_in_plots_df, merged_df_unique, by = "Dor_Box", all.x = TRUE)
nrow(DM_in_plots_df)

DM_in_plots<-vect(DM_in_plots_df, geom=c('Longitude', 'Latitude'))
set.crs(DM_in_plots, "epsg:27700")
writeVector(DM_in_plots, 'DM_in_plots.shp', filetype = "ESRI Shapefile", overwrite=TRUE)

# Check plot
plot.new()
bon_plot<-crop(bon_CHM_2020, bontuchel)
plot(bon_plot, col=viridis(100))

plot(bontuchel, add=TRUE)
plot(DM_in_plots, add=TRUE, col='pink')

#----------------------------------------#
# Plot all dormice records in each year #
#----------------------------------------#

DM_in_plots<-vect('DM_in_plots.shp')

# Get unique years
years <- unique(DM_in_plots$Year_DM_Re)

# Loop through each year and plot
for (yr in years) {
  # Subset dormice observations for the current year
  df_year <- DM_in_plots[DM_in_plots$Year_DM_Re == yr, ]
  
  # Plot Bontuchel outline
  plot(bontuchel, main = paste("Dormouse Observations -", yr), border = 'hotpink')
  
  # Plot dormice observations for the year
  points(df_year, col = 'violet', pch = 19)
}

#--------------#
# DM Abundance #
#--------------#

unique(DM_in_plots_df$No_of_Dormice)

# Replace the specific string with the abundance in the column 'No_of_Dormice'
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "1 adult, 3 live juvenile female, 1 dead juvenile female"] <- 5
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "1 adult female, 1 adult male"] <- 2
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "1 Adult Male"] <- 1
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "1 Sub-Adult Male" ] <- 1
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "1 male"] <- 1
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "1 female, 1 male"] <- 2
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "1 adult male"] <- 1
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "1 adult male, 1 adult female"] <- 2
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "3x Dormice. 1x male juvenile, 12g, active, chipped by LB No.067961. 1x female juvenile, 11g, active, chipped by LB No.067965. 1x male juvenile, 12g, active, chipped by LB No.067977."] <- 3
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "1 x DM Male + nest. Chipped by Lucy Boyett No.077905 Sub-adult Male, 18g, Torpid"] <- 1
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "1 Sub-Adult Male"] <- 1
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "1x Adult DM (not captured) & 4x Juveniles (1x Female 'grey eyes open', 9g, Active, Chipped by Roger Trout No. 074032), (1x Female 'grey eyes open', 9g, Active, Chipped by Roger Trout No. 077837), 1x Female 'grey eyes open', 8g, Active, Chipped by Roger Trout No. 077897) & (1x DEAD Female 'grey eyes open', 8g,) " ] <- 5
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "1 adult, 3 live juvenile female, 1 dead juvenile female"] <- 4
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "1 x dormouse male, chipped No.078185, 27g, active"] <- 1
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "1 male"] <- 1
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "1x dormouse present in very small tight nest 1x Female Adult, 19g, Active (very sleepy). Chip No.078107."] <- 1
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "1 female, 1 male"] <- 2
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "2x DM, 1 x female previously chipped No.078107, 31g, active. 1 x male, chipped No.077832, 29g, active."] <- 2
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "1 adult female, 1 adult male"] <- 2
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "1 Adult Male"] <- 1
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "1 x DM Male + nest. Chipped Previously No.077876 Adult Male, 16.5g, Active"] <- 1
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "2x DM (1x Adult Female, 26g, Active, Chipped this session by Roger Trout No.077856) & (1x Adult Female, 22g, Active, Chipped this session by Roger Trout No.077907)"] <- 2
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "2x DM + Nest"] <- 2
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "1 x Adult DM & at least 3x pink young. Adult chipped by Roger Trout. No.077925. Adult Female, lactating, 20.5g, Active. "] <- 4
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "1 juvenile male"] <- 1
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "2 adult female"] <- 2
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "1 x DM Male + nest. Chipped Previously No.077983 Adult Male, 35g, Active"] <- 1
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "1 Adult Male"] <- 1
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "1 adult male"] <- 1
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "1 adult male, 1 adult female"] <- 2
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "2x DM (1x Adult Male, 18g, Active, Chipped by Roger Trout No.077983) & (1x Adult Female, 16g, Active, Chipped by Roger Trout No.077932)"] <- 2
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "1 DM juv male chipped no.078001, 15g, active"] <- 1
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "1 x DM Male + nest. Chipped by Rhian Hughes 18/05/2022 No.078120. Adult Male, 14g, Active, non-breeding"] <- 1
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "1 adult male (non-breeding)"] <- 1
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "2x DM & Nest (1x Adult Female, 27g, Active, Chipped by Rhian Hughes No.077850) & (1x Adult Male, 33g, Active, Chipped Previously No.077876)"] <- 2
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "1x DM (Juvenile Male, 17g, Active, Chipped this session by ? No.078162)"] <- 1
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "1 adult female (lactating), at least 3 pink young"] <- 4
# DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "x"] <- 0 #check this
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "Shrew"] <- 0 
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "Woodmouse"] <- 0 
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "vacant dormouse nest - old"] <- 0 
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "Possible Dormouse Nest"] <- 0 
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "dormouse nest"] <- 0 
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "Pygmy Shrew"] <- 0 
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "woodmouse nest"] <- 0 
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "Dormouse nest"] <- 0 
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "old woodmouse nest"] <- 0 
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "Old dormouse nest"] <- 0 
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "vacant dormouse"] <- 0 
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "Old Dormouse Nest Material"] <- 0 
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "Potential early dormouse nest. "] <- 0 
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "Vacant Dormouse Nest"] <- 0 
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "woodmouse"] <- 0 
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "empty DM nest"] <- 0 
DM_in_plots_df$No_of_Dormice[DM_in_plots_df$No_of_Dormice == "ecmpty DM nest"] <- 0 

unique(DM_in_plots_df$No_of_Dormice)
# Set X as NAs as have no clarification if this is presence or absence 

# Update the No_of_Dormice column so no strings are left
DM_in_plots_df <- DM_in_plots_df %>%
  mutate(No_of_Dormice = case_when(
    No_of_Dormice %in% c("0","1", "2", "3", "4", "5", "6") ~ No_of_Dormice,
    TRUE ~ NA_character_
  ))

# How many boxes have records?
length(unique(DM_in_plots_df$Dor_Box[!is.na(DM_in_plots_df$No_of_Dormice)])) # 78

# How many plots have records?
length(unique(DM_in_plots_df$ID_No[!is.na(DM_in_plots_df$No_of_Dormice)])) # 30

# Save the dataframe
write.csv(DM_in_plots_df,"DM_in_plots_df.csv",row.names = F)

#------------------------------------------------#
# Plot all positive dormice records in each year #
#------------------------------------------------#

DM_in_plots_positive<- DM_in_plots_df %>%
  dplyr::filter(No_of_Dormice>=0)

# Get unique years
years <- unique(DM_in_plots_positive$Year_DM_Recorded)

# Loop through each year and plot
for (yr in years) {
  # Subset dormice observations for the current year
  df_year <- DM_in_plots_positive[DM_in_plots_positive$Year_DM_Recorded == yr, ]
  
  # Plot Bontuchel outline
  plot(bontuchel, main = paste("Dormouse Observations -", yr), border = 'hotpink')
  
  # Plot dormice observations for the year
  dorm<-vect(df_year, geom=c('Longitude', 'Latitude'))
  
  plot(dorm, add=TRUE, col='violet')
  
}

#-------#
# Merge #
#-------#

# Load
plot_data<- read.csv("Topo_metrics.csv")
DM_in_plots_df<- read.csv("DM_in_plots_df.csv")

# count all non na values in df
sum(!is.na(DM_in_plots_df$No_of_Dormice)) #963
hist(DM_in_plots_df$No_of_Dormice)

# Rename the column 'old_name' to 'new_name'
DM_in_plots_df <- DM_in_plots_df %>% rename(Coppicing_Plot = ID_No)

head(plot_data)
head(DM_in_plots_df)
nrow(plot_data)
nrow(DM_in_plots_df)

length(unique(DM_in_plots_df$Coppicing_Plot))  # Number of unique Coppicing_Plot values in dormice data
length(unique(plot_data$Coppicing_Plot))       # Number of unique Coppicing_Plot values in plot data. 30 of 57 plots had dormice records in.

# Merge
dorm_data <- merge(DM_in_plots_df, plot_data, by = "Coppicing_Plot", all = TRUE)
nrow(dorm_data)

citation('leafR')
# Check names
colnames(dorm_data)
dorm_data$box_yr<-NULL
dorm_data$Area<-NULL

# Remove duplicates
dorm_data <- dorm_data %>%
  distinct()

# Save the dataframe
write.csv(dorm_data,"DM_metrics.csv",row.names = F)

# Load
dorm_data<-read.csv('DM_metrics.csv')
sum(!is.na(dorm_data$No_of_Dormice))
# 744 

# Actual dormice found:
dorm_data_summary <- dorm_data %>%
  group_by(Coppicing_Plot, Year_DM_Recorded, Dor_Box) %>%
  summarise(total_Dormice = sum(No_of_Dormice, na.rm = TRUE), .groups = "drop")

sum(dorm_data_summary$total_Dormice) # 385

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### STEP 8: EARTH TRACK DATA ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#-------------------------------------#
# Load in data and create a dataframe #
#-------------------------------------#

ET_points<-vect('ET_Points.gpkg')
plot(ET_points)
set.crs(ET_points, "epsg:27700")

# Convert the SpatVector to a data frame
ET_points_df <- as.data.frame(ET_points)
colnames(ET_points_df)
length(unique(ET_points_df$ET_ID))
# 71 points in total


ET_in_plots <- terra::intersect(ET_points, bontuchel)
ET_in_plots_df<-as.data.frame(ET_in_plots)

# Give each ET point an ID
ET_in_plots_df$ET_ID <- paste0("ET_", seq_len(nrow(ET_in_plots_df)))
length(unique(ET_in_plots_df$ET_ID))
# 62 of those points were in plots

# Remove unnecessary columns
ET_in_plots_df$Area<-NULL
ET_in_plots_df$DISS<-NULL
ET_in_plots_df$box_yr<-NULL
ET_in_plots_df$Id<-NULL

# Extract the coordinates
coords <- geom(ET_in_plots)[, c("x", "y")]

# Combine the coordinates with the attribute data
ET_in_plots_df <- cbind(ET_in_plots_df, coords)

# Rename the coordinate columns
names(ET_in_plots_df)[names(ET_in_plots_df) == "x"] <- "ET_Latitude"
names(ET_in_plots_df)[names(ET_in_plots_df) == "y"] <- "ET_Longitude"

plot(bontuchel)
plot(ET_points, add=TRUE, col='pink')

# Rename
names(ET_in_plots_df)[names(ET_in_plots_df) == "ID_No"] <- "Coppicing_Plot"
length(unique(ET_in_plots_df$Coppicing_Plot))
# ET points taken in 56 of the plots 

# Remove columns where all is NA
colnames(ET_in_plots_df)
et_data <- ET_in_plots_df %>% select_if(~ !all(is.na(.)))

# Rename
names(et_data)[names(et_data) == "_CREATION_DATE"] <- "ET_Time"
names(et_data)[names(et_data) == "box_yr"] <- "Coppicing_Year"

# Only keep points where accuracy is <=6
et_data<-  et_data %>% filter(locationAccuracy <= 6)

# Remove unnecessary columns
et_data$locationAc<-NULL
et_data$photoSky<-NULL
et_data$photoAhead<-NULL
et_data$photoLeft<-NULL
et_data$photoRight<-NULL
et_data$photoBehin<-NULL
et_data$photoGroun<-NULL
et_data$`_URI`<-NULL
et_data$DISS<-NULL
et_data$Id<-NULL
et_data$radius<-NULL
et_data$`_URI`<-NULL

colnames(et_data)

# Save
write.csv(et_data,"et_data.csv",row.names = F)

length(unique(et_data$ET_ID))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### STEP 9: COMBINE DATASETS ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

dorm_data<-read.csv('DM_metrics.csv')
et_data<-read.csv('et_data.csv')

head(dorm_data, n=1)
head(et_data, n=1)
nrow(dorm_data)
nrow(et_data)

# Check if there are duplicates in the Coppicing_Plot column
sum(duplicated(et_data$Coppicing_Plot))    # Should return 0 if no duplicates
# 6 plots have more than one ET record record in (must be taken into account later)

# Compare unique Coppicing_Plot values from both datasets
setdiff(unique(dorm_data$Coppicing_Plot), unique(et_data$Coppicing_Plot))  # Shows Coppicing_Plot in dorm_data not in et_data
setdiff(unique(et_data$Coppicing_Plot), unique(dorm_data$Coppicing_Plot))  # Shows Coppicing_Plot in et_data not in dorm_data

# Now merge this aggregated data with dorm_data by Coppicing_Plot
final_data <- merge(et_data, dorm_data, by = "Coppicing_Plot", all = TRUE)

# Check the final merged data
head(final_data)
nrow(final_data)
length(unique(final_data$ET_ID))

# Check for any missing values in key columns
sum(is.na(final_data$ET_ID))
sum(is.na(final_data$Coppicing_Plot))

# joining these datasets will introduce duplicates due to the one-to-many relationship between the earthtrack data (et_data) 
# and dormice data (dorm_data). Each Coppicing_Plot in et_data typically has a single record, whereas dorm_data has multiple 
# records per Coppicing_Plot.
# If both dorm_data and et_data contain multiple records for the same Coppicing_Plot, the merge() function will create a Cartesian join (many-to-many merge). 
# This means that each combination of matching records will be represented in the output.

#----------#
# Clean up #
#----------#

# Remove columns where all is NA
final_data <- final_data %>% select_if(~ !all(is.na(.)))

nrow(final_data)
length(unique(et_data$ET_ID))

#---------------#
# Rename layers #
#---------------#

head(final_data)
colnames(final_data)

#	lc_ represent the data from the land cover module.
#	change_ represent the data from the change module.
#	under_ represent the data from the understory module.

names(final_data)[names(final_data) == "lc_cCover"] <- "Layer_0_cover"
names(final_data)[names(final_data) == "lc_cHeight"] <- "Layer_0_height"
names(final_data)[names(final_data) == "lc_phenology"] <- "ET_Phenology"
names(final_data)[names(final_data) == "lc_l3classifcode"] <- "ET_Classification"
names(final_data)[names(final_data) == "ET_Time"] <- "ET_Rec_Date"
names(final_data)[names(final_data) == "lc_leafType"] <- "ET_Leaf_type"
names(final_data)[names(final_data) == "lc_dominantSpecies"] <- "Layer_0_dominantSpecies"
names(final_data)[names(final_data) == "lc_dominantSpeciesPct"] <- "Layer_0_dominantSpeciesPct"
names(final_data)[names(final_data) == "lc_codominantSpecies1"] <- "Layer_0_codominantSpecies"
names(final_data)[names(final_data) == "lc_codominantSpecies1Pct"] <- "Layer_0_codominantSpeciesPct"
final_data$Layer_0_lifeform<-'trees'

# Remove 'lc_' and 'under_under' from column names
colnames(final_data) <- gsub("lc_", "", colnames(final_data))
colnames(final_data) <- gsub("under_under", "", colnames(final_data))

final_data$locationAccuracy<-NULL
final_data$photoBehind<-NULL
final_data$photoGround<-NULL
final_data$waterSeason<-NULL
final_data$comments<-NULL
final_data$change_type<-NULL
final_data$change_impact<-NULL
final_data$change_when<-NULL
final_data$change_pressures<-NULL
final_data$change_pressureUnknown<-NULL
final_data$change_comments<-NULL
final_data$change_legend<- NULL

colnames(final_data)
final_data$indicatorSpecies1<-NULL
length(unique(final_data$ET_ID))

# Keep layer 0 height and cover duplicated and seperate so they are not integrated in the next stage
final_data$ET_Canopy_Height<-final_data$Layer_0_height
final_data$ET_Canopy_Cover<-final_data$Layer_0_cover

min(na.omit(final_data$Layer_0_height)) #10

# Save
write.csv(final_data,"combined_data.csv",row.names = F)

final_data<- read.csv("combined_data.csv")
final_data_shp<- vect("combined_data.csv")
writeVector(final_data_shp, 'final_data.shp', filetype = "ESRI Shapefile", overwrite=TRUE)

head(final_data, n=3)
sum(!is.na(final_data$No_of_Dorm)) 

length(unique(final_data$Dor_Box[!is.na(final_data$No_of_Dormice)]))
# 78 boxes

length(unique(final_data$ET_ID))

#%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### STEP 10: PIVOT DATA ####

##%%%%%%%%%%%%%%%%%%%%%%%%%##

# Any column that begins with 'Layer_' should be put in long format, into a column with all the same suffixes. A new 'Layer column should be made.
# Ie: any value ending with _cover should be in  a new column called Cover
# any value from column ending in _height should be in  a new column called Height
# any value from column ending in _lifeform should be in  a new column called Lifeform
# any value from column ending in _dominantSpecies, _codominantSpecies, _codominantSpecies1, or _codominantSpecies2,  should be in  a new column called Species
# any value from column ending in _Pct should be in  a new column called Proportion of Species
# the new row generated should have all the same information in the other columns as it does originally.

final_data<- read.csv("combined_data.csv")

# Convert all columns starting with "Layer_" to character to ensure consistency
final_data <- final_data %>%
  mutate(across(starts_with("Layer_"), as.character))

# Pivot into long format
long_data <- final_data %>%
  pivot_longer(
    cols = starts_with("Layer_"), 
    names_to = c("Layer", "Attribute"),
    names_pattern = "Layer_(\\d+)_(.*)",
    values_to = "Value"
  ) %>%
  # Map attributes to specific categories
  mutate(
    Attribute = case_when(
      grepl("cover$", Attribute) ~ "Cover",
      grepl("height$", Attribute) ~ "Height",
      grepl("lifeform$", Attribute) ~ "Lifeform",
      grepl("dominantSpecies$", Attribute) ~ "Species",
      grepl("codominantSpecies\\d*$", Attribute) ~ "Species",
      grepl("Pct$", Attribute) ~ "Proportion_of_Species",
      TRUE ~ Attribute
    )
  ) %>%
  # Ensure rows with the same Layer and Attribute are properly separated for species and proportions
  group_by(across(-Value)) %>%
  ungroup() %>%
  pivot_wider(
    names_from = Attribute,
    values_from = Value
  )

# Expand rows so that every species-proportion pair gets its own row
expanded_data <- long_data %>%
  unnest(cols = c(Species, Proportion_of_Species)) %>%
  filter(!is.na(Species) | !is.na(Year_DM_Recorded)) # Remove rows where both are NA

#-----------#
# Tidy data #
#-----------#

colnames(expanded_data)
names(expanded_data)[names(expanded_data) == "indicatorSpecies1"] <- "Indicator_species"
names(expanded_data)[names(expanded_data) == "Layer"] <- "Forest_layer"
names(expanded_data)[names(expanded_data) == "Cover"] <- "ET_total_canopy_cover"
names(expanded_data)[names(expanded_data) == "Height"] <- "ET_layer_height"
names(expanded_data)[names(expanded_data) == "Lifeform"] <- "ET_layer_lifeform"

#-----------------------#
# Years since Coppicing #
#-----------------------#

# Add a column to sp_data for years since coppicing
expanded_data$Yr_since_Coppiced <- NA

# Loop through each row in expanded_data to calculate years since coppicing
for (i in 1:nrow(expanded_data)) {
  coppicing_year <- expanded_data$Coppicing_year[i]
  
  # Check if coppicing_year is NA
  if (is.na(coppicing_year)) {
    expanded_data$Yr_since_Coppiced[i] <- NA
  } else {
    # Calculate years since coppicing up to 2020
    if (coppicing_year > 2020) {
      expanded_data$Yr_since_Coppiced[i] <- 0
    } else {
      expanded_data$Yr_since_Coppiced[i] <- 2020 - coppicing_year
    }
  }
}

#------#
# Save #
#------#

# Unnest list columns in sp_data
sp_data <- expanded_data %>%
  unnest(cols = everything())

length(unique(sp_data$ET_ID))

# Save
write.csv(sp_data,"piv_data.csv", row.names = F)

sp_data<- read.csv('piv_data.csv')

sp_data$ET_Latitude<-as.numeric(sp_data$ET_Latitude)
sp_data$ET_Longitude<-as.numeric(sp_data$ET_Longitude)
sp_data_spat<-vect(sp_data, geom=c('ET_Latitude', 'ET_Longitude'))

# Save
writeVector(sp_data_spat, 'pivoted_dataset.shp', filetype = "ESRI Shapefile", overwrite=TRUE)

forest_layer_range <- sp_data %>%
  group_by(Coppicing_Plot) %>%
  summarise(Forest_layer_range = max(Forest_layer, na.rm = TRUE) - min(Forest_layer, na.rm = TRUE))

# View result
min(forest_layer_range$Forest_layer_range)
max(forest_layer_range$Forest_layer_range)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### STEP 11: CALCULATING CANOPY GAPS ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

sp_data<- read.csv('piv_data.csv')
summary(sp_data$Mean_Canopy_Height_2020)

#------#
# 2020 #
#------#

#-------------------------#
# Determine Gap threshold #
#-------------------------#

# Define thresholds and size
thresholds <- c(6, 8, 10, 12, 14) # Change thresholds as needed
size <- c(2, 20000) # Adjust size limits as required

# Select one plot (e.g., the first plot)
plot_index <- 1
current_id <- bontuchel$ID_No[plot_index]

# Crop CHM to the extent of the selected plot
bon_2020_CHM_crop <- crop(bon_CHM_2020, bontuchel[plot_index, ])

# Create an empty list to store gap layers for different thresholds
gap_layers <- list()

# Loop through thresholds and calculate gaps
for (threshold in thresholds) {
  
  # Get gaps for the selected threshold
  bon_gaps_2020 <- getForestGaps(chm_layer = bon_2020_CHM_crop, threshold = threshold, size = size)
  
  # Store the gap layer for later visualization
  gap_layers[[as.character(threshold)]] <- bon_gaps_2020
  
}

# Open a PNG device
pdf("~/Library/CloudStorage/OneDrive-UniversityofBristol/GitHub/Thesis/SuppInfo/Figures/Chapter4/GapAnalysis.pdf",width = 10, height = 10)

# Plot the original CHM and gaps for each threshold
par(mfrow = c(3, 2)) # Set up a 3x2 plotting grid

# Plot the original CHM
plot(bon_2020_CHM_crop, main = "Original Canopy Height Model", col = viridis::viridis(10))

# Plot gaps for each threshold
for (threshold in thresholds) {
  plot(bon_2020_CHM_crop, main = paste0("Gaps at Threshold ", threshold, "m"), col = viridis::viridis(10))
  plot(gap_layers[[as.character(threshold)]], col = "hotpink", add = TRUE, legend = FALSE)
}

# Empty spaces are where there is no data, not necessarily gaps.

# Turn off dev
dev.off()

#----------------#
# Calculate Gaps #
#----------------#

# Looking at gaps at 10m
threshold2020<- 10
size2020 <- c(2, 20000)

# Initialize gap proportion column
sp_data$Gap_Proportion_2020 <- NA

# Loop through each row in bontuchel
for (i in seq_len(nrow(bontuchel))) {
  
  # Extract plot ID
  current_id <- bontuchel$ID_No[i]
  
  # Crop CHM to the plot extent
  bon_2020_CHM_crop <- crop(bon_CHM_2020, bontuchel[i, ])
  
  # Identify gaps within the cropped area
  bon_gaps_2020 <- getForestGaps(chm_layer = bon_2020_CHM_crop, threshold = threshold2020, size = size2020)
  
  # Calculate total plot area (excluding NA values)
  plot_area <- sum(!is.na(values(bon_2020_CHM_crop))) * res(bon_2020_CHM_crop)[1] * res(bon_2020_CHM_crop)[2]
  
  # Calculate total gap area
  gap_area <- sum(values(bon_gaps_2020) > 0, na.rm = TRUE) * res(bon_gaps_2020)[1] * res(bon_gaps_2020)[2]
  
  # Calculate the proportion of gaps
  proportion_gaps_2020 <- (gap_area / plot_area)*100
  
  # Find the corresponding row in sp_data
  sp_data_index <- which(sp_data$Coppicing_Plot == current_id)
  
  # Update sp_data with the calculated proportions
  sp_data$Gap_Proportion_2020[sp_data_index] <- proportion_gaps_2020
  
  # Print progress
  print(paste0("Completed ", round((i / nrow(bontuchel)) * 100, 2), "%"))
}

length(unique(sp_data$ET_ID))

# Save data
write.csv(sp_data,"gap_data.csv", row.names = F)
sp_data<- read.csv('gap_data.csv')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### STEP 13: CLEAN UP SPECIES NAMES ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

# Load data
sp_data<-read.csv('gap_data.csv')

unique(sp_data$Species)

# Taxonomic groups levels vary, tried for species where possible, but some functional groups remain

# Replace the specific string 
sp_data$Species[sp_data$Species == "Alder (species unknown) (Alnus)"] <- 'Alder'
sp_data$Species[sp_data$Species == "Alder (species unknown)"] <- 'Alder'
sp_data$Species[sp_data$Species == "Ash (species unknown) (Fraxinus)"] <- 'Ash'
sp_data$Species[sp_data$Species == "Ash (species unknown)"] <- 'Ash'
sp_data$Species[sp_data$Species == "Birch (species unknown)"] <- 'Birch'
sp_data$Species[sp_data$Species == "Birches"] <- 'Birch'
sp_data$Species[sp_data$Species == "(Rosaceae)"] <- 'Rosaceae'
sp_data$Species[sp_data$Species == "Oak (Quercus)"] <- 'Oak'
sp_data$Species[sp_data$Species == "Holly (species unknown)"] <- 'Holly'
sp_data$Species[sp_data$Species == "Hazel (species unknown)"] <- 'Hazel'
sp_data$Species[sp_data$Species == " (Ulmus)"] <- 'Elm'
sp_data$Species[sp_data$Species == "Oak (species unknown)"] <- 'Oak'
sp_data$Species[sp_data$Species == "Clovers (species unknown)"] <- 'Clover'
sp_data$Species[sp_data$Species == "Buttercup (species unknown)"] <- 'Buttercup'
sp_data$Species[sp_data$Species == "Meadow grass (species unknown)"] <- 'Grass'
sp_data$Species[sp_data$Species == "Beech (species unknown)"] <- 'Beech'
sp_data$Species[sp_data$Species == "Beech (species unknown) (Fagus)"] <- 'Beech'
sp_data$Species[sp_data$Species == "Oaks"] <- 'Oak'
sp_data$Species[sp_data$Species == " (Rosa canina)"] <- 'Rosa canina'
sp_data$Species[sp_data$Species == "Firs"] <- 'Fir'
sp_data$Species[sp_data$Species == " (Abies)"] <- 'Abies'
sp_data$Species[sp_data$Species == "Elms"] <- 'Elm'
sp_data$Species[sp_data$Species == "Fir (species unknown)"] <- 'Fir'
sp_data$Species[sp_data$Species == " (Ranunculaceae)"] <- 'Ranunculaceae'
sp_data$Species[sp_data$Species == "Larch (species unknown) (Larix)"] <- 'Larch'
sp_data$Species[sp_data$Species == " (Ranunculeae)"] <- 'Ranunculaceae'
sp_data$Species[sp_data$Species == " (Rosaceae)"] <- 'Rosaceae'
sp_data$Species[sp_data$Species == "Geranium (species unknown)"] <- 'Geranium'
sp_data$Species[sp_data$Species == "Willow (species unknown)"] <- 'Willow'
sp_data$Species[sp_data$Species == "Elm (species unknown)"] <- 'Elm'
sp_data$Species[sp_data$Species == "Birch (species unknown) (Betula)"] <- 'Birch'
sp_data$Species[sp_data$Species == "Spruces"] <- 'Spruce'
sp_data$Species[sp_data$Species == "Alders"] <- 'Alder'
sp_data$Species[sp_data$Species == "Ashes"] <- 'Ash'
sp_data$Species[sp_data$Species == "Hollies"] <- 'Holly'
sp_data$Species[sp_data$Species == "Elders"] <- 'Elder'
sp_data$Species[sp_data$Species == "Willows"] <- 'Willow'
sp_data$Species[sp_data$Species == "Rose species"] <- 'Rose'
sp_data$Species[sp_data$Species == "Rosa canina"] <- 'Rose'
sp_data$Species[sp_data$Species == "Hawthorns"] <- 'Hawthorn'
sp_data$Species[sp_data$Species == "Willowherb species"] <- 'Willowherb'
sp_data$Species[sp_data$Species == "Rosebay Willowherb"] <- 'Willowherb'
sp_data$Species[sp_data$Species == "Beeches"] <- 'Beech'
sp_data$Species[sp_data$Species == "Common Ivy"] <- 'Ivy'
sp_data$Species[sp_data$Species == "Common Honeysuckle"] <- 'Honeysuckle'
sp_data$Species[sp_data$Species == "Abies"] <- 'Fir'
sp_data$Species[sp_data$Species == "Common Nettle"] <- 'Nettle'
sp_data$Species[sp_data$Species == "Common Sorrel"] <- 'Sorrel'
sp_data$Species[sp_data$Species == "Western hemlock-spruce"] <- 'Hemlock'
sp_data$Species[sp_data$Species == " (Plantaginaceae)"] <- 'Plantaginaceae'
sp_data$Species[sp_data$Species == " (Saxifragaceae)"] <- 'Saxifragaceae'
sp_data$Species[sp_data$Species == "Willow (Salix)"] <- 'Willow'
sp_data$Species[sp_data$Species == "Hazels"] <- 'Hazel'
sp_data$Species[sp_data$Species == "Beech "] <- 'Beech'
sp_data$Species[sp_data$Species == "Chestnuts"] <- 'Chestnut'
sp_data$Species[sp_data$Species == "Sedge (species unknown)"] <- 'Sedge'
sp_data$Species[sp_data$Species == "Hawthorn (species unknown)"] <- 'Hawthorn'
sp_data$Species[sp_data$Species == "Hazel (species unknown) (Corylus)"] <- 'Hazel'
sp_data$Species[sp_data$Species == "Elder (species unknown)"] <- 'Elder'
sp_data$Species[sp_data$Species == "Hornbeams"] <- 'Hornbeam'
sp_data$Species[sp_data$Species == " (Acer pseudoplatanus)"] <- 'Sycamore'
sp_data <- sp_data %>%
  filter(Species != "ND")

unique(sp_data$Species)
length(unique(sp_data$Species))

# Create a named vector: names are common names, values are scientific names
species_map <- c(
  "Flowering Plant" = "Angiosperm",
  "Alder" = "Alnus",
  "Ash" = "Fraxinus",
  "Hazel" = "Corylus avellana",
  "Willow" = "Salix",
  "Bramble" = "Rubus fruticosus",
  "Grass" = "Poaceae",
  "Fern" = "Pteridophyta",
  "Moss" = "Bryophyta",
  "Birch" = "Betula",
  "Needle-leaf" = 'Needle-leaf',
  "Hawthorn" = "Crataegus",
  "Beech" = "Fagus sylvatica",
  "Elder" = "Sambucus nigra",
  "Ivy" = "Hedera helix",
  "Honeysuckle" = "Lonicera periclymenum",
  "Geranium" = "Geranium robertianum",
  "Oak" = "Quercus",
  "Sycamore" = "Acer pseudoplatanus",
  "Berries (species unknown)" = 'Berries (species unknown)',
  "Holly" = "Ilex aquifolium",
  "Elm" = "Ulmus",
  "Clover" = "Trifolium",
  "Buttercup" = "Ranunculus",
  "Rosaceae" = "Rosaceae",
  "Willowherb" = "Epilobium",
  "Wood Avens" = "Geum urbanum",
  "Ranunculaceae" = "Ranunculaceae",
  "Rowan" = "Sorbus aucuparia",
  "Sedge" = "Carex",
  "Asteraceae" = "Asteraceae",
  "Rose" = "Rosa",
  "Blackthorn" = "Prunus spinosa",
  "Dog's Mercury" = "Mercurialis perennis",
  "Nettle" = "Urtica dioica",
  "Caryophyllaceae" = "Caryophyllaceae",
  "Larch" = "Larix",
  "Selfheal" = "Prunella vulgaris",
  "Spruce" = "Picea",
  "Fir" = "Abies",
  "Hemlock" = "Tsuga",
  "Sorrel" = "Rumex acetosa",
  "Bilberry" = "Vaccinium myrtillus",
  "Enchanter's-nightshade" = "Circaea lutetiana",
  "Saxifragaceae" = "Saxifragaceae",
  "Plantaginaceae" = "Plantaginaceae",
  "Horse chestnut" = "Aesculus hippocastanum",
  "Broad-leaf" = "Broad-leaf (unspecified)",
  "Reed" = "Phragmites australis",
  "Norway spruce" = "Picea abies",
  "Common bluebell" = "Hyacinthoides non-scripta",
  "Chestnut" = "Castanea sativa",
  "Hornbeam" = "Carpinus betulus",
  "Cow Parsley" = "Anthriscus sylvestris",
  "Wild cherry" = "Prunus avium"
)

# Map the species names to scientific names
sp_data$Scientific_Name <- species_map[sp_data$Species]

unique(sp_data$Scientific_Name)

# Find species that are not in the species_map
missing_species <- setdiff(sp_data$Species, names(species_map))
missing_species


# Named vector: common name → taxonomic level
taxonomic_level <- c(
  "Flowering Plant" = "Phylum",
  "Alder" = "Genus",
  "Ash" = "Genus",
  "Hazel" = "Species",
  "Willow" = "Genus",
  "Bramble" = "Species",
  "Grass" = "Group",
  "Fern" = "Group",
  "Moss" = "Group",
  "Birch" = "Genus",
  "Needle-leaf" = "Group",
  "Hawthorn" = "Genus",
  "Beech" = "Species",
  "Elder" = "Species",
  "Ivy" = "Species",
  "Honeysuckle" = "Species",
  "Geranium" = "Genus",
  "Oak" = "Genus",
  "Sycamore" = "Species",
  "Berries (species unknown)" = "Group",
  "Holly" = "Species",
  "Elm" = "Genus",
  "Clover" = "Genus",
  "Buttercup" = "Genus",
  "Rosaceae" = "Family",
  "Willowherb" = "Genus",
  "Wood Avens" = "Species",
  "Ranunculaceae" = "Family",
  "Rowan" = "Species",
  "Sedge" = "Genus",
  "Asteraceae" = "Family",
  "Rose" = "Genus",
  "Blackthorn" = "Species",
  "Dog's Mercury" = "Species",
  "Nettle" = "Species",
  "Caryophyllaceae" = "Family",
  "Larch" = "Genus",
  "Selfheal" = "Species",
  "Spruce" = "Genus",
  "Fir" = "Genus",
  "Hemlock" = "Genus",
  "Sorrel" = "Species",
  "Bilberry" = "Species",
  "Enchanter's-nightshade" = "Species",
  "Saxifragaceae" = "Family",
  "Plantaginaceae" = "Family",
  "Horse chestnut" = "Species",
  "Broad-leaf" = "Group",
  "Reed" = "Species",
  "Norway spruce" = "Species",
  "Common bluebell" = "Species",
  "Chestnut" = "Species",
  "Hornbeam" = "Species",
  "Cow Parsley" = "Species",
  "Wild cherry" = "Species"
)

# Add a new column based on the mapping
sp_data$Taxonomic_Level <- taxonomic_level[sp_data$Species]
unique(sp_data$Taxonomic_Level)

# Save data
write.csv(sp_data,"gap_data_sp.csv", row.names = F)

table(sp_data$Taxonomic_Level)
sum(is.na(sp_data$Scientific_Name))
length(unique(sp_data$Species))
table(sp_data$Species, useNA = "ifany")
length(unique(sp_data$ET_ID))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### STEP 12: DERIVE MEASURE OF BIODIVERSITY- ET ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

sp_data<- read.csv('gap_data_sp.csv')
head(sp_data)
colnames(sp_data)

# Convert Proportion_of_Species from percentage to decimal
sp_data <- sp_data %>%
  mutate(Proportion_of_Species = as.numeric(Proportion_of_Species) / 100)

# Summing proportions for each species within each coppicing plot
species_agg <- sp_data %>%
  group_by(Coppicing_Plot, Species, ET_ID) %>%
  summarise(Total_Proportion = sum(Proportion_of_Species, na.rm = TRUE)) %>%
  ungroup()

# Calculate relative abundance per coppicing plot
species_agg <- species_agg %>%
  group_by(Coppicing_Plot) %>%
  mutate(Relative_Abundance = Total_Proportion / sum(Total_Proportion, na.rm = TRUE)) %>%
  ungroup()

# Calculate Shannon Index for each coppicing plot and keep ET_ID
shannon_index <- species_agg %>%
  group_by(Coppicing_Plot, ET_ID) %>%
  summarise(Shannon_Index = -sum(Relative_Abundance * log(Relative_Abundance), na.rm = TRUE)) %>%
  ungroup()

# Save
write.csv(shannon_index, "shannon_index_data.csv", row.names = FALSE)

# Merge the Shannon Index into the original sp_data dataframe
sp_data <- sp_data %>%
  left_join(shannon_index, by = c("Coppicing_Plot", "ET_ID"))

length(unique(sp_data$ET_ID))
# for plots that have more than one unique ET ID,  i will need to work out the shannon index of both and then get a mean.

# Calculate mean Shannon Index per Coppicing_Plot
mean_shannon <- shannon_index %>%
  group_by(Coppicing_Plot) %>%
  summarise(Mean_Shannon_Index = mean(Shannon_Index, na.rm = TRUE))

# Join mean Shannon back into sp_data
sp_data <- sp_data %>%
  left_join(mean_shannon, by = "Coppicing_Plot")

length(unique(sp_data$ET_ID))

write.csv(sp_data,"final_dataset_1.csv",row.names = F)

##%%%%%%%%%%%%%%%%%%%##

#### Data Analysis ####

##%%%%%%%%%%%%%%%%%%%##

# Set working directory
setwd("~/Library/CloudStorage/OneDrive-UniversityofBristol/Documents/PhD/Project 3_Bontuchel/Data")

select_best_model <- function(data, predictors, response_var,
                              cor_threshold = 0.7,
                              zero_threshold = 0.3,
                              shape = c("linear", "quadratic"),
                              interactions = FALSE) {
  library(MASS); library(pscl); library(car); library(MuMIn)
  shape <- match.arg(shape)
  response <- data[[response_var]]
  is_count <- all(response %% 1 == 0) && all(response >= 0)
  zero_prop <- mean(response == 0)
  
  # 1) Find uncorrelated predictor combos
  if (length(predictors) == 1) {
    valid_combinations <- list(predictors)
  } else {
    cor_matrix <- cor(data[, predictors], use = "pairwise.complete.obs", method = "spearman")
    valid_combinations <- list(); counter <- 1
    for (i in length(predictors):1) {
      combs <- combn(predictors, i, simplify = FALSE)
      for (combo in combs) {
        sub_cor <- cor_matrix[combo, combo]
        if (all(abs(sub_cor[upper.tri(sub_cor)]) < cor_threshold)) {
          valid_combinations[[counter]] <- combo; counter <- counter + 1
        }
      }
      if (length(valid_combinations) > 0) break
    }
    if (!length(valid_combinations)) stop("No uncorrelated predictor combinations found.")
  }
  
  model_results <- list(); fallback_candidates <- list()
  
  for (combo in valid_combinations) {
    linear_terms <- paste(combo, collapse = " + ")
    quadratic_terms <- paste0("I(", combo, "^2)", collapse = " + ")
    if (shape == "quadratic") {
      formula_rhs <- paste(linear_terms, quadratic_terms, sep = " + ")
    } else {
      formula_rhs <- linear_terms
    }
    if (interactions && length(combo) > 1) {
      formula_rhs <- paste0("(", formula_rhs, ")^2")
    }
    formula <- as.formula(paste(response_var, "~", formula_rhs))
    
    # === COUNT MODEL LOGIC (not included here) ===
    if (is_count) {
      # Placeholder: Poisson, NB, ZIP, ZINB etc.
      next
    }
    
    # === LM ===
    lm_model <- tryCatch(lm(formula, data = data), error = function(e) NULL)
    if (!is.null(lm_model)) {
      vif_ok <- tryCatch(all(vif(lm_model) < 5), error = function(e) TRUE)
      shapiro_p <- tryCatch(shapiro.test(residuals(lm_model))$p.value, error = function(e) NA)
      if (vif_ok && !is.na(shapiro_p) && shapiro_p > 0.05) {
        model_results[[length(model_results) + 1]] <- list(
          model = lm_model, type = "Linear Model", aic = AIC(lm_model),
          shapiro_p = shapiro_p, formula = formula,
          reason = "Normal residuals & low collinearity"
        )
      } else {
        # Log-transformed fallback
        response_log <- log(response + 1)
        data_log <- data; data_log[[response_var]] <- response_log
        formula_log <- update(formula, paste(response_var, "~ ."))
        lm_log <- tryCatch(lm(formula_log, data = data_log), error = function(e) NULL)
        if (!is.null(lm_log)) {
          sh_p <- tryCatch(shapiro.test(residuals(lm_log))$p.value, error = function(e) NA)
          fallback_candidates[[length(fallback_candidates) + 1]] <- list(
            model = lm_log, type = "Log-linear (violations)", aic = AIC(lm_log),
            shapiro_p = sh_p, formula = formula_log,
            reason = "Corrected non-normal residuals via log-transform"
          )
        }
      }
    }
    
    # === Gaussian GLM ===
    gauss_glm <- tryCatch(glm(formula, data = data, family = gaussian), error = function(e) NULL)
    if (!is.null(gauss_glm)) {
      residuals_glm <- residuals(gauss_glm, type = "deviance")
      shapiro_p <- tryCatch(shapiro.test(residuals_glm)$p.value, error = function(e) NA)
      if (!is.na(shapiro_p) && shapiro_p > 0.05) {
        model_results[[length(model_results) + 1]] <- list(
          model = gauss_glm, type = "Gaussian GLM", aic = AIC(gauss_glm),
          shapiro_p = shapiro_p, formula = formula,
          reason = "Gaussian GLM with normal residuals"
        )
      } else {
        fallback_candidates[[length(fallback_candidates) + 1]] <- list(
          model = gauss_glm, type = "Gaussian GLM (violations)", aic = AIC(gauss_glm),
          shapiro_p = shapiro_p, formula = formula,
          reason = "Gaussian GLM with some residual non-normality"
        )
      }
    }
    
    # === Gamma GLM (positive-only) ===
    if (all(response > 0)) {
      gamma_model <- tryCatch(glm(formula, data = data, family = Gamma(link = "log")), error = function(e) NULL)
      if (!is.null(gamma_model)) {
        pearson_dispersion <- sum(residuals(gamma_model, type = "pearson")^2) / df.residual(gamma_model)
        if (pearson_dispersion < 2) {
          model_results[[length(model_results) + 1]] <- list(
            model = gamma_model, type = "Gamma GLM", aic = AIC(gamma_model),
            formula = formula,
            reason = paste("Gamma GLM (dispersion =", round(pearson_dispersion, 2), ")")
          )
        } else {
          fallback_candidates[[length(fallback_candidates) + 1]] <- list(
            model = gamma_model, type = "Gamma GLM (dispersion)", aic = AIC(gamma_model),
            formula = formula,
            reason = paste("High dispersion (", round(pearson_dispersion, 2), ") in Gamma GLM")
          )
        }
      }
    }
  }
  
  # === Final selection ===
  if (length(model_results) > 0) {
    best <- model_results[[which.min(sapply(model_results, function(x) x$aic))]]
    cat("\n✅ Best:", best$type,
        "\nFormula:", deparse(best$formula),
        if (!is.null(best$shapiro_p)) paste0("\nShapiro p = ", round(best$shapiro_p, 3)),
        "\nAIC =", round(best$aic, 2),
        "\nReason:", best$reason, "\n\n")
    print(summary(best$model))
    return(best$model)
    
  } else if (length(fallback_candidates) > 0) {
    sorted_fb <- fallback_candidates[order(sapply(fallback_candidates, function(x) x$aic))]
    glm_fb <- Filter(function(x) grepl("GLM", x$type), sorted_fb)
    fb <- if (length(glm_fb)) glm_fb[[1]] else sorted_fb[[1]]
    cat("\n⚠️  Fallback:", fb$type,
        "\nFormula:", deparse(fb$formula),
        if (!is.null(fb$shapiro_p)) paste0("\nShapiro p = ", round(fb$shapiro_p, 3)),
        "\nAIC =", round(fb$aic, 2),
        "\nReason:", fb$reason, "\n\n")
    print(summary(fb$model))
    return(fb$model)
    
  } else {
    stop("❌ No valid models could be fit.")
  }
}

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### STEP 13: CORRELATION OF BIODIVERSITY INDICATORS ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

# Load data
sp_data<-read.csv('final_dataset_1.csv')

# Ensure numeric conversion
sp_data <- sp_data %>%
  mutate(
    No_of_Dormice = as.numeric(No_of_Dormice),
    Shannon_Index = as.numeric(Shannon_Index)
  )

# Remove duplicate entries based on Coppicing_Plot, Year_DM_Recorded, and Dor_Box
sp_data_unique <- sp_data %>%
  distinct(Coppicing_Plot, Year_DM_Recorded, Dor_Box, ET_ID, .keep_all = TRUE)

total_dormice_df <- sp_data_unique %>%
  group_by(Coppicing_Plot) %>%
  summarise(
    Num_Boxes = n_distinct(Dor_Box),  # Unique nest boxes per plot
    Total_Dormice = sum(No_of_Dormice, na.rm = TRUE),  # Total dormice count per plot
    Mean_Dormice = Total_Dormice / Num_Boxes,  # Average dormice per box
    Mean_Dormice = round(Mean_Dormice), #round to int
    Avg_Shannon_Index = mean(Shannon_Index, na.rm = TRUE)
  ) %>%
  ungroup()

# Print final dataset
print(total_dormice_df)

# Add summary data to original sp_data without changing original rows/columns
sp_data <- sp_data %>%
  left_join(total_dormice_df, by = "Coppicing_Plot")

write.csv(sp_data,"final_dataset.csv",row.names = F)

write.csv(total_dormice_df,"total_dormice_df.csv",row.names = F)

# Look at distribution of data
hist(total_dormice_df$Mean_Dormice) # skewed

# Fit a model
fit <- lm(log(Avg_Shannon_Index) ~ Mean_Dormice, data = total_dormice_df)
summary(fit)
shapiro.test(fit$residuals)

# Not normal even with log, skip modelling and do a corr
length(total_dormice_df$Mean_Dormice)

# Correlation
cor.test(total_dormice_df$Mean_Dormice, 
         total_dormice_df$Avg_Shannon_Index, 
         method = "spearman")

# Spearman's rank correlation rho
# data:  total_dormice_df$Mean_Dormice and total_dormice_df$Avg_Shannon_Index
# S = 36011, p-value = 0.0871
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# -0.2307378 

#------#
# Plot #
#------#

# Open a PNG device
pdf("~/Library/CloudStorage/OneDrive-UniversityofBristol/GitHub/Thesis/Chapter4/Figures/DormiceIndicators.pdf", width = 5, height = 5)

ggplot(total_dormice_df, aes(x = Mean_Dormice, y = Avg_Shannon_Index)) +
  geom_point(color = "pink", size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "hotpink", fill = "red", alpha = 0.2) +
  labs(
    x = "Mean Number of Dormice per Panel",
    y = "Biodiversity (Shannon Index)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.line.x = element_line(color = "black", size = 0.25),
    axis.line.y = element_line(color = "black", size = 0.25),
    axis.ticks = element_line(color = "black"),
    panel.border = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    axis.title.x = element_text(margin = ggplot2::margin(15, 0, 0, 0)),
    axis.title.y = element_text(margin = ggplot2::margin(0, 15, 0, 0)),
    axis.text = element_text(size = 12)
  )


# Close the device
dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### STEP 15: GROUND TRUTH LiDAR ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

# Load data
sp_data<-read.csv('final_dataset.csv')
ET_points<-vect('ET_Points.gpkg')
bon_CHM_2020<- rast("bon_CHM_2020_rep.tif")

ET_points<-crop(ET_points, bon_CHM_2020)

# Create a buffer of 15 meters around the ET points
buffer_15m <- buffer(ET_points, width = 15)  # Buffer by 15 meters
buffer_15m$ET_Rec_Date<-buffer_15m$`_CREATION_DATE`

# Check if the extents overlap
ext(bon_CHM_2020)
st_bbox(buffer_15m)

plot(bon_CHM_2020)
plot(bontuchel, border='pink', add=TRUE)
plot(buffer_15m, add=TRUE, col='hotpink')
plot(ET_points, add=TRUE)


#--------------------#
# Mean Canopy Height #
#--------------------#

# For each point in buffer_15m, crop CHM_2020 by it and find the mean height of canopy. Add this to sp_data by ET_Rec_Date.
# Initialize a vector to store mean canopy heights
mean_canopy_heights <- numeric(length(buffer_15m))

# Loop through each point in buffer_15m
for (i in 1:length(buffer_15m)) {
  # Crop CHM_2020 by the buffer around the current point
  cropped_CHM <- crop(bon_CHM_2020, buffer_15m[i])
  
  # Calculate the mean height of the canopy within the cropped area
  mean_canopy_heights[i] <- mean(values(cropped_CHM), na.rm = TRUE)
}

# Normalize both to POSIXct, consistent format and timezone (UTC)
sp_data$ET_Rec_Date <- as.POSIXct(sp_data$ET_Rec_Date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Check if the conversion worked
head(sp_data$ET_Rec_Date)

# For buffer_15m (spatial object with an attribute ET_Rec_Date)
# Extract the 'ET_Rec_Date' attribute as a vector and convert it to POSIXct
buffer_15m$ET_Rec_Date <- as.POSIXct(buffer_15m$ET_Rec_Date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Check if conversion worked
head(buffer_15m$ET_Rec_Date)

# Match up the data
sp_data$LiDAR_Mean_Canopy_Height_2020 <- mean_canopy_heights[match(sp_data$ET_Rec_Date, buffer_15m$ET_Rec_Date)]

summary(as.numeric(sp_data$ET_Canopy_Height))
summary(sp_data$LiDAR_Mean_Canopy_Height_2020)

# Select only the relevant columns for uniqueness check - prevents pseudoreplication 
unique_rows <- !duplicated(sp_data[, c("ET_ID")])

# Keep only unique rows in the dataframe
sp_data_unique <- sp_data[unique_rows, ]

length(sp_data_unique$ET_ID)

#///////////////////#
# Check assumptions #
#///////////////////#

top_layer <- sp_data_unique %>%
  filter(Forest_layer == 0)

min(na.omit(top_layer$ET_Canopy_Height)) #10
min(na.omit(top_layer$LiDAR_Mean_Canopy_Height_2020)) #4.5

hist(top_layer$ET_Canopy_Height)
hist(top_layer$LiDAR_Mean_Canopy_Height_2020)

# Calculate the residuals (difference between the two height measurements)
residuals <- top_layer$ET_Canopy_Height - top_layer$LiDAR_Mean_Canopy_Height_2020
hist(residuals, main = "Residuals (Ground - LiDAR)", xlab = "Residuals", col = "pink")
summary(residuals)
shapiro.test(residuals) 

# Assumptions not met

#/////////////#
# Correlation #
#/////////////#

# Spearmans correlation coefficient
height<-cor.test(top_layer$ET_Canopy_Height, top_layer$LiDAR_Mean_Canopy_Height_2020, method='spearman', use = "complete.obs")

# Spearman's rank correlation rho
# data:  top_layer$ET_Canopy_Height and top_layer$LiDAR_Mean_Canopy_Height_2020
# S = 23913, p-value = 0.001365
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.3978188 

# Assuming top_layer$ET_Canopy_Height and top_layer$LiDAR_Mean_Canopy_Height_2020 exist and have matching rows

# Calculate difference (ground - LiDAR)
diff <- top_layer$ET_Canopy_Height - top_layer$LiDAR_Mean_Canopy_Height_2020

# Calculate average difference, excluding NAs
mean_diff <- mean(diff, na.rm = TRUE)

# Calculate standard deviation of the difference (optional)
sd_diff <- sd(diff, na.rm = TRUE)

mean_diff  # this is the average difference you want
sd_diff    # standard deviation for variability

#--------------------#
# Mean Canopy Cover #
#--------------------#

# When working out canopy cover from LiDAR, we will use the threshold of cover from the min height of the top layer, ie: 10
# We use the min from top layer of et height because lidar will go through all layers of canopy so we can't differentiate.
# Get LiDAR canopy cover
# For each point in buffer_15m, crop CHM_2020 by it and find the mean canopy cover. Add this to sp_data by ET_Rec_Date.

# Initialize a vector to store mean canopy heights
mean_canopy_cover <- numeric(length(buffer_15m))

# Loop through each point in buffer_15m
for (i in 1:length(buffer_15m)) {
  
  # Crop CHM_2020 by the buffer around the current point
  cropped_CHM <- crop(bon_CHM_2020, buffer_15m[i, ])
  Perc_under_8m_2020<-100*(length(which(as.vector(na.omit(as.vector(cropped_CHM)))<=10))/length(as.vector(na.omit(as.vector(cropped_CHM))))) 
  mean_canopy_cover[i]<-100- Perc_under_8m_2020
  
}

# Add the mean canopy heights to sp_data by ET_Rec_Date
sp_data$LiDAR_Mean_Canopy_cover_2020 <- mean_canopy_cover[match(sp_data$ET_Rec_Date, buffer_15m$ET_Rec_Date)]

summary(sp_data$ET_Canopy_Cover)
summary(sp_data$LiDAR_Mean_Canopy_cover_2020)

# Select only the relevant columns for uniqueness check - prevents pseudoreplication 
unique_rows <- !duplicated(sp_data[, c("ET_ID")])

# Keep only unique rows in the dataframe
sp_data_unique <- sp_data[unique_rows, ]

#///////////////////#
# Check assumptions #
#///////////////////#

plot(sp_data_unique$ET_Canopy_Cover, sp_data_unique$LiDAR_Mean_Canopy_cover_2020)
hist(sp_data_unique$ET_Canopy_Cover)
hist(sp_data_unique$LiDAR_Mean_Canopy_cover_2020)

# Calculate the residuals (difference between the two height measurements)
residuals <- sp_data_unique$ET_Canopy_Cover - sp_data_unique$LiDAR_Mean_Canopy_cover_2020
hist(residuals, main = "Residuals (Ground - LiDAR)", xlab = "Residuals", col = "pink")

# Summary of residuals
summary(residuals)
shapiro.test(residuals) 

# Assumption of normality not met

#/////////////#
# Correlation #
#/////////////#

cor.test(sp_data_unique$ET_Canopy_Cover, sp_data_unique$LiDAR_Mean_Canopy_cover_2020, method="spearman", use = "complete.obs")

# Spearman's rank correlation rho
# data:  sp_data_unique$ET_Canopy_Cover and sp_data_unique$LiDAR_Mean_Canopy_cover_2020
# S = 28482, p-value = 0.02596
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.2827615 

# Calculate difference (ground - LiDAR)
diff <- sp_data_unique$ET_Canopy_Cover - sp_data_unique$LiDAR_Mean_Canopy_cover_2020

# Calculate average difference, excluding NAs
mean_diff <- mean(diff, na.rm = TRUE)

# Calculate standard deviation of the difference (optional)
sd_diff <- sd(diff, na.rm = TRUE)

mean_diff  # this is the average difference you want
sd_diff    # standard deviation for variability

#------#
# Plot #
#------#

length(sp_data_unique$ET_Canopy_Cover)
length(sp_data_unique$LiDAR_Mean_Canopy_cover_2020)

# Open a PNG device
pdf("~/Library/CloudStorage/OneDrive-UniversityofBristol/GitHub/Thesis/Chapter4/Figures/LidarVsGround.pdf", width = 4.5, height = 7)

par(mfrow = c(2, 1),         # 2 rows, 1 column
    mar = c(4, 4, 1, 1),     # inner margins (bottom, left, top, right)
    oma = c(0, 0, 0, 0))     # outer margins (bottom, left, top, right)

# First scatter plot
plot(sp_data_unique$LiDAR_Mean_Canopy_Height_2020, sp_data_unique$ET_Canopy_Height,
     xlab = "LiDAR Measured Height", 
     ylab = "Ground Measured Height", 
     xlim = c(0, max(sp_data_unique$ET_Canopy_Height, sp_data_unique$LiDAR_Mean_Canopy_Height_2020)),
     ylim = c(0, max(sp_data_unique$ET_Canopy_Height, sp_data_unique$LiDAR_Mean_Canopy_Height_2020)),
     bty = "l")  # Use 'l' to show only left and bottom box lines (axes)
abline(a = 0, b = 1, col = "hotpink")

# Second scatter plot
plot(sp_data_unique$LiDAR_Mean_Canopy_cover_2020, sp_data_unique$ET_Canopy_Cover,
     xlab = "LiDAR Measured Canopy Cover", 
     ylab = "Ground Measured Canopy Cover", 
     xlim = c(0, max(sp_data_unique$ET_Canopy_Cover, sp_data_unique$LiDAR_Mean_Canopy_cover_2020)),
     ylim = c(0, max(sp_data_unique$ET_Canopy_Cover, sp_data_unique$LiDAR_Mean_Canopy_cover_2020)),
     bty = "l")  # Again, only left and bottom box lines
abline(a = 0, b = 1, col = "hotpink")

# Close PNG
dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### STEP 16: RS METRICS AND BIODIVERSITY ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#GAMS WITH NB ARE JUST GLMS IF NO SMOOTHING REMEBER!

# Load data
sp_data<-read.csv('final_dataset.csv')

#-------------#
# Correlation #
#-------------#

# Determine which metrics may be correlated
colnames(sp_data)
length(unique(sp_data$Coppicing_Plot))

# Select only the relevant columns for uniqueness check - prevents pseudoreplication 
unique_rows <- !duplicated(sp_data[, c("ET_ID")])

# Keep only unique rows in the dataframe
sp_data_unique <- sp_data[unique_rows, ]
colnames(sp_data_unique)

custom_labels <- c(
  Mean_Canopy_Height_2020 = "Mean Canopy Height",
  Perc_under_2m_2020 = "% Vegetation Under 2m",
  Can_cover_2020 = "Canopy Cover",
  Height_cv_2020 = "Canopy Height Variation",
  Max_Canopy_Height_2020 = "Max Canopy Height",
  Gap_Proportion_2020 = "Canopy Gap Proportion",
  Intensity = "LiDAR Intensity",
  VCI = "VCI",
  LAI = "LAI",
  Solar_Radiation = "Solar Radiation",
  FT_5_10m = "First returns (5–10m)",
  FT_1.37m = "First returns  (1.37m)",
  FT_10_20m = "First returns  (10–20m)",
  HALP = "HALP",
  Aspect_Cos = "Aspect (Cosine)",
  TPI = "TPI",
  Roughness = "Roughness",
  TRI = "TRI",
  Elevation = "Elevation",
  Slope = "Slope",
  Mean_Shannon_Index = "Mean Shannon Index",
  Yr_since_Coppiced = "Recovery Period",
  Mean_Curve = "Mean Curvature",
  FlowDir = "Flow Direction",
  Plane_Curve = "Plane Curvature",
  Profile_Curve = "Profile Curvature",
  Canopy_Volume_2020 = "Canopy Volume",
  Mean_Dormice = "Mean Dormice Abundance"
)


cor_matrix <- cor(sp_data_unique[, c("Mean_Canopy_Height_2020", "Perc_under_2m_2020", "Can_cover_2020", 
                                     "Height_cv_2020", "Max_Canopy_Height_2020", "Gap_Proportion_2020", 
                                     "Intensity", "VCI", "LAI", "Solar_Radiation", "FT_5_10m", "FT_1.37m", 
                                     "FT_10_20m", "HALP", "Aspect_Cos", "TPI", "Roughness", "TRI", 
                                     "Elevation", "Slope", "Mean_Shannon_Index", "Yr_since_Coppiced", "Mean_Curve", 
                                     "FlowDir", "Plane_Curve", "Profile_Curve", "Canopy_Volume_2020", "Mean_Dormice")], 
                  use = "pairwise.complete.obs", method='spearman')  

# Apply custom labels
colnames(cor_matrix) <- custom_labels[colnames(cor_matrix)]
rownames(cor_matrix) <- custom_labels[rownames(cor_matrix)]

cor_matrix

#------#
# Plot #
#------#

# Open a pdf device
pdf("~/Library/CloudStorage/OneDrive-UniversityofBristol/Github/Thesis/SuppInfo/Figures/Chapter4/CorrPlot.pdf", width = 10, height = 10)  # in inches, not pixels

# Define the custom color gradient
custom_colors <- colorRampPalette(c("#1A237E", "#C5CAE9", "white", "pink", "hotpink"))

# Set plot parameters
par(mfrow = c(1, 1), mar = c(6, 8, 2, 2), las = 1, xpd = TRUE, cex.axis = 1.4, ps = 10)

# Create the correlation plot
corrplot(cor_matrix, 
         type = 'lower', 
         tl.col = "black", 
         tl.cex = 1, # Make text smaller
         tl.srt = 70,   # Rotate text for better readability
         col = custom_colors(200), 
         cl.cex = 1.5)

# Close the device
dev.off()

#//////////////////////////#

##### Canopy structure #####

#//////////////////////////#

#-------#
# Model #
#-------#

# Load data
sp_data<-read.csv('final_dataset.csv')

# Select only the relevant columns for uniqueness check
unique_rows <- !duplicated(sp_data[, c("Coppicing_Plot")])
sp_data_unique <- sp_data[unique_rows, ]

# Variables of interest (including response)
predictors <- c("Mean_Canopy_Height_2020", "Perc_under_2m_2020", "Can_cover_2020", "Height_cv_2020", 
                "Max_Canopy_Height_2020", "Gap_Proportion_2020", "LAI", "Intensity", "VCI", 
                "Canopy_Volume_2020", "FT_5_10m", "FT_1.37m", "FT_10_20m")

response_var<-c("Mean_Shannon_Index")

best_model <- select_best_model(sp_data_unique, predictors,response_var, cor_threshold = 0.7, shape="linear")
# Also tried with cor thresh of 0.5 but model wasnt as stong

summary(sp_data_unique$Mean_Shannon_Index)
summary(sp_data_unique$VCI)
summary(sp_data_unique$LAI)

#------#
# Plot #
#------#

# Scale predictors so they are comparable
# Identify numeric columns
numeric_vars <- sapply(sp_data_unique, is.numeric)  
sp_data_unique_scaled <- sp_data_unique  # Create a copy of the dataset
sp_data_unique_scaled[numeric_vars] <- scale(sp_data_unique[numeric_vars])  # Scale only numeric columns

model_scaled <- lm(Mean_Shannon_Index ~ Max_Canopy_Height_2020 + Gap_Proportion_2020 +      LAI + Intensity + VCI + Canopy_Volume_2020 , 
                   data = sp_data_unique_scaled)

summary(model_scaled)

# Extract coefficients and confidence intervals
coefs <- summary(model_scaled)$coefficients
coef_names <- rownames(coefs)
estimates <- coefs[, 1]
errors <- coefs[, 2]
lower <- estimates - 1.96 * errors  # 95% CI lower bound
upper <- estimates + 1.96 * errors  # 95% CI upper bound

# Exclude the intercept
coef_names <- coef_names[-1]
estimates <- estimates[-1]
errors <- errors[-1]
lower <- lower[-1]
upper <- upper[-1]

#~~~~~~~~~~~~~~~~~~~#
# Plot Coefficients #
#~~~~~~~~~~~~~~~~~~~#

pdf("~/Library/CloudStorage/OneDrive-UniversityofBristol/GitHub/Thesis/Chapter4/Figures/CoeffsCanopyStructure.pdf", width = 7, height = 5)

par(mar = c(5, 16, 4, 2)) 

# Double text size by setting cex (default is 1)
cex_text <- 1.5

plot(estimates, seq_along(estimates), pch = 16, xlim = c(-1, 1),
     xlab = "Estimate", ylab = "", axes = FALSE, cex = cex_text, cex.lab = cex_text)

axis(1, cex.axis = cex_text)

# Labels on y-axis (custom text)
new_names <- c("Max Canopy Height", "Canopy Gap Proportion", "LAI", "Intensity", "VCI", "Canopy Volume")

axis(2, at = seq_along(coef_names), labels = FALSE, las = 2)
text(x = rep(par("usr")[1]-0.1, length(new_names)), y = seq_along(new_names), 
     labels = new_names, adj = 1, xpd = TRUE, cex = cex_text)

arrows(lower, seq_along(estimates), upper, seq_along(estimates), length = 0.1, angle = 90, code = 3)
abline(v = 0, lty = 2, col = "hotpink")  # Reference line

# Add "a)" in the top-left corner (inside the plot margin)
# Use xpd=NA to allow drawing in the figure margin
text(x = par("usr")[1], y = par("usr")[4] + 0.5, labels = "a)", adj = c(0, 1), cex = cex_text, xpd = NA)

dev.off()


#////////////////////#

##### Topography #####

#////////////////////#

# Load data
sp_data<-read.csv('final_dataset.csv')
colnames(sp_data)

# Select only the relevant columns for uniqueness check
unique_rows <- !duplicated(sp_data[, c("Coppicing_Plot")])
sp_data_unique <- sp_data[unique_rows, ]
colnames(sp_data_unique)
hist(sp_data_unique$Mean_Shannon_Index)

# Variables of interest (including response)
predictors <- c("HALP","Aspect_Cos", "TPI", "Roughness", 
                "TRI", "Elevation", "Slope", "FlowDir","Plane_Curve","Profile_Curve",
                "Mean_Curve","Solar_Radiation")

response_var<-c("Mean_Shannon_Index")

best_model_topo <- select_best_model(sp_data_unique, predictors, response_var,cor_threshold = 0.7)
hist(sp_data_unique$Mean_Shannon_Index)
best_model_topo<- glm(Mean_Shannon_Index~HALP+ Aspect_Cos+ TPI+Slope+ FlowDir+ Profile_Curve, family = Gamma(link = "log"),sp_data_unique)
vif(best_model_topo)
shapiro.test(best_model_topo$residuals)

summary(best_model_topo$model)

plot(fitted(best_model_topo), residuals(best_model_topo, type = "pearson"),
     xlab = "Fitted values", ylab = "Pearson residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red")

res_squared <- residuals(best_model_topo, type = "pearson")^2
fit_vals <- fitted(best_model_topo)
summary(lm(res_squared ~ I(fit_vals^2)))

#------#
# Plot #
#------#

# Make a copy of your data
sp_data_unique_scaled <- sp_data_unique

# Scale only predictor columns (not the response)
sp_data_unique_scaled[predictors] <- scale(sp_data_unique_scaled[predictors])

# Fit Gamma model using unscaled response and scaled predictors
model_scaled <- glm(
  Mean_Shannon_Index ~ HALP + Aspect_Cos + TPI + Slope + FlowDir + Profile_Curve,
  family = Gamma(link = "log"),
  data = sp_data_unique_scaled
)

summary(model_scaled)

# Extract coefficients and confidence intervals
coefs <- summary(model_scaled)$coefficients
coef_names <- rownames(coefs)
estimates <- coefs[, 1]
errors <- coefs[, 2]
lower <- estimates - 1.96 * errors  # 95% CI lower bound
upper <- estimates + 1.96 * errors  # 95% CI upper bound

# Exclude the intercept
coef_names <- coef_names[-1]
estimates <- estimates[-1]
errors <- errors[-1]
lower <- lower[-1]
upper <- upper[-1]

#------#
# Plot #
#------#

# Open a PNG device
pdf("~/Library/CloudStorage/OneDrive-UniversityofBristol/GitHub/Thesis/Chapter4/Figures/CoeffsTopo.pdf", width = 7, height = 5)

par(mar = c(5, 16, 4, 2)) 

# Double text size by setting cex (default is 1)
cex_text <- 1.5

plot(estimates, seq_along(estimates), pch = 16,  xlim = c(- 0.2, 0.2),
     xlab = "Estimate", ylab = "", axes = FALSE, cex = cex_text, cex.lab = cex_text)
axis(1,cex.axis = cex_text)

# New labels
new_names<- c("HALP","Aspect", "TPI", "Slope", "Flow Direction", "Profile Curvature")

# Remove the axis(2) call and replace it with text for custom labels
axis(2, at = seq_along(coef_names), labels = FALSE, las = 2, cex=cex_text)
text(x = rep(par("usr")[1]-0.03, length(new_names)), y = seq_along(new_names), 
     labels = new_names, adj = 1, xpd = TRUE, cex= cex_text)

arrows(lower, seq_along(estimates), upper, seq_along(estimates), length = 0.1, angle = 90, code = 3)
abline(v = 0, lty = 2, col = "hotpink")  # Reference line

# Add "a)" in the top-left corner (inside the plot margin)
# Use xpd=NA to allow drawing in the figure margin
text(x = par("usr")[1], y = par("usr")[4] + 0.5, labels = "c)", adj = c(0, 1), cex = cex_text, xpd = NA)

dev.off()

#////////////////////#

##### Management #####

#////////////////////#

sp_data<-read.csv('final_dataset.csv')

# Get recovery period
sp_data <- sp_data %>%
  mutate(recovery_period = 2020 - Coppicing_year)

# Select only the relevant columns for uniqueness check
unique_rows <- !duplicated(sp_data[, c("Coppicing_Plot")])

# Keep only unique rows in the dataframe
sp_data_unique <- sp_data[unique_rows, ]

# Variables of interest
predictors <- c("recovery_period")

response_var<-c("Mean_Shannon_Index")

model_manag<- select_best_model(sp_data_unique, predictors, response_var,cor_threshold = 0.7, shape="quadratic")






# above i expected a u or n shape, i then tested for flattening (dont inc)
library(mgcv)
gam_model <- gam(Mean_Shannon_Index ~ s(recovery_period, k = 4),
                 data = sp_data_unique, family = Gamma(link = "log"))
summary(gam_model)
plot(gam_model)

plot(gam_model, residuals = TRUE, pch = 16)
shapiro.test(gam_model$residuals)
gam.check(gam_model)
# no flattening 

gam_model2 <- gam(Mean_Shannon_Index ~ s(recovery_period, k=6), 
                  data = sp_data_unique, family=Gamma(link="log"))
summary(gam_model2)
plot(gam_model2, residuals=TRUE)

#------#
# Plot #
#------#

# Open a PNG device
pdf("~/Library/CloudStorage/OneDrive-UniversityofBristol/GitHub/Thesis/Chapter4/Figures/Management.pdf", width = 5, height = 5)

# 1. Create prediction data frame
recovery_seq <- seq(min(sp_data_unique$recovery_period, na.rm = TRUE),
                    max(sp_data_unique$recovery_period, na.rm = TRUE),
                    length.out = 200)

new_data <- data.frame(recovery_period = recovery_seq)
new_data$I.recovery_period.2. <- recovery_seq^2  # Match model term

# 2. Predict fitted values *and* standard errors (on link scale)
pred <- predict(model_manag, newdata = new_data, type = "link", se.fit = TRUE)

# 3. Convert to response scale using inverse link (exp for log link)
new_data$fit <- exp(pred$fit)
new_data$lower <- exp(pred$fit - 1.96 * pred$se.fit)
new_data$upper <- exp(pred$fit + 1.96 * pred$se.fit)

# 4. Plot with confidence ribbon
ggplot() +
  geom_point(data = sp_data_unique, aes(x = recovery_period, y = Mean_Shannon_Index),
             color = "pink", size = 2) +
  geom_ribbon(data = new_data, aes(x = recovery_period, ymin = lower, ymax = upper),
              fill = "red", alpha = 0.2) +
  geom_line(data = new_data, aes(x = recovery_period, y = fit), 
            color = "hotpink", size = 1.2) +
  labs(
    x = "Recovery Period (Years)",
    y = "Mean Shannon Index"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(hjust = 1),
    axis.title = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),  # remove full border
    axis.line = element_line(color = "black", size=0.25)  # draw x and y axes only
  )


dev.off()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### STEP 16: RS METRICS AND DORMICE ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#-------------#
# Correlation #
#-------------#

#//////////////////////////#

##### Canopy structure #####

#//////////////////////////#

#-------#
# Model #
#-------#

# Load data
sp_data<-read.csv('final_dataset.csv')

# Select only the relevant columns for uniqueness check
unique_rows <- !duplicated(sp_data[, c("Coppicing_Plot")])
sp_data_unique <- sp_data[unique_rows, ]
colnames(sp_data_unique)

# Variables of interest (including response)
predictors <- c("Mean_Canopy_Height_2020", "Perc_under_2m_2020", "Can_cover_2020", "Height_cv_2020", 
                "Max_Canopy_Height_2020", "Gap_Proportion_2020", "LAI", "Intensity", "VCI", 
                "Canopy_Volume_2020", "FT_5_10m", "FT_1.37m", "FT_10_20m")

response_var<-c("Mean_Dormice")

best_model <- select_best_model(sp_data_unique, predictors,response_var, cor_threshold = 0.7)
# No model fit

  cor_matrix <- cor(sp_data_unique[, predictors], use = "pairwise.complete.obs", method = "spearman")
  valid_combinations <- list(); counter <- 1
  cor_threshold = 0.7
  
  for (i in (length(predictors):1)) {
    combs <- combn(predictors, i, simplify = FALSE)
    for (combo in combs) {
      sub_cor <- cor_matrix[combo, combo]
      if (all(abs(sub_cor[upper.tri(sub_cor)]) < cor_threshold)) {
        valid_combinations[[counter]] <- combo
        counter <- counter + 1
      }
    }
    if (length(valid_combinations) > 0) break
  }

valid_combinations

# Filter out combinations with more than one FT_ variable
filtered_combinations <- Filter(function(combo) {
  sum(grepl("^FT_", combo)) <= 1
}, valid_combinations)

choose_best_count_glm <- function(response_var, predictors, data) {
  formula <- as.formula(paste(response_var, "~", paste(predictors, collapse = "+")))
  
  # Fit Poisson model
  poisson_model <- glm(formula, data = sp_data_unique, family = poisson())
  
  # Test for overdispersion
  disp_test <- tryCatch({
    dispersiontest(poisson_model)
  }, error = function(e) NULL)
  
  is_overdispersed <- !is.null(disp_test) && disp_test$p.value < 0.05
  
  if (is_overdispersed) {
    message("❗ Overdispersion detected. Using Negative Binomial GLM.")
    nb_model <- glm.nb(formula, data = sp_data_unique)
    return(list(model = nb_model, type = "Negative Binomial", overdispersed = TRUE))
  } else {
    message("✅ No overdispersion detected. Using Poisson GLM.")
    return(list(model = poisson_model, type = "Poisson", overdispersed = FALSE))
  }
}

best_model <- NULL
best_aic <- Inf
best_combo <- NULL

for (combo in filtered_combinations) {
  result <- choose_best_count_glm(response_var = response_var, predictors = combo, data = sp_data_unique)
  current_aic <- AIC(result$model)
  
  if (current_aic < best_aic) {
    best_aic <- current_aic
    best_model <- result
    best_combo <- combo
  }
}

cat("✅ Best model type:", best_model$type, "\n")
cat("📉 AIC:", best_aic, "\n")
cat("🔢 Predictors:", paste(best_combo, collapse = ", "), "\n")
summary(best_model$model)
disp<- deviance(best_model$model) / df.residual(best_model$model)
disp
shapiro.test(best_model$model$residuals)

#dispersion >1.5: use neg bin instead
gam_nb<-gam(Mean_Dormice ~  Max_Canopy_Height_2020+ Gap_Proportion_2020+ LAI+ Intensity+ VCI+ Canopy_Volume_2020 , family = nb(), data = sp_data_unique)
summary(gam_nb)

gam.check(gam_nb)

hist(sp_data_unique$Mean_Dormice)

# alot of zeros so need to zero inflate:
gam_nb<-gam(Mean_Dormice ~  Max_Canopy_Height_2020+ Gap_Proportion_2020+ LAI+ Intensity+ VCI+ Canopy_Volume_2020 , family = nb(), data = sp_data_unique)
summary(gam_nb)

gam.check(gam_nb)

hist(sp_data_unique$Mean_Dormice)
mean(sp_data_unique$Mean_Dormice == 0)

install.packages('mgcViz')
library(DHARMa)
library(mgcViz)
simres <- simulateResiduals(fittedModel = gam_nb, plot = TRUE)
testZeroInflation(simres)

#------#
# Plot #
#------#

# Scale predictors so they are comparable
sp_data_unique_scaled <- sp_data_unique
sp_data_unique_scaled[, predictors] <- scale(sp_data_unique_scaled[, predictors])

model_scaled <- gam(Mean_Dormice ~ Max_Canopy_Height_2020 + Gap_Proportion_2020 + 
                      LAI + Intensity + VCI + Canopy_Volume_2020, 
                    family = nb(), data = sp_data_unique_scaled)

summary(model_scaled)

# Extract coefficients and confidence intervals
coefs <- summary(model_scaled)$p.table  # For parametric coefficients in mgcv::gam
coef_names <- rownames(coefs)
estimates <- coefs[, "Estimate"]
errors <- coefs[, "Std. Error"]
lower <- estimates - 1.96 * errors  # 95% CI lower bound
upper <- estimates + 1.96 * errors  # 95% CI upper bound


# Exclude the intercept
coef_names <- coef_names[-1]
estimates <- estimates[-1]
errors <- errors[-1]
lower <- lower[-1]
upper <- upper[-1]

#~~~~~~~~~~~~~~~~~~~#
# Plot Coefficients #
#~~~~~~~~~~~~~~~~~~~#

# Open a PNG device
pdf("~/Library/CloudStorage/OneDrive-UniversityofBristol/GitHub/Thesis/Chapter4/Figures/CoeffsCanopyStructure_DM.pdf", width = 5, height = 5)

cex_text<-1.5
par(mar = c(5, 5, 4, 2)) 
xlim_range <- range(c(lower, upper)) * 1.1  # Add 10% padding
plot(estimates, seq_along(estimates), pch = 16, xlim = xlim_range,
     xlab = "Estimate", ylab = "", axes = FALSE, cex=cex_text, cex.lab = cex_text)

axis(1, cex.axis = cex_text)

# Remove the axis(2) call and replace it with text for custom labels
axis(2, at = seq_along(coef_names), labels = FALSE, las = 2)

arrows(lower, seq_along(estimates), upper, seq_along(estimates), length = 0.1, angle = 90, code = 3)
abline(v = 0, lty = 2, col = "hotpink")  # Reference line

text(x = par("usr")[1], y = par("usr")[4] + 0.5, labels = "b)", adj = c(0, 1), cex = cex_text, xpd = NA)


dev.off()

#////////////////////#

##### Topography #####

#////////////////////#

# Load data
sp_data<-read.csv('final_dataset.csv')

# Select only the relevant columns for uniqueness check
unique_rows <- !duplicated(sp_data[, c("Coppicing_Plot")])
sp_data_unique <- sp_data[unique_rows, ]

# Variables of interest
predictors <- c("HALP","Aspect_Cos", "TPI", "Roughness", 
                "TRI", "Elevation", "Slope", "FlowDir","Plane_Curve","Profile_Curve",
                "Mean_Curve","Solar_Radiation")

response_var<-c("Mean_Dormice")

best_model_topo <- select_best_model(sp_data_unique, predictors, response_var,cor_threshold = 0.7)

best_model_topo<- lm(log(Mean_Dormice+1) ~HALP+ Aspect_Cos+ TPI+ TRI+FlowDir+ Profile_Curve, sp_data_unique)
summary(best_model_topo)
shapiro.test(best_model_topo$residuals)

cor_matrix <- cor(sp_data_unique[, predictors], use = "pairwise.complete.obs", method = "spearman")
valid_combinations <- list(); counter <- 1
cor_threshold = 0.7

for (i in (length(predictors):1)) {
  combs <- combn(predictors, i, simplify = FALSE)
  for (combo in combs) {
    sub_cor <- cor_matrix[combo, combo]
    if (all(abs(sub_cor[upper.tri(sub_cor)]) < cor_threshold)) {
      valid_combinations[[counter]] <- combo
      counter <- counter + 1
    }
  }
  if (length(valid_combinations) > 0) break
}

valid_combinations

# Filter out combinations with more than one FT_ variable
filtered_combinations <- Filter(function(combo) {
  sum(grepl("^FT_", combo)) <= 1
}, valid_combinations)

choose_best_count_glm <- function(response_var, predictors, data) {
  formula <- as.formula(paste(response_var, "~", paste(predictors, collapse = "+")))
  
  # Fit Poisson model
  poisson_model <- glm(formula, data = sp_data_unique, family = poisson())
  
  # Test for overdispersion
  disp_test <- tryCatch({
    dispersiontest(poisson_model)
  }, error = function(e) NULL)
  
  is_overdispersed <- !is.null(disp_test) && disp_test$p.value < 0.05
  
  if (is_overdispersed) {
    message("❗ Overdispersion detected. Using Negative Binomial GLM.")
    nb_model <- glm.nb(formula, data = sp_data_unique)
    return(list(model = nb_model, type = "Negative Binomial", overdispersed = TRUE))
  } else {
    message("✅ No overdispersion detected. Using Poisson GLM.")
    return(list(model = poisson_model, type = "Poisson", overdispersed = FALSE))
  }
}

best_model <- NULL
best_aic <- Inf
best_combo <- NULL

for (combo in filtered_combinations) {
  result <- choose_best_count_glm(response_var = response_var, predictors = combo, data = sp_data_unique)
  current_aic <- AIC(result$model)
  
  if (current_aic < best_aic) {
    best_aic <- current_aic
    best_model <- result
    best_combo <- combo
  }
}

cat("✅ Best model type:", best_model$type, "\n")
cat("📉 AIC:", best_aic, "\n")
cat("🔢 Predictors:", paste(best_combo, collapse = ", "), "\n")
summary(best_model$model)

disp<- deviance(best_model$model) / df.residual(best_model$model)
disp
shapiro.test((best_model$model$residuals))

#dispersion >1.5: use neg bin instead
gam_nb<-gam(Mean_Dormice ~  HALP+ Aspect_Cos+ TPI+ Slope+ FlowDir+ Profile_Curve, family = nb(), data = sp_data_unique)
summary(gam_nb)
gam.check(gam_nb)

mean(sp_data_unique$Mean_Dormice == 0)

simres <- simulateResiduals(fittedModel = gam_nb, plot = TRUE)
testZeroInflation(simres)

#------#
# Plot #
#------#

# Scale predictors so they are comparable

sp_data_unique_scaled <- sp_data_unique
sp_data_unique_scaled[, predictors] <- scale(sp_data_unique_scaled[, predictors])

model_scaled <- gam(Mean_Dormice ~  HALP+ Aspect_Cos+ TPI+ Slope+ FlowDir+ Profile_Curve, family = nb(), data = sp_data_unique_scaled)

summary(model_scaled)

coefs <- summary(model_scaled)$p.table  # For parametric coefficients in mgcv::gam
coef_names <- rownames(coefs)
estimates <- coefs[, "Estimate"]
errors <- coefs[, "Std. Error"]
lower <- estimates - 1.96 * errors  # 95% CI lower bound
upper <- estimates + 1.96 * errors  # 95% CI upper bound


# Exclude the intercept
coef_names <- coef_names[-1]
estimates <- estimates[-1]
errors <- errors[-1]
lower <- lower[-1]
upper <- upper[-1]


#------#
# Plot #
#------#

# Open a PNG device
pdf("~/Library/CloudStorage/OneDrive-UniversityofBristol/GitHub/Thesis/Chapter4/Figures/CoeffsTopoDM.pdf", width = 5, height = 5)

cex_text<-1.5
par(mar = c(5, 5, 4, 2)) 
xlim_range <- range(c(-1.5, 1.5))
plot(estimates, seq_along(estimates), pch = 16, xlim = xlim_range,
     xlab = "Estimate", ylab = "", axes = FALSE, cex=cex_text,cex.lab=cex_text)
axis(1, at = seq(-1.5, 1.5, by = 0.5), cex.axis=cex_text)  # Adjust ticks


# Remove the axis(2) call and replace it with text for custom labels
axis(2, at = seq_along(coef_names), labels = FALSE, las = 2)

arrows(lower, seq_along(estimates), upper, seq_along(estimates), length = 0.1, angle = 90, code = 3)
abline(v = 0, lty = 2, col = "hotpink")  # Reference line

text(x = par("usr")[1], y = par("usr")[4] + 0.5, labels = "d)", adj = c(0, 1), cex = cex_text, xpd = NA)


dev.off()

#////////////////////#

##### Management #####

#////////////////////#

# Load data
sp_data<-read.csv('final_dataset.csv')

# Select only the relevant columns for uniqueness check
unique_rows <- !duplicated(sp_data[, c("Coppicing_Plot")])
sp_data_unique <- sp_data[unique_rows, ]

# Get recovery period
sp_data_unique <- sp_data_unique %>%
  mutate(recovery_period = 2020 - Coppicing_year)

predictors <- c("recovery_period")
response_var <- "Mean_Dormice"

best_model_manage <- select_best_model(sp_data_unique, predictors, response_var,cor_threshold = 0.7, shape="quadratic")

best_model_manage <-lm(log(Mean_Dormice+1) ~ recovery_period + I(recovery_period^2), data=sp_data_unique)
summary(best_model_manage)
shapiro.test(best_model_manage$residuals)

best_model_manage<- choose_best_count_glm(response_var = response_var, predictors = predictors, data = sp_data_unique)

mod<-glm(Mean_Dormice ~ recovery_period + I(recovery_period^2), data=sp_data_unique, family='poisson')
summary(mod)
dispersion_ratio <- sum(residuals(mod, type = "pearson")^2) / mod$df.residual
dispersion_ratio

observed_zeros <- sum(sp_data_unique$Mean_Dormice == 0)
expected_zeros <- sum(dpois(0, lambda = predict(mod, type = "response")))

cat("Observed zeros:", observed_zeros, "\n")
cat("Expected zeros (Poisson):", round(expected_zeros), "\n")

nb_model <- glm.nb(Mean_Dormice ~ recovery_period + I(recovery_period^2), data = sp_data_unique)

observed_zeros <- sum(nb_model$Mean_Dormice == 0)
expected_zeros <- sum(dpois(0, lambda = predict(mod, type = "response")))

cat("Observed zeros:", observed_zeros, "\n")
cat("Expected zeros (Poisson):", round(expected_zeros), "\n")


# try gam 
gam_count <- gam(Mean_Dormice ~ recovery_period + I(recovery_period^2),
                 family = poisson(link = "log"),
                 data = sp_data_unique)
summary(gam_count)

dispersion <- deviance(gam_count) / df.residual(gam_count)
dispersion #overdispersed

library(mgcv)
gam_nb <- gam(Mean_Dormice ~ recovery_period+ I(recovery_period^2), family = nb(), data = sp_data_unique)
summary(gam_nb)


dispersion <- deviance(gam_nb) / df.residual(gam_nb)
dispersion #overdispersed
observed_zeros <- sum(gam_nb$Mean_Dormice == 0)
expected_zeros <- sum(dpois(0, lambda = predict(gam_nb, type = "response")))

cat("Observed zeros:", observed_zeros, "\n")
cat("Expected zeros (Poisson):", round(expected_zeros), "\n")

#------#
# Plot #
#------#

# Open a PNG device
pdf("~/Library/CloudStorage/OneDrive-UniversityofBristol/GitHub/Thesis/Chapter4/Figures/ManagementDM.pdf", width = 5, height = 5)

ggplot(sp_data_unique, aes(x = recovery_period, y = Mean_Dormice)) +
  geom_point(color = "pink", size = 2) +  # Data points
  geom_smooth(
    method = "glm",
    method.args = list(family = "gaussian"),
    formula = y ~ poly(x, 2),  # Quadratic GLM fit
    color = "hotpink",
    fill = "red",
    alpha = 0.2
  ) +
  labs(
    x = "Recovery Period (Years)",
    y = "Mean Dormice Abundance"
  ) +
  scale_y_continuous(breaks = seq(0, 8, by = 2))+  # Or higher if needed
  theme_minimal() +
  theme(
    axis.text.x = element_text(hjust = 1),
    axis.title = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.25),  # Thin x and y axes,
    plot.margin = unit(c(20, 10, 10, 10), "pt"),
    axis.text = element_text(size = 12)
    
  )

dev.off()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### STEP 16: FULL MODELS ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#--------------#
# Mean Dormice #
#--------------#

# Load data
sp_data<-read.csv('final_dataset.csv')

# Select only the relevant columns for uniqueness check
unique_rows <- !duplicated(sp_data[, c("Coppicing_Plot")])
sp_data_unique <- sp_data[unique_rows, ]

# Get recovery period
sp_data_unique <- sp_data_unique %>%
  mutate(recovery_period = 2020 - Coppicing_year)

colnames(sp_data_unique)

predictors <- c("Mean_Canopy_Height_2020", "Perc_under_2m_2020", "Can_cover_2020", "Height_cv_2020", 
                "Max_Canopy_Height_2020", "Gap_Proportion_2020", "LAI", "Intensity", "VCI", 
                "Canopy_Volume_2020", "FT_5_10m", "FT_1.37m", "FT_10_20m", "HALP","Aspect_Cos", "TPI", "Roughness", 
                "TRI", "Elevation", "Slope", "FlowDir","Plane_Curve","Profile_Curve",
                "Mean_Curve","Solar_Radiation", "recovery_period", "Mean_Shannon_Index")

response_var<- "Mean_Dormice"

#best_full_model<- select_best_model(sp_data_unique, predictors, response_var,cor_threshold = 0.7)

cor_matrix <- cor(sp_data_unique[, predictors], use = "pairwise.complete.obs", method = "spearman")
valid_combinations <- list(); counter <- 1
cor_threshold = 0.7

for (i in (length(predictors):1)) {
  combs <- combn(predictors, i, simplify = FALSE)
  for (combo in combs) {
    sub_cor <- cor_matrix[combo, combo]
    if (all(abs(sub_cor[upper.tri(sub_cor)]) < cor_threshold)) {
      valid_combinations[[counter]] <- combo
      counter <- counter + 1
    }
  }
  if (length(valid_combinations) > 0) break
}

valid_combinations

# Filter out combinations with more than one FT_ variable
filtered_combinations <- Filter(function(combo) {
  sum(grepl("^FT_", combo)) <= 1
}, valid_combinations)

choose_best_count_glm <- function(response_var, predictors, data) {
  formula <- as.formula(paste(response_var, "~", paste(predictors, collapse = "+")))
  
  # Fit Poisson model
  poisson_model <- glm(formula, data = sp_data_unique, family = poisson())
  
  # Test for overdispersion
  disp_test <- tryCatch({
    dispersiontest(poisson_model)
  }, error = function(e) NULL)
  
  is_overdispersed <- !is.null(disp_test) && disp_test$p.value < 0.05
  
  if (is_overdispersed) {
    message("❗ Overdispersion detected. Using Negative Binomial GLM.")
    nb_model <- glm.nb(formula, data = sp_data_unique)
    return(list(model = nb_model, type = "Negative Binomial", overdispersed = TRUE))
  } else {
    message("✅ No overdispersion detected. Using Poisson GLM.")
    return(list(model = poisson_model, type = "Poisson", overdispersed = FALSE))
  }
}

best_model <- NULL
best_aic <- Inf
best_combo <- NULL

for (combo in filtered_combinations) {
  result <- choose_best_count_glm(response_var = response_var, predictors = combo, data = sp_data_unique)
  current_aic <- AIC(result$model)
  
  if (current_aic < best_aic) {
    best_aic <- current_aic
    best_model <- result
    best_combo <- combo
  }
}

cat("✅ Best model type:", best_model$type, "\n")
cat("📉 AIC:", best_aic, "\n")
cat("🔢 Predictors:", paste(best_combo, collapse = ", "), "\n")
summary(best_model$model)

mod<-best_model$model

shapiro.test(mod$residuals)
# Calculate dispersion
dispersion <- deviance(mod) / df.residual(mod)
dispersion

#not over dispersed

# cat("✅ Best model type:", best_model$type, "\n")
# ✅ Best model type: Poisson 
# > cat("📉 AIC:", best_aic, "\n")
# 📉 AIC: 146.2439 
# > cat("🔢 Predictors:", paste(best_combo, collapse = ", "), "\n")
# 🔢 Predictors: Max_Canopy_Height_2020, Gap_Proportion_2020, LAI, Intensity, VCI, Canopy_Volume_2020, Aspect_Cos, TRI, Elevation, FlowDir, Profile_Curve, Mean_Curve, recovery_period, Mean_Shannon_Index 
# > summary(best_model$model)
# 
# Call:
#   glm(formula = formula, family = poisson(), data = sp_data_unique)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)            -3.336e+01  1.180e+01  -2.827 0.004691 ** 
#   Max_Canopy_Height_2020  2.051e-01  1.149e-01   1.786 0.074106 .  
# Gap_Proportion_2020     6.600e-02  1.823e-02   3.621 0.000294 ***
#   LAI                     1.749e+00  4.937e-01   3.543 0.000395 ***
#   Intensity               3.936e-03  1.778e-03   2.214 0.026828 *  
#   VCI                     3.375e+01  9.932e+00   3.398 0.000680 ***
#   Canopy_Volume_2020     -8.473e-06  1.022e-05  -0.829 0.407114    
# Aspect_Cos             -2.129e+00  8.998e-01  -2.366 0.017968 *  
#   TRI                     3.170e+01  1.620e+01   1.957 0.050296 .  
# Elevation              -7.332e-02  2.541e-02  -2.885 0.003908 ** 
#   FlowDir                -7.601e-02  5.680e-02  -1.338 0.180866    
# Profile_Curve          -2.486e+02  1.167e+02  -2.131 0.033093 *  
#   Mean_Curve              3.709e+02  1.190e+02   3.118 0.001824 ** 
#   recovery_period         1.920e-02  4.027e-02   0.477 0.633486    
# Mean_Shannon_Index     -1.584e+00  6.234e-01  -2.542 0.011036 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 126.96  on 55  degrees of freedom
# Residual deviance:  55.86  on 41  degrees of freedom
# AIC: 146.24
# 
# Number of Fisher Scoring iterations: 6
range(sp_data_unique$LAI)
#------#
# Plot #
#------#

# Scale predictors so they are comparable

sp_data_unique_scaled <- sp_data_unique
sp_data_unique_scaled[, predictors] <- scale(sp_data_unique_scaled[, predictors])

model_scaled <- glm(Mean_Dormice ~ Max_Canopy_Height_2020 + Gap_Proportion_2020 +      LAI + Intensity + VCI + Canopy_Volume_2020 + Aspect_Cos +      Roughness + Elevation + FlowDir + Plane_Curve + Profile_Curve +      recovery_period + Mean_Shannon_Index, family='poisson', sp_data_unique_scaled)

summary(model_scaled)

# Extract coefficients and confidence intervals
coefs <- summary(model_scaled)$coefficients
coef_names <- rownames(coefs)
estimates <- coefs[, 1]
errors <- coefs[, 2]
lower <- estimates - 1.96 * errors  # 95% CI lower bound
upper <- estimates + 1.96 * errors  # 95% CI upper bound

# Exclude the intercept
coef_names <- coef_names[-1]
estimates <- estimates[-1]
errors <- errors[-1]
lower <- lower[-1]
upper <- upper[-1]

#------#
# Plot #
#------#

# Open a PNG device
png("~/Library/CloudStorage/OneDrive-UniversityofBristol/Documents/PhD/Project 3_Bontuchel/Figures/FullModelDM.png", width = 900, height = 900, res=150)

par(mar = c(8, 10,2, 1)) 
xlim_range <- range(c(-1.5, 2.5))
plot(estimates, seq_along(estimates), pch = 16, xlim = xlim_range,
     xlab = "Estimate", ylab = "", axes = FALSE)
axis(1, at = seq(-1.5, 2.5, by = 0.5))  # Adjust ticks

# New labels
new_names<- c("Max Canopy Height","Gap proportion","LAI", "Intensity","VCI", "Canopy Volume", "Aspect", "Roughness", "Elevation", "Flow Direction", "Plane Curvature", "Profile Curvature", "Recovery Period", "Mean Shannon Index")

# Remove the axis(2) call and replace it with text for custom labels
axis(2, at = seq_along(coef_names), labels = FALSE, las = 2)
text(x = rep(par("usr")[1]-0.03, length(new_names)), y = seq_along(new_names), 
     labels = new_names, adj = 1.1, xpd = TRUE)

arrows(lower, seq_along(estimates), upper, seq_along(estimates), length = 0.1, angle = 90, code = 3)
abline(v = 0, lty = 2, col = "hotpink")  # Reference line

dev.off()



#---------------#
# sjannon #
#---------------#

# Load data
sp_data<-read.csv('final_dataset.csv')

# Select only the relevant columns for uniqueness check
unique_rows <- !duplicated(sp_data[, c("Coppicing_Plot")])
sp_data_unique <- sp_data[unique_rows, ]


# Get recovery period
sp_data_unique <- sp_data_unique %>%
  mutate(recovery_period = 2020 - Coppicing_year)

predictors <- c("Mean_Canopy_Height_2020", "Perc_under_2m_2020", "Can_cover_2020", "Height_cv_2020", 
                "Max_Canopy_Height_2020", "Gap_Proportion_2020", "LAI", "Intensity", "VCI", 
                "Canopy_Volume_2020", "FT_5_10m", "FT_1.37m", "FT_10_20m", "HALP","Aspect_Cos", "TPI", "Roughness", 
                "TRI", "Elevation", "Slope", "FlowDir","Plane_Curve","Profile_Curve",
                "Mean_Curve","Solar_Radiation", "recovery_period", "Mean_Dormice")

response_var<- "Mean_Shannon_Index"

best_full_model_DM<- select_best_model(sp_data_unique, predictors, response_var,cor_threshold = 0.7)

# ✅ Best: Gaussian GLM 
# Formula: Mean_Shannon_Index ~ Perc_under_2m_2020 + Max_Canopy_Height_2020 +      Gap_Proportion_2020 + Intensity + VCI + Canopy_Volume_2020 +      HALP + Aspect_Cos + Slope + FlowDir + Plane_Curve + Profile_Curve +      recovery_period + Mean_Dormice 
# Shapiro p = 0.075 
# AIC = 46.07 
# Reason: Gaussian GLM with normal residuals 
# 
# 
# Call:
#   glm(formula = formula, family = gaussian, data = data)
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            -3.049e+00  1.609e+00  -1.896 0.065057 .  
# Perc_under_2m_2020     -2.695e-02  1.963e-02  -1.373 0.177194    
# Max_Canopy_Height_2020  2.105e-02  2.235e-02   0.942 0.351787    
# Gap_Proportion_2020     1.103e-02  4.468e-03   2.470 0.017778 *  
#   Intensity               5.892e-04  6.458e-04   0.912 0.366910    
# VCI                     4.857e+00  1.360e+00   3.572 0.000923 ***
#   Canopy_Volume_2020     -4.771e-10  3.034e-06   0.000 0.999875    
# HALP                   -1.287e-02  5.704e-03  -2.256 0.029496 *  
#   Aspect_Cos             -9.955e-02  2.181e-01  -0.456 0.650524    
# Slope                  -7.804e-01  9.137e-01  -0.854 0.398022    
# FlowDir                 9.295e-03  1.342e-02   0.693 0.492506    
# Plane_Curve             2.454e+01  1.419e+01   1.729 0.091367 .  
# Profile_Curve          -4.920e+01  2.496e+01  -1.971 0.055475 .  
# recovery_period         2.396e-04  1.093e-02   0.022 0.982610    
# Mean_Dormice           -7.498e-02  3.123e-02  -2.401 0.020977 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 0.1028129)
# 
# Null deviance: 7.6618  on 55  degrees of freedom
# Residual deviance: 4.2153  on 41  degrees of freedom
# AIC: 46.07
# 
# Number of Fisher Scoring iterations: 2

shapiro.test(best_full_model_DM$residuals)

#------#
# Plot #
#------#

# Scale predictors so they are comparable

sp_data_unique_scaled <- sp_data_unique
sp_data_unique_scaled[, predictors] <- scale(sp_data_unique_scaled[, predictors])

model_scaled <- glm(Mean_Shannon_Index ~ Perc_under_2m_2020 + Max_Canopy_Height_2020 +      Gap_Proportion_2020 + Intensity + VCI + Canopy_Volume_2020 +      HALP + Aspect_Cos + Slope + FlowDir + Plane_Curve + Profile_Curve +      recovery_period + Mean_Dormice , family='gaussian', data= sp_data_unique_scaled)

summary(model_scaled)

# Extract coefficients and confidence intervals
coefs <- summary(model_scaled)$coefficients
coef_names <- rownames(coefs)
estimates <- coefs[, 1]
errors <- coefs[, 2]
lower <- estimates - 1.96 * errors  # 95% CI lower bound
upper <- estimates + 1.96 * errors  # 95% CI upper bound

# Exclude the intercept
coef_names <- coef_names[-1]
estimates <- estimates[-1]
errors <- errors[-1]
lower <- lower[-1]
upper <- upper[-1]

#------#
# Plot #
#------#

# Open a PNG device
png("~/Library/CloudStorage/OneDrive-UniversityofBristol/Documents/PhD/Project 3_Bontuchel/Figures/FullModel.png", width = 900, height = 900, res=150)

par(mar = c(8, 14,2, 1)) 
xlim_range <- range(c(-0.5, 0.5))
plot(estimates, seq_along(estimates), pch = 16, xlim = xlim_range,
     xlab = "Estimate", ylab = "", axes = FALSE)
axis(1, at = seq(-0.5, 0.5, by = 0.25))  # Adjust ticks


# New labels
new_names<- c("Proportion of cover below 2m","Max Canopy Height","Gap proportion", "Intensity","VCI", "Flow Canopy Volume", "HALP", "Aspect", "Slope", "Flow Direction", "Plane Curvature", "Profile Curvature", "Recovery Period", "Mean Dormouse Abundance")

# Remove the axis(2) call and replace it with text for custom labels
axis(2, at = seq_along(coef_names), labels = FALSE, las = 2)
text(x = rep(par("usr")[1]-0.03, length(new_names)), y = seq_along(new_names), 
     labels = new_names, adj = 1.1, xpd = TRUE)

arrows(lower, seq_along(estimates), upper, seq_along(estimates), length = 0.1, angle = 90, code = 3)
abline(v = 0, lty = 2, col = "hotpink")  # Reference line

dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#### STEP 17: FOREST VERTICAL DIVERITY ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

library(dplyr)
library(readr)

setwd("~/Library/CloudStorage/OneDrive-UniversityofBristol/Documents/PhD/Project 3_Bontuchel/Data")

# Load data
sp_data <- read.csv("final_dataset.csv")

# Step 1: Remove duplicate species entries per plot × layer
sp_data_unique <- sp_data %>%
  distinct(Coppicing_Plot, Forest_layer, Species, .keep_all = TRUE)

# Step 3: Define strata based on ET_layer_height
sp_data_unique$Stratum <- cut(
  sp_data_unique$ET_layer_height,
  breaks = c(-Inf, 1.37, 5, 10, 20, Inf),
  labels = c("Understorey", "Shrub", "Lower_Canopy", "Upper_Canopy", "Emergent")
)

# Step 2: Scale proportions so they sum to 1 within each plot × layer
sp_data_corrected <- sp_data_unique %>%
  group_by(Coppicing_Plot, Stratum) %>%
  mutate(Sum_Proportion = sum(Proportion_of_Species, na.rm = TRUE)) %>%
  mutate(Proportion_of_Species = ifelse(Sum_Proportion > 1, Proportion_of_Species / Sum_Proportion, Proportion_of_Species)) %>%
  ungroup()

# Step 4: Calculate Shannon Index per ET_ID
shannon_per_point <- sp_data_corrected %>%
  group_by(Coppicing_Plot, Stratum, ET_ID) %>%
  mutate(
    p = Proportion_of_Species / sum(Proportion_of_Species, na.rm = TRUE),
    p_log_p = -p * log(p)
  ) %>%
  summarise(Shannon_ET = sum(p_log_p, na.rm = TRUE), .groups = "drop")

# Step 5: Average Shannon Index across ET_IDs per plot × stratum
shannon_index_normalised <- shannon_per_point %>%
  group_by(Coppicing_Plot, Stratum) %>%
  summarise(Shannon_Index_strata = mean(Shannon_ET, na.rm = TRUE), .groups = "drop")

# Step 6: Merge back with species proportions (1 row per species per stratum)
model_data_species_wide <- sp_data_corrected %>%
  left_join(shannon_index_normalised, by = c("Coppicing_Plot", "Stratum")) %>%
  dplyr::select(Coppicing_Plot, Stratum, Species, Proportion_of_Species, Shannon_Index_strata)

# Step 7: Add dormice data
dorm <- read.csv("dormice_data_Panels.csv")

total_dormice_sum <- dorm %>%
  group_by(Coppicing_Plot) %>%
  summarise(Total_dormice_sum = sum(as.numeric(Total_dormice), na.rm = TRUE))

# Step 8: Merge dormice total and mean
model_data_species_wide <- model_data_species_wide %>%
  left_join(total_dormice_sum, by = "Coppicing_Plot") %>%
  left_join(
    sp_data %>% dplyr::select(Coppicing_Plot, Mean_Dormice) %>% distinct(),
    by = "Coppicing_Plot"
  )

model_data_species_wide %>%
  group_by(Coppicing_Plot, Stratum) %>%
  summarise(total_prop = sum(Proportion_of_Species), .groups = "drop") %>%
  filter(total_prop > 1)

#filter for Emergent stratum
emergent <- model_data_species_wide %>%
  filter(Stratum == "Emergent")

# Variables of interest
predictors <- c("Shannon_Index_strata")
response_var<-c("Mean_Dormice")

mod<-glm(Mean_Dormice ~ Shannon_Index_strata, data=emergent, family="poisson")
summary(mod)
#Overdisp
mod_nb <- glm.nb(Mean_Dormice ~ Shannon_Index_strata, data = emergent)
summary(mod_nb)
sum(emergent$Mean_Dormice == 0)  # actual
preds <- predict(mod_nb, type = "response")
sum(dnbinom(0, size = mod_nb$theta, mu = preds))  # expected


Upper <- model_data_species_wide %>%
  filter(Stratum == "Upper_Canopy")
mod<-glm(Mean_Dormice ~ Shannon_Index_strata, data=Upper, family="poisson")
summary(mod)
pearson_resid <- residuals(mod, type = "pearson")

# Calculate overdispersion ratio
overdispersion_ratio <- sum(pearson_resid^2) / df.residual(mod)

print(overdispersion_ratio)
mod<-glm.nb(Mean_Dormice ~ Shannon_Index_strata, data=Upper)
summary(mod)


Lower <- model_data_species_wide %>%
  filter(Stratum == "Lower_Canopy")
modellow<-glm(Mean_Dormice ~ Shannon_Index_strata, data=Lower, family="poisson")
summary(modellow)
overdispersion_ratio <- sum(pearson_resid^2) / df.residual(modellow)
print(overdispersion_ratio)
mod<-glm.nb(Mean_Dormice ~ Shannon_Index_strata, data=Lower)
summary(mod)

Shrub <- model_data_species_wide %>%
  filter(Stratum == "Shrub")
mod<-glm(Mean_Dormice ~ Shannon_Index_strata, data=Shrub, family="poisson")
overdispersion_ratio <- sum(pearson_resid^2) / df.residual(mod)
print(overdispersion_ratio)
mod<-glm.nb(Mean_Dormice ~ Shannon_Index_strata, data=Shrub)
summary(mod)

Understorey <- model_data_species_wide %>%
  filter(Stratum == "Understorey")
mod<-glm(Mean_Dormice ~ Shannon_Index_strata, data=Understorey, family="poisson")
summary(mod)
overdispersion_ratio <- sum(pearson_resid^2) / df.residual(mod)
print(overdispersion_ratio)
mod<-glm.nb(Mean_Dormice ~ Shannon_Index_strata, data=Understorey)
summary(mod)
shapiro.test(mod$residuals)

# Diversity in the upper canopy reflects no of dormice found, lower diversity= more mice

#------#
# PLOT #
#------#

library(dplyr)
library(ggplot2)
library(readr)

# Set working directory
setwd("~/Library/CloudStorage/OneDrive-UniversityofBristol/Documents/PhD/Project 3_Bontuchel")

# Load data
sp_data <- read.csv("Data/final_dataset.csv")

# Step 1: Remove duplicate species entries per plot × layer
sp_data_unique <- sp_data %>%
  distinct(Coppicing_Plot, Forest_layer, Species, .keep_all = TRUE)

# Step 2: Define strata based on ET_layer_height
sp_data_unique$Stratum <- cut(
  sp_data_unique$ET_layer_height,
  breaks = c(-Inf, 1.37, 5, 10, 20, Inf),
  labels = c("Understorey", "Shrub", "Lower_Canopy", "Upper_Canopy", "Emergent")
)

# Step 3: Scale proportions
sp_data_corrected <- sp_data_unique %>%
  group_by(Coppicing_Plot, Stratum) %>%
  mutate(Sum_Proportion = sum(Proportion_of_Species, na.rm = TRUE)) %>%
  mutate(Proportion_of_Species = ifelse(Sum_Proportion > 1, Proportion_of_Species / Sum_Proportion, Proportion_of_Species)) %>%
  ungroup()

# Step 4: Shannon Index per ET_ID
shannon_per_point <- sp_data_corrected %>%
  group_by(Coppicing_Plot, Stratum, ET_ID) %>%
  mutate(
    p = Proportion_of_Species / sum(Proportion_of_Species, na.rm = TRUE),
    p_log_p = -p * log(p)
  ) %>%
  summarise(Shannon_ET = sum(p_log_p, na.rm = TRUE), .groups = "drop")

# Step 5: Mean Shannon Index per plot × stratum
shannon_index_normalised <- shannon_per_point %>%
  group_by(Coppicing_Plot, Stratum) %>%
  summarise(Shannon_Index_strata = mean(Shannon_ET, na.rm = TRUE), .groups = "drop")

# Step 6: Merge back
model_data_species_wide <- sp_data_corrected %>%
  left_join(shannon_index_normalised, by = c("Coppicing_Plot", "Stratum")) %>%
  dplyr::select(Coppicing_Plot, Stratum, Species, Proportion_of_Species, Shannon_Index_strata)

# Step 7: Add dormice data
dorm <- read.csv("Data/dormice_data_Panels.csv")

total_dormice_sum <- dorm %>%
  group_by(Coppicing_Plot) %>%
  summarise(Total_dormice_sum = sum(as.numeric(Total_dormice), na.rm = TRUE))

model_data_species_wide <- model_data_species_wide %>%
  left_join(total_dormice_sum, by = "Coppicing_Plot") %>%
  left_join(
    sp_data %>% dplyr::select(Coppicing_Plot, Mean_Dormice) %>% distinct(),
    by = "Coppicing_Plot"
  )

# Step 8: Add Coppicing_Year
model_data_species_wide <- model_data_species_wide %>%
  left_join(
    sp_data %>% dplyr::select(Coppicing_Plot, Coppicing_year) %>% distinct(),
    by = "Coppicing_Plot"
  )

# Step 9: Define consistent color palette (blues, purples, pinks)
all_species <- sort(unique(model_data_species_wide$Species))
num_species <- length(all_species)

# Use cohesive palette (cool tones)
cool_palette <- c(
  "#2e3b8f", "#5c6bc0", "#7986cb", "#9fa8da", # blues
  "#8e24aa", "#ab47bc", "#ba68c8", "#ce93d8", # purples
  "#c2185b", "#ec407a", "#f48fb1", "#f8bbd0"  # pinks
)
# Extend palette if more species than colors
if (num_species > length(cool_palette)) {
  cool_palette <- colorRampPalette(cool_palette)(num_species)
}
species_color_map <- setNames(cool_palette[1:num_species], all_species)

# Step 10: Plot by Coppicing Year
years <- unique(model_data_species_wide$Coppicing_year)

for (year in years) {
  data_year <- model_data_species_wide %>%
    filter(Coppicing_year == year) %>%
    mutate(
      Species = factor(Species, levels = all_species),
      Proportion_of_Species = Proportion_of_Species * 100
    )
  
  p <- ggplot(data_year, aes(x = Stratum, y = Proportion_of_Species, fill = Species)) +
    geom_col(position = "stack", width = 0.8) +
    facet_wrap(~ Coppicing_Plot, scales = "fixed") +
    coord_flip() +
    scale_fill_manual(values = species_color_map, drop = TRUE) +
    scale_x_discrete(labels = c(
      "Upper_Canopy" = "Upper canopy",
      "Lower_Canopy" = "Lower canopy",
      "Shrub" = "Shrub layer",
      "Understorey" = "Understorey",
      "Emergent" = "Emergent"
    )) +
    labs(
      title = paste("Coppicing Year:", year),
      x = "", y = "Proportional Cover (%)"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      axis.title.x = element_text(size = 16, margin = ggplot2::margin(t = 15)),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      strip.text = element_text(size = 16)
    )
  
  # Save plot
  file_name <- paste0("~/Library/CloudStorage/OneDrive-UniversityofBristol/Github/Thesis/SuppInfo/Figures/Chapter4/Species_Cover_CoppicingYear_", year, ".pdf")
  ggsave(file_name, plot = p, width = 12, height = 8, dpi = 300)
  
  print(p)
}

# Save the dataframe
write.csv(model_data_species_wide,"sp_data_shan.csv",row.names = F)


# Overal Plot
# Diversity of each strata on x, (strata different colours), dormice abundance on y

library(ggplot2)
library(dplyr)

# Assuming your processed data is in `model_data_species_wide`

# For clarity, let's just keep unique Shannon index values per plot and stratum
shan_strata <- model_data_species_wide %>%
  distinct(Coppicing_Plot, Stratum, Shannon_Index_strata) %>%
  filter(!is.na(Shannon_Index_strata))  # remove any NA

# Reorder Stratum factor for meaningful order on y-axis (optional)
shan_strata$Stratum <- factor(shan_strata$Stratum,
                              levels = c("Understorey", "Shrub", "Lower_Canopy", "Upper_Canopy", "Emergent"))

# Plot mean Shannon Index on x, stratum on y
p <- ggplot(shan_strata, aes(x = Shannon_Index_strata, y = Stratum)) +
  geom_jitter(height = 0.15, size = 3, alpha = 0.7, color = "#2e3b8f") +  # jitter to avoid overlap
  stat_summary(fun = mean, geom = "point", shape = 18, size = 5, color = "red") +  # mean per stratum
  labs(
    title = "Mean Shannon Index across Strata",
    x = "Mean Shannon Index",
    y = "Stratum"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5)
  )

print(p)



library(ggplot2)
library(dplyr)

# First, get mean dormice abundance per stratum and plot (or just per stratum if you want aggregated)
# Since dormice abundance is per plot, and Shannon Index is per plot × stratum, we can group by both.

# Prepare summary data: average dormice and mean Shannon index by stratum
summary_data <- model_data_species_wide %>%
  group_by(Stratum) %>%
  summarise(
    Mean_Shannon = mean(Shannon_Index_strata, na.rm = TRUE),
    Mean_Dormice = mean(Mean_Dormice, na.rm = TRUE)
  ) %>%
  filter(!is.na(Stratum))

# If you want lines per stratum *across* Coppicing_Plots (e.g. plotting each Coppicing_Plot as x), 
# and lines connecting points across plots ordered by Shannon_Index, you can do:

line_data <- model_data_species_wide %>%
  filter(!is.na(Shannon_Index_strata), !is.na(Mean_Dormice)) %>%
  group_by(Coppicing_Plot, Stratum) %>%
  summarise(
    Mean_Shannon = mean(Shannon_Index_strata, na.rm = TRUE),
    Mean_Dormice = mean(Mean_Dormice, na.rm = TRUE),
    .groups = "drop"
  )

# Order Coppicing_Plot by increasing Shannon to get sensible lines:
line_data <- line_data %>%
  arrange(Coppicing_Plot, Stratum) %>%
  mutate(Coppicing_Plot = factor(Coppicing_Plot, levels = unique(Coppicing_Plot[order(Mean_Shannon)])))

# Plotting
p <- ggplot(line_data, aes(x = Mean_Shannon, y = Mean_Dormice, color = Stratum, group = Stratum)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  labs(
    title = "Mean Dormice Abundance vs Mean Shannon Index by Stratum",
    x = "Mean Shannon Index",
    y = "Mean Dormice Abundance",
    color = "Stratum"
  ) +
  theme_minimal(base_size = 14) +
  scale_color_brewer(palette = "Set1") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

print(p)

library(ggplot2)
library(dplyr)

# Prepare plot data: one row per Coppicing_Plot × Stratum with Mean Shannon and Dormice
plot_data <- model_data_species_wide %>%
  filter(!is.na(Mean_Dormice), !is.na(Shannon_Index_strata)) %>%
  group_by(Coppicing_Plot, Stratum) %>%
  summarise(
    Mean_Shannon = mean(Shannon_Index_strata, na.rm = TRUE),
    Mean_Dormice = mean(Mean_Dormice, na.rm = TRUE),
    .groups = "drop"
  )

# Define nicer labels for strata
stratum_labels <- c(
  "Understorey" = "Understorey",
  "Shrub" = "Shrub layer",
  "Lower_Canopy" = "Lower canopy",
  "Upper_Canopy" = "Upper canopy",
  "Emergent" = "Emergent"
)

# Recode Stratum factor levels in plot_data
plot_data <- plot_data %>%
  mutate(
    Stratum = factor(Stratum, levels = names(stratum_labels), labels = stratum_labels)
  )


# Custom palette inspired by your example but a bit different:
custom_plot_colors <- c(
  "#C9E8FF",
  "#9fa8da",  # soft periwinkle (lighter blue)
  "#9b70d1",  # medium lavender purple
  "#d86b9c",  # light pink
  "#880e4f"   # deep magenta/dark pink
)

# If you have more than 5 plots, extend the palette smoothly
num_plots <- length(unique(plot_data$Stratum))
if (num_plots > length(custom_plot_colors)) {
  custom_plot_colors <- colorRampPalette(custom_plot_colors)(num_plots)
}

# Plot with updated colors:
p <- ggplot(plot_data, aes(x = Mean_Shannon, y = Mean_Dormice)) +
  geom_line(aes(group = Stratum, color = Stratum), size = 1) +
  geom_point(aes(color = Stratum), size = 2) +
  facet_wrap(~ Stratum, scales = "free_x") +
  scale_color_manual(values = custom_plot_colors) +
  labs(
    x = "Mean Shannon Index",
    y = "Mean Dormice Abundance"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(p)
ggsave("~/Library/CloudStorage/OneDrive-UniversityofBristol/Github/Thesis/SuppInfo/Figures/Chapter4/DMbyStratum.pdf", plot = p, width = 8, height = 6)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#### STEP 17: DORMICE PREFERENCE ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Set working directory
setwd("~/Library/CloudStorage/OneDrive-UniversityofBristol/Documents/PhD/Project 3_Bontuchel")

# Load data
model_data<-read.csv('Data/final_dataset.csv')
head(model_data)

# Summarise prop of each species
model_data_species <- model_data %>%
  group_by(Coppicing_Plot, Species) %>%
  summarise(
    Proportion_of_Species = sum(Proportion_of_Species, na.rm = TRUE),
    Mean_Dormice = first(Mean_Dormice),
    Coppicing_year = first(Coppicing_year),
    .groups = "drop"
  )

head(model_data_species)

#rescale so sum of prop of all species across each coppicing panel is never above 1
model_data_species_wide_scaled <- model_data_species %>%
  group_by(Coppicing_Plot) %>%
  mutate(
    total_prop = sum(Proportion_of_Species, na.rm = TRUE),
    Proportion_of_Species = Proportion_of_Species / total_prop
  ) %>%
  dplyr::select(-total_prop) %>%
  ungroup()

model_data<-model_data_species_wide_scaled
# Filter species:
# 1. Hazel (Corylus avellana): One of the most important food sources (hazelnuts) and good for nesting material.
# 2. Honeysuckle (Lonicera periclymenum): Provides nesting material (bark stripping), and nectar in summer.
# 3. Bramble (Rubus fruticosus agg.): Offers fruit in late summer and dense cover for protection and nesting.
# 4. Oak (Quercus robur, Quercus petraea): Supports many invertebrates (an indirect food source) and provides structure for arboreal movement.
# 5. Hawthorn (Crataegus monogyna) and Blackthorn (Prunus spinosa): Flowers in spring (nectar), fruits in autumn, and good structural habitat.
# 6. Sycamore (Acer pseudoplatanus): Often less preferred and can shade out understorey — sometimes considered negatively associated with dormice if dominant.

predictors= c("Species")
response_var=c("Mean_Dormice")

select_best_model(model_data, predictors,response_var,cor_threshold = 0.7,)

# Keep only required species:
key_species <- c(
  "Hazel",
  "Honeysuckle",
  "Bramble",
  "Oak",
  "Sycamore"
)

model_data_filtered <- model_data %>%
  filter(Species %in% key_species)

predictors= c("Species")
response_var=c("Mean_Dormice")

select_best_model(model_data_filtered, predictors,response_var,cor_threshold = 0.7,)

modlog<-glm(Mean_Dormice~Species, data=model_data_filtered, family="poisson")
summary(modlog)
shapiro.test(modlog$residuals)
# Calculate dispersion
dispersion <- deviance(modlog) / df.residual(modlog)
dispersion

modlog<-gam(Mean_Dormice~Species, data=model_data_filtered, family=nb())
summary(modlog)
shapiro.test(modlog$residuals)

# Multiply proportion by 100 if you want percentages (optional)
plot_data <- model_data_filtered %>%
  mutate(Proportion_of_Species = Proportion_of_Species * 100) %>%
  filter(!is.na(Mean_Dormice), !is.na(Proportion_of_Species)) %>%
  mutate(Species = factor(Species, levels = names(species_color_map)))


library(ggplot2)
library(dplyr)

# Prepare plot data
plot_data <- model_data_filtered %>%
  mutate(
    Proportion_of_Species = Proportion_of_Species * 100,
    Species = factor(Species, levels = names(species_color_map))
  ) %>%
  filter(!is.na(Mean_Dormice), !is.na(Proportion_of_Species))

# Define custom color palette for 7 species
custom_colors <- c(
  "#C5CAE9",  # light blue
  "#ce93d8",  # lavender
  "pink",     # pink
  "hotpink",  # hot pink
  "#c2185b" # dark pink 
)

# Create facetted dot-to-dot line plot
p2 <- ggplot(plot_data, aes(x = Proportion_of_Species, y = Mean_Dormice)) +
  geom_line(aes(group = Species, color = Species)) +
  geom_point(aes(color = Species), size = 2) +
  scale_color_manual(values = custom_colors) +
  labs(
    x = "Proportion of Species (%)",
    y = "Mean Dormouse Abundance"
  ) +
  facet_wrap(~ Species, scales = "free_x") +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",  # Remove legend since each panel is labeled
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

# Save as PDF
file_name <- "~/Library/CloudStorage/OneDrive-UniversityofBristol/Github/Thesis/SuppInfo/Figures/Chapter4/Key_species.pdf"
ggsave(file_name, plot = p2, width = 8, height = 5)




summary(sp_data$Avg_Shannon_Index)




# 1. Keep unique combinations
unique_combos <- model_data |> 
  dplyr::select(Species, Forest_layer, Coppicing_Plot, Taxonomic_Level) |> 
  dplyr::distinct()

# 2. Count number of rows by Taxonomic_Level
level_counts <- unique_combos |> 
  dplyr::count(Taxonomic_Level, name = "Count")

# 3. Add proportion column
level_counts <- level_counts |> 
  dplyr::mutate(Proportion = Count / sum(Count))

# 4. View result
print(level_counts)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#### STEP 17: management and structure ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Set working directory and load data
setwd("~/Library/CloudStorage/OneDrive-UniversityofBristol/Documents/PhD/Project 3_Bontuchel")
sp_data <- read.csv('Data/final_dataset.csv')

# Select only the relevant columns for uniqueness check
unique_rows <- !duplicated(sp_data[, c("Coppicing_Plot")])
sp_data_unique <- sp_data[unique_rows, ]

# Get recovery period
sp_data_unique <- sp_data_unique %>%
  mutate(recovery_period = 2020 - Coppicing_year)

predictors <- c("recovery_period")
response_var <- c("Mean_Canopy_Height_2020", "Perc_under_2m_2020", "Can_cover_2020", "Height_cv_2020", 
                    "Max_Canopy_Height_2020", "Gap_Proportion_2020", "LAI", "Intensity", "VCI", 
                    "Canopy_Volume_2020", "FT_5_10m", "FT_1.37m", "FT_10_20m")

data <- sp_data_unique

# cor_matrix <- cor(sp_data_unique[, response_var], use = "pairwise.complete.obs", method = "spearman")
# valid_combinations <- list(); counter <- 1
# cor_threshold = 0.7
# 
# for (i in (length(response_var):1)) {
#   combs <- combn(response_var, i, simplify = FALSE)
#   for (combo in combs) {
#     sub_cor <- cor_matrix[combo, combo]
#     if (all(abs(sub_cor[upper.tri(sub_cor)]) < cor_threshold)) {
#       valid_combinations[[counter]] <- combo
#       counter <- counter + 1
#     }
#   }
#   if (length(valid_combinations) > 0) break
# }
# 
# # Filter out combinations with more than one FT_ variable
# filtered_combinations <- Filter(function(combo) {
#   sum(grepl("^FT_", combo)) <= 1
# }, valid_combinations)
# 
# library(vegan)
# library(dplyr)
# 
# # Prepare a dataframe to store results
# results <- data.frame(
#   combo_num = integer(),
#   variables = character(),
#   pseudo_F = numeric(),
#   p_value = numeric(),
#   R2 = numeric(),
#   stringsAsFactors = FALSE
# )
# 
# for (i in seq_along(filtered_combinations)) {
#   vars <- filtered_combinations[[i]]
#   
#   # Prepare data subset and remove NA rows
#   temp_data <- sp_data_unique[, c(vars, "recovery_period")]
#   temp_data <- na.omit(temp_data)
#   
#   # Scale the response variables before distance calculation
#   response_scaled <- scale(temp_data[, vars])
#   
#   # Calculate Bray-Curtis distance on scaled data
#   euc_dist <- vegdist(response_scaled, method = "euclidean")
#   predictor <- temp_data$recovery_period
#   
#   # Run PERMANOVA
#   adonis_res <- adonis2(euc_dist ~ predictor, permutations = 999)
#   
#   # Extract summary stats
#   pseudo_F <- adonis_res$F[1]        # F-value for recovery_period
#   p_val <- adonis_res$`Pr(>F)`[1]   # p-value for recovery_period
#   R2 <- adonis_res$R2[1]             # R² for recovery_period
#   
#   # Append to results
#   results <- rbind(results, data.frame(
#     combo_num = i,
#     variables = paste(vars, collapse = ", "),
#     pseudo_F = pseudo_F,
#     p_value = p_val,
#     R2 = R2,
#     stringsAsFactors = FALSE
#   ))
# }
# 
# # Order results by R² (descending) or by p-value (ascending)
# results_sorted <- results %>% arrange(desc(R2))
# print(results_sorted)
# 
# for (i in seq_along(filtered_combinations)) {
#   vars <- filtered_combinations[[i]]
#   
#   # Prepare data subset and remove NA rows
#   temp_data <- sp_data_unique[, c(vars, "recovery_period")]
#   temp_data <- na.omit(temp_data)
#   
#   response_matrix <- temp_data[, vars]
#   predictor <- temp_data$recovery_period
#   
#   # Run PERMANOVA
#   adonis_res <- adonis2(response_matrix ~ predictor, permutations = 999, method = "euclidean")
#   
#   # Extract summary stats
#   pseudo_F <- adonis_res$F[1]        # F-value for recovery_period
#   p_val <- adonis_res$`Pr(>F)`[1]   # p-value for recovery_period
#   R2 <- adonis_res$R2[1]             # R² for recovery_period
#   
#   # Append to results
#   results <- rbind(results, data.frame(
#     combo_num = i,
#     variables = paste(vars, collapse = ", "),
#     pseudo_F = pseudo_F,
#     p_value = p_val,
#     R2 = R2,
#     stringsAsFactors = FALSE
#   ))
# }
# 
# # Order results by R² (descending) or by p-value (ascending)
# results_sorted <- results %>% arrange(desc(R2))
# print(results_sorted)
# #MANOVAs were a poor fit, assumptions violated
# 
# 
# # PERMANOVA #
# 
# 
# # Set working directory and load data
# setwd("~/Library/CloudStorage/OneDrive-UniversityofBristol/Documents/PhD/Project 3_Bontuchel")
# sp_data <- read.csv('Data/final_dataset.csv')
# 
# # Select only the relevant columns for uniqueness check
# unique_rows <- !duplicated(sp_data[, c("Coppicing_Plot")])
# sp_data_unique <- sp_data[unique_rows, ]
# 
# # Get recovery period
# sp_data_unique <- sp_data_unique %>%
#   mutate(recovery_period = 2020 - Coppicing_year)
# 
# predictors <- c("recovery_period")
# response_var <- c("Mean_Canopy_Height_2020", "Perc_under_2m_2020", "Can_cover_2020", "Height_cv_2020", 
#                   "Max_Canopy_Height_2020", "Gap_Proportion_2020", "LAI", "Intensity", "VCI", 
#                   "Canopy_Volume_2020", "FT_5_10m", "FT_1.37m", "FT_10_20m")
# 
# data <- sp_data_unique
# 
# cor_matrix <- cor(sp_data_unique[, response_var], use = "pairwise.complete.obs", method = "spearman")
# valid_combinations <- list(); counter <- 1
# cor_threshold = 0.7
# 
# for (i in (length(response_var):1)) {
#   combs <- combn(response_var, i, simplify = FALSE)
#   for (combo in combs) {
#     sub_cor <- cor_matrix[combo, combo]
#     if (all(abs(sub_cor[upper.tri(sub_cor)]) < cor_threshold)) {
#       valid_combinations[[counter]] <- combo
#       counter <- counter + 1
#     }
#   }
#   if (length(valid_combinations) > 0) break
# }
# 
# # Filter out combinations with more than one FT_ variable
# filtered_combinations <- Filter(function(combo) {
#   sum(grepl("^FT_", combo)) <= 1
# }, valid_combinations)
# 
# 
# # Prepare a dataframe to store results
# results <- data.frame(
#   combo_num = integer(),
#   variables = character(),
#   pseudo_F = numeric(),
#   p_value = numeric(),
#   R2 = numeric(),
#   stringsAsFactors = FALSE
# )
# 
# for (i in seq_along(filtered_combinations)) {
#   vars <- filtered_combinations[[i]]
#   
#   temp_data <- sp_data_unique[, c(vars, "recovery_period")]
#   temp_data <- na.omit(temp_data)
#   
#   # Scale the response variables
#   response_scaled <- scale(temp_data[, vars])
#   
#   # Compute Euclidean distance matrix
#   dist_matrix <- dist(response_scaled, method = "euclidean")
#   
#   # Run PERMANOVA with distance matrix
#   adonis_res <- adonis2(dist_matrix ~ recovery_period, data = temp_data, permutations = 999)
#   
#   pseudo_F <- adonis_res$F[1]
#   p_val <- adonis_res$`Pr(>F)`[1]
#   R2 <- adonis_res$R2[1]
#   
# library(vegan)
# 
#   results <- data.frame(
#     combo_num = integer(),
#     variables = character(),
#     pseudo_F = numeric(),
#     p_value = numeric(),
#     R2 = numeric(),
#     stringsAsFactors = FALSE
#   )
#   
#   for (i in seq_along(filtered_combinations)) {
#     vars <- filtered_combinations[[i]]
#     
#     temp_data <- sp_data_unique[, c(vars, "recovery_period")]
#     temp_data <- na.omit(temp_data)
#     
#     response_scaled <- scale(temp_data[, vars])
#     dist_matrix <- dist(response_scaled, method = "euclidean")
#     
#     adonis_res <- adonis2(dist_matrix ~ recovery_period, data = temp_data, permutations = 999)
#     
#     pseudo_F <- adonis_res$F[1]
#     p_val <- adonis_res$`Pr(>F)`[1]
#     R2 <- adonis_res$R2[1]
#     
#     results <- rbind(results, data.frame(
#       combo_num = i,
#       variables = paste(vars, collapse = ", "),
#       pseudo_F = pseudo_F,
#       p_value = p_val,
#       R2 = R2,
#       stringsAsFactors = FALSE
#     ))
#   }
# }
#   # Now order and print once
#   results_sorted <- results %>% arrange(desc(R2))
#   print(results_sorted)
#   
#   for (i in seq_along(filtered_combinations)) {
#     vars <- filtered_combinations[[i]]
#     temp_data <- sp_data_unique[, c(vars, "recovery_period")]
#     temp_data <- na.omit(temp_data)
#     response_scaled <- scale(temp_data[, vars])
#     dist_matrix <- dist(response_scaled, method = "euclidean")
#     
#     # Test homogeneity of dispersions
#     disp <- betadisper(dist_matrix, temp_data$recovery_period)
#     disp_test <- permutest(disp, permutations = 999)
#     
#     cat("Combo", i, "variables:", paste(vars, collapse = ", "), "\n")
#     print(disp_test)
#   }
#   
# #All combos have non-significant p-values here, so the assumption is not violated for any.
# 
# # choose combo 8 to match others: Combo 1: "Perc_under_2m_2020, Max_Canopy_Height_2020, Gap_Proportion_2020, Intensity, VCI, Canopy_Volume_2020"
# # Strongest R² (~0.08), significant PERMANOVA p = 0.005
# # Homogeneity test p = 0.38 (no dispersion issue)
# 
#   # Interpretation:
#   #   Overall effect of Combo 1 variables on multivariate response:
#   #   The PERMANOVA test p-value is 0.005, which is significant at the 0.05 level.
#   # This means the environmental variables in Combo 1 collectively explain a statistically significant amount of variation in the multivariate community or response data you are analyzing.
#   # The null hypothesis (that groups do not differ based on these variables) is rejected, indicating differences between groups.
#   # Effect size (R² = 0.08):
#   #   The R² value of 0.08 means that about 8% of the variation in the multivariate response is explained by this combination of predictors.
#   # While 8% might seem modest, ecological data are often complex and noisy, so explaining this proportion can still be meaningful.
#   # Assumption check: homogeneity of dispersions
#   # The permutation test for homogeneity of multivariate dispersions gave a p-value of 0.38, which is not significant.
#   # This means the spread or variability within groups is similar across groups.
#   # Therefore, the significant PERMANOVA result is unlikely to be due to differences in group variances (dispersion), but more likely due to actual differences in the centroid locations (group means in multivariate space).
#   # Biological/Practical implications:
#   #   Variables like Percentage of vegetation under 2m, maximum canopy height, gap proportion, canopy volume, and indices such as Intensity and VCI (vegetation condition index?) appear to influence community composition or your multivariate response.
#   # These structural and condition-related canopy variables might be important drivers or indicators of ecological variation in your sites/groups.
#   # This suggests management or monitoring could focus on these variables to understand or predict changes in community or ecosystem structure.
#   # 
#   
#   library(tidyverse)
#   library(scales)  # for rescale()
# 
#   # Select the needed variables plus recovery_period
#   df <- sp_data_unique %>%
#     dplyr::select(recovery_period, Can_cover_2020, Max_Canopy_Height_2020, Gap_Proportion_2020, Intensity, VCI, Canopy_Volume_2020)
#   
#   # Standardize (rescale to 0-1) each response variable
#   df_std <- df %>%
#     mutate(across(c(Can_cover_2020, Max_Canopy_Height_2020, Gap_Proportion_2020, Intensity, VCI, Canopy_Volume_2020), 
#                   ~ rescale(.x, to = c(0, 1)))) 
#   
#   # Pivot longer for ggplot
#   df_long <- df_std %>%
#     pivot_longer(cols = -recovery_period, names_to = "Variable", values_to = "Standardized_Value")
#   
#   # Plot
#   ggplot(df_long, aes(x = recovery_period, y = Standardized_Value, color = Variable)) +
#     geom_line(stat = "summary", fun = mean, size = 1) +  # mean line per variable
#     geom_point(stat = "summary", fun = mean) +           # mean points
#     labs(x = "Recovery Period", y = "Standardized Metric (0-1)", color = "Variable",
#          title = "Standardized Response Variables Over Recovery Period") +
#     theme_minimal() +
#     theme(legend.position = "bottom")
#   
# #one per plot
#   
#   # Assume df is your data frame, standardize the columns first
#   min_max_scale <- function(x) {
#     (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
#   }
#   
#   metrics <- c("Can_cover_2020", "Max_Canopy_Height_2020", "Gap_Proportion_2020", 
#                "Intensity", "VCI", "Canopy_Volume_2020")
#   
#   df_std <- df %>%
#     mutate(across(all_of(metrics), min_max_scale))
#   
#   # Now plot one graph per metric
#   plots <- list()
#   
#   for (metric in metrics) {
#     p <- ggplot(df_std, aes_string(x = "recovery_period", y = metric)) +
#       geom_point() +
#       geom_smooth(method = "loess", se = TRUE) +
#       labs(title = paste("Recovery Period vs", metric),
#            x = "Recovery Period",
#            y = paste("Standardized", metric)) +
#       theme_minimal()
#     
#     plots[[metric]] <- p
#   }
#   
#   # To print plots one by one:
#   for (metric in metrics) {
#     print(plots[[metric]])
#   }
# 
# 
# ### manova wouldnt work as a hump relationship is expected: 
#   
# #try pca with recovery: 
#   # 
#   # PCA reduces multivariate canopy structure into a few interpretable gradients.
#   # Fitting quadratic models to PCs lets you see how those gradients respond to recovery time.
#   # If PC1 captures overall "forest maturity", a hump in PC1 shows that maturity peaks at intermediate ages.
#   # 
#   
#   # Set working directory and load data
#   setwd("~/Library/CloudStorage/OneDrive-UniversityofBristol/Documents/PhD/Project 3_Bontuchel")
#   sp_data <- read.csv('Data/final_dataset.csv')
#   
#   # Select only the relevant columns for uniqueness check
#   unique_rows <- !duplicated(sp_data[, c("Coppicing_Plot")])
#   sp_data_unique <- sp_data[unique_rows, ]
#   
#   # Get recovery period
#   sp_data_unique <- sp_data_unique %>%
#     mutate(recovery_period = 2020 - Coppicing_year)
#   
#   predictors <- c("recovery_period")
#   response_var <- c("Mean_Canopy_Height_2020", "Perc_under_2m_2020", "Can_cover_2020", "Height_cv_2020", 
#                     "Max_Canopy_Height_2020", "Gap_Proportion_2020", "LAI", "Intensity", "VCI", 
#                     "Canopy_Volume_2020", "FT_5_10m", "FT_1.37m", "FT_10_20m")
#   
#   sp_data_unique[,response_var] <- scale(sp_data_unique[,response_var])
#   
#   pca <- prcomp(sp_data_unique[, response_var], center = FALSE, scale. = FALSE)
#   summary(pca)  # To see variance explained
#   
#   scores <- as.data.frame(pca$x)
#   scores$recovery_period <- sp_data_unique$recovery_period
#   
#   lm_PC1 <- lm(PC1 ~ recovery_period + I(recovery_period^2), data = scores)
#   lm_PC2 <- lm(PC2 ~ recovery_period + I(recovery_period^2), data = scores)
#   
#   summary(lm_PC1)
#   summary(lm_PC2)
#   
#   loadings <- pca$rotation[, 1:2]
#   round(loadings, 2)
#   
#   ggplot(scores, aes(x = recovery_period, y = PC1)) +
#     geom_point() +
#     geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
#     theme_minimal()
#   
# #obscured due to opposing hump shapes, only option is to do each variable seperately:
#   
#   results <- lapply(response_var, function(var) {
#     model <- lm(as.formula(paste0(var, " ~ recovery_period + I(recovery_period^2)")),
#                 data = sp_data_unique)
#     summary(model)
#   })
#   names(results) <- response_var
#   
#   results
#   
  
  # GAMS
  
  # Load packages
  library(mgcv)
  library(ggplot2)
  library(dplyr)
  
  # Set working directory and load data
  setwd("~/Library/CloudStorage/OneDrive-UniversityofBristol/Documents/PhD/Project 3_Bontuchel")
  sp_data <- read.csv('Data/final_dataset.csv')
  
  # Keep one row per Coppicing_Plot
  sp_data_unique <- sp_data[!duplicated(sp_data$Coppicing_Plot), ]
  
  # Add recovery period
  sp_data_unique <- sp_data_unique %>%
    mutate(recovery_period = 2020 - Coppicing_year)
  
  # Variables to model
  response_vars <- c("Max_Canopy_Height_2020", "Gap_Proportion_2020", 
                     "LAI", "Intensity", "VCI", "Canopy_Volume_2020")

  
  # Function to fit GAM and check assumptions
  fit_and_check_gam <- function(response) {
    formula <- as.formula(paste0(response, " ~ s(recovery_period, k = 4)"))
    model <- gam(formula, data = sp_data_unique, method = "REML")
    
    # Justifications for using k = 4:
    #   Sample size and model parsimony: You have only 56 observations (plots), so a smaller k (like 4) reduces risk of overfitting.
    # Biological or ecological reasoning: If the recovery process is expected to follow a smooth, gently curving pattern (e.g., linear or gently nonlinear), a lower k suffices.
    # Diagnostics: If basis dimension checks (e.g., k-index test) show no evidence that k is too low, and smooth terms don’t show complex shapes, k = 4 is reasonable.
    # Interpretability: Smoother curves are easier to interpret ecologically (less chance of fitting noise).
    
    cat("\n\n=== GAM for", response, "===\n")
    print(summary(model))
    
    vci_gamma<-gam(VCI ~ s(recovery_period, k=4), family=Gamma(link="log"), data = sp_data_unique)
    summary(vci_gamma)
    gam.check(vci_gamma)
    
    
    library(mgcv)
    library(gratia)
    library(ggplot2)
    library(patchwork)  # for combining plots
    
    # Fit models
    gam_models <- list(
      max_height = gam(Max_Canopy_Height_2020 ~ s(recovery_period, k = 4), data = sp_data_unique, method = "REML"),
      lai         = gam(LAI ~ s(recovery_period, k = 4), data = sp_data_unique, method = "REML"),
      vci         = gam(VCI ~ s(recovery_period, k = 4), family = Gamma(link = "log"), data = sp_data_unique, method = "REML"),
      intensity   = gam(Intensity ~ s(recovery_period, k = 4), data = sp_data_unique, method = "REML"),
      gap         = gam(Gap_Proportion_2020 ~ s(recovery_period, k = 4), data = sp_data_unique, method = "REML"),
      volume      = gam(Canopy_Volume_2020 ~ s(recovery_period, k = 4), data = sp_data_unique, method = "REML")
    )
    summary(gam_models$max_height)
    summary(gam_models$lai)
    summary(gam_models$vci)
    summary(gam_models$intensity)
    summary(gam_models$gap)
    summary(gam_models$volume)
    
    
    
    summary(gam_models$vci)
    gratia::smooth_estimates(gam_models$vci) %>% pull(.smooth)
    gam_models$gap <- mgcv::gam(gam_models$gap ~ s(recovery_period, k=10), data=sp_data_unique, family='gaussian')
    plot(sp_data_unique$recovery_period, sp_data_unique$gap)
    lines(sm$recovery_period, sm$fitted, col = "red", lwd = 2)
    
    
    
    plot_style <- function(model, y_label) {
      sm <- gratia::smooth_estimates(model)
      
      sm <- sm %>%
        filter(.smooth == "s(recovery_period)") %>%
        rename(
          recovery_period = recovery_period,
          fitted = .estimate,
          se = .se
        ) %>%
        mutate(
          lower = fitted - 2 * se,
          upper = fitted + 2 * se
        )
      
      ggplot(sm, aes(x = recovery_period)) +
        geom_ribbon(aes(ymin = lower, ymax = upper), fill = "red", alpha = 0.15) +
        geom_line(aes(y = fitted), color = "hotpink", size = 1.2) +
        labs(x = "Time since coppicing (years)", y = y_label) +
        theme_minimal() +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", size = 0.25),
          axis.text.x = element_text(hjust = 1),
          axis.title = element_text(size = 12)
        )
    }
    
    
    
    # Generate plots
    p1 <- plot_style(gam_models$max_height, "Max Canopy Height (m)")
    p2 <- plot_style(gam_models$lai, "Leaf Area Index")
    p3 <- plot_style(gam_models$vci, "VCI (log scale)")
    p4 <- plot_style(gam_models$intensity, "Intensity")
    p5 <- plot_style(gam_models$gap, "Gap Proportion")
    p6 <- plot_style(gam_models$volume, "Canopy Volume (m³)")
    
    # Combine into 2-column panel
    final_plot <- (p1 | p2) / (p3 | p4) / (p5 | p6)
    
    # Print
    print(final_plot)
    
    pdf("~/Library/CloudStorage/OneDrive-UniversityofBristol/GitHub/Thesis/Chapter4/Figures/Structure.pdf", width = 7, height = 7)
    print(final_plot)
    dev.off()
    
    par(mar = c(5, 16, 4, 2)) 
    
    
    
    
    
    
    
    
    
    
    # Check residuals (normality, patterns)
    par(mfrow = c(2, 2))
    plot(model, residuals = TRUE, pch = 16, main = response)
    gam.check(model)
    shapiro_test <- shapiro.test(residuals(model))
    cat("Shapiro-Wilk test p-value:", shapiro_test$p.value, "\n")
    
    # Reset plotting window
    par(mfrow = c(1, 1))
    
    return(model)
  }
  
  # Run for all variables
  gam_models <- lapply(response_vars, fit_and_check_gam)
  names(gam_models) <- response_vars
  #lai vc and gap prop are sig with shap <0.001
  
  plot(gam_models$LAI, residuals=TRUE, pch=16, main="LAI ~ s(recovery_period)")
  
  gam_models[1]
  
  hist(sp_data_unique$Max_Canopy_Height_2020)
  hist(sp_data_unique$Gap_Proportion_2020)
  hist(sp_data_unique$Intensity)
  hist(sp_data_unique$VCI)
  hist(sp_data_unique$Canopy_Volume_2020)
  
  # For Gap_Proportion_2020
  gam_gap <- gam(Gap_Proportion_2020 ~ s(recovery_period, k=4), 
                 family = Gamma(link = "log"), 
                 data = sp_data_unique, method = "REML")
  
  # For Intensity
  gam_intensity <- gam(Intensity ~ s(recovery_period, k=4), 
                       family = Gamma(link = "log"), 
                       data = sp_data_unique, method = "REML")
  
  # For VCI
  gam_vci <- gam(VCI ~ s(recovery_period, k=4), 
                 family = Gamma(link = "log"), 
                 data = sp_data_unique, method = "REML")
  #do gamma ones, see if it helps because residuals were off
  summary(gam_gap)
  summary(gam_intensity)
  summary(gam_vci)
  
  # Plot residuals to check fit
  par(mfrow=c(3,2))
  
  plot(gam_gap, residuals=TRUE, pch=16, cex=0.5)
  gam.check(gam_gap)
  
  plot(gam_intensity, residuals=TRUE, pch=16, cex=0.5)
  gam.check(gam_intensity)
  
  plot(gam_vci, residuals=TRUE, pch=16, cex=0.5)
  gam.check(gam_vci)
  
  #plots
  # Keep a copy of original responses before scaling
  sp_data_scaled <- sp_data_unique  # for plotting only
  # Scale response vars only in the plotting dataset
  sp_data_scaled[response_vars] <- scale(sp_data_scaled[response_vars])
  
  
  get_gam_plot_data_with_partial_resids <- function(model, response_name) {
    y_raw <- sp_data_unique[[response_name]]
    y_scaled <- sp_data_scaled[[response_name]]
    
    # Predict smooth fit and SE on response scale
    pred <- predict(model, se.fit = TRUE)
    fitted <- pred$fit
    se <- pred$se.fit
    
    y_mean <- mean(y_raw, na.rm = TRUE)
    y_sd <- sd(y_raw, na.rm = TRUE)
    
    standardised_fitted <- (fitted - y_mean) / y_sd
    standardised_lower <- (fitted - 2 * se - y_mean) / y_sd
    standardised_upper <- (fitted + 2 * se - y_mean) / y_sd
    
    # Partial residuals for the smooth term
    term_pred <- predict(model, type = "terms", se.fit = TRUE)
    smooth_term <- term_pred$fit[, "s(recovery_period)"]
    partial_residuals <- smooth_term + residuals(model)
    
    standardised_partial_resids <- (partial_residuals - y_mean) / y_sd
    
    tibble(
      recovery_period = sp_data_unique$recovery_period,
      response = y_scaled,
      fitted = standardised_fitted,
      lower = standardised_lower,
      upper = standardised_upper,
      partial_resid = standardised_partial_resids
    )
  }
  
  plot_data_lai <- get_gam_plot_data_with_partial_resids(gam_models$LAI, "LAI")
  plot_data_vci <- get_gam_plot_data_with_partial_resids(gam_models$VCI, "VCI")
  plot_data_gapprop <- get_gam_plot_data_with_partial_resids(gam_models$Gap_Proportion_2020, "Gap_Proportion_2020")
  
  plot_style <- function(data, y_label) {
    ggplot(data, aes(x = recovery_period, y = response)) +
      geom_point(color = "pink", size = 2, alpha = 0.6) +
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "red", alpha = 0.15) +
      geom_line(aes(y = fitted), color = "hotpink", size = 1.2) +
      theme_minimal() +
      labs(x = "Recovery Period (years)", y = y_label) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", size = 0.25),
        axis.text.x = element_text(hjust = 1),
        axis.title = element_text(size = 12)
      )
  }
  
  plot_lai <- plot_style(plot_data_lai, "Standardised LAI")
  plot_vci <- plot_style(plot_data_vci, "Standardised VCI")
  plot_gapprop <- plot_style(plot_data_gapprop, "Standardised Gap Proportion")
  
  # Print or arrange them as you want:
  plot_lai
  plot_vci
  plot_gapprop
  
  
  # all in one:
  
  library(dplyr)
  library(ggplot2)
  
  # Add a column identifying the response variable
  plot_data_lai$variable <- "LAI"
  plot_data_vci$variable <- "VCI"
  plot_data_gapprop$variable <- "Gap Proportion"
  
  # Combine into one dataframe
  combined_data <- bind_rows(plot_data_lai, plot_data_vci, plot_data_gapprop)
  
  # Define colors
  colors <- c("LAI" = "blue", "VCI" = "pink", "Gap Proportion" = "purple")
  
  # Plot
  ggplot(combined_data, aes(x = recovery_period, y = response, color = variable, fill = variable)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
    geom_line(aes(y = fitted), size = 1.2) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    labs(
      x = "Recovery Period (years)",
      y = "Standardised Response",
      color = "Variable",
      fill = "Variable"
    ) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", size = 0.25),
      axis.text.x = element_text(hjust = 1),
      axis.title = element_text(size = 12)
    )
  
  
  
  ### COMBINED PLOT ####
  
  library(magick)
  
  # Set working directory if needed
  setwd("~/Library/CloudStorage/OneDrive-UniversityofBristol/Github/Thesis/Chapter4")
  
  # Correct relative file paths (no leading slash)
  pdf_files <- c("Figures/CoeffsCanopyStructure.pdf", 
                 "Figures/CoeffsCanopyStructure_DM.pdf", 
                 "Figures/CoeffsTopo.pdf", 
                 "Figures/CoeffsTopoDM.pdf")
  
  # Read first page of each PDF
  images <- lapply(pdf_files, function(file) {
    image_read_pdf(file, density = 150)[1]
  })
  
  # OPTIONAL: Crop with margins preserved (adjust geometry as needed)
  # Format: "widthxheight+x_offset+y_offset"
  images <- lapply(images, function(img) {
    image_crop(img, geometry = "980x980+40+40")# Keep tight control
  })
  # Combine images in 2 rows, 2 columns (2x2 grid)
  row1 <- image_append(c(images[[1]], images[[2]]), stack = FALSE)  # horizontally append first two
  row2 <- image_append(c(images[[3]], images[[4]]), stack = FALSE)  # horizontally append next two
  
  # Stack the two rows vertically
  final_image <- image_append(c(row1, row2), stack = TRUE)
  
  # Save the combined image as a PDF
  image_write(final_image, path = "Figures/Combined_2x2.pdf", format = "pdf")
  #tidy this plot up, tried 10x10 and 5x5 maybe something in between
  
  
  
  
  
  
  
  # Set working directory if needed
  setwd("~/Library/CloudStorage/OneDrive-UniversityofBristol/Github/Thesis/Chapter4")
  
  # Correct relative file paths (no leading slash)
  pdf_files <- c("Figures/Management.pdf", 
                 "Figures/ManagementDM.pdf")
  
  # Read first page of each PDF
  images <- lapply(pdf_files, function(file) {
    image_read_pdf(file, density = 150)[1]
  })
  
  # OPTIONAL: Crop with margins preserved (adjust geometry as needed)
  # Format: "widthxheight+x_offset+y_offset"

  # Combine images in 2 rows, 2 columns (2x2 grid)
  row1 <- image_append(images[[1]])  # horizontally append first two
  row2 <- image_append(images[[2]])  # horizontally append next two
  
  # Stack the two rows vertically
  final_image <- image_append(c(row1, row2), stack = TRUE)
  
  # Save the combined image as a PDF
  image_write(final_image, path = "Figures/ManagementCombo.pdf", format = "pdf")
  #tidy this plot up, tried 10x10 and 5x5 maybe something in between
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  
  #### STEP 17: BIAS ####
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
  
  
  # how well did each sample represent each plot
  
  # Load data
  sp_data<-read.csv('final_dataset.csv')
  ET_points<-vect('ET_Points.gpkg')
  bon_CHM_2020<- rast("bon_CHM_2020_rep.tif")
  bontuchel<-vect('Bontuchel_Coppice_Panels.shp')
  
  ET_points<-crop(ET_points, bontuchel)
  
  # Create a buffer of 15 meters around the ET points
  buffer_15m <- buffer(ET_points, width = 15)  # Buffer by 15 meters
  buffer_15m$ET_Rec_Date<-buffer_15m$`_CREATION_DATE`
  
  # Check if the extents overlap
  ext(bon_CHM_2020)
  st_bbox(buffer_15m)
  
  plot(bon_CHM_2020)
  plot(bontuchel, border='pink', add=TRUE)
  plot(buffer_15m, add=TRUE, col='hotpink')
  plot(ET_points, add=TRUE)
  
  # 1. Normalize datetime format in both ET_points and buffer_15m
  ET_points$ET_Rec_Date <- as.POSIXct(ET_points$`_CREATION_DATE`, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  buffer_15m$ET_Rec_Date <- as.POSIXct(buffer_15m$'_CREATION_DATE', format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  
  # 2. Ensure buffer_15m and ET_points match by ET_Rec_Date
  # Safe inner join based on ET_Rec_Date to get the same order
  buffer_15m <- buffer_15m[order(buffer_15m$ET_Rec_Date), ]
  ET_points  <- ET_points[order(ET_points$ET_Rec_Date), ]
  
  # Confirm matching
  stopifnot(all(buffer_15m$ET_Rec_Date == ET_points$ET_Rec_Date))
  
  # 3. Now do the spatial join to get coppicing plot IDs
  matches <- terra::extract(bontuchel["ID_No"], ET_points)
  
  # 4. Assign matched IDs to buffer_15m
  buffer_15m$Coppicing_Plot <- matches$ID_No
  
  # Check if the conversion worked
  head(sp_data$ET_Rec_Date)
  
  summary(as.numeric(sp_data$ET_Canopy_Height))
  summary(sp_data$Mean_Canopy_Height_2020)
  
  # Select only the relevant columns for uniqueness check - prevents pseudoreplication 
  unique_rows <- !duplicated(sp_data[, c("ET_ID")])
  
  # Keep only unique rows in the dataframe
  sp_data_unique <- sp_data[unique_rows, ]
  
  length(sp_data_unique$ET_ID)
  
  # Check structure
  names(sp_data)
  names(buffer_15m)
  names(bontuchel)

  bon_Slope <- rast('bon_Slope.tif')
  bon_TPI<- rast('bon_TPI.tif')
  bon_TRI<- rast('bon_TRI.tif')
  bon_Aspect<- rast('bon_Aspect.tif')
  bon_Roughness<- rast('bon_Roughness.tif')
  bon_flowdir<- rast('bon_Flowdir.tif')
  bon_planc<- rast('bon_Planc.tif')
  bon_profc<- rast('bon_Profc.tif')
  bon_DTM_2020_hillshade<- rast('bon_Hillshade.tif')
  bon_meanc<- rast('bon_Meanc.tif')
  bon_TWI<- rast('bon_TWI.tif')
  bon_solar<- rast('bon_Solar.tif')
  
  
  # Create data frame to store data
  plot_data <- data.frame(matrix(NA, nrow=dim(buffer_15m)[1], ncol=26))
  names(plot_data)<-c("Coppicing_Plot","radius_plotvalue_ID","ET_Rec_Date", "radius_Mean_Canopy_Height_2020", "radius_Perc_under_2m_2020", 
                      "radius_Can_cover_2020", "radius_Height_cv_2020",
                      "radius_Elevation","radius_Slope","radius_TRI","radius_Max_Canopy_Height_2020", "radius_FT_5_10m","radius_FT_1.37m",
                      "radius_FT_10_20m","radius_Intensity","radius_Aspect","radius_TPI", 'radius_Roughness', "radius_LAI", "radius_VCI", "radius_Canopy_Volume_2020",
                      "radius_FlowDir", "radius_Plane_Curve", "radius_Profile_Curve", "radius_Mean_Curve", "radius_Solar_Radiation")
  
  plot_data$ET_Rec_Date <- buffer_15m$ET_Rec_Date
  
  # Run loop across each plot and extract data
  for (i in 1:dim(plot_data)[1]){

    plot_data$Coppicing_Plot[i] <- as.character(buffer_15m$Coppicing_Plot[i])

    #//////////#
    # CHM 2020 #
    #//////////#
    
    # Crop
    bon_2020_CHM_crop<-crop(bon_CHM_2020, buffer_15m[i,])
    
    #--------------------------------#
    # Mean Top of Canopy Height 2020 #
    #--------------------------------#
    plot_data[i,"radius_Mean_Canopy_Height_2020"]<-global(bon_2020_CHM_crop, fun='mean', na.rm=TRUE) 
    
    #-----------------------------------------#
    # Coefficient of variation in Height 2020 #
    #-----------------------------------------#
    # Coefficient of variation 2020
    plot_data[i,"radius_Height_cv_2020"]<-global(bon_2020_CHM_crop, fun='sd', na.rm=TRUE)/global(bon_2020_CHM_crop, fun='mean', na.rm=TRUE) 
    
    #----------------------------------#
    # % pixels below 2m threshold 2020 #
    #----------------------------------#
    plot_data[i,"radius_Perc_under_2m_2020"]<-100*(length(which(as.vector(na.omit(as.vector(bon_2020_CHM_crop)))<2))/length(as.vector(na.omit(as.vector(bon_2020_CHM_crop))))) 
    
    #----------------------------------#
    # % pixels above 2m threshold 2020 #
    #----------------------------------#
    plot_data[i,"radius_Can_cover_2020"]<-100- plot_data[i,"radius_Perc_under_2m_2020"]
    
    #------------------------#
    # Max Canopy Height 2020 #
    #------------------------#
    values_MCH_crop_2020 <- values(bon_2020_CHM_crop, na.rm = TRUE)
    Max_Canopy_Height_2020<-stats::quantile(values_MCH_crop_2020, probs=0.99, na.rm = TRUE)
    plot_data[i,"radius_Max_Canopy_Height_2020"]<-Max_Canopy_Height_2020
    
    #---------------#
    # Canopy volume #
    #---------------#
    plot_data[i,"radius_Canopy_Volume_2020"]<-canopy_volume(bon_2020_CHM_crop)
    
    #/////////////#
    # Point Cloud #
    #/////////////#
    
    # Convert the shapefile to a spatial object compatible with lidR
    polygon_ <- as(buffer_15m[i,], "Spatial")
    target <- st_crs(27700)  # Define the target CRS using sf
    
    # Convert the target CRS to sp's CRS format
    target_sp <- CRS(st_as_text(target))
    
    # Transform the polygon_lidar using spTransform()
    polygon_ <- spTransform(polygon_, target_sp)
    plot(polygon_,add=TRUE, col='yellow')
    # Crop point cloud
    fr_pt_crop<-clip_roi(first_returns, polygon_)
    
    #-------------------------------------------------------------------#
    # Calculate the percentage of first returns between 5 and 10 meters #
    #-------------------------------------------------------------------#
    
    # Filter
    first_returns_5_10m <- filter_poi(fr_pt_crop,(Z >= 5 & Z <= 10))
    
    # Get the total number of first returns and those between 5 and 10 meters
    total_first_returns <- nrow(fr_pt_crop)
    first_returns_5_10m_count <- nrow(first_returns_5_10m)
    
    # Calculate the percentage
    percent_5_10m <- (first_returns_5_10m_count / total_first_returns) * 100
    plot_data[i,"radius_FT_5_10m"]<-percent_5_10m
    
    #-------------------------------------------------------------#
    # Calculate the percentage of first returns above 1.37 meters #
    #-------------------------------------------------------------#
    
    # Filter
    first_returns_1.37m <- filter_poi(fr_pt_crop,(Z >= 1.37))
    
    # Get the total number of first returns and those above 1.37 meters
    first_returns_1.37m_count <- nrow(first_returns_1.37m)
    
    # Calculate the percentage
    percent_1.37m <- (first_returns_1.37m_count / total_first_returns) * 100
    plot_data[i,"radius_FT_1.37m"]<-percent_1.37m
    
    #-------------------------------------------------------------#
    # Calculate the percentage of first returns above 10-20m meters #
    #-------------------------------------------------------------#
    
    # Calculate the percentage of first returns 10-20m meters
    first_returns_10_20m <- filter_poi(fr_pt_crop,(Z >= 10 & Z <= 20))
    
    # Get the total number of first returns and those above 1.37 meters
    first_returns_10_20m_count <- nrow(first_returns_10_20m)
    
    # Calculate the percentage
    percent_10_20m <- (first_returns_10_20m_count / total_first_returns) * 100
    plot_data[i,"radius_FT_10_20m"]<-percent_10_20m
    
    #-----------#
    # Intensity #
    #-----------#
    
    # Crop
    pt_crop<-clip_roi(pt_cloud, polygon_)
    
    # Mean intensity
    plot_data[i,"radius_Intensity"] <- mean(pt_crop$Intensity)
    
    #-----#
    # LAI #
    #-----#
    
    # Save the cropped point cloud to a temporary file
    temp_file <- tempfile(fileext = ".laz")
    writeLAS(pt_crop, temp_file)
    
    # Calculate Leaf Area Density (LAD) from voxelization
    VOXELS_LAD <- lad.voxels(temp_file, grain.size = 2)
    
    # Calculate the LAD profile
    lad_profile <- lad.profile(VOXELS_LAD)
    
    # Calculate LAI derived from LAD profile
    lidar_lai <- lai(lad_profile)
    
    # Add LAI to dataframe
    plot_data[i, "radius_LAI"] <- lidar_lai
    
    #-----#
    # VCI #
    #-----#
    
    
    # Calculate Vertical Complexity Index (VCI)
    vci <- VCI(pt_crop@data$Z, by = 1, zmax = 20)
    
    # Add VCI to dataframe
    plot_data[i, "radius_VCI"] <- vci
    
    #/////////#
    # Terrain #
    #/////////#
    
    #-----------#
    # Elevation #
    #-----------#
    
    bon_DTM_2020_crop<-crop(bon_DTM_2020,buffer_15m[i,])
    plot_data[i,"radius_Elevation"]<-global(bon_DTM_2020_crop, fun='mean', na.rm=TRUE)
    
    #--------#
    # Aspect #
    #--------#
    
    Aspect_crop<-crop(bon_Aspect,buffer_15m[i,])
    plot_data[i,"radius_Aspect"]<-global(Aspect_crop, fun='mean', na.rm=TRUE)
    
    #-----#
    # TPI #
    #-----#
    
    TPI_crop<-crop(bon_TPI,buffer_15m[i,])
    plot_data[i,"radius_TPI"]<-global(TPI_crop, fun='mean', na.rm=TRUE)
    
    #-------#
    # Slope #
    #-------#
    
    Slope_crop<-crop(bon_Slope,buffer_15m[i,])
    plot_data[i,"radius_Slope"]<-global(Slope_crop, fun='mean', na.rm=TRUE)
    
    #-----#
    # TRI #
    #-----#
    TRI_crop<-crop(bon_TRI,buffer_15m[i,])
    plot_data[i,"radius_TRI"]<-global(TRI_crop, fun='mean', na.rm=TRUE)
    
    #-----------#
    # Roughness #
    #-----------#
    
    roughness_crop<-crop(bon_Roughness,buffer_15m[i,])
    plot_data[i,"radius_Roughness"]<-global(roughness_crop, fun='mean', na.rm=TRUE)
    
    #---------#
    # FlowDir #
    #---------#
    
    flowdir_crop<-crop(bon_flowdir,buffer_15m[i,])
    plot_data[i,"radius_FlowDir"]<-global(flowdir_crop, fun='mean', na.rm=TRUE)
    
    #-----------------#
    # Plane Curvature #
    #-----------------#
    
    planc_crop<-crop(bon_planc,buffer_15m[i,])
    plot_data[i,"radius_Plane_Curve"]<-global(planc_crop, fun='mean', na.rm=TRUE)
    
    #-------------------#
    # Profile Curvature #
    #-------------------#
    
    profc_crop<-crop(bon_profc,buffer_15m[i,])
    plot_data[i,"radius_Profile_Curve"]<-global(profc_crop, fun='mean', na.rm=TRUE)
    
    #----------------#
    # Mean Curvature #
    #----------------#
    
    meanc_crop<-crop(bon_meanc,buffer_15m[i,])
    plot_data[i,"radius_Mean_Curve"]<-global(meanc_crop, fun='mean', na.rm=TRUE)
    
    #-----------------#
    # Solar Radiation #
    #-----------------#
    
    solar_crop<-crop(bon_solar,buffer_15m[i,])
    plot_data[i,"radius_Solar_Radiation"]<-global(solar_crop, fun='mean', na.rm=TRUE)
    
    # Give each plotvalue a unique ID
    plot_data[i, "radius_plotvalue_ID"] <- paste0("pv_", i)
    
    # Progress (% of plots completed)
    print((i/dim(plot_data)[1])*100)
    
  }
  
  # Minimum DTM Value (for Height above lowest point, HALP)
  min_DTM <- min(values(bon_DTM_2020), na.rm = TRUE)
  plot_data$radius_HALP<-plot_data$radius_Elevation-min_DTM
  
  # Aspect (cos transformed)
  plot_data$radius_Aspect_Cos<-cos(plot_data$radius_Aspect)
  
  #------------#
  # Check data #
  #------------#
  
  summary(plot_data)
  
  # Now join together by coppicing plot and ET_Rec_Date
  
  # Ensure ET_Rec_Date is in POSIXct format in both
  sp_data$ET_Rec_Date <- as.POSIXct(sp_data$ET_Rec_Date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  plot_data$ET_Rec_Date <- as.POSIXct(plot_data$ET_Rec_Date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  
  # Inner join by Coppicing_plot and ET_Rec_Date
  combined_data <- inner_join(sp_data, plot_data, by = c("Coppicing_Plot", "ET_Rec_Date"))
  
  sum(is.na(sp_data$Coppicing_Plot))      # Should be 0
  sum(is.na(sp_data$ET_Rec_Date))         # Should be 0
  sum(is.na(plot_data$Coppicing_Plot))    # Should be 0
  sum(is.na(plot_data$ET_Rec_Date))       # Should be 0
  
  class(sp_data$ET_Rec_Date)
  class(plot_data$ET_Rec_Date)
  
  #--------------------#
  # Mean Canopy Height #
  #--------------------#

  # Select only the relevant columns for uniqueness check
  unique_rows <- !duplicated(combined_data[, c("ET_ID")])
  sp_data_unique <- combined_data[unique_rows, ]
  
  
  # Calculate the residuals (difference between the two height measurements)
  residuals <- sp_data_unique$Mean_Canopy_Height_2020 - sp_data_unique$radius_Mean_Canopy_Height_2020
  hist(residuals, main = "Residuals (Ground - LiDAR)", xlab = "Residuals", col = "pink")
  summary(residuals)
  shapiro.test(residuals) 
  
  # Assumptions  met
  
  #/////////////#
  # Correlation #
  #/////////////#
  
  # Spearmans correlation coefficient
  
  combined_data$radius_Can_cover_2020<-as.numeric(combined_data$radius_Can_cover_2020)
  height<-cor.test(combined_data$Can_cover_2020, combined_data$radius_Can_cover_2020, use = "complete.obs")
  height
# so mean canopy height for points sampled and plots overall are correlated
  
  #--------------------#
  # Mean Canopy Cover #
  #--------------------#

  hist(combined_data$radius_Can_cover_2020)
  hist(combined_data$Can_cover_2020)
  
  # Calculate the residuals (difference between the two height measurements)
  residuals <- combined_data$Can_cover_2020 - combined_data$radius_Can_cover_2020
  hist(residuals, main = "Residuals (Ground - LiDAR)", xlab = "Residuals", col = "pink")
  
  # Summary of residuals
  summary(residuals)
  shapiro.test(residuals) 
  
  # Assumptions notmet
  
  cor.test(combined_data$Can_cover_2020, combined_data$radius_Can_cover_2020,method='spearman', use = "complete.obs")
  
  # so mean canopy height for points sampled and plots overall are correlated
  
  #--------------------#
  # VCI  #
  #--------------------#
  
  hist(combined_data$radius_VCI)
  hist(combined_data$VCI)
  
  # Calculate the residuals (difference between the two height measurements)
  residuals <- combined_data$VCI - combined_data$radius_VCI
  hist(residuals, main = "Residuals (Ground - LiDAR)", xlab = "Residuals", col = "pink")
  
  # Summary of residuals
  summary(residuals)
  shapiro.test(residuals) 
  
  # Assumptions not met
  
  cor.test(combined_data$VCI, combined_data$radius_VCI,method='spearman', use = "complete.obs")
  
  # so vci for points sampled and plots overall are correlated
  
  #--------------------#
  # Intensity #
  #--------------------#
  
  hist(combined_data$Intensity)
  hist(combined_data$radius_Intensity)
  
  # Calculate the residuals (difference between the two height measurements)
  residuals <- combined_data$Intensity - combined_data$radius_Intensity
  hist(residuals, main = "Residuals (Ground - LiDAR)", xlab = "Residuals", col = "pink")
  
  # Summary of residuals
  summary(residuals)
  shapiro.test(residuals) 
  
  # Assumptions not met
  
  cor.test(combined_data$Intensity, combined_data$radius_Intensity,method='spearman', use = "complete.obs")
  
  # so intensity for points sampled and plots overall are correlated
  
  #--------------------#
  # gap Proportion #
  #--------------------#
  
  # Looking at gaps at 10m
  threshold2020<- 10
  size2020 <- c(2, 20000)
  
  # Initialize gap proportion column
  combined_data$radius_Gap_Proportion_2020 <- NA
  
  # Loop through each row in bontuchel
  for (i in seq_len(nrow(buffer_15m))) {
    
    # Extract plot ID
    current_id <- buffer_15m$Coppicing_Plot[i]
    
    # Crop CHM to the plot extent
    bon_2020_CHM_crop <- crop(bon_CHM_2020, buffer_15m[i, ])
    
    # Identify gaps within the cropped area
    bon_gaps_2020 <- getForestGaps(chm_layer = bon_2020_CHM_crop, threshold = threshold2020, size = size2020)
    
    # Calculate total plot area (excluding NA values)
    plot_area <- sum(!is.na(values(bon_2020_CHM_crop))) * res(bon_2020_CHM_crop)[1] * res(bon_2020_CHM_crop)[2]
    
    # Calculate total gap area
    gap_area <- sum(values(bon_gaps_2020) > 0, na.rm = TRUE) * res(bon_gaps_2020)[1] * res(bon_gaps_2020)[2]
    
    # Calculate the proportion of gaps
    proportion_gaps_2020 <- (gap_area / plot_area)*100
    
    # Find the corresponding row in sp_data
    sp_data_index <- which(combined_data$Coppicing_Plot == current_id)
    
    # Update sp_data with the calculated proportions
    combined_data$radius_Gap_Proportion_2020[sp_data_index] <- proportion_gaps_2020
    
    # Print progress
    print(paste0("Completed ", round((i / nrow(bontuchel)) * 100, 2), "%"))
  }
  
  hist(combined_data$Gap_Proportion_2020)
  hist(combined_data$radius_Gap_Proportion_2020)
  
  # Calculate the residuals (difference between the two height measurements)
  residuals <- combined_data$Gap_Proportion_2020 - combined_data$radius_Gap_Proportion_2020
  hist(residuals, main = "Residuals (Ground - LiDAR)", xlab = "Residuals", col = "pink")
  
  # Summary of residuals
  summary(residuals)
  shapiro.test(residuals) 
  
  # Assumptions not met
  
  cor.test(combined_data$Gap_Proportion_2020, combined_data$radius_Gap_Proportion_2020,method='spearman', use = "complete.obs")
  
  # so gap prop for points sampled and plots overall are correlated
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #------#
  # Plot #
  #------#
  
  length(sp_data_unique$ET_Canopy_Cover)
  length(sp_data_unique$LiDAR_Mean_Canopy_cover_2020)
  
  # Open a PNG device
  pdf("~/Library/CloudStorage/OneDrive-UniversityofBristol/GitHub/Thesis/Chapter4/Figures/LidarVsGround.pdf", width = 4.5, height = 7)
  
  par(mfrow = c(2, 1),         # 2 rows, 1 column
      mar = c(4, 4, 1, 1),     # inner margins (bottom, left, top, right)
      oma = c(0, 0, 0, 0))     # outer margins (bottom, left, top, right)
  
  # First scatter plot
  plot(sp_data_unique$LiDAR_Mean_Canopy_Height_2020, sp_data_unique$ET_Canopy_Height,
       xlab = "LiDAR Measured Height", 
       ylab = "Ground Measured Height", 
       xlim = c(0, max(sp_data_unique$ET_Canopy_Height, sp_data_unique$LiDAR_Mean_Canopy_Height_2020)),
       ylim = c(0, max(sp_data_unique$ET_Canopy_Height, sp_data_unique$LiDAR_Mean_Canopy_Height_2020)),
       bty = "l")  # Use 'l' to show only left and bottom box lines (axes)
  abline(a = 0, b = 1, col = "hotpink")
  
  # Second scatter plot
  plot(sp_data_unique$LiDAR_Mean_Canopy_cover_2020, sp_data_unique$ET_Canopy_Cover,
       xlab = "LiDAR Measured Canopy Cover", 
       ylab = "Ground Measured Canopy Cover", 
       xlim = c(0, max(sp_data_unique$ET_Canopy_Cover, sp_data_unique$LiDAR_Mean_Canopy_cover_2020)),
       ylim = c(0, max(sp_data_unique$ET_Canopy_Cover, sp_data_unique$LiDAR_Mean_Canopy_cover_2020)),
       bty = "l")  # Again, only left and bottom box lines
  abline(a = 0, b = 1, col = "hotpink")
  
  # Close PNG
  dev.off()
  
  #--------------------#
  # Canopy_Volume_2020 #
  #--------------------#
  
  hist(combined_data$Canopy_Volume_2020)
  hist(combined_data$radius_Canopy_Volume_2020)
  
  # Calculate the residuals (difference between the two height measurements)
  residuals <- combined_data$Canopy_Volume_2020 - combined_data$radius_Canopy_Volume_2020
  hist(residuals, main = "Residuals (Ground - LiDAR)", xlab = "Residuals", col = "pink")
  
  # Summary of residuals
  summary(residuals)
  shapiro.test(residuals) 
  
  # Assumptions not met
  
  cor.test(combined_data$Canopy_Volume_2020, combined_data$radius_Canopy_Volume_2020,method='spearman', use = "complete.obs")
  
  # so Canopy_Volume_2020 for points sampled and plots overall are correlated
  
  #--------------------#
  # Perc_under_2m_2020 #
  #--------------------#
  
  hist(combined_data$Perc_under_2m_2020)
  hist(combined_data$radius_Perc_under_2m_2020)
  
  # Calculate the residuals (difference between the two height measurements)
  residuals <- combined_data$Perc_under_2m_2020 - combined_data$radius_Perc_under_2m_2020
  hist(residuals, main = "Residuals (Ground - LiDAR)", xlab = "Residuals", col = "pink")
  
  # Summary of residuals
  summary(residuals)
  shapiro.test(residuals) 
  
  # Assumptions not met
  
  cor.test(combined_data$Perc_under_2m_2020, combined_data$radius_Perc_under_2m_2020,method='spearman', use = "complete.obs")
  
  # so Perc_under_2m_2020 for points sampled and plots overall are correlated
  
  #--------------------#
  # Height_cv_2020 #
  #--------------------#
  
  hist(combined_data$Height_cv_2020)
  hist(combined_data$radius_Height_cv_2020)
  
  # Calculate the residuals (difference between the two height measurements)
  residuals <- combined_data$Height_cv_2020 - combined_data$radius_Height_cv_2020
  hist(residuals, main = "Residuals (Ground - LiDAR)", xlab = "Residuals", col = "pink")
  
  # Summary of residuals
  summary(residuals)
  shapiro.test(residuals) 
  
  # Assumptions not met
  
  cor.test(combined_data$Height_cv_2020, combined_data$radius_Height_cv_2020,method='spearman', use = "complete.obs")
  
  # so Height_cv_2020 for points sampled and plots overall are correlated
  
  #--------------------#
  # Max_Canopy_Height_2020 #
  #--------------------#
  
  hist(combined_data$Max_Canopy_Height_2020)
  hist(combined_data$radius_Max_Canopy_Height_2020)
  
  # Calculate the residuals (difference between the two height measurements)
  residuals <- combined_data$Max_Canopy_Height_2020 - combined_data$radius_Max_Canopy_Height_2020
  hist(residuals, main = "Residuals (Ground - LiDAR)", xlab = "Residuals", col = "pink")
  
  # Summary of residuals
  summary(residuals)
  shapiro.test(residuals) 
  
  # Assumptions not met
  
  cor.test(combined_data$Max_Canopy_Height_2020, combined_data$radius_Max_Canopy_Height_2020,method='spearman', use = "complete.obs")
  
  # so Max_Canopy_Height_2020_cv_2020 for points sampled and plots overall are correlated
  
  #--------------------#
  # LAI #
  #--------------------#
  
  hist(combined_data$LAI)
  hist(combined_data$radius_LAI)
  
  # Calculate the residuals (difference between the two height measurements)
  residuals <- combined_data$LAI - combined_data$radius_LAI
  hist(residuals, main = "Residuals (Ground - LiDAR)", xlab = "Residuals", col = "pink")
  
  # Summary of residuals
  summary(residuals)
  shapiro.test(residuals) 
  
  # Assumptions not met
  
  cor.test(combined_data$LAI, combined_data$radius_LAI,method='spearman', use = "complete.obs")
  
  # so LAI for points sampled and plots overall are correlated
  
  
  
  # Load required package
  install.packages('e1071')
  library(e1071)
  # Load required package
  library(e1071)
  
  # Create a data frame to store results
  cor_results <- data.frame(
    Metric = character(),
    Skewness = numeric(),
    Correlation = numeric(),
    p_value = numeric(),
    Method = character(),
    stringsAsFactors = FALSE
  )
  
  # Get the list of metrics
  metrics <- gsub("radius_", "", names(plot_data)[grepl("^radius_", names(plot_data))])
  metrics <- metrics[metrics %in% names(combined_data)]  # Make sure both versions exist
  
  # Loop through metrics
  for (metric in metrics) {
    radius_col <- paste0("radius_", metric)
    
    x <- combined_data[[metric]]
    y <- combined_data[[radius_col]]
    
    # Skip non-numeric columns
    if (!is.numeric(x) || !is.numeric(y)) next
    
    # Remove rows with missing data
    valid_rows <- complete.cases(x, y)
    x_clean <- x[valid_rows]
    y_clean <- y[valid_rows]
    
    # Skip if not enough data
    if (length(x_clean) < 10) next
    
    # Calculate skewness
    skew_x <- skewness(x_clean, na.rm = TRUE)
    skew_y <- skewness(y_clean, na.rm = TRUE)
    
    # Choose method
    method <- ifelse(abs(skew_x) > 1 | abs(skew_y) > 1, "spearman", "pearson")
    
    # Run correlation test
    cor_test <- cor.test(x_clean, y_clean, method = method)
    
    # Store results
    cor_results <- rbind(cor_results, data.frame(
      Metric = metric,
      Skewness = round(c(skew_x, skew_y), 3) |> paste(collapse = ", "),
      Correlation = unname(cor_test$estimate),
      p_value = cor_test$p.value,
      Method = method,
      stringsAsFactors = FALSE
    ))
  }
  cor_results 
  
  
  plot(combined_data$Canopy_Volume_2020, combined_data$radius_Canopy_Volume_2020)
  
  
# Ranges #
range(sp_data$Elevation)
range(sp_data$Slope)
range(sp_data$Aspect_Cos)
range(sp_data$TRI)
range(sp_data$TPI)
range(sp_data$HALP)
range(sp_data$Roughness)
range(sp_data$FlowDir)
range(sp_data$Solar_Radiation)
range(sp_data$Mean_Curve)
range(sp_data$Profile_Curve)
range(sp_data$Plane_Curve)

# Variables to analyze
topo_vars <- c("Elevation", "Slope", "Aspect_Cos", "TRI", "TPI", "HALP", "Roughness", "Solar_Radiation", 
               "Mean_Curve", "Profile_Curve", "Plane_Curve")

# Function to calculate stats for each variable
summary_stats <- function(x) {
  rng <- diff(range(x, na.rm = TRUE))
  sd_val <- sd(x, na.rm = TRUE)
  mean_val <- mean(x, na.rm = TRUE)
  cv <- ifelse(mean_val != 0, sd_val / abs(mean_val), NA)  # avoid division by zero
  c(Range = rng, SD = sd_val, CV = cv)
}

# Apply function to each variable and transpose for readability
stats <- t(sapply(sp_data[topo_vars], summary_stats))

# Show results nicely
print(stats)

hist(sp_data$HALP)
hist(sp_data$Elevation)
