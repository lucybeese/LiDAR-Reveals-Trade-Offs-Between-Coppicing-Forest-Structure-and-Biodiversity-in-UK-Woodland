##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

### Bontuchel Biodiversity Analysis ###

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

##%%%%%%%%%%%%%%%%%%%%%##

### Metric Extraction ###

##%%%%%%%%%%%%%%%%%%%%%##

#//////////////////////#
#### LOAD LIBRARIES ####
#//////////////////////#

library(terra)
library(sf)
library(MultiscaleDTM)
library(lidR)
library(sp)
library(leafR)
library(tidyr)
library(ggplot2)

#/////////////////#
#### LOAD DATA ####
#/////////////////#

# Load LAS file 
pt_cloud <- readLAS("Data/Dormouse_Normalised.laz")

# Set the CRS to EPSG:27700 (British National Grid)
projection(pt_cloud) <- "EPSG:27700"

bontuchel<-vect('Data/Bontuchel_Coppice_Panels.shp')
bon_CHM_2020<- rast("Data/bon_CHM_2020_rep.tif")
bon_DTM_2020<- rast("Data/bon_DTM_2020_rep.tif")

#------------------#
# Load the Rasters #
#------------------#

bon_Slope <- rast('Data/bon_Slope.tif')
bon_TPI<- rast('Data/bon_TPI.tif')
bon_TRI<- rast('Data/bon_TRI.tif')
bon_Aspect<- rast('Data/bon_Aspect.tif')
bon_Roughness<- rast('Data/bon_Roughness.tif')
bon_flowdir<- rast('Data/bon_Flowdir.tif')
bon_planc<- rast('Data/bon_Planc.tif')
bon_profc<- rast('Data/bon_Profc.tif')
bon_DTM_2020_hillshade<- rast('Data/bon_Hillshade.tif')
bon_meanc<- rast('Data/bon_Meanc.tif')
bon_TWI<- rast('Data/bon_TWI.tif')
bon_solar<- rast('Data/bon_Solar.tif')

#////////////////////////////////////////#
#### EXTRACT METRICS FROM POINT CLOUD ####
#////////////////////////////////////////#

#---------------#
# First Returns #
#---------------#

# Summarise and plot
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

#------#
# Save #
#------#

writeLAS(pt_cloud, 'Data/pt_cloud.laz', index = FALSE)

#//////////////////////#
#### LOAD FUNCTIONS ####
#//////////////////////#

#---------------#
# Canopy Volume #
#---------------#

canopy_volume <- function(chm) {
  volume <- sum(values(chm), na.rm = TRUE) * res(chm)[1] * res(chm)[2]
  return(volume)
}

#/////////////////////////////////////#
#### EXTRACT METRICS FOR EACH PLOT ####
#/////////////////////////////////////#

#---------------------------------------------------#
# Extract canopy and topographic data for each plot #
#---------------------------------------------------#

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
  target_sp <- sp::CRS(st_as_text(target))
  
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
  
  #---------------------------------------------------------------#
  # Calculate the percentage of first returns above 10-20m meters #
  #---------------------------------------------------------------#
  
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

#------------------------#
# Set columns to numeric #
#------------------------#

# List of columns to be numeric
cols_to_numeric <- c("Mean_Canopy_Height_2020", "Perc_under_2m_2020", "Can_cover_2020", 
                     "Height_cv_2020", "Elevation", "Slope", "TRI", "Max_Canopy_Height_2020",
                     "FT_5_10m", "FT_1.37m", "FT_10_20m", "Intensity", "Aspect",
                     "TPI", "Roughness", "LAI", "VCI", "Canopy_Volume_2020", "FlowDir",
                     "Plane_Curve", "Profile_Curve", "Mean_Curve", "Solar_Radiation",
                     "HALP", "Aspect_Cos")

# Convert these columns to numeric
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
write.csv(plot_data,"Data/Topo_metrics.csv",row.names = F)
