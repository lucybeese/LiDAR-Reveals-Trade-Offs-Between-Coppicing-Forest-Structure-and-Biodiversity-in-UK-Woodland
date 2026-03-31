##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

### Bontuchel Biodiversity Analysis ###

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

##%%%%%%%%%%%%%%%%%%%%%%%##

### Supplementary stats ###

##%%%%%%%%%%%%%%%%%%%%%%%##

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
library(dplyr)
library(corrplot)
library(mgcv)
library(DHARMa)
library(MASS)
library(MuMIn)
library(colorspace)
library(patchwork)
library(TeachingDemos)
library(e1071)

#/////////////////#
#### LOAD DATA ####
#/////////////////#

sp_data<-read.csv('Data/final_data_dedup.csv')
ET_points<-vect('Data/ET_Points.gpkg')
bon_CHM_2020<- rast("Data/bon_CHM_2020_rep.tif")
bontuchel<-vect('Data/Bontuchel_Coppice_Panels.shp')
bon_DTM_2020<- rast("Data/bon_DTM_2020_rep.tif")
boxes<-vect('Data/Boxcoordinates.shp')
ET_points<-crop(ET_points, bontuchel)

#/////////////////#
#### PLOT BIAS ####
#/////////////////#

# Here we want to look at how well did each sample represent each plot
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

# Normalize datetime format in both ET_points and buffer_15m
ET_points$ET_Rec_Date <- as.POSIXct(ET_points$`_CREATION_DATE`, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
buffer_15m$ET_Rec_Date <- as.POSIXct(buffer_15m$'_CREATION_DATE', format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Ensure buffer_15m and ET_points match by ET_Rec_Date
# Safe inner join based on ET_Rec_Date to get the same order
buffer_15m <- buffer_15m[order(buffer_15m$ET_Rec_Date), ]
ET_points  <- ET_points[order(ET_points$ET_Rec_Date), ]

# Confirm matching
stopifnot(all(buffer_15m$ET_Rec_Date == ET_points$ET_Rec_Date))

# Spatial join to get coppicing plot IDs
matches <- terra::extract(bontuchel["ID_No"], ET_points)

# Assign matched IDs to buffer_15m
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

bon_Slope <- rast('Data/bon_Slope.tif')
bon_Aspect<- rast('Data/bon_Aspect.tif')
bon_flowdir<- rast('Data/bon_Flowdir.tif')
bon_planc<- rast('Data/bon_Planc.tif')
bon_profc<- rast('Data/bon_Profc.tif')

# Create data frame to store data
plot_data <- data.frame(matrix(NA, nrow=dim(buffer_15m)[1], ncol=16))
names(plot_data)<-c("Coppicing_Plot","radius_plotvalue_ID","ET_Rec_Date", 
                    "radius_Can_cover_2020", "radius_Height_cv_2020"
                    ,"radius_Slope","radius_Max_Canopy_Height_2020", 
                    "radius_FT_10_20m","radius_Intensity","radius_Aspect", "radius_LAI", "radius_VCI", "radius_Canopy_Volume_2020",
                    "radius_FlowDir", "radius_Plane_Curve", "radius_Profile_Curve")

plot_data$ET_Rec_Date <- buffer_15m$ET_Rec_Date

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
  return(proportion)}

#---------------#
# First Returns #
#---------------#

# Load LAS file 
pt_cloud <- readLAS("Data/Dormouse_Normalised.laz")

# Set the CRS to EPSG:27700 (British National Grid)
projection(pt_cloud) <- "EPSG:27700"

# Summarise and plot
summary(pt_cloud$Z)
plot(pt_cloud)

# Set any first returns below 0 to zero
pt_cloud@data$Z[pt_cloud@data$Z<0]<-0

# Filter the first returns (Class 1 corresponds to first returns)
first_returns <- filter_poi(pt_cloud, ReturnNumber == 1)

min(first_returns$Z) #0
max(first_returns$Z) #31.41
mean(first_returns$Z) #12.2
hist(first_returns$Z)

#-----------#
# Intensity #
#-----------#

# Filter the for Intensity
pt_intensity <- pt_cloud$Intensity

#--------------#
# Extract Data #
#--------------#

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
  
  total_first_returns <- nrow(fr_pt_crop)
  
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
  
  
  #-------#
  # Slope #
  #-------#
  
  Slope_crop<-crop(bon_Slope,buffer_15m[i,])
  plot_data[i,"radius_Slope"]<-global(Slope_crop, fun='mean', na.rm=TRUE)
  

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
  
  
  # Give each plotvalue a unique ID
  plot_data[i, "radius_plotvalue_ID"] <- paste0("pv_", i)
  
  # Progress (% of plots completed)
  print((i/dim(plot_data)[1])*100)
  
}

# Minimum DTM Value (for Height above lowest point, HALP)
min_DTM <- min(values(bon_DTM_2020), na.rm = TRUE) #121.38
plot_data$radius_HALP<-plot_data$radius_Elevation-min_DTM

# Aspect (cos transformed)
plot_data$radius_Aspect_Cos<-cos(plot_data$radius_Aspect)

#------------#
# Check data #
#------------#

summary(plot_data)

# Now join data together by coppicing plot and ET_Rec_Date

# Ensure ET_Rec_Date is in POSIXct format in both
sp_data$ET_Rec_Date <- as.POSIXct(sp_data$ET_Rec_Date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
plot_data$ET_Rec_Date <- as.POSIXct(plot_data$ET_Rec_Date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Inner join by Coppicing_plot and ET_Rec_Date
combined_data <- inner_join(sp_data, plot_data, by = c("Coppicing_Plot", "ET_Rec_Date"))

sum(is.na(sp_data$Coppicing_Plot))      # Should be 0
sum(is.na(sp_data$ET_Rec_Date))         # Should be 0
sum(is.na(plot_data$Coppicing_Plot))    # Should be 0
sum(is.na(plot_data$ET_Rec_Date))       # Should be 0

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

# See results
cor_results 

#/////////////////////////////////#
#### TOPOGRAPHIC HETEROGENEITY ####
#/////////////////////////////////#

# Ranges 
range(sp_data$Slope)
range(sp_data$HALP)
range(sp_data$FlowDir)
range(sp_data$Aspect_Cos)
range(sp_data$Profile_Curve)
range(sp_data$Plane_Curve)

# Variables to analyse
topo_vars <- c( "Slope", "Aspect_Cos", "HALP",
                "Profile_Curve", "Plane_Curve")

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

# Range          SD        CV
# Slope         14.96062522 4.166706983 0.2722807
# Aspect_Cos     1.99992589 0.637912822 5.1469204
# HALP          43.59685220 9.736999615 0.3823697
# Profile_Curve  0.01177338 0.002420335 1.9046668
# Plane_Curve    0.02398893 0.005689547 1.1541215

#////////////////////////////#
#### CANOPY HETEROGENEITY ####
#////////////////////////////#

# Ranges 
range(sp_data$Mean_Canopy_Height_2020)
range(sp_data$Perc_under_2m_2020)
range(sp_data$Can_cover_2020)
range(sp_data$Height_cv_2020)
range(sp_data$Max_Canopy_Height_2020)
range(sp_data$FT_5_10m)
range(sp_data$FT_1.37m)
range(sp_data$FT_10_20m)
range(sp_data$Intensity)
range(sp_data$LAI)
range(sp_data$VCI)
range(sp_data$Canopy_Volume_2020)

# Variables to analyze
cano_vars <- c( "Perc_under_2m_2020", "Can_cover_2020", "Max_Canopy_Height_2020", "FT_10_20m", 
               "Intensity", "LAI", "VCI", "Canopy_Volume_2020")

# Function to calculate stats for each variable
summary_stats <- function(x) {
  rng <- diff(range(x, na.rm = TRUE))
  sd_val <- sd(x, na.rm = TRUE)
  mean_val <- mean(x, na.rm = TRUE)
  cv <- ifelse(mean_val != 0, sd_val / abs(mean_val), NA)  # avoid division by zero
  c(Range = rng, SD = sd_val, CV = cv)
}

# Apply function to each variable and transpose for readability
stats <- t(sapply(sp_data[cano_vars], summary_stats))

# Show results nicely
print(stats)

# Range           SD         CV
# Perc_under_2m_2020     1.651625e+01 3.626420e+00 0.93261063
# Can_cover_2020         1.651625e+01 3.626420e+00 0.03773137
# Max_Canopy_Height_2020 2.109890e+01 2.715270e+00 0.12076349
# FT_10_20m              8.523832e+01 1.141598e+01 0.17427770
# Intensity              4.044492e+02 9.165108e+01 0.17737146
# LAI                    3.158546e+00 6.905009e-01 0.18487828
# VCI                    2.417597e-01 2.617911e-02 0.02741910
# Canopy_Volume_2020     7.855137e+04 1.890948e+04 0.39882911

summary(sp_data$Mean_Canopy_Height_2020)

#///////////////////#
#### CORRELATION ####
#///////////////////#

final_data <- sp_data

# Select only the relevant columns for uniqueness check - prevents pseudoreplication 
unique_rows <- !duplicated(final_data[, c("ET_ID")])

# Keep only unique rows in the dataframe
sp_data_unique <- final_data[unique_rows, ]
colnames(sp_data_unique)

custom_labels <- c(
  Max_Canopy_Height_2020 = "Max Canopy Height",
  Gap_Proportion_2020 = "Canopy Gap Proportion",
  Intensity = "LiDAR Intensity",
  VCI = "VCI",
  LAI = "LAI",
  Solar_Radiation = "Solar Radiation",
  FT_10_20m = "First returns  (10–20m)",
  HALP = "HALP",
  Aspect_Cos = "Aspect (Cosine)",
  Elevation = "Elevation",
  Slope = "Slope",
  Avg_Shannon_Index = "Mean Shannon Index",
  recovery_period = "Time since Coppicing",
  Mean_Curve = "Mean Curvature",
  FlowDir = "Flow Direction",
  Plane_Curve = "Plane Curvature",
  Profile_Curve = "Profile Curvature",
  Canopy_Volume_2020 = "Canopy Volume",
  Mean_Dormice = "Mean Dormice Abundance",
  num_unique_records ="No. Dormice Records",
  num_ET_points =" No. ET Records"
)

colnames(sp_data_unique)

cor_matrix <- cor(sp_data_unique[, c("Max_Canopy_Height_2020", "Gap_Proportion_2020", 
                                     "Intensity", "VCI", "LAI", "Solar_Radiation",
                                     "FT_10_20m", "HALP", "Aspect_Cos", "num_ET_points", "num_unique_records",
                                     "Elevation", "Slope", "Avg_Shannon_Index", "recovery_period", "Mean_Curve", 
                                     "FlowDir", "Plane_Curve", "Profile_Curve", "Canopy_Volume_2020", "Mean_Dormice")], 
                  use = "pairwise.complete.obs", method='spearman')  

# Apply custom labels
colnames(cor_matrix) <- custom_labels[colnames(cor_matrix)]
rownames(cor_matrix) <- custom_labels[rownames(cor_matrix)]

cor_matrix

#%%%%%%%%%%%%%%%%%%%%#
###### Plotting ######
#%%%%%%%%%%%%%%%%%%%%#

# Open PDF device
cairo_pdf("Figures/CorrPlot.pdf",
          width = 12,
          height = 12,
          family = "Cambria",
          bg = "white")

par(mfrow = c(1,1), mar = c(4,4,1,1), las = 1, family = "Cambria")

# Define the custom color gradient
# Diverging colour palette
neg <- colorRampPalette(c("#D34A68", "pink"))(100)
pos <- colorRampPalette(c("#F7E7B6", "#F5B856"))(100)
custom_colors <- c(neg, pos)

# Set plot parameters
par(mfrow = c(1, 1), mar = c(6, 8, 2, 2), las = 1, xpd = TRUE, cex.axis = 1.4, ps = 10)

# Create the correlation plot
corrplot(cor_matrix, 
         type = 'lower', 
         tl.col = "black", 
         tl.cex = 1, # Make text smaller
         tl.srt = 70,   # Rotate text for better readability
         col = custom_colors, 
         cl.cex = 1.5)

# Close the device
dev.off()
