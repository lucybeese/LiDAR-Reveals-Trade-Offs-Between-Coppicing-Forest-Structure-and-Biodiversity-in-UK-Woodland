##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

### Bontuchel Biodiversity Analysis ###

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#%%%%%%%%%%%%%%%%%%%%##

### EarthTrack Data ###

##%%%%%%%%%%%%%%%%%%%##

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

#/////////////////#
#### LOAD DATA ####
#/////////////////#

bontuchel<-vect('Data/Bontuchel_Coppice_Panels.shp')
ET_points<-vect('Data/ET_Points.gpkg')

#//////////////////////////#
#### CREATE A DATAFRAME ####
#//////////////////////////#

plot(ET_points)
set.crs(ET_points, "epsg:27700")

# Convert the SpatVector to a data frame
ET_points_df <- as.data.frame(ET_points)
colnames(ET_points_df)
length(ET_points_df$radius) # 71 points in total

ET_in_plots <- terra::intersect(ET_points, bontuchel)
ET_in_plots_df<-as.data.frame(ET_in_plots)

# Give each ET point an ID
ET_in_plots_df$ET_ID <- paste0("ET_", seq_len(nrow(ET_in_plots_df)))
length(unique(ET_in_plots_df$ET_ID)) # 62 of those points were in plots

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
length(unique(ET_in_plots_df$Coppicing_Plot)) # ET points taken in 56 of the plots 

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
write.csv(et_data,"Data/et_data.csv",row.names = F)
