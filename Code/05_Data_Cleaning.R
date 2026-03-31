##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

### Bontuchel Biodiversity Analysis ###

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

##%%%%%%%%%%%%%%%%%%%##

### Data Cleaning ###

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

dorm_data<-read.csv('Data/DM_metrics.csv')
et_data<-read.csv('Data/et_data.csv')

head(dorm_data, n=1)
head(et_data, n=1)
nrow(dorm_data)
nrow(et_data)

#//////////////////#
#### CLEAN DATA ####
#//////////////////#

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

# Joining these datasets will introduce duplicates due to the one-to-many relationship between the earthtrack data (et_data) 
# and dormice data (dorm_data). Each Coppicing_Plot in et_data typically has a single record, whereas dorm_data has multiple 
# records per Coppicing_Plot.
# If both dorm_data and et_data contain multiple records for the same Coppicing_Plot, the merge() function will create a Cartesian join (many-to-many merge). 
# This means that each combination of matching records will be represented in the output.

#------------------#
##### Clean up #####
#------------------#

# Remove columns where all is NA
final_data <- final_data %>% select_if(~ !all(is.na(.)))

nrow(final_data)
length(unique(et_data$ET_ID))

#-----------------------#
##### Rename layers #####
#-----------------------#

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

min(na.omit(final_data$Layer_0_height)) # 10

# Save
write.csv(final_data,"Data/combined_data.csv",row.names = F)

# Read as vector
final_data_shp<- vect("Data/combined_data.csv")

# Write as vector
writeVector(final_data_shp, 'Data/final_data.shp', filetype = "ESRI Shapefile", overwrite=TRUE)
sum(!is.na(final_data$No_of_Dorm)) 

length(unique(final_data$Dor_Box[!is.na(final_data$No_of_Dormice)]))
# 78 boxes

length(unique(final_data$ET_ID))

#//////////////////#
#### PIVOT DATA ####
#//////////////////#

#---------------#
##### Pivot #####
#---------------#

# Any column that begins with 'Layer_' should be put in long format, into a column with all the same suffixes. A new 'Layer column should be made.
# Ie: any value ending with _cover should be in  a new column called Cover
# any value from column ending in _height should be in  a new column called Height
# any value from column ending in _lifeform should be in  a new column called Lifeform
# any value from column ending in _dominantSpecies, _codominantSpecies, _codominantSpecies1, or _codominantSpecies2,  should be in  a new column called Species
# any value from column ending in _Pct should be in  a new column called Proportion of Species
# the new row generated should have all the same information in the other columns as it does originally.

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

#-------------------#
##### Tidy data #####
#-------------------#

colnames(expanded_data)
names(expanded_data)[names(expanded_data) == "indicatorSpecies1"] <- "Indicator_species"
names(expanded_data)[names(expanded_data) == "Layer"] <- "Forest_layer"
names(expanded_data)[names(expanded_data) == "Cover"] <- "ET_total_canopy_cover"
names(expanded_data)[names(expanded_data) == "Height"] <- "ET_layer_height"
names(expanded_data)[names(expanded_data) == "Lifeform"] <- "ET_layer_lifeform"

#-------------------------------#
##### Years since Coppicing #####
#-------------------------------#

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
write.csv(sp_data,"Data/piv_data.csv", row.names = F)

sp_data$ET_Latitude<-as.numeric(sp_data$ET_Latitude)
sp_data$ET_Longitude<-as.numeric(sp_data$ET_Longitude)
sp_data$Forest_layer <- as.numeric(as.character(sp_data$Forest_layer))
sp_data_spat<-vect(sp_data, geom=c('ET_Latitude', 'ET_Longitude'))

# Save
writeVector(sp_data_spat, 'Data/pivoted_dataset.shp', filetype = "ESRI Shapefile", overwrite=TRUE)

forest_layer_range <- sp_data %>%
  group_by(Coppicing_Plot) %>%
  summarise(Forest_layer_range = max(Forest_layer, na.rm = TRUE) - min(Forest_layer, na.rm = TRUE))

# View result
min(forest_layer_range$Forest_layer_range) 
max(forest_layer_range$Forest_layer_range)

