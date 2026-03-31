##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

####Bontuchel Biodiversity Analysis ###

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

##%%%%%%%%%%%%%%%%##

### Dormice Data ###

##%%%%%%%%%%%%%%%%##

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

# Load in data and create a dataframe 
dormice<-vect('Data/Individual Species - Muscardinus avellanarius_points.shp')
dor_df<- as.data.frame(dormice)

dormice2<-vect('Data/Coed Fron Wyllt, Bontuchel - Dormouse Monitoring Results May21 to Oct24.xlsx', layer = "Records")
dor_df_2<- as.data.frame(dormice2)

dormice3<-vect('Data/Bontuchel NDMP.xlsx', layer='NDMP_box')
dor_df_3<- as.data.frame(dormice3)

dormice4<-vect('Data/E130908_new export inc custom fields.xlsx', layer='Records')
dor_df_4<- as.data.frame(dormice4)

boxes_df<-read.csv('Data/boxes_df.csv')

bontuchel<-vect('Data/Bontuchel_Coppice_Panels.shp')

bon_CHM_2020<- rast("Data/bon_CHM_2020_rep.tif")

plot_data<- read.csv("Data/Topo_metrics.csv")

#//////////////////#
#### FORMATTING ####
#//////////////////#

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

#////////////////#
##### dor_df #####
#////////////////#

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

#//////////////////#
##### dor_df_2 #####
#//////////////////#

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

#//////////////////#
##### dor_df_3 #####
#//////////////////#

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

#//////////////////#
##### dor_df_4 #####
#//////////////////#

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

#//////////////////////#
#### MERGE DATASETS ####
#//////////////////////#

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

# Drop undated absence-style legacy records while retaining undated presence records.
dormice_df <- dormice_df %>%
  filter(!(is.na(Date_DM_Recorded) & is.na(Negative_Rec) & No_of_Dormice == "0"))

# Check the final number of rows
nrow(dormice_df) # 4408

#/////////////////////#
#### TIDY DATASETS ####
#/////////////////////#

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
writeVector(merged, 'Data/DM_and_boxes.shp', filetype = "ESRI Shapefile", overwrite=TRUE)

# See which Dormice records intersect which plots
DM_in_plots <- terra::intersect(merged, bontuchel)
DM_in_plots_df<-as.data.frame(DM_in_plots)
head(DM_in_plots_df)
nrow(DM_in_plots_df) # Only 1180 dormice records out of 4869 were in our plots

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
writeVector(DM_in_plots, 'Data/DM_in_plots.shp', filetype = "ESRI Shapefile", overwrite=TRUE)

# Check plot
plot.new()
bon_plot<-crop(bon_CHM_2020, bontuchel)
plot(bon_plot)

plot(bontuchel, add=TRUE)
plot(DM_in_plots, add=TRUE, col='pink')

#---------------------------------------#
# Plot all dormice records in each year #
#---------------------------------------#

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

#------------------#
# Fix DM Abundance #
#------------------#

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
write.csv(DM_in_plots_df,"Data/DM_in_plots_df_V2.csv",row.names = F)

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

# count all non na values in df
sum(!is.na(DM_in_plots_df$No_of_Dormice)) # 599
DM_in_plots_df$No_of_Dormice<-as.numeric(DM_in_plots_df$No_of_Dormice)
hist(DM_in_plots_df$No_of_Dormice)

head(plot_data)
head(DM_in_plots_df)
nrow(plot_data)
nrow(DM_in_plots_df)

# Rename ID_No to coppicing plot
DM_in_plots_df <- DM_in_plots_df %>%
  rename(Coppicing_Plot = ID_No)

length(unique(DM_in_plots_df$Coppicing_Plot))  # Number of unique Coppicing_Plot values in dormice data
length(unique(plot_data$Coppicing_Plot))       # Number of unique Coppicing_Plot values in plot data. 30 of 57 plots had dormice records in.

# Merge
dorm_data <- merge(DM_in_plots_df, plot_data, by = "Coppicing_Plot", all = TRUE)
nrow(dorm_data)

# Check names
colnames(dorm_data)
dorm_data$box_yr<-NULL
dorm_data$Area<-NULL

# Remove duplicates
dorm_data <- dorm_data %>%
  distinct()

# Save the dataframe
write.csv(dorm_data,"Data/DM_metrics.csv",row.names = F)
sum(!is.na(dorm_data$No_of_Dormice))
# 588

# Actual dormice found:
dorm_data_summary <- dorm_data %>%
  group_by(Coppicing_Plot, Year_DM_Recorded, Dor_Box) %>%
  summarise(total_Dormice = sum(No_of_Dormice, na.rm = TRUE), .groups = "drop")

sum(dorm_data_summary$total_Dormice) # 385
