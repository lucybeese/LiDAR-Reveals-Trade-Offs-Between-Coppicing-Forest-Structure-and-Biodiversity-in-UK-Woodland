##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

### Bontuchel Biodiversity Analysis ###

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

##%%%%%%%%%%%%%%%%%%%%%%%%%%%##

### Calculating Canopy Gaps ###

##%%%%%%%%%%%%%%%%%%%%%%%%%%%##

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

sp_data<- read.csv('Data/piv_data.csv')
bontuchel<-vect('Data/Bontuchel_Coppice_Panels.shp')
bon_CHM_2020<- rast("Data/bon_CHM_2020_rep.tif")

#//////////////////////#
#### LOAD FUNCTIONS ####
#//////////////////////#

# Create and load updated ForestGapR Functions to work with terra

#---------------#
# GetForestGaps #
#---------------#

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

#///////////////////////////////#
#### CALCULATING CANOPY GAPS ####
#///////////////////////////////#

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

#%%%%%%%%%%%%%%%%%%%%#
###### Plotting ######
#%%%%%%%%%%%%%%%%%%%%#

# Plot figure
cairo_pdf("Figures/GapAnalysis.pdf",
          width = 12,
          height = 12,
          family = "Cambria",
          bg = "white")
par(family = "Cambria")
par(mfrow=c(3,2),mar=c(6,8,2,2),las=1,xpd=T)

# Plot the original CHM and gaps for each threshold
par(mfrow = c(3, 2)) # Set up a 3x2 plotting grid

# Plot the original CHM
plot(bon_2020_CHM_crop, main = "Original Canopy Height Model", col = viridis::viridis(10))

# Plot gaps for each threshold
for (threshold in thresholds) {
  plot(bon_2020_CHM_crop, main = paste0("Gaps at Threshold ", threshold, "m"), col = viridis::viridis(10))
  plot(gap_layers[[as.character(threshold)]], col = "#D34A68", add = TRUE, legend = FALSE)
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

#-----------#
# Save data #
#-----------#

write.csv(sp_data,"Data/gap_data.csv", row.names = F)
sp_data<- read.csv('Data/gap_data.csv')

#//////////////////////////////#
#### CLEAN UP SPECIES NAMES ####
#//////////////////////////////#

unique(sp_data$Species)

# Taxonomic groups levels vary, tried for species where possible, but some functional groups remain

#-----------------------------#
# Replace the specific string #
#-----------------------------#

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

#----------------------------------------------#
# Find species that are not in the species_map #
#----------------------------------------------#

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

#-----------#
# Save data #
#-----------#

write.csv(sp_data,"Data/gap_data_sp.csv", row.names = F)

table(sp_data$Taxonomic_Level)
sum(is.na(sp_data$Scientific_Name))
length(unique(sp_data$Species))
table(sp_data$Species, useNA = "ifany")
length(unique(sp_data$ET_ID))
