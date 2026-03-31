##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

### Bontuchel Biodiversity Analysis ###

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

##%%%%%%%%%%%%%%%%##

### Biodiversity ###

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

sp_data<- read.csv('Data/gap_data_sp.csv')

#//////////////////////////////////////#
#### DERIVE MEASURE OF BIODIVERSITY ####
#//////////////////////////////////////#

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

#------#
# Save #
#------#

write.csv(shannon_index, "Data/shannon_index_data.csv", row.names = FALSE)

# Merge the Shannon Index into the original sp_data dataframe
sp_data <- sp_data %>%
  left_join(shannon_index, by = c("Coppicing_Plot", "ET_ID"))

length(unique(sp_data$ET_ID))
# for plots that have more than one unique ET ID,  need to work out the shannon index of both and then get a mean.

# Calculate mean Shannon Index per Coppicing_Plot
mean_shannon <- shannon_index %>%
  group_by(Coppicing_Plot) %>%
  summarise(Mean_Shannon_Index = mean(Shannon_Index, na.rm = TRUE))

# Join mean Shannon back into sp_data
sp_data <- sp_data %>%
  left_join(mean_shannon, by = "Coppicing_Plot")

length(unique(sp_data$ET_ID))

#------#
# Save #
#------#
write.csv(sp_data,"Data/sp_data.csv",row.names = F)
