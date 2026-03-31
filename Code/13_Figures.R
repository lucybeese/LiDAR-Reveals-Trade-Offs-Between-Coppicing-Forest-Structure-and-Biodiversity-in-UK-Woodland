##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

### Bontuchel Biodiversity Analysis ###

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

library(terra)
library(ggplot2)
library(magick)
library(MASS)
library(mgcv)
library(colorspace)
library(dplyr)
library(TeachingDemos)
library(patchwork)
library(corrplot)
library(RColorBrewer)
library(paletteer)

##%%%%%%%%%%%##

### Figures ###

##%%%%%%%%%%%##

##%%%%%%%%%%%%%%%%%%%%%%%%##

#### 05 Data Processing ####

##%%%%%%%%%%%%%%%%%%%%%%%%##

#///////////////////////////////#
###### CHM 2020 Bontuchel #######
#///////////////////////////////#

# Load data
bon_CHM_2020<-rast('Data/bon_CHM_2020_rep.tif')
bontuchel<-vect('Data/Bontuchel_Coppice_Panels.shp')

# Open a pdf device
pdf("Figures/CHM_2020.pdf")

# Plot
plot(bon_CHM_2020)
plot(bontuchel, add=TRUE, border='#810847', lwd= 2)

# Close the device
dev.off()

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### 06 Calculating Canopy Gaps ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#/////////////////////////#
###### Gap Analysis #######
#/////////////////////////#

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

#------#
# Plot #
#------#

cairo_pdf("Figures/GapAnalysis.pdf",
          width = 14,
          height = 14,
          family = "Cambria",
          bg = "white")

par(
  family = "Cambria",
  mfrow = c(3, 2),
  mar = c(3, 3, 3, 2),
  las = 1,
  xpd = TRUE,
  cex = 1.5,
  cex.main = 1.6,
  cex.axis = 1.2
)

# Plot original CHM
plot(bon_2020_CHM_crop,
     main = "Original Canopy Height Model")

# Plot gaps
for (threshold in thresholds) {
  plot(bon_2020_CHM_crop,
       main = paste0("Gaps at Threshold ", threshold, " m"))
  
  plot(gap_layers[[as.character(threshold)]],
       col = "#D59AAE",
       add = TRUE,
       legend = FALSE)
}

dev.off()

##%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### 08 Preliminary Stats ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%##

#///////////////////////////////#
###### Dormice Indicators #######
#///////////////////////////////#

#-----------#
# Load Data #
#-----------#

total_dormice_df <- read.csv("Data/total_dormice_df.csv")

#------#
# Plot #
#------#

cairo_pdf("Figures/DormiceIndicators.pdf",
          width = 6,
          height = 6,
          family = "Cambria",
          bg = "white")
par(family = "Cambria")
par(mfrow=c(1,1),mar=c(6,8,2,2),las=1,xpd=T)

ggplot(total_dormice_df, aes(x =Avg_Shannon_Index , y = Mean_Dormice)) +
  geom_smooth(method = "lm", se = TRUE, color = "#f48fb1", fill = "#f48fb1", alpha = 0.18) +
  labs(
    x = "Mean Shannon Index",
    y = "Mean Number of Dormice per Panel"
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
  
  ) +
geom_point(color = adjustcolor("#f48fb1", alpha.f = 0.6), size = 2)


# Close the device
dev.off()

#////////////////////////////#
###### Lidar Vs Ground #######
#////////////////////////////#

#-----------#
# Load Data #
#-----------#

bontuchel<-vect('Data/Bontuchel_Coppice_Panels.shp')
bon_CHM_2020<-rast('Data/bon_CHM_2020_rep.tif')
final_dataset <- read.csv("Data/final_dataset.csv")
ET_points<-vect('Data/ET_Points.gpkg')

#------#
# Crop #
#------#

ET_points<-crop(ET_points, bon_CHM_2020)

# Create a buffer of 15 meters around the ET points
buffer_15m <- buffer(ET_points, width = 15)  # Buffer by 15 meters
buffer_15m$ET_Rec_Date<-buffer_15m$`_CREATION_DATE`

#--------------------#
# Mean Canopy Height #
#--------------------#

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
final_dataset$ET_Rec_Date <- as.POSIXct(final_dataset$ET_Rec_Date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Check if the conversion worked
head(final_dataset$ET_Rec_Date)

# For buffer_15m (spatial object with an attribute ET_Rec_Date)
# Extract the 'ET_Rec_Date' attribute as a vector and convert it to POSIXct
buffer_15m$ET_Rec_Date <- as.POSIXct(buffer_15m$ET_Rec_Date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Check if conversion worked
head(buffer_15m$ET_Rec_Date)

#-------------------#
# Match up the data #
#-------------------#

final_dataset$LiDAR_Mean_Canopy_Height_2020 <- mean_canopy_heights[match(final_dataset$ET_Rec_Date, buffer_15m$ET_Rec_Date)]

summary(as.numeric(final_dataset$ET_Canopy_Height))
summary(final_dataset$LiDAR_Mean_Canopy_Height_2020)

# Select only the relevant columns for uniqueness check - prevents pseudoreplication 
unique_rows <- !duplicated(final_dataset[, c("ET_ID")])

# Keep only unique rows in the dataframe
sp_data_unique <- final_dataset[unique_rows, ]

length(sp_data_unique$ET_ID)  # 62

#-------------------#
# Mean Canopy Cover #
#-------------------#

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
final_dataset$LiDAR_Mean_Canopy_cover_2020 <- mean_canopy_cover[match(final_dataset$ET_Rec_Date, buffer_15m$ET_Rec_Date)]

summary(final_dataset$ET_Canopy_Cover)
summary(final_dataset$LiDAR_Mean_Canopy_cover_2020)

# Select only the relevant columns for uniqueness check - prevents pseudoreplication 
unique_rows <- !duplicated(final_dataset[, c("ET_ID")])

# Keep only unique rows in the dataframe
sp_data_unique <- final_dataset[unique_rows, ]

#------#
# Plot #
#------#

# Open PDF device
cairo_pdf("Figures/LidarVsGround.pdf",
          width = 6,
          height = 12,
          family = "Cambria",
          bg = "white")

par(mfrow = c(2,1), mar = c(4,4,1,1), las = 1, family = "Cambria")

lims1 <- c(0, max(c(sp_data_unique$ET_Canopy_Height,
                    sp_data_unique$LiDAR_Mean_Canopy_Height_2020), na.rm = TRUE))

plot(sp_data_unique$LiDAR_Mean_Canopy_Height_2020,
     sp_data_unique$ET_Canopy_Height,
     xlab = "LiDAR Measured Height",
     ylab = "Ground Measured Height",
     xlim = lims1,
     ylim = lims1,
     pch  = 16,
     cex  = 1.5,
     col  = adjustcolor("#f48fb1", alpha.f = 0.6),
     bty  = "l")

abline(0, 1, col = "#f48fb1", lwd = 3)

lims2 <- c(0, max(c(sp_data_unique$ET_Canopy_Cover,
                    sp_data_unique$LiDAR_Mean_Canopy_cover_2020), na.rm = TRUE))

plot(sp_data_unique$LiDAR_Mean_Canopy_cover_2020,
     sp_data_unique$ET_Canopy_Cover,
     xlab = "LiDAR Measured Canopy Cover",
     ylab = "Ground Measured Canopy Cover",
     xlim = lims2,
     ylim = lims2,
     pch  = 16,
     cex  = 1.5,
     col  = adjustcolor("#f48fb1", alpha.f = 0.6),
     bty  = "l")

abline(0, 1, col = "#f48fb1", lwd = 3)

# Close PDF
dev.off()

##%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### 09 Hypothesis Testing ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#-----------#
# Load Data #
#-----------#

final_data <- read.csv("Data/final_data_dedup.csv")

# Select only the relevant columns for uniqueness check - prevents pseudoreplication 
unique_rows <- !duplicated(final_data[, c("ET_ID")])

# Keep only unique rows in the dataframe
sp_data_unique <- final_data[unique_rows, ]

#//////////////////////////#
##### Helper functions #####
#//////////////////////////#

extract_model_coefs <- function(model, label_map) {
  coefs <- summary(model)$coefficients
  
  coef_names <- rownames(coefs)[-1]
  estimates <- coefs[-1, 1]
  errors <- coefs[-1, 2]
  lower <- estimates - 1.96 * errors
  upper <- estimates + 1.96 * errors
  
  ord <- order(abs(estimates), decreasing = FALSE)
  
  list(
    coef_names = coef_names[ord],
    estimates = estimates[ord],
    errors = errors[ord],
    lower = lower[ord],
    upper = upper[ord],
    labels = label_map[coef_names[ord]]
  )
}

plot_coef_panel <- function(coef_obj, file_name, panel_label,
                            xlim_vals = c(-1, 1),
                            left_text_offset = 0.2,
                            mar_vals = c(5, 16, 4, 2)) {
  
  cairo_pdf(file_name,
            width = 7,
            height = 5,
            family = "Cambria",
            bg = "white")
  
  par(mfrow = c(1, 1), mar = mar_vals, las = 1, family = "Cambria")
  cex_text <- 1.3
  
  plot(coef_obj$estimates, seq_along(coef_obj$estimates),
       pch = 16,
       xlim = xlim_vals,
       xlab = "Estimate",
       ylab = "",
       axes = FALSE,
       cex = cex_text,
       cex.lab = cex_text)
  
  axis(1, cex.axis = cex_text)
  axis(2, at = seq_along(coef_obj$coef_names), labels = FALSE, las = 2)
  
  text(x = rep(par("usr")[1] - left_text_offset, length(coef_obj$labels)),
       y = seq_along(coef_obj$labels),
       labels = coef_obj$labels,
       adj = 1,
       xpd = TRUE,
       cex = cex_text)
  
  arrows(coef_obj$lower, seq_along(coef_obj$estimates),
         coef_obj$upper, seq_along(coef_obj$estimates),
         length = 0.1, angle = 90, code = 3)
  
  abline(v = 0, lty = 2, col = "#BA436B")
  
  text(x = par("usr")[1],
       y = par("usr")[4] + 0.5,
       labels = panel_label,
       adj = c(0, 1),
       cex = cex_text,
       xpd = NA)
  
  dev.off()
}

#/////////////////#
##### SHANNON #####
#/////////////////#

#-----------#
# Load Data #
#-----------#

final_data <- read.csv("Data/final_data_dedup.csv")
plot_data <- read.csv("Data/plot_data.csv")

# Select only the relevant columns for uniqueness check - prevents pseudoreplication 
unique_rows <- !duplicated(final_data[, c("ET_ID")])

# Keep only unique rows in the dataframe
sp_data_unique <- final_data[unique_rows, ]

#////////////////////////////////#
###### Label maps for plots ######
#////////////////////////////////#

label_map_struc <- c(
  "Max_Canopy_Height_2020" = "Max Canopy Height",
  "LAI" = "LAI",
  "Intensity" = "Intensity",
  "VCI" = "VCI",
  "Canopy_Volume_2020" = "Canopy Volume",
  "FT_10_20m" = "First returns 10–20 m",
  "num_ET_points" = "No. of ET points",
  "num_unique_records" = "No. of dormice records"
)

label_map_topo <- c(
  "HALP" = "HALP",
  "Aspect_Cos" = "Aspect (cos)",
  "Slope" = "Slope",
  "FlowDir" = "Flow direction",
  "Plane_Curve" = "Plane curvature",
  "Profile_Curve" = "Profile curvature",
  "num_ET_points" = "No. of ET points",
  "num_unique_records" = "No. of dormice records"
)

#////////////////////////////////#
###### Q4: Canopy structure ######
#////////////////////////////////#

shannon_struc_data <- sp_data_unique
shannon_struc_numeric <- sapply(shannon_struc_data, is.numeric)
shannon_struc_data[shannon_struc_numeric] <- scale(shannon_struc_data[shannon_struc_numeric])

shannon_q4_struc_model <- lm(
  Mean_Shannon_Index ~ Max_Canopy_Height_2020 + LAI + Intensity + VCI +
    Canopy_Volume_2020 + FT_10_20m + num_ET_points,
  data = shannon_struc_data
)

shannon_q4_struc_coefs <- extract_model_coefs(shannon_q4_struc_model, label_map_struc)

#------#
# Plot #
#------#

plot_coef_panel(
  coef_obj = shannon_q4_struc_coefs,
  file_name = "Figures/CoeffsStruc.pdf",
  panel_label = "a)",
  xlim_vals = c(-1, 1),
  left_text_offset = 0.2,
  mar_vals = c(5, 16, 4, 2)
)

#//////////////////////////#
###### Q4: Topography ######
#//////////////////////////#

shannon_topo_data <- sp_data_unique

shannon_topo_predictors <- c(
  "HALP", "Aspect_Cos", "Slope", "FlowDir",
  "Plane_Curve", "Profile_Curve", "num_ET_points"
)

shannon_topo_data[shannon_topo_predictors] <- scale(shannon_topo_data[shannon_topo_predictors])

shannon_q4_topo_model <- glm(
  Mean_Shannon_Index ~ HALP + Aspect_Cos + Slope + FlowDir +
    Plane_Curve + Profile_Curve + num_ET_points,
  family = Gamma(link = "log"),
  data = shannon_topo_data
)

shannon_q4_topo_coefs <- extract_model_coefs(shannon_q4_topo_model, label_map_topo)

#------#
# Plot #
#------#

plot_coef_panel(
  coef_obj = shannon_q4_topo_coefs,
  file_name = "Figures/CoeffsTopo.pdf",
  panel_label = "c)",
  xlim_vals = c(-0.2, 0.2),
  left_text_offset = 0.03,
  mar_vals = c(5, 16, 4, 2)
)

#//////////////////////////#
###### Q2: Management ######
#//////////////////////////#

plot_data$recovery_c <- plot_data$recovery_period - mean(plot_data$recovery_period, na.rm = TRUE)

shannon_q2_management_model <- gam(
  Mean_Shannon_Index ~ recovery_c + I(recovery_c^2) + num_ET_points,
  family = Gamma(link = "log"),
  data = plot_data
)

#------#
# Plot #
#------#

cairo_pdf("Figures/Management.pdf",
          width = 6,
          height = 6,
          family = "Cambria",
          bg = "white")

par(family = "Cambria")

ggplot(plot_data, aes(x = recovery_period, y = Mean_Shannon_Index)) +
  geom_smooth(
    method = "lm",
    formula = y ~ poly(x, 2),
    color = "#f48fb1",
    fill = "#f48fb1",
    alpha = 0.18
  ) +
  geom_point(color = "#f48fb1", size = 2, alpha = 0.6) +
  labs(
    x = "Time since coppicing (Years)",
    y = "Mean Shannon Index"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.25),
    axis.text.x = element_text(hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 12)
  )

dev.off()

#/////////////////#
##### DORMICE #####
#/////////////////#

#////////////////////////////////#
###### Q4: Canopy structure ######
#////////////////////////////////#

dormice_struc_data <- sp_data_unique

dormice_struc_predictors <- c(
  "Max_Canopy_Height_2020",
  "FT_10_20m",
  "LAI",
  "Intensity",
  "VCI",
  "Canopy_Volume_2020",
  "num_unique_records"
)

dormice_struc_data[dormice_struc_predictors] <- scale(dormice_struc_data[dormice_struc_predictors])

dormice_q4_struc_model <- glm(
  Mean_Dormice ~ Max_Canopy_Height_2020 + FT_10_20m + LAI +
    Intensity + VCI + Canopy_Volume_2020 + num_unique_records,
  family = poisson(link = "log"),
  data = dormice_struc_data
)

dormice_q4_struc_coefs <- extract_model_coefs(dormice_q4_struc_model, label_map_struc)

#------#
# Plot #
#------#

plot_coef_panel(
  coef_obj = dormice_q4_struc_coefs,
  file_name = "Figures/CoeffsStrucDM.pdf",
  panel_label = "b)",
  xlim_vals = c(-1, 1.5),
  left_text_offset = 0.2,
  mar_vals = c(5, 14.5, 4, 2)
)

#//////////////////////////#
###### Q4: Topography ######
#//////////////////////////#

dormice_topo_data <- sp_data_unique

dormice_topo_predictors <- c(
  "HALP", "Aspect_Cos", "Slope", "FlowDir",
  "Plane_Curve", "Profile_Curve", "num_unique_records"
)

dormice_topo_data[dormice_topo_predictors] <- scale(dormice_topo_data[dormice_topo_predictors])

dormice_q4_topo_model <- glm(
  Mean_Dormice ~ HALP + Aspect_Cos + FlowDir + Profile_Curve +
    Plane_Curve + Slope + num_unique_records,
  family = poisson(link = "log"),
  data = dormice_topo_data
)

dormice_q4_topo_coefs <- extract_model_coefs(dormice_q4_topo_model, label_map_topo)

#------#
# Plot #
#------#

plot_coef_panel(
  coef_obj = dormice_q4_topo_coefs,
  file_name = "Figures/CoeffsTopoDM.pdf",
  panel_label = "d)",
  xlim_vals = c(-1, 1.5),
  left_text_offset = 0.2,
  mar_vals = c(5, 14.5, 4, 2)
)

#//////////////////////////#
###### Q3: Management ######
#//////////////////////////#

dormice_q3_management_model <- glm(
  Mean_Dormice ~ recovery_c + I(recovery_c^2) + num_ET_points,
  family = gaussian(),
  data = plot_data
)

#------#
# Plot #
#------#

cairo_pdf("Figures/ManagementDM.pdf",
          width = 6,
          height = 6,
          family = "Cambria",
          bg = "white")

par(family = "Cambria")

ggplot(plot_data, aes(x = recovery_period, y = Mean_Dormice)) +
  geom_smooth(
    method = "glm",
    method.args = list(family = "gaussian"),
    formula = y ~ poly(x, 2),
    color = "#f48fb1",
    fill = "#f48fb1",
    alpha = 0.18
  ) +
  geom_point(color = "#f48fb1", size = 2, alpha = 0.6) +
  labs(
    x = "Time since coppicing (Years)",
    y = "Mean Dormice Abundance"
  ) +
  scale_y_continuous(breaks = seq(0, 8, by = 2)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.25),
    axis.text.x = element_text(hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 12)
  )

dev.off()

# Joined:

img1 <- image_read_pdf("Figures/Management.pdf", density = 300)
img2 <- image_read_pdf("Figures/ManagementDM.pdf", density = 300)

combined <- image_append(c(img1, img2))

image_write(combined, path = "Figures/combined_management.pdf", format = "pdf")

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### 10 Hypothesis Q2 and Q3 ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#-----------#
# Load Data #
#-----------#

model_data_species_wide <- read.csv("Data/sp_data_shan.csv")

#---------#
# Palette #
#---------#

# Define consistent color palette (blues, purples, pinks)
all_species <- sort(unique(model_data_species_wide$Species))
num_species <- length(all_species)

cool_palette <- c(
  "#E6C6DC", "pink",
  "#BA436B", "#843873", "#7A0016",
  "#4C032A", "#D66F90", "#EAC2CA",
  "#6E1E3A","#B5DA88","#276419", "#E2F3CB",
  "#337357", "#6D9F71", "#949e90","#d4a9cf"
  ,"#ed98ca", "#e85a95")

# Create darker and lighter versions
cool_palette_dark  <- darken(cool_palette, 0.25)
cool_palette_light <- lighten(cool_palette, 0.25)

# Combine
cool_palette <- c(cool_palette, cool_palette_dark, cool_palette_light)

# Extend palette 
if (num_species > length(cool_palette)) {
  cool_palette <- colorRampPalette(cool_palette)(num_species)
}
species_color_map <- setNames(cool_palette[1:num_species], all_species)

#----------------------#
##### Plot Species #####
#----------------------#

# Plot by Coppicing Year
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
    facet_wrap(~Coppicing_Plot, scales = "fixed") +
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
      x = "",
      y = "Proportional Cover (%)"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      text = element_text(family = "Cambria"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      axis.title.x = element_text(size = 16, margin = ggplot2::margin(t = 15)),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      strip.text = element_text(size = 16)
    )
  
  file_name <- paste0(
    "Figures/Species_Cover_CoppicingYear_",
    year, ".pdf"
  )
  
  cairo_pdf(file_name,
            width = 12,
            height = 8,
            family = "Cambria",
            bg = "white")
  
  print(p)
  
  dev.off()
}

# Now just 2015

# Your fixed colours
species_color_map <- c(
  "Honeysuckle" = "#DC8B7F",
  "Hazel" = "#96BDD9",
  "Oak" = "#f48fb1",
  "Bramble" = "#8986C8",
  "Sycamore" = "#A3C57B",
  "Ash" = "#E6C6DC",
  "Beech" = "#9A3268",
  "Birch" = "#65317B",
  "Clover" = "#F5B856",
  "Dog's Mercury" = "#4C032A",
  "Fern" = "#D34A68",
  "Grass" = "#9FA5B0",
  "Holly" = "#BA436B",
  "Ivy" = "#C18FCC",
  "Moss" = "#00497A"
)

#----------------------#
##### Plot Species #####
#----------------------#

data_year <- model_data_species_wide %>%
  filter(Coppicing_year == year) %>%
  mutate(
    Species = factor(Species, levels = all_species),
    Proportion_of_Species = Proportion_of_Species * 100
  )

p <- ggplot(data_year, aes(x = Stratum, y = Proportion_of_Species, fill = Species)) +
  geom_col(position = "stack", width = 0.8) +
  facet_wrap(~Coppicing_Plot, scales = "fixed") +
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
    x = "",
    y = "Proportional Cover (%)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    text = element_text(family = "Cambria"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.title.x = element_text(size = 16, margin = ggplot2::margin(t = 15)),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    strip.text = element_text(size = 16)
  )

cairo_pdf("Figures/Species_Cover_CoppicingYear_2015.pdf",
          width = 12,
          height = 8,
          family = "Cambria",
          bg = "white")

print(p)
dev.off()

#------------------------------------------#
##### Plot Shannon Index across Strata #####
#------------------------------------------#

#-----------#
# Load Data #
#-----------#

model_data_strata <- read.csv("Data/model_data_strata.csv")

# Prepare dataset
shan_strata <- model_data_strata %>%
  distinct(Coppicing_Plot, Stratum, Shannon_Index_strata) %>%
  filter(!is.na(Shannon_Index_strata))

min(shan_strata$Shannon_Index_strata) # 0
max(shan_strata$Shannon_Index_strata) # 2.2

# Ensure consistent order
shan_strata$Stratum <- factor(shan_strata$Stratum,
                              levels = c("Understorey","Shrub","Lower_Canopy","Upper_Canopy","Emergent"))

# Load violin plot function
source(file="Scripts/wvioplot.r")

# Open PDF device
cairo_pdf("Figures/SI_Strata_violin.pdf",
          width = 6,
          height = 6,
          family = "Cambria",
          bg = "white")
par(family = "Cambria")
par(mfrow=c(1,1),mar=c(6,6,2,2),las=1,xpd=T)

plot(1,bty="l",xlab="",ylab="Mean Shannon Index", type="n",yaxt="n",ylim=c(0,2.5),xlim=c(0.2,5.5),cex.lab=1.2,axes=F)
clipplot(abline(h=0,col="grey75",lty=2),xlim=c(0.55,par("usr")[2]))
wvioplot(shan_strata$Shannon_Index_strata[shan_strata$Stratum=="Understorey"],at=1,add=T,
         col=alpha("#C8E8C4",0.9),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(shan_strata$Shannon_Index_strata[shan_strata$Stratum=="Shrub"],at=2,add=T,
         col=alpha("#3F828D",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(shan_strata$Shannon_Index_strata[shan_strata$Stratum=="Lower_Canopy"],at=3,add=T,
         col=alpha("#E6C6DC",0.9),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(shan_strata$Shannon_Index_strata[shan_strata$Stratum=="Upper_Canopy"],at=4,add=T,
         col=alpha("#D34A68",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(shan_strata$Shannon_Index_strata[shan_strata$Stratum=="Emergent"],at=5,add=T,
         col=alpha("#8986C8",0.9),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
abline(h=par("usr")[4],col="white",lwd=4);abline(h=par("usr")[3],col="white",lwd=4)
abline(v=par("usr")[1],col="white",lwd=4);abline(v=par("usr")[2],col="white",lwd=4)
axis(2,at=seq(0,2.5,0.5),cex.axis=1,las=1,col="grey20",col.axis="grey20",lwd=1.5,line=-1)
axis(1, at = seq(1, 5, 1), labels = FALSE, col = "grey20", col.axis = "grey20", lwd = 1.5, line = 0.5)
mtext(side = 1, text = c("Under-\nstorey","\nShrub","Lower\ncanopy","Upper\ncanopy","\nEmergent"),
      at = seq(1, 5, 1), padj = 1.5, cex = 0.9, col = "grey20")

dev.off()

#-----------------------------------------#
##### Plot DM abundance across Strata #####
#-----------------------------------------#

#-----------#
# Load Data #
#-----------#

model_data_strata <- read.csv("Data/model_data_strata.csv")

#-------------------#
# Prepare plot data #
#-------------------#

plot_data <- model_data_strata %>%
  filter(!is.na(Mean_Dormice), !is.na(Shannon_Index_strata)) %>%
  group_by(Coppicing_Plot, Stratum) %>%
  summarise(
    Mean_Shannon = mean(Shannon_Index_strata),
    Mean_Dormice = mean(Mean_Dormice),
    .groups = "drop"
  )

# Labels
stratum_labels <- c(
  "Understorey" = "Understorey",
  "Shrub" = "Shrub layer",
  "Lower_Canopy" = "Lower canopy",
  "Upper_Canopy" = "Upper canopy",
  "Emergent" = "Emergent"
)

plot_data$Stratum <- factor(
  plot_data$Stratum,
  levels = names(stratum_labels),
  labels = stratum_labels
)

# Colour palette
custom_plot_colors <- c("#C8E8C4","#3F828D","#E6C6DC","#D34A68","#8986C8")

#--------------#
# Store models #
#--------------#

emergent <- model_data_strata %>%
  filter(Stratum == "Emergent")
mod_em_nb <- glm.nb(Mean_Dormice ~ Shannon_Index_strata + num_unique_records,data = emergent)

Upper <- model_data_strata %>%
  filter(Stratum == "Upper_Canopy")
mod_up<-glm(Mean_Dormice ~ Shannon_Index_strata+ num_unique_records, data=Upper, family="poisson")

Lower <- model_data_strata %>%
  filter(Stratum == "Lower_Canopy")
mod_low_nb<-glm.nb(Mean_Dormice ~ Shannon_Index_strata + num_unique_records, data = Lower)

Shrub <- model_data_strata %>%
  filter(Stratum == "Shrub")
mod_shrub<-glm(Mean_Dormice ~ Shannon_Index_strata+ num_unique_records, data=Shrub, family="poisson")

Understorey <- model_data_strata %>%
  filter(Stratum == "Understorey")
mod_un<-glm(Mean_Dormice ~ Shannon_Index_strata+ num_unique_records, data=Understorey, family="poisson")

mods <- list(
  "Emergent" = mod_em_nb,
  "Upper canopy" = mod_up,
  "Lower canopy" = mod_low_nb,
  "Shrub layer" = mod_shrub,
  "Understorey" = mod_un
)

#----------------------#
# Prediction dataframe #
#----------------------#

pred_data <- lapply(names(mods), function(s){
  
  df <- plot_data %>% filter(Stratum == s)
  
  newdat <- data.frame(
    Shannon_Index_strata = seq(
      min(df$Mean_Shannon, na.rm = TRUE),
      max(df$Mean_Shannon, na.rm = TRUE),
      length.out = 100
    ),
    num_unique_records = median(model_data_strata$num_unique_records, na.rm = TRUE)
  )
  
  pr <- predict(mods[[s]], newdat, type = "link", se.fit = TRUE)
  
  data.frame(
    Stratum = s,
    Shannon_Index_strata = newdat$Shannon_Index_strata,
    fit = exp(pr$fit),
    lwr = exp(pr$fit - 1.96 * pr$se.fit),
    upr = exp(pr$fit + 1.96 * pr$se.fit)
  )
  
}) %>%
  bind_rows() %>%
  mutate(Stratum = factor(Stratum, levels = levels(plot_data$Stratum))) %>%
  arrange(Stratum, Shannon_Index_strata)

#---------------#
# Create figure #
#---------------#

cairo_pdf(
  "Figures/Q3DMbyStratum.pdf",
  width = 12,
  height = 10,
  family = "Cambria",
  bg = "white"
)

p <- ggplot(plot_data,
            aes(Mean_Shannon, Mean_Dormice)) +
  
  geom_point(
    aes(color = Stratum),
    size = 2,
    alpha = 0.6
  ) +
  
  geom_ribbon(
    data = pred_data,
    aes(
      x = Shannon_Index_strata,
      ymin = lwr,
      ymax = upr,
      fill = Stratum,
      group = Stratum
    ),
    alpha = 0.25,
    inherit.aes = FALSE
  ) +
  
  geom_line(
    data = pred_data,
    aes(
      x = Shannon_Index_strata,
      y = fit,
      color = Stratum,
      group = Stratum
    ),
    size = 1,
    inherit.aes = FALSE
  ) +
  
  facet_wrap(~ Stratum, scales = "free_x") +
  
  scale_color_manual(values = custom_plot_colors) +
  scale_fill_manual(values = custom_plot_colors) +
  
  labs(
    x = "Mean Shannon Index",
    y = "Mean Dormouse Abundance"
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(p)
dev.off()

#//////////////////////////////////////////////#
###### Q3 DM abundance across key species ######
#//////////////////////////////////////////////#

model_data_species_wide<- read.csv('Data/sp_data_shan.csv')

# Key species included in the analysis
key_species <- c(
  "Hazel",
  "Honeysuckle",
  "Bramble",
  "Oak",
  "Sycamore"
)

# Filter dataset
model_data_filtered <- model_data_species_wide %>%
  filter(Species %in% key_species) %>%
  mutate(Species = factor(Species, levels = key_species))

#--------------------------#
# Prepare data for plotting
#--------------------------#

plot_data <- model_data_filtered %>%
  filter(!is.na(Mean_Dormice), !is.na(Proportion_of_Species))

custom_colors <- c(
  "Honeysuckle" = "#D16159",
  "Hazel" = "#96BDD9",
  "Oak" = "#f48fb1",
  "Bramble" = "#8986C8",
  "Sycamore" = "#A3C57B"
)

#------#
# Plot #
#------#

cairo_pdf(
  "Figures/Key_species.pdf",
  width = 12,
  height = 10,
  family = "Cambria",
  bg = "white"
)

p2 <- ggplot(
  plot_data,
  aes(x = Proportion_of_Species * 100,
      y = Mean_Dormice,
      color = Species,
      fill = Species)
) +
  
  geom_point(
    size = 2,
    alpha = 0.7
  ) +
  
  stat_smooth(
    method = "glm",
    method.args = list(family = poisson(link = "log")),
    se = TRUE,
    linewidth = 1,
    alpha = 0.25
  ) +
  
  facet_wrap(~ Species, scales = "free_x") +
  
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  
  labs(
    x = "Proportion of Species (%)",
    y = "Mean Dormouse Abundance"
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(p2)
dev.off()

##%%%%%%%%%%%%%%%%%%%%%%%##

#### 11 Hypothesis Q1  ####

##%%%%%%%%%%%%%%%%%%%%%%%##

#-----------#
# Load data #
#-----------#

sp_data <- read.csv("Data/final_data_dedup.csv")

# Keep one row per coppicing plot
sp_data_unique <- sp_data[!duplicated(sp_data$Coppicing_Plot), ]

# Response variables
response_vars <- c(
  "Max_Canopy_Height_2020",
  "Gap_Proportion_2020",
  "LAI",
  "Intensity",
  "VCI",
  "Canopy_Volume_2020"
)

#---------------------#
# Fit Gaussian models #
#---------------------#

family_map <- list(
  Max_Canopy_Height_2020 = gaussian(),
  Gap_Proportion_2020    = gaussian(),
  LAI                    = gaussian(),
  Intensity              = gaussian(),
  VCI                    = gaussian(),
  Canopy_Volume_2020     = gaussian()
)

# Function to fit the final retained GAM
fit_final_gam <- function(response_var, data = sp_data_unique) {
  
  form <- as.formula(
    paste0(response_var, " ~ s(recovery_period, k = 4)")
  )
  
  fam <- family_map[[response_var]]
  
  model <- gam(
    formula = form,
    data = data,
    family = fam,
    method = "REML"
  )
  
  cat("\n\n==============================\n")
  cat("Final GAM for:", response_var, "\n")
  cat("Family:", fam$family, "\n")
  cat("==============================\n")
  print(summary(model))
  
  return(model)
}

# Fit final retained models
gam_models <- setNames(
  lapply(response_vars, fit_final_gam),
  response_vars
)

# Named colours for each response variable
panel_colors <- c(
  Max_Canopy_Height_2020 = "#f48fb1",
  Gap_Proportion_2020    = "#f48fb1",
  LAI                    = "#f48fb1",
  Intensity              = "#f48fb1",
  VCI                    = "#f48fb1",
  Canopy_Volume_2020     = "#f48fb1"
)

# Pretty labels for y-axis
response_labels <- c(
  Max_Canopy_Height_2020 = "Maximum canopy height (m)",
  Gap_Proportion_2020    = "Gap proportion (%)",
  LAI                    = "LAI",
  Intensity              = "Intensity",
  VCI                    = "VCI",
  Canopy_Volume_2020     = "Canopy volume (m³)"
)

plot_gam_smooth <- function(model, data, response_var, y_label, col, show_x = TRUE) {
  
  newdat <- data.frame(
    recovery_period = seq(
      min(data$recovery_period, na.rm = TRUE),
      max(data$recovery_period, na.rm = TRUE),
      length.out = 200
    )
  )
  
  # Predict on link scale, then transform back
  pred_link <- predict(model, newdata = newdat, type = "link", se.fit = TRUE)
  
  fit_link   <- pred_link$fit
  se_link    <- pred_link$se.fit
  lower_link <- fit_link - 1.96 * se_link
  upper_link <- fit_link + 1.96 * se_link
  
  linkinv <- model$family$linkinv
  
  sm <- newdat %>%
    dplyr::mutate(
      fitted = linkinv(fit_link),
      lower  = linkinv(lower_link),
      upper  = linkinv(upper_link)
    )
  
  obs <- data %>%
    dplyr::select(recovery_period, dplyr::all_of(response_var)) %>%
    dplyr::rename(response = dplyr::all_of(response_var))
  
  # Rescale canopy volume for plotting only
  scale_factor <- if (response_var == "Canopy_Volume_2020") 1000 else 1
  
  obs <- obs %>%
    dplyr::mutate(response = response / scale_factor)
  
  sm <- sm %>%
    dplyr::mutate(
      fitted = fitted / scale_factor,
      lower  = lower / scale_factor,
      upper  = upper / scale_factor
    )
  
  p <- ggplot() +
    geom_point(
      data = obs,
      aes(x = recovery_period, y = response),
      size = 2,
      alpha = 0.6,
      color = col
    ) +
    geom_ribbon(
      data = sm,
      aes(x = recovery_period, ymin = lower, ymax = upper),
      fill = col,
      alpha = 0.18
    ) +
    geom_line(
      data = sm,
      aes(x = recovery_period, y = fitted),
      color = col,
      linewidth = 1.2
    ) +
    labs(
      x = NULL,
      y = y_label
    ) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.25),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      axis.title.y = element_text(size = 12)
    )
  
  if (!show_x) {
    p <- p + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  }
  
  return(p)
}

# Generate plots
plots <- lapply(seq_along(names(gam_models)), function(i) {
  rv <- names(gam_models)[i]
  
  plot_gam_smooth(
    model = gam_models[[rv]],
    data = sp_data_unique,
    response_var = rv,
    y_label = response_labels[[rv]],
    col = panel_colors[[rv]],
    show_x = i %in% c(5, 6)   # only bottom row shows x-axis text/ticks
  )
})

# Combine into 2-column layout
final_plot <- (plots[[1]] | plots[[2]]) /
  (plots[[3]] | plots[[4]]) /
  (plots[[5]] | plots[[6]])

# Add shared bottom x-axis label
final_plot <- final_plot + patchwork::plot_annotation(
  theme = theme(
    plot.margin = margin(5.5, 5.5, 20, 5.5)
  )
)

# Save
cairo_pdf(
  "Figures/Structure.pdf",
  width = 9,
  height = 12,
  family = "Cambria",
  bg = "white"
)

xlab_plot <- ggplot() +
  annotate("text", x = 0.5, y = 0.8,
           label = "Time since coppicing (years)", size = 5) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_void()

print(
  patchwork::wrap_elements(final_plot) /
    xlab_plot +
    plot_layout(heights = c(20, 1))
)

dev.off()

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
  num_unique_records ="No. of dormice records",
  num_ET_points =" No. of ET records"
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

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### 12 Supplementary Stats  ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#-----------#
# Load Data #
#-----------#

final_data<-read.csv('Data/final_data_dedup.csv')

#-----------#
# Data Prep #
#-----------#

# Select only the relevant columns for uniqueness check - prevents pseudoreplication 
unique_rows <- !duplicated(final_data[, c("ET_ID")])

# Keep only unique rows in the dataframe
sp_data_unique <- final_data[unique_rows, ]

custom_labels <- c(
  Max_Canopy_Height_2020 = "Max Canopy Height",
  Gap_Proportion_2020 = "Canopy Gap Proportion",
  Intensity = "LiDAR Intensity",
  VCI = "VCI",
  LAI = "LAI",
  FT_10_20m = "First returns  (10–20m)",
  HALP = "HALP",
  Aspect_Cos = "Aspect (Cosine)",
  Slope = "Slope",
  Avg_Shannon_Index = "Mean Shannon Index",
  recovery_period = "Time since Coppicing",
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
                                     "Intensity", "VCI", "LAI",
                                     "FT_10_20m", "HALP", "Aspect_Cos", "num_ET_points", "num_unique_records",
                                      "Slope", "Avg_Shannon_Index", "recovery_period", 
                                     "FlowDir", "Plane_Curve", "Profile_Curve", "Canopy_Volume_2020", "Mean_Dormice")], 
                  use = "pairwise.complete.obs", method='spearman')  

# Apply custom labels
colnames(cor_matrix) <- custom_labels[colnames(cor_matrix)]
rownames(cor_matrix) <- custom_labels[rownames(cor_matrix)]

cor_matrix

#------#
# Plot #
#------#

# Open PDF device
cairo_pdf("Figures/CorrPlot.pdf",
          width = 12,
          height = 12,
          family = "Cambria",
          bg = "white")

par(mfrow = c(1,1), mar = c(4,4,1,1), las = 1, family = "Cambria")

# Define the custom color gradient
custom_colors <- colorRampPalette(brewer.pal(11, "PiYG"))(200)

# Set plot parameters
par(mfrow = c(1, 1), mar = c(6, 8, 2, 2), las = 1, xpd = TRUE, cex.axis = 1.4, ps = 10)

# Create the correlation plot
corrplot(cor_matrix, 
         type = 'lower', 
         tl.col = "black", 
         tl.cex = 1,
         tl.srt = 70,
         col = custom_colors,
         cl.cex = 1.5,
         cl.lim = c(-1, 1))

# Close the device
dev.off()


##%%%%%%%%%%%%%%%%%%%%%%%##

#### Combine Q4 Plots  ####

##%%%%%%%%%%%%%%%%%%%%%%%##

img1 <- image_read_pdf("Figures/CoeffsStruc.pdf", density = 300)
img2 <- image_read_pdf("Figures/CoeffsStrucDM.pdf", density = 300)
img3 <- image_read_pdf("Figures/CoeffsTopo.pdf", density = 300)
img4 <- image_read_pdf("Figures/CoeffsTopoDM.pdf", density = 300)

# Top row: 1 and 2
top_row <- image_append(c(img1, img2))

# Bottom row: 3 and 4
bottom_row <- image_append(c(img3, img4))

# Combine rows vertically
combined <- image_append(c(top_row, bottom_row), stack = TRUE)

image_write(combined, path = "Figures/Combined_2x2.pdf", format = "pdf")
  