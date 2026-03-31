##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

### Bontuchel Biodiversity Analysis ###

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

##%%%%%%%%%%%%%%%%%%%%%##

### Preliminary Stats ###

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
library(dplyr)

#/////////////////#
#### LOAD DATA ####
#/////////////////#

sp_data<-read.csv('Data/sp_data.csv')
ET_points<-vect('Data/ET_Points.gpkg')
bon_CHM_2020<- rast("Data/bon_CHM_2020_rep.tif")
bontuchel<-vect('Data/Bontuchel_Coppice_Panels.shp')

#//////////////////////////////////////////////#
#### CORRELATION OF BIODIVERSITY INDICATORS ####
#//////////////////////////////////////////////#

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
    Avg_Shannon_Index = mean(Shannon_Index, na.rm = TRUE)
  ) %>%
  ungroup()

# Print final dataset
print(total_dormice_df)
total_dormice_df$Mean_Dormice

# Add summary data to original sp_data without changing original rows/columns
final_dataset <- sp_data %>%
  left_join(total_dormice_df, by = "Coppicing_Plot")

#------#
# Save #
#------#

write.csv(final_dataset,"Data/final_dataset.csv",row.names = F)
write.csv(total_dormice_df,"Data/total_dormice_df.csv",row.names = F)

# Look at distribution of data
hist(total_dormice_df$Mean_Dormice) # skewed

#-------------#
# Fit a model #
#-------------#

fit <- lm(log(Avg_Shannon_Index) ~ Mean_Dormice, data = total_dormice_df)
summary(fit)
 
# Call:
#   lm(formula = log(Avg_Shannon_Index) ~ Mean_Dormice, data = total_dormice_df)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.81397 -0.04142  0.05983  0.13946  0.28986 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   0.69996    0.03511  19.936   <2e-16 ***
#   Mean_Dormice -0.02092    0.01928  -1.085    0.283    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2182 on 54 degrees of freedom
# Multiple R-squared:  0.02132,	Adjusted R-squared:  0.003194 
# F-statistic: 1.176 on 1 and 54 DF,  p-value: 0.2829

shapiro.test(fit$residuals)

# Shapiro-Wilk normality test
# 
# data:  fit$residuals
# W = 0.84006, p-value = 3.025e-06

# Not normal even with log, skip modelling and do a corr
# Correlation
cor.test(total_dormice_df$Mean_Dormice, 
         total_dormice_df$Avg_Shannon_Index, 
         method = "spearman")

# Spearman's rank correlation rho
# 
# data:  total_dormice_df$Mean_Dormice and total_dormice_df$Avg_Shannon_Index
# S = 30077, p-value = 0.8381
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#         rho 
# -0.02792462 

#%%%%%%%%%%%%%%%%%%%#
###### Plotting #####
#%%%%%%%%%%%%%%%%%%%#

# Plot figure
cairo_pdf("Figures/DormiceIndicators.pdf",
          width = 6,
          height = 6,
          family = "Cambria",
          bg = "white")
par(family = "Cambria")
par(mfrow=c(1,1),mar=c(6,8,2,2),las=1,xpd=T)

ggplot(total_dormice_df, aes(x = Mean_Dormice, y = Avg_Shannon_Index)) +
  geom_point(color = "pink", size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "hotpink", fill = "hotpink", alpha = 0.2) +
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

#//////////////////////////#
#### GROUND TRUTH LiDAR ####
#//////////////////////////#

#------#
# Crop #
#------#

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

#----------------------------#
##### Mean Canopy Height #####
#----------------------------#

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

length(sp_data_unique$ET_ID) #62

#-------------------#
# Check assumptions #
#-------------------#

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

# Shapiro-Wilk normality test
# 
# data:  residuals
# W = 0.95592, p-value = 0.02595

# Assumptions not met, do correlation

#-------------#
# Correlation #
#-------------#

# Spearmans correlation coefficient
height<-cor.test(top_layer$ET_Canopy_Height, top_layer$LiDAR_Mean_Canopy_Height_2020, method='spearman', use = "complete.obs")
height

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

mean_diff  # this is the average difference (9.98)
sd_diff    # standard deviation for variability (5.20)

#---------------------------#
##### Mean Canopy Cover #####
#---------------------------#

# When working out canopy cover from LiDAR, we will use the threshold of cover from the min height of the top layer, ie: 10
# We use the min from top layer of et height because lidar will go through all layers of canopy so we can't differentiate.
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
final_dataset$LiDAR_Mean_Canopy_cover_2020 <- mean_canopy_cover[match(final_dataset$ET_Rec_Date, buffer_15m$ET_Rec_Date)]

summary(final_dataset$ET_Canopy_Cover)
summary(final_dataset$LiDAR_Mean_Canopy_cover_2020)

# Select only the relevant columns for uniqueness check - prevents pseudoreplication 
unique_rows <- !duplicated(final_dataset[, c("ET_ID")])

# Keep only unique rows in the dataframe
sp_data_unique <- final_dataset[unique_rows, ]

#-------------------#
# Check assumptions #
#-------------------#

plot(sp_data_unique$ET_Canopy_Cover, sp_data_unique$LiDAR_Mean_Canopy_cover_2020)
hist(sp_data_unique$ET_Canopy_Cover)
hist(sp_data_unique$LiDAR_Mean_Canopy_cover_2020)

# Calculate the residuals (difference between the two height measurements)
residuals <- sp_data_unique$ET_Canopy_Cover - sp_data_unique$LiDAR_Mean_Canopy_cover_2020
hist(residuals, main = "Residuals (Ground - LiDAR)", xlab = "Residuals", col = "pink")

# Summary of residuals
summary(residuals)
shapiro.test(residuals) 

# Shapiro-Wilk normality test
# 
# data:  residuals
# W = 0.95297, p-value = 0.01859

# Assumption of normality not met

#-------------#
# Correlation #
#-------------#

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

mean_diff  # this is the average difference (-9.15)
sd_diff    # standard deviation for variability (31.64)

#%%%%%%%%%%%%%%%%%%%#
###### Plotting #####
#%%%%%%%%%%%%%%%%%%%#

length(sp_data_unique$ET_Canopy_Cover)
length(sp_data_unique$LiDAR_Mean_Canopy_cover_2020)

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
     col  = adjustcolor("#7A0016", alpha.f = 0.6),
     bty  = "l")

abline(0, 1, col = "#D66F90", lwd = 3)

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
     col  = adjustcolor("#7A0016", alpha.f = 0.6),
     bty  = "l")

abline(0, 1, col = "#D66F90", lwd = 3)

# Close PDF
dev.off()
