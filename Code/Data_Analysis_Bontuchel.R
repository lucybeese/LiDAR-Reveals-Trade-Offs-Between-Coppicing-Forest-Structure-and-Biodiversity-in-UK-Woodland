##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### Bontuchel Biodiversity Analysis ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

##%%%%%%%%%%%%%%%%%%%##

#### Data Analysis ####

##%%%%%%%%%%%%%%%%%%%##

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
library(DHARMa)
library(mgcViz)
library(readr)
library(gratia)
library(magick)
library(patchwork)
library(e1071)

##%%%%%%%%%%%%%%%%%%%%%%%##

#### STEP 1: LOAD DATA ####

##%%%%%%%%%%%%%%%%%%%%%%%##

#-----------------------#
# Set working directory #
#-----------------------#

setwd("~/Library/CloudStorage/OneDrive-UniversityofBristol/Github/Bontuchel_Chapter")

#-----------------------#
##### Load all data #####
#-----------------------#

bon_DTM_2020<- rast("Data/bon_DTM_2020_rep.tif")
bon_CHM_2020<- rast("Data/bon_CHM_2020_rep.tif")
bon_df<-read.csv('Data/bon_df.csv')
boxes_df<-read.csv('Data/boxes_df.csv')
bontuchel<-vect('Data/Bontuchel_Coppice_Panels.shp')
boxes<-vect('Data/Boxcoordinates.shp')
sp_data1<-read.csv('Data/final_dataset_1.csv')
ET_points<-vect('Data/ET_Points.gpkg')
bon_CHM_2020<- rast("Data/bon_CHM_2020_rep.tif")

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### STEP 2: CORRELATION OF BIODIVERSITY INDICATORS ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

# Load data
sp_data1<-read.csv('Data/final_dataset_1.csv')

# Ensure numeric conversion
sp_data1 <- sp_data1 %>%
  mutate(
    No_of_Dormice = as.numeric(No_of_Dormice),
    Shannon_Index = as.numeric(Shannon_Index)
  )

# Remove duplicate entries based on Coppicing_Plot, Year_DM_Recorded, and Dor_Box
sp_data_unique1 <- sp_data1 %>%
  distinct(Coppicing_Plot, Year_DM_Recorded, Dor_Box, ET_ID, .keep_all = TRUE)

total_dormice_df <- sp_data_unique1 %>%
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
sp_data <- sp_data1 %>%
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

# Open a PDF device
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### STEP 3: GROUND TRUTH LiDAR ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

# Load data
sp_data<-read.csv('Data/final_dataset.csv')

# Crop
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

mean_diff  # this is the average difference
sd_diff    # standard deviation for variability

#--------------------#
# Mean Canopy Cover #
#--------------------#

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

mean_diff  # this is the average difference 
sd_diff    # standard deviation for variability

#------#
# Plot #
#------#

length(sp_data_unique$ET_Canopy_Cover)
length(sp_data_unique$LiDAR_Mean_Canopy_cover_2020)

# Open a PDF device
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

# Close PDF
dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### STEP 4: RS METRICS AND BIODIVERSITY ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

# GAMS WITH NB ARE JUST GLMS IF NO SMOOTHING REMEMBER!

# Load data
sp_data<-read.csv('final_dataset.csv')

#----------------------------------------------------#
# Function to help select the most appropriate model #
#----------------------------------------------------#

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

#//////////////////////////////#

##### Q4: Canopy structure #####

#//////////////////////////////#

#-------#
# Model #
#-------#

# Load data
sp_data<-read.csv('final_dataset.csv')

# Count unique ET_ID per Coppicing_Plot
sp_data_count <- sp_data %>%
  group_by(Coppicing_Plot) %>%
  summarise(num_ET_points = n_distinct(ET_ID), .groups = "drop")

# Join back to sp_data_unique
sp_data <- sp_data_count %>%
  left_join(sp_data, by = "Coppicing_Plot")

# Select only the relevant columns for uniqueness check
unique_rows <- !duplicated(sp_data[, c("Coppicing_Plot")])
sp_data_unique <- sp_data[unique_rows, ]

# Variables of interest (including response)
predictors <- c("Mean_Canopy_Height_2020", "Perc_under_2m_2020", "Can_cover_2020", "Height_cv_2020", 
                "Max_Canopy_Height_2020", "Gap_Proportion_2020", "LAI", "Intensity", "VCI", 
                "Canopy_Volume_2020", "FT_5_10m", "FT_1.37m", "FT_10_20m", "num_ET_points")

response_var<-c("Mean_Shannon_Index")

best_model <- select_best_model(sp_data_unique, predictors,response_var, cor_threshold = 0.7, shape="linear")
# Also tried with cor threshold of 0.5 but model wasnt as stong

summary(sp_data$Mean_Shannon_Index)
summary(sp_data$Intensity)

#------#
# Plot #
#------#

# Scale predictors so they are comparable
# Identify numeric columns
numeric_vars <- sapply(sp_data_unique, is.numeric)  
sp_data_unique_scaled <- sp_data_unique  # Create a copy of the dataset
sp_data_unique_scaled[numeric_vars] <- scale(sp_data_unique[numeric_vars])  # Scale only numeric columns

model_scaled <- lm(Mean_Shannon_Index ~ Max_Canopy_Height_2020 + LAI + Intensity + VCI + Canopy_Volume_2020 + FT_10_20m + num_ET_points, 
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
new_names <- c("Max Canopy Height", "LAI", "Intensity", "VCI", "Canopy Volume", "First returns 10-20m", "Number of ET points")

axis(2, at = seq_along(coef_names), labels = FALSE, las = 2)
text(x = rep(par("usr")[1]-0.2, length(new_names)), y = seq_along(new_names), 
     labels = new_names, adj = 1, xpd = TRUE, cex = cex_text)

arrows(lower, seq_along(estimates), upper, seq_along(estimates), length = 0.1, angle = 90, code = 3)
abline(v = 0, lty = 2, col = "hotpink")  # Reference line

# Add "a)" in the top-left corner (inside the plot margin)
# Use xpd=NA to allow drawing in the figure margin
text(x = par("usr")[1], y = par("usr")[4] + 0.5, labels = "a)", adj = c(0, 1), cex = cex_text, xpd = NA)

dev.off()

#////////////////////////#

##### Q4: Topography #####

#////////////////////////#

# Load data
sp_data<-read.csv('final_dataset.csv')
colnames(sp_data)

# Count unique ET_ID per Coppicing_Plot
sp_data_count <- sp_data %>%
  group_by(Coppicing_Plot) %>%
  summarise(num_ET_points = n_distinct(ET_ID), .groups = "drop")

# Join back to sp_data_unique
sp_data <- sp_data_count %>%
  left_join(sp_data, by = "Coppicing_Plot")

# Select only the relevant columns for uniqueness check
unique_rows <- !duplicated(sp_data[, c("Coppicing_Plot")])
sp_data_unique <- sp_data[unique_rows, ]
colnames(sp_data_unique)
hist(sp_data_unique$Mean_Shannon_Index)

# Variables of interest (including response)
predictors <- c("HALP","Aspect_Cos", "TPI", "Roughness", 
                "TRI", "Elevation", "Slope", "FlowDir","Plane_Curve","Profile_Curve",
                "Mean_Curve","Solar_Radiation", "num_ET_points")

response_var<-c("Mean_Shannon_Index")

best_model_topo <- select_best_model(sp_data_unique, predictors, response_var,cor_threshold = 0.7)
hist(sp_data_unique$Mean_Shannon_Index)
best_model_topo<- glm(Mean_Shannon_Index~HALP+ Aspect_Cos+Slope+ FlowDir+ Plane_Curve+  Profile_Curve + num_ET_points , family = Gamma(link = "log"),sp_data_unique)
vif(best_model_topo)

summary(best_model_topo)

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
  Mean_Shannon_Index~HALP+ Aspect_Cos+Slope+ FlowDir+ Plane_Curve+  Profile_Curve + num_ET_points,
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

# Open a PDF device
pdf("~/Library/CloudStorage/OneDrive-UniversityofBristol/GitHub/Thesis/Chapter4/Figures/CoeffsTopo.pdf", width = 7, height = 5)

par(mar = c(5, 16, 4, 2)) 

# Double text size by setting cex (default is 1)
cex_text <- 1.5

plot(estimates, seq_along(estimates), pch = 16,  xlim = c(- 0.2, 0.2),
     xlab = "Estimate", ylab = "", axes = FALSE, cex = cex_text, cex.lab = cex_text)
axis(1,cex.axis = cex_text)

# New labels
new_names<- c("HALP","Aspect", "Slope", "Flow Direction", "Plane Curvature","Profile Curvature", "Number of ET points")

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

#////////////////////////#

##### Q2: Management #####

#////////////////////////#

sp_data<-read.csv('final_dataset.csv')

# Get recovery period
sp_data <- sp_data %>%
  mutate(recovery_period = 2020 - Coppicing_year)

# Count unique ET_ID per Coppicing_Plot
et_counts <- sp_data %>%
  group_by(Coppicing_Plot) %>%
  summarise(num_ET_points = n_distinct(ET_ID), .groups = "drop")

# Join back to sp_data_unique
sp_data <- sp_data %>%
  left_join(et_counts, by = "Coppicing_Plot")

# Select only the relevant columns for uniqueness check
unique_rows <- !duplicated(sp_data[, c("Coppicing_Plot")])

# Keep only unique rows in the dataframe
sp_data_unique <- sp_data[unique_rows, ]

# Variables of interest
predictors <- c("recovery_period")

response_var<-c("Mean_Shannon_Index")

model_manag<- select_best_model(sp_data_unique, predictors, response_var,cor_threshold = 0.7, shape="quadratic")

# Above I expected an n or u-shape, I then tested for flattening (don't include:)
# gam_model <- gam(Mean_Shannon_Index ~ s(recovery_period, k = 4),
#                  data = sp_data_unique, family = Gamma(link = "log"))
# summary(gam_model)
# plot(gam_model)
# plot(gam_model, residuals = TRUE, pch = 16)
# shapiro.test(gam_model$residuals)
# gam.check(gam_model)
# no flattening 

gam_model2 <- gam(Mean_Shannon_Index ~ s(recovery_period, k=6), 
                  data = sp_data_unique, family=Gamma(link="log"))
summary(gam_model2)
plot(gam_model2, residuals=TRUE)

# Add the number of samples per plot as a covariate
gam_model3 <- gam(
  Mean_Shannon_Index ~ s(recovery_period, k = 6) + num_ET_points,
  data = sp_data_unique,
  family = Gamma(link = "log")
)
summary(gam_model3)
plot(gam_model3, residuals=TRUE)

dev_resid <- deviance(gam_model3)
df_resid <- df.residual(gam_model3)

dispersion <- dev_resid / df_resid


# Add the number of samples per plot as a covariate
gam_model <- gam(
  Mean_Shannon_Index ~ recovery_period + I(recovery_period^2)  + num_ET_points,
  data = sp_data_unique,
  family = Gamma(link = "log")
)
summary(gam_model)

#------#
# Plot #
#------#

# Open a PDF device
pdf("~/Library/CloudStorage/OneDrive-UniversityofBristol/GitHub/Thesis/Chapter4/Figures/Management.pdf", width = 3.5, height = 4)

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
    x = "Time since coppicing (Years)",
    y = "Mean Shannon index"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.25),
    axis.text.x = element_text(hjust = 1, size=11),
    axis.text.y = element_text(size=11),
    axis.title = element_text(size = 12),
  )


dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### STEP 5: RS METRICS AND DORMICE ####

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

# Helper function: Fit Poisson or NB GLM depending on overdispersion
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

#//////////////////////////#

##### Q4:Canopy structure #####

#//////////////////////////#

#-------#
# Model #
#-------#

# Load data
sp_data<-read.csv('final_dataset.csv')

# Count unique combinations per Coppicing_Plot
plot_summary <- sp_data %>%
  distinct(Coppicing_Plot, Dor_Box, Year_DM_Recorded, No_of_Dormice) %>%
  group_by(Coppicing_Plot) %>%
  summarise(
    num_unique_records = n()
  )

# Join back to unique plot data
sp_data <- sp_data %>%
  distinct(Coppicing_Plot, .keep_all = TRUE) %>%
  filter(Coppicing_Plot %in% plot_summary$Coppicing_Plot) %>%
  left_join(plot_summary, by = "Coppicing_Plot")


# Select only the relevant columns for uniqueness check
unique_rows <- !duplicated(sp_data[, c("Coppicing_Plot")])
sp_data_unique <- sp_data[unique_rows, ]
colnames(sp_data_unique)

# Variables of interest
predictors <- c("Mean_Canopy_Height_2020", "Perc_under_2m_2020", "Can_cover_2020", "Height_cv_2020", 
                "Max_Canopy_Height_2020", "Gap_Proportion_2020", "LAI", "Intensity", "VCI", 
                "Canopy_Volume_2020", "FT_5_10m", "FT_1.37m", "FT_10_20m", "num_unique_records" )

response_var<-c("Mean_Dormice")

# Create Spearman correlation matrix
cor_matrix <- cor(sp_data_unique[, predictors], use = "pairwise.complete.obs", method = "spearman")
valid_combinations <- list(); counter <- 1
cor_threshold = 0.7

# Identify predictor sets with all pairwise correlations < threshold (0.7)
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

# Loop over all valid combinations to find best (lowest AIC) model
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

# Check dispersion ratio manually
disp<- deviance(best_model$model) / df.residual(best_model$model)
disp
# Residual normality check (not crucial for GLMs, but included anyway)
shapiro.test(best_model$model$residuals)

# dispersion >1.5: use negative binomial model instead:
gam_nb<-glm.nb(Mean_Dormice ~  Max_Canopy_Height_2020+ FT_5_10m+ LAI+ Intensity+ VCI+ Canopy_Volume_2020 +num_unique_records, data = sp_data_unique)
summary(gam_nb)
disp<- deviance(gam_nb) / df.residual(gam_nb)
disp

# Explore distribution
hist(sp_data_unique$Mean_Dormice) # alot of zeros so need to check for zero inflation:
simres <- simulateResiduals(fittedModel = gam_nb, plot = TRUE)
testZeroInflation(simres)
# p-value = 0.68 → This is a high p-value, which means there is no evidence of a statistically significant difference 
# between the number of observed and expected zeros.

#------#
# Plot #
#------#

# Scale predictors so they are comparable
sp_data_unique_scaled <- sp_data_unique
sp_data_unique_scaled[, predictors] <- scale(sp_data_unique_scaled[, predictors])

model_scaled <- glm.nb(Mean_Dormice ~  Max_Canopy_Height_2020+ FT_5_10m+ LAI+ Intensity+ VCI+ Canopy_Volume_2020 +num_unique_records,  data = sp_data_unique_scaled)

summary(model_scaled)

# Extract coefficients and confidence intervals
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

# Open a PDF device
pdf("~/Library/CloudStorage/OneDrive-UniversityofBristol/GitHub/Thesis/Chapter4/Figures/CoeffsCanopyStructure_DM.pdf", width = 7, height = 5)

cex_text<-1.5

par(mar = c(5, 16, 4, 2)) 
xlim_range <- range(c(lower, upper)) * 1.1  # Add 10% padding
plot(estimates, seq_along(estimates), pch = 16, xlim = c(-1,1.5),
     xlab = "Estimate", ylab = "", axes = FALSE, cex=cex_text, cex.lab = cex_text)

axis(1, cex.axis = cex_text)

# Remove the axis(2) call and replace it with text for custom labels
axis(2, at = seq_along(coef_names), labels = FALSE, las = 2)

# Labels on y-axis (custom text)
new_names <- c("Max Canopy Height","First returns 5-10m","LAI", "Intensity", "VCI", "Canopy Volume",  "No. dormice records")

axis(2, at = seq_along(coef_names), labels = FALSE, las = 2)
text(x = rep(par("usr")[1]-0.2, length(new_names)), y = seq_along(new_names), 
     labels = new_names, adj = 1, xpd = TRUE, cex = cex_text)


arrows(lower, seq_along(estimates), upper, seq_along(estimates), length = 0.1, angle = 90, code = 3)
abline(v = 0, lty = 2, col = "hotpink")  # Reference line

text(x = par("usr")[1], y = par("usr")[4] + 0.5, labels = "b)", adj = c(0, 1), cex = cex_text, xpd = NA)


dev.off()

#////////////////////#

##### Q4: Topography #####

#////////////////////#

# Load data
sp_data<-read.csv('final_dataset.csv')

# Count unique combinations per Coppicing_Plot
plot_summary <- sp_data %>%
  distinct(Coppicing_Plot, Dor_Box, Year_DM_Recorded, No_of_Dormice) %>%
  group_by(Coppicing_Plot) %>%
  summarise(
    num_unique_records = n()
  )

# Join back to unique plot data
sp_data <- sp_data %>%
  distinct(Coppicing_Plot, .keep_all = TRUE) %>%
  filter(Coppicing_Plot %in% plot_summary$Coppicing_Plot) %>%
  left_join(plot_summary, by = "Coppicing_Plot")

# Select only the relevant columns for uniqueness check
unique_rows <- !duplicated(sp_data[, c("Coppicing_Plot")])
sp_data_unique <- sp_data[unique_rows, ]

# Variables of interest
predictors <- c("HALP","Aspect_Cos", "TPI", "Roughness", 
                "TRI", "Elevation", "Slope", "FlowDir","Plane_Curve","Profile_Curve",
                "Mean_Curve","Solar_Radiation", "num_unique_records")

response_var<-c("Mean_Dormice")

# lm with log-transformed response
best_model_topo<- lm(log(Mean_Dormice+1) ~HALP+ Aspect_Cos+ TPI+ TRI+FlowDir+ Profile_Curve, sp_data_unique)
summary(best_model_topo)
shapiro.test(best_model_topo$residuals)

# Calculate correlation matrix
cor_matrix <- cor(sp_data_unique[, predictors], use = "pairwise.complete.obs", method = "spearman")

# Find valid predictor subsets with all pairwise Spearman correlations < 0.7
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

# Fit all valid models and choose the best
best_model <- NULL
best_aic <- Inf
best_combo <- NULL

for (combo in valid_combinations) {
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

# dispersion >1.5: use neg bin instead
gam_nb<-glm.nb(Mean_Dormice ~  HALP+ Aspect_Cos+ FlowDir+ Profile_Curve+ Mean_Curve+ Solar_Radiation+num_unique_records, data = sp_data_unique)
summary(gam_nb)
gam.check(gam_nb)

# Proportion of observed zeros
mean(sp_data_unique$Mean_Dormice == 0)

# Check for zero inflation
simres <- simulateResiduals(fittedModel = gam_nb, plot = TRUE)
testZeroInflation(simres)
# p-value = 0.64 → This is a high p-value, which means there is no evidence of a statistically significant difference 
# between the number of observed and expected zeros.

#------#
# Plot #
#------#

# Scale predictors so they are comparable
sp_data_unique_scaled <- sp_data_unique
sp_data_unique_scaled[, predictors] <- scale(sp_data_unique_scaled[, predictors])

model_scaled <- glm.nb(Mean_Dormice ~  HALP+ Aspect_Cos+ FlowDir+ Profile_Curve+ Mean_Curve+ Solar_Radiation+num_unique_records, data = sp_data_unique_scaled)

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

# Open a PDF device
pdf("~/Library/CloudStorage/OneDrive-UniversityofBristol/GitHub/Thesis/Chapter4/Figures/CoeffsTopoDM.pdf", width = 7, height = 5)

cex_text<-1.5
par(mar = c(5, 16, 4, 2)) 

xlim_range <- range(c(lower, upper)) * 1.1  # Add 10% padding
plot(estimates, seq_along(estimates), pch = 16, xlim = c(-1,1.5),
     xlab = "Estimate", ylab = "", axes = FALSE, cex=cex_text, cex.lab = cex_text)

axis(1, at = seq(-2, 2, by = 1), cex.axis=cex_text)  # Adjust ticks


# Remove the axis(2) call and replace it with text for custom labels
axis(2, at = seq_along(coef_names), labels = FALSE, las = 2)

# Labels on y-axis (custom text)
new_names <- c("HALP","Aspect","Flow direction", "Profile curvature", "Mean curvature", "Solar radiation",  "No. dormice records")

axis(2, at = seq_along(coef_names), labels = FALSE, las = 2)
text(x = rep(par("usr")[1]-0.2, length(new_names)), y = seq_along(new_names), 
     labels = new_names, adj = 1, xpd = TRUE, cex = cex_text)

arrows(lower, seq_along(estimates), upper, seq_along(estimates), length = 0.1, angle = 90, code = 3)
abline(v = 0, lty = 2, col = "hotpink")  # Reference line

text(x = par("usr")[1], y = par("usr")[4] + 0.5, labels = "d)", adj = c(0, 1), cex = cex_text, xpd = NA)


dev.off()

#////////////////////////#

##### Q2: Management #####

#////////////////////////#

# Load data
sp_data<-read.csv('final_dataset.csv')

# Filter out rows with any NAs in the key columns
sp_data_clean <- sp_data %>%
  filter(!is.na(Dor_Box), !is.na(Year_DM_Recorded), !is.na(No_of_Dormice))

# Count unique combinations per Coppicing_Plot
plot_summary <- sp_data_clean %>%
  distinct(Coppicing_Plot, Dor_Box, Year_DM_Recorded, No_of_Dormice) %>%
  group_by(Coppicing_Plot) %>%
  summarise(
    num_unique_records = n()
  )

# Join back to unique plot data
sp_data_unique <- sp_data %>%
  distinct(Coppicing_Plot, .keep_all = TRUE) %>%
  filter(Coppicing_Plot %in% plot_summary$Coppicing_Plot) %>%
  left_join(plot_summary, by = "Coppicing_Plot")

# Get recovery period
sp_data_unique <- sp_data_unique %>%
  mutate(recovery_period = 2020 - Coppicing_year)

predictors <- c("recovery_period")
response_var <- "Mean_Dormice"

# lm with quadratic term
best_model_manage <-lm(log(Mean_Dormice+1) ~ recovery_period + I(recovery_period^2), data=sp_data_unique)
summary(best_model_manage)

# Normality check
shapiro.test(best_model_manage$residuals)

# Find best count model
best_model_manage<- choose_best_count_glm(response_var = response_var, predictors = predictors, data = sp_data_unique)

# Poisson model with quadratic
mod<-glm(Mean_Dormice ~ recovery_period + I(recovery_period^2), data=sp_data_unique, family='poisson')
summary(mod)

# Check for overdispersion
dispersion_ratio <- sum(residuals(mod, type = "pearson")^2) / mod$df.residual
dispersion_ratio # overdispersed

# Fit NB GLM with quadratic term
nb_model <- glm.nb(Mean_Dormice ~ recovery_period + I(recovery_period^2), data = sp_data_unique)
summary(nb_model)

# Zero inflation check
simres <- simulateResiduals(fittedModel = nb_model, plot = TRUE)
testZeroInflation(simres) # no zero inflation

# Try a poisson gam 
gam_count <- gam(Mean_Dormice ~ recovery_period + I(recovery_period^2),
                 family = poisson(link = "log"),
                 data = sp_data_unique)
summary(gam_count)

# Check for overdispersion
dispersion <- deviance(gam_count) / df.residual(gam_count)
dispersion # overdispersed

# Fit neg binomial gam
gam_nb <- gam(Mean_Dormice ~ recovery_period+ I(recovery_period^2), family = nb(), data = sp_data_unique)
summary(gam_nb)

# Check for overdispersion
dispersion <- deviance(gam_nb) / df.residual(gam_nb)
dispersion # not overdispersed

# Zero inflation check
simres <- simulateResiduals(fittedModel = gam_nb, plot = TRUE)
testZeroInflation(simres) # no zero inflation

# Use the neg binoial GLM , statistically sound and simpler than GAM:
# summary(nb_model)
# Call:
#   glm.nb(formula = Mean_Dormice ~ recovery_period + I(recovery_period^2), 
#          data = sp_data_unique, init.theta = 3.143384603, link = log)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)          -0.051956   0.658536  -0.079    0.937
# recovery_period       0.224357   0.159802   1.404    0.160
# I(recovery_period^2) -0.012282   0.008363  -1.469    0.142
# 
# (Dispersion parameter for Negative Binomial(3.1434) family taken to be 1)
# 
# Null deviance: 33.674  on 28  degrees of freedom
# Residual deviance: 31.456  on 26  degrees of freedom
# AIC: 112.71
# 
# Number of Fisher Scoring iterations: 1

# Try with covariate:

# Fit neg binomial
nb_modelcov <- glm.nb(Mean_Dormice ~ recovery_period + I(recovery_period^2) + num_unique_records, data = sp_data_unique)
summary(nb_model)

# Check for overdispersion
dispersion <- deviance(gam_nb) / df.residual(gam_nb)
dispersion # not overdispersed

# Zero inflation check
simres <- simulateResiduals(fittedModel = gam_nb, plot = TRUE)
testZeroInflation(simres) # no zero inflation

#------#
# Plot #
#------#

# Open a PDF device
pdf("~/Library/CloudStorage/OneDrive-UniversityofBristol/GitHub/Thesis/Chapter4/Figures/ManagementDM.pdf", width = 3.5, height = 4)

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
    x = "Time since coppicing (Years)",
    y = "Mean Dormice Abundance"
  ) +
  scale_y_continuous(breaks = seq(0, 8, by = 2))+  # Or higher if needed
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.25),
    axis.text.x = element_text(hjust = 1, size=11),
    axis.text.y = element_text(size=11),
    axis.title = element_text(size = 12),
  )

dev.off()

#---------------------------#
##### Combining Figures #####
#---------------------------#

#---------------------------------#
# Canopy Structure and Topography #
#---------------------------------#

# Specify file paths
pdf_files <- c("~/Library/CloudStorage/OneDrive-UniversityofBristol/Github/Thesis/Chapter4/Figures/CoeffsCanopyStructure.pdf", 
               "~/Library/CloudStorage/OneDrive-UniversityofBristol/Github/Thesis/Chapter4/Figures/CoeffsCanopyStructure_DM.pdf", 
               "~/Library/CloudStorage/OneDrive-UniversityofBristol/Github/Thesis/Chapter4/Figures/CoeffsTopo.pdf", 
               "~/Library/CloudStorage/OneDrive-UniversityofBristol/Github/Thesis/Chapter4/Figures/CoeffsTopoDM.pdf")

# Read first page of each PDF
images <- lapply(pdf_files, function(file) {
  image_read_pdf(file, density = 150)[1]
})

images <- lapply(images, function(img) {
  image_crop(img, geometry = "980x980+40+40")# Keep tight control
})
# Combine images in 2 rows, 2 columns (2x2 grid)
row1 <- image_append(c(images[[1]], images[[2]]), stack = FALSE)  # horizontally append first two
row2 <- image_append(c(images[[3]], images[[4]]), stack = FALSE)  # horizontally append next two

# Stack the two rows vertically
final_image <- image_append(c(row1, row2), stack = TRUE)

# Save the combined image as a PDF
image_write(final_image, path = "~/Library/CloudStorage/OneDrive-UniversityofBristol/Github/Thesis/Chapter4/Figures/Combined_2x2.pdf", format = "pdf")

#------------#
# Management #
#------------#

# Specify file paths
pdf_files <- c("~/Library/CloudStorage/OneDrive-UniversityofBristol/Github/Thesis/Chapter4/Figures/Management.pdf", 
               "~/Library/CloudStorage/OneDrive-UniversityofBristol/Github/Thesis/Chapter4/Figures/ManagementDM.pdf")

# Read first page of each PDF
images <- lapply(pdf_files, function(file) {
  image_read_pdf(file, density = 600)[1]
})

# Combine images in 2 rows, 2 columns (2x2 grid)
row1 <- image_append(c(images[[1]], images[[2]]), stack = FALSE)  # horizontally append first two

# Stack the two rows vertically
final_image <- image_append(c(row1), stack = TRUE)

# Save the combined image as a PDF
image_write(final_image, path = "~/Library/CloudStorage/OneDrive-UniversityofBristol/Github/Thesis/Chapter4/Figures/ManagementCombo.pdf", format = "pdf")

#%%%%%%%%%%%%%%%%%%%%%%%%%%##

#### STEP 6: FULL MODELS ####

##%%%%%%%%%%%%%%%%%%%%%%%%%##

#----------------------#
##### Mean Dormice #####
#----------------------#

# Load data
sp_data<-read.csv('final_dataset.csv')

# Count unique combinations per Coppicing_Plot
plot_summary <- sp_data %>%
  distinct(Coppicing_Plot, Dor_Box, Year_DM_Recorded, No_of_Dormice) %>%
  group_by(Coppicing_Plot) %>%
  summarise(
    num_unique_records = n()
  )

# Join back to unique plot data
sp_data <- sp_data %>%
  distinct(Coppicing_Plot, .keep_all = TRUE) %>%
  filter(Coppicing_Plot %in% plot_summary$Coppicing_Plot) %>%
  left_join(plot_summary, by = "Coppicing_Plot")

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
                "Mean_Curve","Solar_Radiation", "recovery_period", "Mean_Shannon_Index", "num_unique_records")

response_var<- "Mean_Dormice"

# Correlation matrix
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
dispersion #not over dispersed

# > mod<-best_model$model
# > cat("✅ Best model type:", best_model$type, "\n")
# ✅ Best model type: Poisson 
# > cat("📉 AIC:", best_aic, "\n")
# 📉 AIC: 128.5043 
# > cat("🔢 Predictors:", paste(best_combo, collapse = ", "), "\n")
# 🔢 Predictors: Max_Canopy_Height_2020, LAI, Intensity, VCI, Canopy_Volume_2020, FT_10_20m, HALP, Aspect_Cos, FlowDir, Profile_Curve, Mean_Curve, Solar_Radiation, recovery_period, Mean_Shannon_Index, num_unique_records 
# > summary(best_model$model)

# Call:
#   glm(formula = formula, family = poisson(), data = sp_data_unique)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)            -4.485e+01  1.193e+01  -3.759 0.000170 ***
#   Max_Canopy_Height_2020  3.699e-01  1.174e-01   3.151 0.001625 ** 
#   LAI                     2.071e+00  6.270e-01   3.303 0.000958 ***
#   Intensity              -1.261e-03  1.968e-03  -0.641 0.521809    
# VCI                     3.945e+01  1.177e+01   3.351 0.000806 ***
#   Canopy_Volume_2020     -3.720e-05  1.365e-05  -2.725 0.006422 ** 
#   FT_10_20m              -1.106e-01  2.996e-02  -3.693 0.000222 ***
#   HALP                   -9.466e-02  3.207e-02  -2.952 0.003162 ** 
#   Aspect_Cos             -1.399e+00  4.154e-01  -3.367 0.000759 ***
#   FlowDir                 7.140e-02  6.408e-02   1.114 0.265119    
# Profile_Curve          -1.706e+02  1.284e+02  -1.329 0.183846    
# Mean_Curve              8.002e+02  1.898e+02   4.215 2.50e-05 ***
#   Solar_Radiation         7.228e+00  4.691e+00   1.541 0.123381    
# recovery_period        -1.399e-02  4.762e-02  -0.294 0.768931    
# Mean_Shannon_Index     -2.309e+00  6.791e-01  -3.399 0.000675 ***
#   num_unique_records      7.625e-02  1.702e-02   4.481 7.44e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 126.96  on 55  degrees of freedom
# Residual deviance:  36.12  on 40  degrees of freedom
# AIC: 128.5
# 
# Number of Fisher Scoring iterations: 6

#------#
# Plot #
#------#

# Scale predictors so they are comparable

sp_data_unique_scaled <- sp_data_unique
sp_data_unique_scaled[, predictors] <- scale(sp_data_unique_scaled[, predictors])

model_scaled <- glm(Mean_Dormice ~ Max_Canopy_Height_2020  +      LAI + Intensity + VCI + Canopy_Volume_2020 + FT_10_20m+ HALP+ Aspect_Cos + FlowDir + Mean_Curve + Profile_Curve  + Mean_Shannon_Index + Solar_Radiation + recovery_period +num_unique_records, family='poisson', sp_data_unique_scaled)

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

# Open a PDF device
pdf("~/Library/CloudStorage/OneDrive-UniversityofBristol/Github/Thesis/Chapter4/Figures/FullModelDM.pdf", width = 900, height = 900)

par(mar = c(8, 10,2, 1)) 
xlim_range <- range(c(-4, 4))
plot(estimates, seq_along(estimates), pch = 16, xlim = xlim_range,
     xlab = "Estimate", ylab = "", axes = FALSE)
axis(1, at = seq(-4, 4, by = 1))  # Adjust ticks

# New labels
new_names<- c("Max Canopy Height","LAI", "Intensity","VCI","Canopy volume", "First returns 10-20m", "HALP" , "Aspect", "Flow Direction", "Mean Curvature", "Profile Curvature", "Mean Shannon Index", "Solar radiation", "Time since coppicing", "No. dormice records")

# Remove the axis(2) call and replace it with text for custom labels
axis(2, at = seq_along(coef_names), labels = FALSE, las = 2)
text(x = rep(par("usr")[1]-0.2, length(new_names)), y = seq_along(new_names), 
     labels = new_names, adj = 1.1, xpd = TRUE)

arrows(lower, seq_along(estimates), upper, seq_along(estimates), length = 0.1, angle = 90, code = 3)
abline(v = 0, lty = 2, col = "hotpink")  # Reference line

dev.off()

#----------------------#
##### Mean Shannon #####
#----------------------#

# Load data
sp_data<-read.csv('final_dataset.csv')

# Count unique ET_ID per Coppicing_Plot
sp_data_count <- sp_data %>%
  group_by(Coppicing_Plot) %>%
  summarise(num_ET_points = n_distinct(ET_ID), .groups = "drop")

# Join back to sp_data_unique
sp_data <- sp_data_count %>%
  left_join(sp_data, by = "Coppicing_Plot")

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
                "Mean_Curve","Solar_Radiation", "recovery_period", "Mean_Dormice", "num_ET_points")

response_var<- "Mean_Shannon_Index"

best_full_model_DM<- select_best_model(sp_data_unique, predictors, response_var,cor_threshold = 0.7)

# ✅ Best: Gaussian GLM 
# Formula: Mean_Shannon_Index ~ Max_Canopy_Height_2020 + Gap_Proportion_2020 +      LAI + Intensity + VCI + Canopy_Volume_2020 + Aspect_Cos +      Elevation + Slope + FlowDir + Plane_Curve + Profile_Curve +      recovery_period + Mean_Dormice + num_ET_points 
# Shapiro p = 0.174 
# AIC = 12.06 
# Reason: Gaussian GLM with normal residuals 
# 
# 
# Call:
#   glm(formula = formula, family = gaussian, data = data)
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             2.658e+00  1.542e+00   1.724   0.0924 .  
# Max_Canopy_Height_2020  3.348e-03  1.695e-02   0.197   0.8445    
# Gap_Proportion_2020     3.234e-03  3.874e-03   0.835   0.4088    
# LAI                    -1.144e-01  1.080e-01  -1.058   0.2962    
# Intensity              -7.507e-05  4.281e-04  -0.175   0.8617    
# VCI                     1.919e+00  1.131e+00   1.696   0.0976 .  
# Canopy_Volume_2020      1.865e-06  2.366e-06   0.788   0.4352    
# Aspect_Cos             -7.485e-02  5.147e-02  -1.454   0.1537    
# Elevation              -8.783e-03  4.547e-03  -1.932   0.0605 .  
# Slope                  -1.837e-02  1.067e-02  -1.722   0.0928 .  
# FlowDir                 2.842e-03  9.952e-03   0.286   0.7767    
# Plane_Curve             1.720e+01  9.710e+00   1.772   0.0841 .  
# Profile_Curve          -2.117e+00  1.917e+01  -0.110   0.9126    
# recovery_period         1.568e-02  9.243e-03   1.696   0.0976 .  
# Mean_Dormice           -5.399e-02  2.338e-02  -2.309   0.0262 *  
#   num_ET_points          -6.700e-01  1.176e-01  -5.697 1.26e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 0.05539536)
# 
# Null deviance: 7.6618  on 55  degrees of freedom
# Residual deviance: 2.2158  on 40  degrees of freedom
# AIC: 12.056
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

# Open a PDF device
pdf("~/Library/CloudStorage/OneDrive-UniversityofBristol/Github/Thesis/Chapter4/Figures/FullModelShan.pdf", width = 900, height = 900, res=150)

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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#### STEP 7: Q3 FOREST VERTICAL DIVERITY COVAR ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Load data
sp_data <- read.csv("Data/final_dataset.csv")

# Filter out rows with any NAs in the key columns
sp_data_clean <- sp_data %>%
  filter(!is.na(Dor_Box), !is.na(Year_DM_Recorded), !is.na(No_of_Dormice))

# Count unique combinations per Coppicing_Plot
plot_summary <- sp_data_clean %>%
  distinct(Coppicing_Plot, Dor_Box, Year_DM_Recorded, No_of_Dormice) %>%
  group_by(Coppicing_Plot) %>%
  summarise(
    num_unique_records = n()
  )

sp_data_unique <- sp_data_clean %>%
  filter(Coppicing_Plot %in% plot_summary$Coppicing_Plot) %>%
  left_join(plot_summary, by = "Coppicing_Plot")

# Remove duplicate species entries per plot × layer
sp_data_unique <- sp_data_unique %>%
  distinct(Coppicing_Plot, Forest_layer, Species, .keep_all = TRUE)

# Define strata based on ET_layer_height
sp_data_unique$Stratum <- cut(
  sp_data_unique$ET_layer_height,
  breaks = c(-Inf, 1.37, 5, 10, 20, Inf),
  labels = c("Understorey", "Shrub", "Lower_Canopy", "Upper_Canopy", "Emergent")
)

# Scale proportions so they sum to 1 within each plot × layer
sp_data_corrected <- sp_data_unique %>%
  group_by(Coppicing_Plot, Stratum) %>%
  mutate(Sum_Proportion = sum(Proportion_of_Species, na.rm = TRUE)) %>%
  mutate(Proportion_of_Species = ifelse(Sum_Proportion > 1, Proportion_of_Species / Sum_Proportion, Proportion_of_Species)) %>%
  ungroup()

# Calculate Shannon Index per ET_ID
shannon_per_point <- sp_data_corrected %>%
  group_by(Coppicing_Plot, Stratum, ET_ID) %>%
  mutate(
    p = Proportion_of_Species / sum(Proportion_of_Species, na.rm = TRUE),
    p_log_p = -p * log(p)
  ) %>%
  summarise(Shannon_ET = sum(p_log_p, na.rm = TRUE), .groups = "drop")

# Average Shannon Index across ET_IDs per plot × stratum
shannon_index_normalised <- shannon_per_point %>%
  group_by(Coppicing_Plot, Stratum) %>%
  summarise(Shannon_Index_strata = mean(Shannon_ET, na.rm = TRUE), .groups = "drop")

# Merge back with species proportions (1 row per species per stratum)
model_data_species_wide <- sp_data_corrected %>%
  left_join(shannon_index_normalised, by = c("Coppicing_Plot", "Stratum")) %>%
  dplyr::select(Coppicing_Plot, Stratum, Species, Proportion_of_Species, Shannon_Index_strata)

# Add dormice data
dorm <- read.csv("Data/dormice_data_Panels.csv")

total_dormice_sum <- dorm %>%
  group_by(Coppicing_Plot) %>%
  summarise(Total_dormice_sum = sum(as.numeric(Total_dormice), na.rm = TRUE))

# Merge dormice total and mean
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

model_data_species_wide <- model_data_species_wide %>%
  left_join(plot_summary, by = "Coppicing_Plot")

# Function for model choice:
strata_list <- c("Emergent", "Upper_Canopy", "Lower_Canopy", "Shrub", "Understorey")

for (stratum in strata_list) {
  cat("Stratum:", stratum, "\n")
  data_sub <- model_data_species_wide %>% filter(Stratum == stratum)
  
  # Poisson model
  mod_pois <- glm(Mean_Dormice ~ Shannon_Index_strata + num_unique_records, data = data_sub, family = "poisson")
  aic_pois <- AIC(mod_pois)
  
  # Negative binomial model
  mod_nb <- try(glm.nb(Mean_Dormice ~ Shannon_Index_strata + num_unique_records, data = data_sub), silent = TRUE)
  
  if (inherits(mod_nb, "try-error")) {
    cat("NB model failed to converge for", stratum, "\n\n")
    next
  }
  
  aic_nb <- AIC(mod_nb)
  
  cat("Poisson AIC: ", aic_pois, "\n")
  cat("Neg Binomial AIC: ", aic_nb, "\n")
  
  if (aic_nb < aic_pois) {
    cat("Negative Binomial model preferred (lower AIC)\n\n")
  } else {
    cat("Poisson model preferred (lower AIC)\n\n")
  }
}

#------------------#
##### Emergent #####
#------------------#

# Filter for Emergent stratum
emergent <- model_data_species_wide %>%
  filter(Stratum == "Emergent")

mod<-glm(Mean_Dormice ~ Shannon_Index_strata + num_unique_records, data=emergent, family="poisson")
summary(mod)

# Check dispersion ratio manually
disp<- deviance(mod) / df.residual(mod)
disp # not  Overdispersed

#---------------#
##### Upper #####
#---------------#

Upper <- model_data_species_wide %>%
  filter(Stratum == "Upper_Canopy")
mod<-glm(Mean_Dormice ~ Shannon_Index_strata+ num_unique_records, data=Upper, family="poisson")
summary(mod)
shapiro.test(mod$residuals)

# Check dispersion ratio manually
disp<- deviance(mod) / df.residual(mod)
disp # Not Overdispersed

#---------------#
##### Lower #####
#---------------#

Lower <- model_data_species_wide %>%
  filter(Stratum == "Lower_Canopy")
modellow<-glm(Mean_Dormice ~ Shannon_Index_strata+ num_unique_records, data=Lower, family="poisson")
summary(modellow)
shapiro.test(mod$residuals)
# Check dispersion ratio manually

disp<- deviance(modellow) / df.residual(modellow)
disp # not Overdispersed

#nb
nb<-glm.nb(Mean_Dormice ~ Shannon_Index_strata + num_unique_records, data = Lower)
summary(nb)

#---------------#
##### Shrub #####
#---------------#

Shrub <- model_data_species_wide %>%
  filter(Stratum == "Shrub")
mod<-glm(Mean_Dormice ~ Shannon_Index_strata+ num_unique_records, data=Shrub, family="poisson")
shapiro.test(mod$residuals)

# Check dispersion ratio manually
disp<- deviance(mod) / df.residual(mod)
disp # not Overdispersed

#nb
nb<-glm.nb(Mean_Dormice ~ Shannon_Index_strata + num_unique_records, data = Shrub)
summary(nb)

#---------------------#
##### Understorey #####
#---------------------#

Understorey <- model_data_species_wide %>%
  filter(Stratum == "Understorey")
mod<-glm(Mean_Dormice ~ Shannon_Index_strata+ num_unique_records, data=Understorey, family="poisson")
summary(mod)
shapiro.test(mod$residuals)


# Check dispersion ratio manually
disp<- deviance(mod) / df.residual(mod)
disp # not Overdispersed

#nb
nb<-glm.nb(Mean_Dormice ~ Shannon_Index_strata + num_unique_records, data = Understorey)
summary(nb)


# Diversity in the upper canopy reflects no of dormice found, lower diversity= more mice



#------#
# PLOT #
#------#

# Add Coppicing_Year to wide data
model_data_species_wide <- model_data_species_wide %>%
  left_join(
    sp_data %>% dplyr::select(Coppicing_Plot, Coppicing_year) %>% distinct(),
    by = "Coppicing_Plot"
  )

# Define consistent color palette (blues, purples, pinks)
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
write.csv(model_data_species_wide,"Data/sp_data_shan.csv",row.names = F)

#-----------------------------#
# Shannon Index across Strata #
#-----------------------------#

# For clarity, keep unique Shannon index values per plot and stratum
shan_strata <- model_data_species_wide %>%
  distinct(Coppicing_Plot, Stratum, Shannon_Index_strata) %>%
  filter(!is.na(Shannon_Index_strata))  # remove any NA

# Reorder Stratum factor for meaningful order on y-axis (optional)
shan_strata$Stratum <- factor(shan_strata$Stratum,
                              levels = c("Understorey", "Shrub", "Lower_Canopy", "Upper_Canopy", "Emergent"))

# Plot mean Shannon Index on x, stratum on y
p <- ggplot(shan_strata, aes(x = Shannon_Index_strata, y = Stratum)) +
  geom_jitter(height = 0.15, size = 3, alpha = 0.7, color = "pink") +  # jitter to avoid overlap
  stat_summary(fun = mean, geom = "point", shape = 18, size = 5, color = "hotpink") +  # mean per stratum
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

#----------------------------#
# DM abundance across Strata #
#----------------------------#

# Prepare plot data: one row per Coppicing_Plot × Stratum with Mean Shannon and Dormice
plot_data <- model_data_species_wide %>%
  filter(!is.na(Mean_Dormice), !is.na(Shannon_Index_strata)) %>%
  group_by(Coppicing_Plot, Stratum) %>%
  summarise(
    Mean_Shannon = mean(Shannon_Index_strata, na.rm = TRUE),
    Mean_Dormice = mean(Mean_Dormice, na.rm = TRUE),
    .groups = "drop"
  )

# Define labels for strata
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#### STEP 8: Q3: DORMICE PREFERENCE ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Load data
model_data<-read.csv('Data/final_dataset.csv')
head(model_data)

# Count unique ET_ID per Coppicing_Plot
model_data_count <- model_data %>%
  group_by(Coppicing_Plot) %>%
  summarise(num_ET_points = n_distinct(ET_ID), .groups = "drop")

# Join back to sp_data_unique
model_data <- model_data_count %>%
  left_join(model_data, by = "Coppicing_Plot")

model_data$num_ET_points

# Summarise prop of each species
model_data_species <- model_data %>%
  group_by(Coppicing_Plot, Species) %>%
  summarise(
    Proportion_of_Species = sum(Proportion_of_Species, na.rm = TRUE),
    Mean_Dormice = first(Mean_Dormice),
    Coppicing_year = first(Coppicing_year),
    num_ET_points = first(num_ET_points),  # keep this column
    .groups = "drop"
  )

head(model_data_species)
model_data_species$num_ET_points

# Rescale so sum of prop of all species across each coppicing panel is never above 1
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

# Poisson GLM
modlog<-glm(Mean_Dormice~Species, data=model_data_filtered, family="poisson")
summary(modlog)
shapiro.test(modlog$residuals)

# Calculate dispersion
dispersion <- deviance(modlog) / df.residual(modlog)
dispersion # Overdispersed

# Negative binomial 
modlog<-glm.nb(Mean_Dormice~Species, data=model_data_filtered)
summary(modlog)
shapiro.test(modlog$residuals)

# With covariate: 
# Poisson GLM
modlog<-glm(Mean_Dormice~Species+num_ET_points, data=model_data_filtered, family="poisson")
summary(modlog)
shapiro.test(modlog$residuals)
dispersion <- deviance(modlog) / df.residual(modlog)
dispersion # Overdispersed

# Negative binomial 
modlog<-gam(Mean_Dormice~Species+num_ET_points, data=model_data_filtered, family=nb())
summary(modlog)
shapiro.test(modlog$residuals)

mod_nb <- glm.nb(Mean_Dormice ~ Species + num_ET_points, data = model_data_filtered)
summary(mod_nb)

# Multiply proportion by 100 for percentages
plot_data <- model_data_filtered %>%
  mutate(Proportion_of_Species = Proportion_of_Species * 100) %>%
  filter(!is.na(Mean_Dormice), !is.na(Proportion_of_Species)) %>%
  mutate(Species = factor(Species, levels = names(species_color_map)))

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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#### STEP 9: Q1: MANAGEMENT AND STRUCTURE ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Load data
sp_data <- read.csv('Data/final_dataset.csv')

# Keep one row per Coppicing_Plot
sp_data_unique <- sp_data[!duplicated(sp_data$Coppicing_Plot), ]

# Add recovery period (years since coppicing)
sp_data_unique <- sp_data_unique %>%
  mutate(recovery_period = 2020 - Coppicing_year)

# Variables to model
response_vars <- c(
  "Max_Canopy_Height_2020", "Gap_Proportion_2020", 
  "LAI", "Intensity", "VCI", "Canopy_Volume_2020"
)

# Function to fit a GAM for a single response variable
fit_gam <- function(response_var, data = sp_data_unique) {
  formula <- as.formula(paste0(response_var, " ~ s(recovery_period, k = 4)"))
  
  # Fit GAM with REML method for smoothness estimation
  model <- gam(formula, data = data, method = "REML")
  
  cat("\n\n=== GAM for", response_var, "===\n")
  print(summary(model))
  
  return(model)
}

# Fit GAM models for all response variables and store in a named list
gam_models <- setNames(lapply(response_vars, fit_gam), response_vars)

# Check residuals and diagnostics for each model individually
for (rv in names(gam_models)) {
  model <- gam_models[[rv]]
  cat("\n\nResidual diagnostics for:", rv, "\n")
  
  # Plot residuals vs fitted values and QQ plots
  par(mfrow = c(2, 2))
  plot(model, residuals = TRUE, pch = 16, main = rv)
  gam.check(model)
  
  # Shapiro-Wilk normality test on residuals
  shapiro_test <- shapiro.test(residuals(model))
  cat("Shapiro-Wilk test p-value:", shapiro_test$p.value, "\n")
}

# Try Gamma family with log link for VCI, because residuals were off under Gaussian.

# For VCI
gam_vci <- gam(VCI ~ s(recovery_period, k=4), 
               family = Gamma(link = "log"), 
               data = sp_data_unique, method = "REML")

summary(gam_vci)

# Plot residuals to check fit

plot(gam_vci, residuals=TRUE, pch=16, cex=0.5)
gam.check(gam_vci)

# VCI was the only model with a residual distribution that clearly violated normality under your pre-defined conservative threshold (p < 0.001).
# All the other models, while not perfect, had p-values well above 0.001, so did not meet my own criterion for transformation.

#----------#
# Plotting #
#----------#

summary(sp_data_unique$Max_Canopy_Height_2020)

# Function to create ggplot for GAM smooths with predictions on the response scale
plot_gam_smooth <- function(model, y_label) {
  # Generate a sequence of recovery_period values for prediction
  newdat <- data.frame(
    recovery_period = seq(min(model$model$recovery_period),
                          max(model$model$recovery_period),
                          length.out = 200)
  )
  
  # Get predictions on the response scale
  pred <- predict(model, newdata = newdat, type = "response", se.fit = TRUE)
  
  sm <- newdat %>%
    mutate(
      fitted = pred$fit,
      se     = pred$se.fit,
      lower  = fitted - 2 * se,
      upper  = fitted + 2 * se
    )
  
  # Build ggplot
  p <- ggplot(sm, aes(x = recovery_period)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "red", alpha = 0.15) +
    geom_line(aes(y = fitted), color = "hotpink", size = 1.2) +
    labs(x = "Time since coppicing (years)", y = y_label) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", size = 0.25),
      axis.text.x = element_text(hjust = 1, size=11),
      axis.text.y = element_text(size=11),
      axis.title = element_text(size = 12)
    )
  
  return(p)
}


# Fit GAM models for all structural variables
gam_models <- list(
  Max_Canopy_Height_2020 = gam(Max_Canopy_Height_2020 ~ s(recovery_period, k = 4),
                               data = sp_data_unique, method = "REML"),
  Gap_Proportion_2020 = gam(Gap_Proportion_2020 ~ s(recovery_period, k = 4),
                            data = sp_data_unique, method = "REML"),
  LAI = gam(LAI ~ s(recovery_period, k = 4),
            data = sp_data_unique, method = "REML"),
  Intensity = gam(Intensity ~ s(recovery_period, k = 4),
                  data = sp_data_unique, method = "REML"),
  VCI = gam(VCI ~ s(recovery_period, k = 4),
            family = Gamma(link = "log"),
            data = sp_data_unique, method = "REML"),
  Canopy_Volume_2020 = gam(Canopy_Volume_2020 ~ s(recovery_period, k = 4),
                           data = sp_data_unique, method = "REML")
)

# Generate plots for all GAMs
plots <- lapply(names(gam_models), function(rv) {
  plot_gam_smooth(gam_models[[rv]], y_label = rv)
})

# Pretty labels for y-axis
  response_labels <- list(
    "Maximum Canopy Height (m)",
    "Gap Proportion (%)",
    "LAI",
    "Intensity",
    "VCI",
    expression("Canopy Volume (m"^3*")")
  )
  
# Generate plots with pretty labels
plots <- mapply(function(model, label) {
  plot_gam_smooth(model, y_label = label)
}, model = gam_models, label = response_labels, SIMPLIFY = FALSE)


# Combine into 2-column layout
final_plot <- (plots[[1]] | plots[[2]]) / (plots[[3]] | plots[[4]]) / (plots[[5]] | plots[[6]])

# Print and save
print(final_plot)
ggsave(filename = "~/Library/CloudStorage/OneDrive-UniversityofBristol/GitHub/Thesis/Chapter4/Figures/Structure.pdf",
       plot = final_plot, width = 7, height = 8)

#%%%%%%%%%%%%%%%%%%%%%#

#### STEP 10: BIAS ####

#%%%%%%%%%%%%%%%%%%%%%#

# Here we want to look at how well did each sample represent each plot

# Load data
sp_data<-read.csv('Data/final_dataset.csv')
ET_points<-vect('Data/ET_Points.gpkg')
bon_CHM_2020<- rast("Data//bon_CHM_2020_rep.tif")
bontuchel<-vect('Data/Bontuchel_Coppice_Panels.shp')
bon_DTM_2020<- rast("Data/bon_DTM_2020_rep.tif")
boxes<-vect('Data/Boxcoordinates.shp')

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

# Create data frame to store data
plot_data <- data.frame(matrix(NA, nrow=dim(buffer_15m)[1], ncol=26))
names(plot_data)<-c("Coppicing_Plot","radius_plotvalue_ID","ET_Rec_Date", "radius_Mean_Canopy_Height_2020", "radius_Perc_under_2m_2020", 
                    "radius_Can_cover_2020", "radius_Height_cv_2020",
                    "radius_Elevation","radius_Slope","radius_TRI","radius_Max_Canopy_Height_2020", "radius_FT_5_10m","radius_FT_1.37m",
                    "radius_FT_10_20m","radius_Intensity","radius_Aspect","radius_TPI", 'radius_Roughness', "radius_LAI", "radius_VCI", "radius_Canopy_Volume_2020",
                    "radius_FlowDir", "radius_Plane_Curve", "radius_Profile_Curve", "radius_Mean_Curve", "radius_Solar_Radiation")

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
  
  min(first_returns$Z)
  max(first_returns$Z)
  mean(first_returns$Z)
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#### STEP 11: TOPOGRAPHIC HETEROGENEITY ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Ranges 
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#### STEP 12: CANOPY HETEROGENEITY ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

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
cano_vars <- c("Mean_Canopy_Height_2020", "Perc_under_2m_2020", "Can_cover_2020", "Height_cv_2020", "Max_Canopy_Height_2020", "FT_5_10m", "FT_1.37m", "FT_10_20m", 
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

summary(sp_data$Mean_Canopy_Height_2020)

#//////////////////////////////////////////#

#### STEP 12: HYPOTHESIS RESULT SUMMARY ####

#//////////////////////////////////////////#

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

##### Q1: MANAGEMENT AND STRUCTURE #####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Load data
sp_data <- read.csv('Data/final_dataset.csv')

# Keep one row per Coppicing_Plot
sp_data_unique <- sp_data[!duplicated(sp_data$Coppicing_Plot), ]

# Add recovery period (years since coppicing)
sp_data_unique <- sp_data_unique %>%
  mutate(recovery_period = 2020 - Coppicing_year)

# Variables to model
response_vars <- c(
  "Max_Canopy_Height_2020", "Gap_Proportion_2020", 
  "LAI", "Intensity", "VCI", "Canopy_Volume_2020"
)

# Function to fit a GAM for a single response variable
fit_gam <- function(response_var, data = sp_data_unique) {
  formula <- as.formula(paste0(response_var, " ~ s(recovery_period, k = 4)"))
  
  # Fit GAM with REML method for smoothness estimation
  model <- gam(formula, data = data, method = "REML")
  
  cat("\n\n=== GAM for", response_var, "===\n")
  print(summary(model))
  
  return(model)
}

# Fit GAM models for all response variables and store in a named list
gam_models <- setNames(lapply(response_vars, fit_gam), response_vars)

# Check residuals and diagnostics for each model individually
for (rv in names(gam_models)) {
  model <- gam_models[[rv]]
  cat("\n\nResidual diagnostics for:", rv, "\n")
  
  # Plot residuals vs fitted values and QQ plots
  par(mfrow = c(2, 2))
  plot(model, residuals = TRUE, pch = 16, main = rv)
  gam.check(model)
  
  # Shapiro-Wilk normality test on residuals
  shapiro_test <- shapiro.test(residuals(model))
  cat("Shapiro-Wilk test p-value:", shapiro_test$p.value, "\n")
}

# Gamma family with log link for VCI, because residuals were off under Gaussian.

# For VCI
gam_vci <- gam(VCI ~ s(recovery_period, k=4), 
               family = Gamma(link = "log"), 
               data = sp_data_unique, method = "REML")

summary(gam_vci)

# 5 GAMS, all log link, VCI is gamma instead of gaussian:

# === GAM for Max_Canopy_Height_2020 ===
#   
#   Family: gaussian 
# Link function: identity 
# 
# Formula:
#   Max_Canopy_Height_2020 ~ s(recovery_period, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  21.7657     0.4663   46.68   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value  
# s(recovery_period)   1      1 3.094  0.0842 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0367   Deviance explained = 5.42%
# -REML = 148.13  Scale est. = 12.176    n = 56
# 
# 
# === GAM for Gap_Proportion_2020 ===
#   
#   Family: gaussian 
# Link function: identity 
# 
# Formula:
#   Gap_Proportion_2020 ~ s(recovery_period, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   34.150      2.951   11.57 3.04e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value  
# s(recovery_period)   1      1 5.978  0.0178 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.083   Deviance explained = 9.97%
# -REML = 247.77  Scale est. = 487.73    n = 56
# 
# 
# === GAM for LAI ===
#   
#   Family: gaussian 
# Link function: identity 
# 
# Formula:
#   LAI ~ s(recovery_period, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  3.66246    0.08186   44.74   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value   
# s(recovery_period) 2.457  2.789 5.003 0.00367 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.211   Deviance explained = 24.7%
# -REML = 55.772  Scale est. = 0.37524   n = 56
# 
# 
# === GAM for Intensity ===
#   
#   Family: gaussian 
# Link function: identity 
# 
# Formula:
#   Intensity ~ s(recovery_period, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   490.89      12.83   38.26   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value
# s(recovery_period)   1  1.001 1.447   0.234
# 
# R-sq.(adj) =  0.00807   Deviance explained = 2.61%
# -REML = 327.13  Scale est. = 9217.3    n = 56
# 
# 
# === GAM for VCI ===
#   
#   Family: gaussian 
# Link function: identity 
# 
# Formula:
#   VCI ~ s(recovery_period, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.930891   0.006793     137   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value   
# s(recovery_period)   1      1 7.573 0.00805 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.107   Deviance explained = 12.3%
# -REML = -80.231  Scale est. = 0.0025838  n = 56
# 
# 
# === GAM for Canopy_Volume_2020 ===
#   
#   Family: gaussian 
# Link function: identity 
# 
# Formula:
#   Canopy_Volume_2020 ~ s(recovery_period, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    42081       2706   15.55   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value
# s(recovery_period) 1.019  1.037 0.046   0.854
# 
# R-sq.(adj) =  -0.0171   Deviance explained = 0.17%
# -REML = 616.12  Scale est. = 4.1005e+08  n = 56
#
# === GAM for VCI ===
#   
# VCI ~ s(recovery_period, k = 4)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.071815   0.007365   -9.75 1.66e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value   
# s(recovery_period)   1      1 7.393 0.00878 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.107   Deviance explained = 11.4%
# -REML = -78.106  Scale est. = 0.0030379  n = 56

#/////////////////////////////////#

##### Q2: Management/ Shannon #####

#/////////////////////////////////#

sp_data<-read.csv('final_dataset.csv')

# Get recovery period
sp_data <- sp_data %>%
  mutate(recovery_period = 2020 - Coppicing_year)

# Count unique ET_ID per Coppicing_Plot
et_counts <- sp_data %>%
  group_by(Coppicing_Plot) %>%
  summarise(num_ET_points = n_distinct(ET_ID), .groups = "drop")

# Join back to sp_data_unique
sp_data <- sp_data %>%
  left_join(et_counts, by = "Coppicing_Plot")

# Select only the relevant columns for uniqueness check
unique_rows <- !duplicated(sp_data[, c("Coppicing_Plot")])

# Keep only unique rows in the dataframe
sp_data_unique <- sp_data[unique_rows, ]

# Add the number of samples per plot as a covariate
gam_model <- gam(
  Mean_Shannon_Index ~ recovery_period + I(recovery_period^2)  + num_ET_points,
  data = sp_data_unique,
  family = Gamma(link = "log")
)
summary(gam_model)

# Family: Gamma 
# Link function: log 
# 
# Formula:
#   Mean_Shannon_Index ~ recovery_period + I(recovery_period^2) + 
#   num_ET_points
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           1.112e+00  8.731e-02  12.740  < 2e-16 ***
#   recovery_period       5.149e-03  1.539e-02   0.335    0.739    
# I(recovery_period^2)  4.759e-05  8.388e-04   0.057    0.955    
# num_ET_points        -4.253e-01  6.308e-02  -6.741 1.28e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# 
# R-sq.(adj) =  0.355   Deviance explained = 39.1%
# GCV = 0.029446  Scale est. = 0.020976  n = 56

#/////////////////////////////////#

##### Q3: Management/ Dormice #####

#/////////////////////////////////#

# Load data
sp_data<-read.csv('final_dataset.csv')

# Filter out rows with any NAs in the key columns
sp_data_clean <- sp_data %>%
  filter(!is.na(Dor_Box), !is.na(Year_DM_Recorded), !is.na(No_of_Dormice))

# Count unique combinations per Coppicing_Plot
plot_summary <- sp_data_clean %>%
  distinct(Coppicing_Plot, Dor_Box, Year_DM_Recorded, No_of_Dormice) %>%
  group_by(Coppicing_Plot) %>%
  summarise(
    num_unique_records = n()
  )

# Join back to unique plot data
sp_data_unique <- sp_data %>%
  distinct(Coppicing_Plot, .keep_all = TRUE) %>%
  filter(Coppicing_Plot %in% plot_summary$Coppicing_Plot) %>%
  left_join(plot_summary, by = "Coppicing_Plot")

# Get recovery period
sp_data_unique <- sp_data_unique %>%
  mutate(recovery_period = 2020 - Coppicing_year)

# Fit neg binomial
nb_modelcov <- glm.nb(Mean_Dormice ~ recovery_period + I(recovery_period^2) + num_unique_records, data = sp_data_unique)
summary(nb_modelcov)

# Call:
#   glm.nb(formula = Mean_Dormice ~ recovery_period + I(recovery_period^2) + 
#            num_unique_records, data = sp_data_unique, init.theta = 4.155253067, 
#          link = log)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)          -0.659015   0.737403  -0.894   0.3715  
# recovery_period       0.255494   0.153681   1.663   0.0964 .
# I(recovery_period^2) -0.014091   0.008164  -1.726   0.0843 .
# num_unique_records    0.028814   0.018348   1.570   0.1163  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Negative Binomial(4.1553) family taken to be 1)
# 
# Null deviance: 36.543  on 28  degrees of freedom
# Residual deviance: 31.554  on 25  degrees of freedom
# AIC: 112.29
# 
# Number of Fisher Scoring iterations: 1
# 
# 
# Theta:  4.16 
# Std. Err.:  3.54 
# 
# 2 x log-likelihood:  -102.291 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#####  Q3 FOREST VERTICAL DIVERITY #####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Load data
sp_data <- read.csv("Data/final_dataset.csv")

# Filter out rows with any NAs in the key columns
sp_data_clean <- sp_data %>%
  filter(!is.na(Dor_Box), !is.na(Year_DM_Recorded), !is.na(No_of_Dormice))

# Count unique combinations per Coppicing_Plot
plot_summary <- sp_data_clean %>%
  distinct(Coppicing_Plot, Dor_Box, Year_DM_Recorded, No_of_Dormice) %>%
  group_by(Coppicing_Plot) %>%
  summarise(
    num_unique_records = n()
  )

sp_data_unique <- sp_data_clean %>%
  filter(Coppicing_Plot %in% plot_summary$Coppicing_Plot) %>%
  left_join(plot_summary, by = "Coppicing_Plot")

# Remove duplicate species entries per plot × layer
sp_data_unique <- sp_data_unique %>%
  distinct(Coppicing_Plot, Forest_layer, Species, .keep_all = TRUE)

# Define strata based on ET_layer_height
sp_data_unique$Stratum <- cut(
  sp_data_unique$ET_layer_height,
  breaks = c(-Inf, 1.37, 5, 10, 20, Inf),
  labels = c("Understorey", "Shrub", "Lower_Canopy", "Upper_Canopy", "Emergent")
)

# Scale proportions so they sum to 1 within each plot × layer
sp_data_corrected <- sp_data_unique %>%
  group_by(Coppicing_Plot, Stratum) %>%
  mutate(Sum_Proportion = sum(Proportion_of_Species, na.rm = TRUE)) %>%
  mutate(Proportion_of_Species = ifelse(Sum_Proportion > 1, Proportion_of_Species / Sum_Proportion, Proportion_of_Species)) %>%
  ungroup()

# Calculate Shannon Index per ET_ID
shannon_per_point <- sp_data_corrected %>%
  group_by(Coppicing_Plot, Stratum, ET_ID) %>%
  mutate(
    p = Proportion_of_Species / sum(Proportion_of_Species, na.rm = TRUE),
    p_log_p = -p * log(p)
  ) %>%
  summarise(Shannon_ET = sum(p_log_p, na.rm = TRUE), .groups = "drop")

# Average Shannon Index across ET_IDs per plot × stratum
shannon_index_normalised <- shannon_per_point %>%
  group_by(Coppicing_Plot, Stratum) %>%
  summarise(Shannon_Index_strata = mean(Shannon_ET, na.rm = TRUE), .groups = "drop")

# Merge back with species proportions (1 row per species per stratum)
model_data_species_wide <- sp_data_corrected %>%
  left_join(shannon_index_normalised, by = c("Coppicing_Plot", "Stratum")) %>%
  dplyr::select(Coppicing_Plot, Stratum, Species, Proportion_of_Species, Shannon_Index_strata)

# Add dormice data
dorm <- read.csv("Data/dormice_data_Panels.csv")

total_dormice_sum <- dorm %>%
  group_by(Coppicing_Plot) %>%
  summarise(Total_dormice_sum = sum(as.numeric(Total_dormice), na.rm = TRUE))

# Merge dormice total and mean
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

model_data_species_wide <- model_data_species_wide %>%
  left_join(plot_summary, by = "Coppicing_Plot")

#------------------#
##### Emergent #####
#------------------#

# Filter for Emergent stratum
emergent <- model_data_species_wide %>%
  filter(Stratum == "Emergent")

mod<-glm(Mean_Dormice ~ Shannon_Index_strata + num_unique_records, data=emergent, family="poisson")
summary(mod)

# Call:
#   glm(formula = Mean_Dormice ~ Shannon_Index_strata + num_unique_records, 
#       family = "poisson", data = emergent)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)           0.43524    0.26723   1.629 0.103382    
# Shannon_Index_strata -0.84600    0.38839  -2.178 0.029391 *  
#   num_unique_records    0.04395    0.01144   3.843 0.000122 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 64.157  on 35  degrees of freedom
# Residual deviance: 46.410  on 33  degrees of freedom
# AIC: 134.86
# 
# Number of Fisher Scoring iterations: 5

#---------------#
##### Upper #####
#---------------#

Upper <- model_data_species_wide %>%
  filter(Stratum == "Upper_Canopy")
mod<-glm(Mean_Dormice ~ Shannon_Index_strata+ num_unique_records, data=Upper, family="poisson")
summary(mod)

# Call:
#   glm(formula = Mean_Dormice ~ Shannon_Index_strata + num_unique_records, 
#       family = "poisson", data = Upper)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)           1.264939   0.287117   4.406 1.05e-05 ***
#   Shannon_Index_strata -1.220624   0.242258  -5.039 4.69e-07 ***
#   num_unique_records    0.003216   0.012269   0.262    0.793    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 85.66  on 53  degrees of freedom
# Residual deviance: 58.23  on 51  degrees of freedom
# AIC: 166.4
# 
# Number of Fisher Scoring iterations: 5


#---------------#
##### Lower #####
#---------------#

Lower <- model_data_species_wide %>%
  filter(Stratum == "Lower_Canopy")
#nb
nb<-glm.nb(Mean_Dormice ~ Shannon_Index_strata + num_unique_records, data = Lower)
summary(nb)

# glm.nb(formula = Mean_Dormice ~ Shannon_Index_strata + num_unique_records, 
#        data = Lower, init.theta = 6.386004066, link = log)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)           0.81414    0.30218   2.694  0.00706 **
#   Shannon_Index_strata  0.29207    0.33423   0.874  0.38220   
# num_unique_records   -0.03324    0.01549  -2.146  0.03190 * 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Negative Binomial(6.386) family taken to be 1)
# 
# Null deviance: 53.449  on 49  degrees of freedom
# Residual deviance: 48.830  on 47  degrees of freedom
# AIC: 171.92
# 
# Number of Fisher Scoring iterations: 1
# 
# 
# Theta:  6.39 
# Std. Err.:  5.34 
# 
# 2 x log-likelihood:  -163.918 

#---------------#
##### Shrub #####
#---------------#

Shrub <- model_data_species_wide %>%
  filter(Stratum == "Shrub")
#nb
nb<-glm.nb(Mean_Dormice ~ Shannon_Index_strata + num_unique_records, data = Shrub)
summary(nb)

# glm.nb(formula = Mean_Dormice ~ Shannon_Index_strata + num_unique_records, 
#        data = Shrub, init.theta = 3.035206923, link = log)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)           0.14305    0.51732   0.277    0.782    
# Shannon_Index_strata -1.57794    0.37636  -4.193 2.76e-05 ***
#   num_unique_records    0.07816    0.01967   3.973 7.11e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Negative Binomial(3.0352) family taken to be 1)
# 
# Null deviance: 51.378  on 35  degrees of freedom
# Residual deviance: 28.065  on 33  degrees of freedom
# AIC: 103.86
# 
# Number of Fisher Scoring iterations: 1
# 
# 
# Theta:  3.04 
# Std. Err.:  1.62 
# Warning while fitting theta: alternation limit reached 
# 
# 2 x log-likelihood:  -95.864 

#---------------------#
##### Understorey #####
#---------------------#

Understorey <- model_data_species_wide %>%
  filter(Stratum == "Understorey")
#nb
nb<-glm.nb(Mean_Dormice ~ Shannon_Index_strata + num_unique_records, data = Understorey)
summary(nb)
# Call:
#   glm.nb(formula = Mean_Dormice ~ Shannon_Index_strata + num_unique_records, 
#          data = Understorey, init.theta = 5.424206458, link = log)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)           1.472922   0.304789   4.833 1.35e-06 ***
#   Shannon_Index_strata -0.738131   0.162868  -4.532 5.84e-06 ***
#   num_unique_records    0.013758   0.007611   1.808   0.0707 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Negative Binomial(5.4242) family taken to be 1)
# 
# Null deviance: 210.56  on 176  degrees of freedom
# Residual deviance: 188.14  on 174  degrees of freedom
# AIC: 591.5
# 
# Number of Fisher Scoring iterations: 1
# 
# 
# Theta:  5.42 
# Std. Err.:  2.31 
# 
# 2 x log-likelihood:  -583.497 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#### Q3: DORMICE PREFERENCE ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Load data
model_data<-read.csv('Data/final_dataset.csv')
head(model_data)

# Count unique ET_ID per Coppicing_Plot
model_data_count <- model_data %>%
  group_by(Coppicing_Plot) %>%
  summarise(num_ET_points = n_distinct(ET_ID), .groups = "drop")

# Join back to sp_data_unique
model_data <- model_data_count %>%
  left_join(model_data, by = "Coppicing_Plot")

model_data$num_ET_points

# Summarise prop of each species
model_data_species <- model_data %>%
  group_by(Coppicing_Plot, Species) %>%
  summarise(
    Proportion_of_Species = sum(Proportion_of_Species, na.rm = TRUE),
    Mean_Dormice = first(Mean_Dormice),
    Coppicing_year = first(Coppicing_year),
    num_ET_points = first(num_ET_points),  # keep this column
    .groups = "drop"
  )

head(model_data_species)
model_data_species$num_ET_points

# Rescale so sum of prop of all species across each coppicing panel is never above 1
model_data_species_wide_scaled <- model_data_species %>%
  group_by(Coppicing_Plot) %>%
  mutate(
    total_prop = sum(Proportion_of_Species, na.rm = TRUE),
    Proportion_of_Species = Proportion_of_Species / total_prop
  ) %>%
  dplyr::select(-total_prop) %>%
  ungroup()

model_data<-model_data_species_wide_scaled

predictors= c("Species")
response_var=c("Mean_Dormice")

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

# With covariate: 
mod_nb <- glm.nb(Mean_Dormice ~ Species + num_ET_points, data = model_data_filtered)
summary(mod_nb)
# Call:
#   glm.nb(formula = Mean_Dormice ~ Species + num_ET_points, data = model_data_filtered, 
#          init.theta = 0.5886994045, link = log)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)         0.33117    0.51880   0.638    0.523
# SpeciesHazel        0.18377    0.34531   0.532    0.595
# SpeciesHoneysuckle -0.27776    0.43941  -0.632    0.527
# SpeciesOak         -0.05533    0.43874  -0.126    0.900
# SpeciesSycamore     0.39694    0.46709   0.850    0.395
# num_ET_points      -0.48422    0.42255  -1.146    0.252
# 
# (Dispersion parameter for Negative Binomial(0.5887) family taken to be 1)
# 
# Null deviance: 142.67  on 159  degrees of freedom
# Residual deviance: 138.87  on 154  degrees of freedom
# AIC: 416.69
# 
# Number of Fisher Scoring iterations: 1
# 
# 
# Theta:  0.589 
# Std. Err.:  0.151 
# 
# 2 x log-likelihood:  -402.692 

#//////////////////////////#

##### Q4: Canopy structure Shannon #####

#//////////////////////////#

#-------#
# Model #
#-------#

# Load data
sp_data<-read.csv('final_dataset.csv')

# Count unique ET_ID per Coppicing_Plot
sp_data_count <- sp_data %>%
  group_by(Coppicing_Plot) %>%
  summarise(num_ET_points = n_distinct(ET_ID), .groups = "drop")

# Join back to sp_data_unique
sp_data <- sp_data_count %>%
  left_join(sp_data, by = "Coppicing_Plot")

# Select only the relevant columns for uniqueness check
unique_rows <- !duplicated(sp_data[, c("Coppicing_Plot")])
sp_data_unique <- sp_data[unique_rows, ]

# Variables of interest (including response)
predictors <- c("Mean_Canopy_Height_2020", "Perc_under_2m_2020", "Can_cover_2020", "Height_cv_2020", 
                "Max_Canopy_Height_2020", "Gap_Proportion_2020", "LAI", "Intensity", "VCI", 
                "Canopy_Volume_2020", "FT_5_10m", "FT_1.37m", "FT_10_20m", "num_ET_points")

response_var<-c("Mean_Shannon_Index")

best_model <- select_best_model(sp_data_unique, predictors,response_var, cor_threshold = 0.7, shape="linear")

# Call:
#   lm(formula = formula, data = data)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.72752 -0.12631 -0.00567  0.15192  0.48955 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             2.802e+00  9.962e-01   2.812   0.0071 ** 
#   Max_Canopy_Height_2020  9.527e-03  1.382e-02   0.690   0.4938    
# LAI                    -1.935e-01  7.341e-02  -2.637   0.0112 *  
#   Intensity              -8.732e-04  4.207e-04  -2.075   0.0433 *  
#   VCI                     8.607e-01  1.088e+00   0.791   0.4329    
# Canopy_Volume_2020     -7.250e-07  2.236e-06  -0.324   0.7471    
# FT_10_20m               1.559e-03  2.889e-03   0.540   0.5920    
# num_ET_points          -6.526e-01  1.234e-01  -5.287 3.01e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2712 on 48 degrees of freedom
# Multiple R-squared:  0.5394,	Adjusted R-squared:  0.4722 
# F-statistic: 8.029 on 7 and 48 DF,  p-value: 1.962e-06

#//////////////////////////#

##### Q4:Canopy structure DM #####

#//////////////////////////#

#-------#
# Model #
#-------#

# Load data
sp_data<-read.csv('final_dataset.csv')

# Count unique combinations per Coppicing_Plot
plot_summary <- sp_data %>%
  distinct(Coppicing_Plot, Dor_Box, Year_DM_Recorded, No_of_Dormice) %>%
  group_by(Coppicing_Plot) %>%
  summarise(
    num_unique_records = n()
  )

# Join back to unique plot data
sp_data <- sp_data %>%
  distinct(Coppicing_Plot, .keep_all = TRUE) %>%
  filter(Coppicing_Plot %in% plot_summary$Coppicing_Plot) %>%
  left_join(plot_summary, by = "Coppicing_Plot")


# Select only the relevant columns for uniqueness check
unique_rows <- !duplicated(sp_data[, c("Coppicing_Plot")])
sp_data_unique <- sp_data[unique_rows, ]
colnames(sp_data_unique)

# dispersion >1.5: use negative binomial model instead:
gam_nb<-glm.nb(Mean_Dormice ~  Max_Canopy_Height_2020+ FT_5_10m+ LAI+ Intensity+ VCI+ Canopy_Volume_2020 +num_unique_records, data = sp_data_unique)
summary(gam_nb)

# Call:
#   glm.nb(formula = Mean_Dormice ~ Max_Canopy_Height_2020 + FT_5_10m + 
#            LAI + Intensity + VCI + Canopy_Volume_2020 + num_unique_records, 
#          data = sp_data_unique, init.theta = 1.794482978, link = log)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)            -1.611e+01  7.140e+00  -2.256 0.024039 *  
#   Max_Canopy_Height_2020  1.460e-02  8.966e-02   0.163 0.870681    
# FT_5_10m                2.941e-02  2.469e-02   1.191 0.233605    
# LAI                     9.498e-01  3.602e-01   2.637 0.008373 ** 
#   Intensity               3.813e-03  2.044e-03   1.866 0.062056 .  
# VCI                     9.970e+00  6.299e+00   1.583 0.113490    
# Canopy_Volume_2020     -1.769e-05  1.304e-05  -1.356 0.175037    
# num_unique_records      4.435e-02  1.155e-02   3.841 0.000122 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Negative Binomial(1.7945) family taken to be 1)
# 
# Null deviance: 82.681  on 55  degrees of freedom
# Residual deviance: 45.169  on 48  degrees of freedom
# AIC: 142.95
# 
# Number of Fisher Scoring iterations: 1
# 
# 
# Theta:  1.79 
# Std. Err.:  1.06 
# 
# 2 x log-likelihood:  -124.946 

# Load data
sp_data<-read.csv('final_dataset.csv')
colnames(sp_data)

# Count unique ET_ID per Coppicing_Plot
sp_data_count <- sp_data %>%
  group_by(Coppicing_Plot) %>%
  summarise(num_ET_points = n_distinct(ET_ID), .groups = "drop")

# Join back to sp_data_unique
sp_data <- sp_data_count %>%
  left_join(sp_data, by = "Coppicing_Plot")

# Select only the relevant columns for uniqueness check
unique_rows <- !duplicated(sp_data[, c("Coppicing_Plot")])
sp_data_unique <- sp_data[unique_rows, ]
colnames(sp_data_unique)
hist(sp_data_unique$Mean_Shannon_Index)

best_model_topo<- glm(Mean_Shannon_Index~HALP+ Aspect_Cos+Slope+ FlowDir+ Plane_Curve+  Profile_Curve + num_ET_points , family = Gamma(link = "log"),sp_data_unique)
summary(best_model_topo)

# Call:
#   glm(formula = Mean_Shannon_Index ~ HALP + Aspect_Cos + Slope + 
#         FlowDir + Plane_Curve + Profile_Curve + num_ET_points, family = Gamma(link = "log"), 
#       data = sp_data_unique)
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    1.544230   0.205001   7.533 1.13e-09 ***
#   HALP          -0.006578   0.002214  -2.971  0.00462 ** 
#   Aspect_Cos    -0.039222   0.028533  -1.375  0.17564    
# Slope         -0.010737   0.005515  -1.947  0.05745 .  
# FlowDir       -0.001683   0.004545  -0.370  0.71276    
# Plane_Curve    4.283007   4.644367   0.922  0.36104    
# Profile_Curve  8.259708   9.559564   0.864  0.39187    
# num_ET_points -0.416442   0.063395  -6.569 3.37e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Gamma family taken to be 0.01923946)
# 
# Null deviance: 2.3329  on 55  degrees of freedom
# Residual deviance: 1.1026  on 48  degrees of freedom
# AIC: 33.172

#////////////////////#

##### Q4: Topography DM #####

#////////////////////#

# Load data
sp_data<-read.csv('final_dataset.csv')

# Count unique combinations per Coppicing_Plot
plot_summary <- sp_data %>%
  distinct(Coppicing_Plot, Dor_Box, Year_DM_Recorded, No_of_Dormice) %>%
  group_by(Coppicing_Plot) %>%
  summarise(
    num_unique_records = n()
  )

# Join back to unique plot data
sp_data <- sp_data %>%
  distinct(Coppicing_Plot, .keep_all = TRUE) %>%
  filter(Coppicing_Plot %in% plot_summary$Coppicing_Plot) %>%
  left_join(plot_summary, by = "Coppicing_Plot")

# Select only the relevant columns for uniqueness check
unique_rows <- !duplicated(sp_data[, c("Coppicing_Plot")])
sp_data_unique <- sp_data[unique_rows, ]

# dispersion >1.5: use neg bin instead
gam_nb<-glm.nb(Mean_Dormice ~  HALP+ Aspect_Cos+ FlowDir+ Profile_Curve+ Mean_Curve+ Solar_Radiation+num_unique_records, data = sp_data_unique)
summary(gam_nb)

# Call:
#   glm.nb(formula = Mean_Dormice ~ HALP + Aspect_Cos + FlowDir + 
#            Profile_Curve + Mean_Curve + Solar_Radiation + num_unique_records, 
#          data = sp_data_unique, init.theta = 1.351654205, link = log)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)        -8.751e+00  4.251e+00  -2.059   0.0395 *  
#   HALP                7.641e-04  2.369e-02   0.032   0.9743    
# Aspect_Cos         -2.027e-01  3.174e-01  -0.639   0.5231    
# FlowDir             6.571e-02  6.268e-02   1.048   0.2945    
# Profile_Curve      -8.014e+01  1.251e+02  -0.641   0.5217    
# Mean_Curve          1.783e+02  1.155e+02   1.544   0.1226    
# Solar_Radiation     8.071e+00  4.556e+00   1.772   0.0764 .  
# num_unique_records  5.385e-02  1.212e-02   4.443 8.88e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Negative Binomial(1.3517) family taken to be 1)
# 
# Null deviance: 75.525  on 55  degrees of freedom
# Residual deviance: 42.662  on 48  degrees of freedom
# AIC: 144.57
# 
# Number of Fisher Scoring iterations: 1
# 
# 
# Theta:  1.352 
# Std. Err.:  0.652 
# Warning while fitting theta: alternation limit reached 
# 
# 2 x log-likelihood:  -126.570 

# Count unique records per year
year_summary <- sp_data_clean %>%
  distinct(Coppicing_Plot, Dor_Box, Year_DM_Recorded, No_of_Dormice) %>%
  group_by(Year_DM_Recorded) %>%
  summarise(
    num_unique_records = n(),
    .groups = "drop"
  )

year_summary
