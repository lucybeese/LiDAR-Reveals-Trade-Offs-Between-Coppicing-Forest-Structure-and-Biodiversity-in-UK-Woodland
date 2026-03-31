##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

### Bontuchel Biodiversity Analysis ###

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

### Hypothesis testing Q2, Q3 and Q4 ###

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

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
library(lmtest)
library(car)
library(sandwich)
library(AICcmodavg)

#/////////////////#
#### LOAD DATA ####
#/////////////////#

final_data<-read.csv('Data/final_dataset.csv')
colnames(final_data)

#/////////////////////////#
#### DE-DUPLICATE DATA ####
#/////////////////////////#

# Ensure numeric conversion
final_data <- final_data %>%
  mutate(
    No_of_Dormice = as.numeric(No_of_Dormice),
    Shannon_Index = as.numeric(Shannon_Index)
  )

# Get recovery period
final_data <- final_data %>%
  mutate(recovery_period = 2020 - Coppicing_year)

# Count unique ET_ID per Coppicing_Plot
plot_data_count <- final_data %>%
  group_by(Coppicing_Plot) %>%
  summarise(num_ET_points = n_distinct(ET_ID), .groups = "drop")

# Join back
final_data <- plot_data_count %>%
  left_join(final_data, by = "Coppicing_Plot")

# Count unique records per Coppicing_Plot
plot_dormouse_count <- final_data %>%
  distinct(Coppicing_Plot, Dor_Box, Year_DM_Recorded, No_of_Dormice) %>%
  group_by(Coppicing_Plot) %>%
  summarise(num_unique_records = n(), .groups = "drop")

# Join back
final_data <- plot_dormouse_count %>%
  left_join(final_data, by = "Coppicing_Plot")

#------#
# Save #
#------#

write.csv(final_data,"Data/data_allcols.csv",row.names = F)

# Remove duplicate entries based on dormouse observation structure
dorm_data <- final_data %>%
  distinct(Coppicing_Plot, Year_DM_Recorded, Dor_Box, ET_ID, .keep_all = TRUE)

vars_to_check <- c(
  "Mean_Canopy_Height_2020",
  "Perc_under_2m_2020",
  "Can_cover_2020",
  "Height_cv_2020",
  "Max_Canopy_Height_2020",
  "Gap_Proportion_2020",
  "Intensity",
  "VCI",
  "LAI",
  "Solar_Radiation",
  "FT_5_10m",
  "FT_1.37m",
  "FT_10_20m",
  "HALP",
  "Aspect_Cos",
  "TPI",
  "Roughness",
  "TRI",
  "Elevation",
  "Slope",
  "Mean_Shannon_Index",
  "Yr_since_Coppiced",
  "Mean_Curve",
  "FlowDir",
  "Plane_Curve",
  "Profile_Curve",
  "Canopy_Volume_2020",
  "Mean_Dormice",
  "num_unique_records"
)

#-------------------------------------------------#
# Check all these variables are the same per plot #
#-------------------------------------------------#

check_variation <- final_data %>%
  group_by(Coppicing_Plot) %>%
  summarise(across(all_of(vars_to_check), ~ n_distinct(.x, na.rm = TRUE))) %>%
  ungroup()

check_variation %>%
  summarise(across(-Coppicing_Plot, max)) %>%
  pivot_longer(everything(),
               names_to = "variable",
               values_to = "max_unique_values") %>%
  filter(max_unique_values > 1)

# One row per plot for modelling Shannon Index
plot_data <- final_data %>%
  distinct(Coppicing_Plot, .keep_all = TRUE)

#------#
# Save #
#------#

write.csv(plot_data,"Data/plot_data.csv",row.names = F)
write.csv(final_data,"Data/final_data_dedup.csv",row.names = F)

#------------------------#
# Look at DM Distrbution #
#------------------------#

hist(plot_data$Mean_Dormice, breaks = 20)
mean(plot_data$Mean_Dormice == 0)
unique_vals <- unique(plot_data$Mean_Dormice)
head(unique_vals)

# Select only the relevant columns for uniqueness check - prevents pseudoreplication 
unique_rows <- !duplicated(final_data[, c("ET_ID")])

# Keep only unique rows in the dataframe
sp_data_unique <- final_data[unique_rows, ]
colnames(sp_data_unique)

#/////////////////////#
#### SHANNON INDEX ####
#/////////////////////#

#------------------------------#
##### Q4: Canopy structure #####
#------------------------------#

#-------#
# Model #
#-------#

colnames(plot_data)
Q4_CS <- lm(Mean_Shannon_Index~Max_Canopy_Height_2020+LAI+Intensity+VCI+Canopy_Volume_2020+FT_10_20m+num_ET_points, data=plot_data)
summary(Q4_CS)

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

#-------------------#
# Check Assumptions #
#-------------------#

#plot(Q4_CS)

# Residuals
res <- residuals(Q4_CS)
fit <- fitted(Q4_CS)

# Normality
shapiro.test(res) # p=0.139>0.05 there is no statistical evidence to reject normality

# Homoscedasticity
bptest(Q4_CS) # p=0.038<0.05 there is statistical evidence of heteroscedascity, use heteroscedasticity-robust standard errors
robust_q4 <- coeftest(Q4_CS, vcov = vcovHC(Q4_CS, type = "HC3"))
print(robust_q4)

# Multicollinearity
vif(Q4_CS) # All values are below 3

# Influential points
cooks.distance(Q4_CS) # obs. 38

# Refit without 38
Q4_CS_no38 <- lm(
  Mean_Shannon_Index ~ Max_Canopy_Height_2020 + LAI + Intensity +
    VCI + Canopy_Volume_2020 + FT_10_20m + num_ET_points,
  data = plot_data[-38, ]
)

summary(Q4_CS_no38)

AIC(Q4_CS, Q4_CS_no38)

models <- list(Q4_CS, Q4_CS_no38)
model_names <- c("Full", "No38")

aictab(cand.set = models, modnames = model_names)


summary(Q4_CS)$adj.r.squared
summary(Q4_CS_no38)$adj.r.squared
# 38 did have influence, but didn't change direction of relationship

#------#
# Plot #
#------#

label_map <- c(
  "Max_Canopy_Height_2020" = "Max Canopy Height",
  "LAI" = "LAI",
  "Intensity" = "Intensity",
  "VCI" = "VCI",
  "Canopy_Volume_2020" = "Canopy Volume",
  "FT_10_20m" = "First returns 10–20 m",
  "num_ET_points" = "Number of ET points",
  "num_unique_records" = "No. of DM records"
)

# Scale predictors so they are comparable
# Identify numeric columns
numeric_vars <- sapply(sp_data_unique, is.numeric)  
sp_data_unique_scaled <- sp_data_unique  # Create a copy of the dataset
sp_data_unique_scaled[numeric_vars] <- scale(sp_data_unique[numeric_vars])  # Scale only numeric columns

Q4_CS_scaled <- lm(Mean_Shannon_Index ~ Max_Canopy_Height_2020 + LAI + Intensity + VCI + Canopy_Volume_2020 + FT_10_20m + num_ET_points, 
                   data = sp_data_unique_scaled)

# Extract coefficients and confidence intervals
coefs <- summary(Q4_CS_scaled)$coefficients
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

# Order predictors by absolute coefficient size 
ord <- order(abs(estimates), decreasing = FALSE)

estimates <- estimates[ord]
lower <- lower[ord]
upper <- upper[ord]
coef_names <- coef_names[ord]
new_names <- label_map[coef_names]

#%%%%%%%%%%%%%%%%%%%%%%#
###### Plotting  #######
#%%%%%%%%%%%%%%%%%%%%%%#

# Open PDF device
cairo_pdf("Figures/CoeffsStruc.pdf",
          width = 7,
          height = 5,
          family = "Cambria",
          bg = "white")

par(mfrow = c(1,1), mar = c(5, 16, 4, 2), las = 1, family = "Cambria")
cex_text <- 1.3

plot(estimates, seq_along(estimates), pch = 16, xlim = c(-1, 1),
     xlab = "Estimate", ylab = "", axes = FALSE, cex = cex_text, cex.lab = cex_text)

axis(1, cex.axis = cex_text)

axis(2, at = seq_along(coef_names), labels = FALSE, las = 2)
text(x = rep(par("usr")[1]-0.2, length(new_names)), y = seq_along(new_names), 
     labels = new_names, adj = 1, xpd = TRUE, cex = cex_text)

arrows(lower, seq_along(estimates), upper, seq_along(estimates), length = 0.1, angle = 90, code = 3)
abline(v = 0, lty = 2, col = "hotpink")  # Reference line

# Add "a)" in the top-left corner (inside the plot margin)
# Use xpd=NA to allow drawing in the figure margin
text(x = par("usr")[1], y = par("usr")[4] + 0.5, labels = "a)", adj = c(0, 1), cex = cex_text, xpd = NA)

dev.off()

#------------------------#
##### Q4: Topography #####
#------------------------#

#-------#
# Model #
#-------#

# lm
Q4_Topo_lm <- lm(
  Mean_Shannon_Index ~ HALP + Aspect_Cos + Slope + FlowDir +
    Plane_Curve + Profile_Curve + num_ET_points,
  data = plot_data
)

#-------------------#
# Check Assumptions #
#-------------------#

# plot(Q4_Topo_lm)

# Residuals
res <- residuals(Q4_Topo_lm)
fit <- fitted(Q4_Topo_lm)

# Normality
shapiro.test(res) # p=0.004095 reject normality

hist(plot_data$Mean_Shannon_Index) 

#-------#
# Model #
#-------#

# Gamma GLM

Q4_Topo <- glm(
  Mean_Shannon_Index~HALP+ Aspect_Cos+Slope+ FlowDir+ Plane_Curve+  Profile_Curve + num_ET_points,
  family = Gamma(link = "log"),
  data = plot_data
)

summary(Q4_Topo)

# Call:
#   glm(formula = Mean_Shannon_Index ~ HALP + Aspect_Cos + Slope + 
#         FlowDir + Plane_Curve + Profile_Curve + num_ET_points, family = Gamma(link = "log"), 
#       data = plot_data)
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    1.546108   0.204950   7.544 1.08e-09 ***
#   HALP          -0.006592   0.002213  -2.979  0.00453 ** 
#   Aspect_Cos    -0.039251   0.028533  -1.376  0.17532    
# Slope         -0.010750   0.005514  -1.950  0.05705 .  
# FlowDir       -0.001738   0.004549  -0.382  0.70416    
# Plane_Curve    4.272204   4.644901   0.920  0.36230    
# Profile_Curve  8.277633   9.558893   0.866  0.39082    
# num_ET_points -0.416493   0.063395  -6.570 3.36e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Gamma family taken to be 0.0192383)
# 
# Null deviance: 2.3329  on 55  degrees of freedom
# Residual deviance: 1.1024  on 48  degrees of freedom
# AIC: 33.163
# 
# Number of Fisher Scoring iterations: 5

#-------------------#
# Check assumptions #
#-------------------#

# Fitted vs resid
plot(fitted(Q4_Topo), residuals(Q4_Topo, type = "pearson"),
     xlab = "Fitted values", ylab = "Pearson residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red")

res_squared <- residuals(Q4_Topo, type = "pearson")^2
fit_vals <- fitted(Q4_Topo)
summary(lm(res_squared ~ I(fit_vals^2)))

# Dispersion
summary(Q4_Topo)$dispersion

# Outlier
plot(cooks.distance(Q4_Topo), type = "h")
abline(h = 4 / nrow(plot_data), col = "red") # 38 again

# Run again without 38

Q4_Topo_no38 <- glm(
  Mean_Shannon_Index ~ HALP + Aspect_Cos + Slope + FlowDir +
    Plane_Curve + Profile_Curve + num_ET_points,
  family = Gamma(link = "log"),
  data = plot_data[-38, ]
)

summary(Q4_Topo_no38) # 38 did have influence, but didn't change direction of relationship

#------#
# Plot #
#------#

label_map_topo <- c(
  "HALP" = "HALP",
  "Aspect_Cos" = "Aspect (cos)",
  "Slope" = "Slope",
  "FlowDir" = "Flow direction",
  "Plane_Curve" = "Plane curvature",
  "Profile_Curve" = "Profile curvature",
  "num_ET_points" = "Number of ET points",
  "num_unique_records" = "No. of DM records"
)

# Make a copy of your data
sp_data_unique_scaled <- sp_data_unique

predictors <- c("HALP","Aspect_Cos",
                "Slope", "FlowDir","Plane_Curve","Profile_Curve",
               "num_ET_points")

# Scale only predictor columns (not the response)
sp_data_unique_scaled[predictors] <- scale(sp_data_unique_scaled[predictors])

# Fit Gamma model using unscaled response and scaled predictors
model_scaled <- glm(
  Mean_Shannon_Index~HALP+ Aspect_Cos+Slope+ FlowDir+ Plane_Curve+  Profile_Curve + num_ET_points,
  family = Gamma(link = "log"),
  data = sp_data_unique_scaled
)

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

# Order predictors
ord <- order(abs(estimates), decreasing = FALSE)

estimates <- estimates[ord]
errors <- errors[ord]
lower <- lower[ord]
upper <- upper[ord]
coef_names <- coef_names[ord]
new_names <- new_names[ord]
new_names <- label_map_topo[coef_names]

#%%%%%%%%%%%%%%%%%%%%%%#
###### Plotting  #######
#%%%%%%%%%%%%%%%%%%%%%%#

# Open PDF device
cairo_pdf("Figures/CoeffsTopo.pdf",
          width = 7,
          height = 5,
          family = "Cambria",
          bg = "white")

par(mfrow = c(1,1), mar = c(5, 16, 4, 2), las = 1, family = "Cambria")
cex_text <- 1.3

plot(estimates, seq_along(estimates), pch = 16,  xlim = c(- 0.2, 0.2),
     xlab = "Estimate", ylab = "", axes = FALSE, cex = cex_text, cex.lab = cex_text)
axis(1,cex.axis = cex_text)

# Remove the axis(2) call and replace it with text for custom labels
axis(2, at = seq_along(coef_names), labels = FALSE, las = 2, cex=cex_text)
text(x = rep(par("usr")[1]-0.03, length(new_names)), y = seq_along(new_names), 
     labels = new_names, adj = 1, xpd = TRUE, cex= cex_text)

arrows(lower, seq_along(estimates), upper, seq_along(estimates), length = 0.1, angle = 90, code = 3)
abline(v = 0, lty = 2, col = "hotpink")  # Reference line

# Add "c)" in the top-left corner (inside the plot margin)
# Use xpd=NA to allow drawing in the figure margin
text(x = par("usr")[1], y = par("usr")[4] + 0.5, labels = "c)", adj = c(0, 1), cex = cex_text, xpd = NA)

dev.off()

#------------------------#
##### Q2: Management #####
#------------------------#

#-------#
# Model #
#-------#

# Add the number of samples per plot as a covariate
# Center due to polnomial
plot_data$recovery_c <- plot_data$recovery_period - mean(plot_data$recovery_period, na.rm = TRUE)

model_manag <- gam(
  Mean_Shannon_Index ~ recovery_c + I(recovery_c^2) + num_ET_points,
  family = Gamma(link = "log"),
  data = plot_data
)
summary(model_manag)

# Family: Gamma 
# Link function: log 
# 
# Formula:
#   Mean_Shannon_Index ~ recovery_c + I(recovery_c^2) + num_ET_points
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      1.161e+00  7.764e-02  14.954  < 2e-16 ***
#   recovery_c       5.982e-03  3.924e-03   1.524    0.133    
# I(recovery_c^2)  4.759e-05  8.388e-04   0.057    0.955    
# num_ET_points   -4.253e-01  6.308e-02  -6.741 1.28e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# 
# R-sq.(adj) =  0.355   Deviance explained = 39.1%
# GCV = 0.029446  Scale est. = 0.020976  n = 56

#-------------------#
# Check assumptions #
#-------------------#

# Fitted vs resid
plot(fitted(model_manag), residuals(model_manag, type = "pearson"),
     xlab = "Fitted values", ylab = "Pearson residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red")

res_squared <- residuals(model_manag, type = "pearson")^2
fit_vals <- fitted(model_manag)
summary(lm(res_squared ~ I(fit_vals^2)))

# Dispersion
summary(model_manag)$dispersion

# Outlier
plot(cooks.distance(model_manag), type = "h")
abline(h = 4 / nrow(plot_data), col = "red") # 38 again

# Run again without 38
model_manag_no38 <- gam(
  Mean_Shannon_Index ~ recovery_period + I(recovery_period^2)  + num_ET_points,
  data = plot_data[-38, ],
  family = Gamma(link = "log")
)
summary(model_manag_no38) # 38 did have influence, but didn't change direction of relationship

#------#
# Plot #
#------#

# Sequence of recovery years
recovery_seq <- seq(min(plot_data$recovery_period, na.rm = TRUE),
                    max(plot_data$recovery_period, na.rm = TRUE),
                    length.out = 200)

# Convert to centred variable used in the model
recovery_c_seq <- recovery_seq - mean(plot_data$recovery_period, na.rm = TRUE)

new_data <- data.frame(
  recovery_c = recovery_c_seq,
  num_ET_points = mean(plot_data$num_ET_points, na.rm = TRUE)
)

# Add quadratic term
new_data$`I(recovery_c^2)` <- recovery_c_seq^2

# Predict on link scale
pred <- predict(model_manag, newdata = new_data, type = "link", se.fit = TRUE)

# Back-transform to response scale
new_data$fit <- exp(pred$fit)
new_data$lower <- exp(pred$fit - 1.96 * pred$se.fit)
new_data$upper <- exp(pred$fit + 1.96 * pred$se.fit)

# Add original recovery period back for plotting
new_data$recovery_period <- recovery_seq

#%%%%%%%%%%%%%%%%%%%%%%#
###### Plotting  #######
#%%%%%%%%%%%%%%%%%%%%%%#

# Open PDF device
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
    color = "#D66F90",
    fill = "pink",
    alpha = 0.5
  ) +
  geom_point(color = "#7A0016", size = 2, alpha = 0.6) +
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

#///////////////#
#### DORMICE ####
#///////////////#

#------------------------------#
##### Q4: Canopy structure #####
#------------------------------#

#-------#
# Model #
#-------#

# Poisson GLM
pois_DM <- glm(
  Mean_Dormice ~ Max_Canopy_Height_2020 + FT_10_20m + LAI +
    Intensity + VCI + Canopy_Volume_2020 + num_unique_records,
  family = poisson(link = "log"),
  data = plot_data
)

summary(pois_DM)

# Call:
#   glm(formula = Mean_Dormice ~ Max_Canopy_Height_2020 + FT_10_20m + 
#         LAI + Intensity + VCI + Canopy_Volume_2020 + num_unique_records, 
#       family = poisson(link = "log"), data = plot_data)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)            -1.147e+01  5.205e+00  -2.204   0.0275 *  
#   Max_Canopy_Height_2020 -1.394e-03  6.028e-02  -0.023   0.9815    
# FT_10_20m              -7.641e-03  1.511e-02  -0.506   0.6130    
# LAI                     5.741e-01  2.759e-01   2.081   0.0375 *  
#   Intensity               2.466e-03  1.359e-03   1.814   0.0697 .  
# VCI                     8.460e+00  5.400e+00   1.567   0.1172    
# Canopy_Volume_2020     -8.793e-06  9.298e-06  -0.946   0.3443    
# num_unique_records      3.999e-02  8.270e-03   4.836 1.33e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 111.99  on 55  degrees of freedom
# Residual deviance:  59.70  on 48  degrees of freedom
# AIC: Inf
# 
# Number of Fisher Scoring iterations: 6

#-------------------#
# Check assumptions #
#-------------------#

# Check dispersion
disp <- deviance(pois_DM) / df.residual(pois_DM)
disp

# Zero inflation
simres <- simulateResiduals(pois_DM)
plot(simres)
testZeroInflation(simres)

#VIF
car::vif(pois_DM)

# Outliers
plot(cooks.distance(pois_DM), type = "h")
abline(h = 4 / nrow(plot_data), col = 2, lty = 2)
testOutliers(simres)
testUniformity(simres)

#------#
# Plot #
#------#

predictors_DM <- c(
  "Max_Canopy_Height_2020",
  "FT_10_20m",
  "LAI",
  "Intensity",
  "VCI",
  "Canopy_Volume_2020",
  "num_unique_records"
)

# Scale predictors so they are comparable
sp_data_unique_scaled <- sp_data_unique
sp_data_unique_scaled[predictors_DM] <- scale(sp_data_unique_scaled[predictors_DM])

model_scaled <- glm(
  Mean_Dormice ~ Max_Canopy_Height_2020 + FT_10_20m + LAI +
    Intensity + VCI + Canopy_Volume_2020 + num_unique_records,
  family = poisson(link = "log"), data = sp_data_unique_scaled)

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

# Order predictors
ord <- order(abs(estimates), decreasing = FALSE)

estimates <- estimates[ord]
errors <- errors[ord]
lower <- lower[ord]
upper <- upper[ord]
coef_names <- coef_names[ord]
new_names <- new_names[ord]
new_names <- label_map[coef_names]

#%%%%%%%%%%%%%%%%%%%%%%#
###### Plotting  #######
#%%%%%%%%%%%%%%%%%%%%%%#

# Open PDF device
cairo_pdf("Figures/CoeffsStrucDM.pdf",
          width = 7,
          height = 5,
          family = "Cambria",
          bg = "white")

par(mfrow = c(1,1), mar = c(5, 14.5, 4, 2), las = 1, family = "Cambria")
cex_text<-1.3

plot(estimates, seq_along(estimates), pch = 16,  xlim = c(- 1, 1.5),
     xlab = "Estimate", ylab = "", axes = FALSE, cex = cex_text, cex.lab = cex_text)
axis(1,cex.axis = cex_text)

axis(2, at = seq_along(coef_names), labels = FALSE, las = 2)
text(x = rep(par("usr")[1]-0.2, length(new_names)), y = seq_along(new_names), 
     labels = new_names, adj = 1, xpd = TRUE, cex = cex_text)

arrows(lower, seq_along(estimates), upper, seq_along(estimates), length = 0.1, angle = 90, code = 3)
abline(v = 0, lty = 2, col = "hotpink")  # Reference line

text(x = par("usr")[1], y = par("usr")[4] + 0.5, labels = "b)", adj = c(0, 1), cex = cex_text, xpd = NA)

dev.off()

#------------------------#
##### Q4: Topography #####
#------------------------#

#-------#
# Model #
#-------#

# Poisson GLM
glm_pois <- glm(
  Mean_Dormice ~ HALP + Aspect_Cos + Slope + FlowDir +
    Plane_Curve + Profile_Curve + num_unique_records,
  family = poisson(link = "log"),
  data = sp_data_unique
)

summary(glm_pois)

# # Call:
# glm(formula = Mean_Dormice ~ HALP + Aspect_Cos + Slope + FlowDir + 
#       Plane_Curve + Profile_Curve + num_unique_records, family = poisson(link = "log"), 
#     data = sp_data_unique)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)         -0.936132   1.558074  -0.601    0.548    
# HALP                -0.001981   0.016876  -0.117    0.907    
# Aspect_Cos           0.007767   0.218589   0.036    0.972    
# Slope                0.016573   0.040499   0.409    0.682    
# FlowDir             -0.010969   0.043619  -0.251    0.801    
# Plane_Curve         31.198269  37.074280   0.842    0.400    
# Profile_Curve      -64.630890  84.261647  -0.767    0.443    
# num_unique_records   0.044495   0.007265   6.124  9.1e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 122.613  on 61  degrees of freedom
# Residual deviance:  69.572  on 54  degrees of freedom
# AIC: Inf
# 
# Number of Fisher Scoring iterations: 6

#-------------------#
# Check assumptions #
#-------------------#

# Check dispersion
disp <- deviance(glm_pois) / df.residual(glm_pois)
disp

# Zero inflation
simres <- simulateResiduals(glm_pois)
plot(simres)
testZeroInflation(simres)

#VIF
car::vif(glm_pois)

# Outliers
plot(cooks.distance(glm_pois), type = "h")
abline(h = 4 / nrow(plot_data), col = 2, lty = 2)
testOutliers(simres)
testUniformity(simres)

#------#
# Plot #
#------#

# Make a copy of your data
sp_data_unique_scaled <- sp_data_unique

predictors <- c("HALP","Aspect_Cos",
                "Slope", "FlowDir","Plane_Curve","Profile_Curve",
                "num_unique_records")

sp_data_unique_scaled[, predictors] <- scale(sp_data_unique_scaled[, predictors])

model_scaled <- glm(
  Mean_Dormice ~ HALP + Aspect_Cos + FlowDir + Profile_Curve +
    Plane_Curve + Slope + num_unique_records,
  family = poisson(link = "log"),
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

coef_names <- coef_names[ord]
new_names <- new_names[ord]
new_names <- label_map_topo[coef_names]

#%%%%%%%%%%%%%%%%%%%%%%#
###### Plotting  #######
#%%%%%%%%%%%%%%%%%%%%%%#

# Open PDF device
cairo_pdf("Figures/CoeffsTopoDM.pdf",
          width = 7,
          height = 5,
          family = "Cambria",
          bg = "white")

par(mfrow = c(1,1), mar = c(5, 14.5, 4, 2), las = 1, family = "Cambria")
cex_text<-1.3

plot(estimates, seq_along(estimates), pch = 16,  xlim = c(- 1, 1.5),
     xlab = "Estimate", ylab = "", axes = FALSE, cex = cex_text, cex.lab = cex_text)
axis(1,cex.axis = cex_text)

axis(2, at = seq_along(coef_names), labels = FALSE, las = 2)
text(x = rep(par("usr")[1]-0.2, length(new_names)), y = seq_along(new_names), 
     labels = new_names, adj = 1, xpd = TRUE, cex = cex_text)

arrows(lower, seq_along(estimates), upper, seq_along(estimates), length = 0.1, angle = 90, code = 3)
abline(v = 0, lty = 2, col = "hotpink")  # Reference line

text(x = par("usr")[1], y = par("usr")[4] + 0.5, labels = "d)", adj = c(0, 1), cex = cex_text, xpd = NA)

dev.off()

#------------------------#
##### Q3: Management #####
#------------------------#

#-------#
# Model #
#-------#

#  Gaussian LM (log response)

# Used because Mean_Dormice is continuous and right-skewed
model_manage <- lm(
  log(Mean_Dormice + 1) ~ recovery_period + I(recovery_period^2)+ num_unique_records,
  data = plot_data
)

#-------------------#
# Check assumptions #
#-------------------#

summary(model_manage)
shapiro.test(residuals(model_manage)) # Residuals showed deviation from normality,

#-------#
# Model #
#-------#

#  Poisson GLM

# Poisson models are commonly used for count data
mod <- glm(
  Mean_Dormice ~ recovery_period + I(recovery_period^2)+ num_unique_records,
  family = poisson(link = "log"),
  data = plot_data
)
summary(mod)

#-------------------#
# Check assumptions #
#-------------------#

# Check for overdispersion
dispersion_ratio <- sum(residuals(mod, type = "pearson")^2) / mod$df.residual
dispersion_ratio
# Dispersion ratio > 2 indicates strong overdispersion,

#-------#
# Model #
#-------#

# Negative binomial

# The number of survey records per plot is included as a covariate.
nb_model_final <- glm.nb(
  Mean_Dormice ~ recovery_period + I(recovery_period^2) + num_unique_records,
  data = plot_data
)
summary(nb_model_final)

#-------------------#
# Check assumptions #
#-------------------#

# Check for overdispersion
dispersion <- deviance(nb_model_final) / df.residual(nb_model_final)
dispersion # not overdispersed

# Zero inflation check
simres <- simulateResiduals(fittedModel = nb_model_final, plot = TRUE)
testZeroInflation(simres) # no zero inflation

# VIF
car::vif(nb_model_final) # high because of polynomial

# Outliers
plot(cooks.distance(nb_model_final), type = "h")
abline(h = 4 / nrow(plot_data), col = 2, lty = 2)
testOutliers(simres)
testUniformity(simres)

plot_data$recovery_c <- plot_data$recovery_period - mean(plot_data$recovery_period, na.rm = TRUE)

#-------#
# Model #
#-------#
# centered to improve vif

nb_model_final <- glm.nb(
  Mean_Dormice ~ recovery_c + I(recovery_c^2) + num_unique_records,
  data = plot_data
)
summary(nb_model_final)

# Call:
#   glm.nb(formula = Mean_Dormice ~ recovery_c + I(recovery_c^2) + 
#            num_unique_records, data = plot_data, init.theta = 2.277152915, 
#          link = log)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)        -0.980526   0.328634  -2.984  0.00285 ** 
#   recovery_c          0.016153   0.036736   0.440  0.66016    
# I(recovery_c^2)    -0.013617   0.008685  -1.568  0.11689    
# num_unique_records  0.054431   0.009312   5.845 5.05e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Negative Binomial(2.2772) family taken to be 1)
# 
# Null deviance: 80.421  on 55  degrees of freedom
# Residual deviance: 45.906  on 52  degrees of freedom
# AIC: 137
# 
# Number of Fisher Scoring iterations: 1
# 
# 
# Theta:  2.28 
# Std. Err.:  1.32 
# Warning while fitting theta: alternation limit reached 
# 
# 2 x log-likelihood:  -126.999 

#%%%%%%%%%%%%%%%%%%%%%%#
###### Plotting  #######
#%%%%%%%%%%%%%%%%%%%%%%#

# Open PDF device
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
    color = "#D66F90",
    fill = "pink",
    alpha=0.5
  ) +
  geom_point(color = "#7A0016", size = 2, alpha = 0.6) +  # Transparent points
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
