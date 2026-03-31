##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

### Bontuchel Biodiversity Analysis ###

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

### Hypothesis testing Q2 and Q4 ###

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

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
library(TeachingDemos)

#/////////////////#
#### LOAD DATA ####
#/////////////////#

final_data<-read.csv("Data/data_allcols.csv")

#////////////////////#
#### Process data ####
#////////////////////#

#---------------------------------------------------#
# Remove duplicate species entries per plot × layer #
#---------------------------------------------------#

sp_data_unique <- final_data %>%
  distinct(Coppicing_Plot, Forest_layer, Species, .keep_all = TRUE)

#------------------------#
# Define vertical strata #
#------------------------#

sp_data_unique <- sp_data_unique %>%
  mutate(
    Stratum = cut(
      ET_layer_height,
      breaks = c(-Inf, 1.37, 5, 10, 20, Inf),
      labels = c("Understorey", "Shrub", "Lower_Canopy", "Upper_Canopy", "Emergent")
    )
  )

#--------------------------------------------------------#
# Ensure species proportions sum to 1 per plot × stratum #
#--------------------------------------------------------#

sp_data_corrected <- sp_data_unique %>%
  group_by(Coppicing_Plot, Stratum) %>%
  mutate(
    Sum_Proportion = sum(Proportion_of_Species, na.rm = TRUE),
    Proportion_of_Species = ifelse(
      Sum_Proportion > 1,
      Proportion_of_Species / Sum_Proportion,
      Proportion_of_Species
    )
  ) %>%
  ungroup()

summary(sp_data_corrected$Proportion_of_Species)

#------------------------------------------#
# Calculate Shannon diversity per ET point #
#------------------------------------------#

shannon_per_point <- sp_data_corrected %>%
  group_by(Coppicing_Plot, Stratum, ET_ID) %>%
  mutate(
    p = Proportion_of_Species / sum(Proportion_of_Species, na.rm = TRUE),
    p_log_p = -p * log(p)
  ) %>%
  summarise(
    Shannon_ET = sum(p_log_p, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    sp_data_corrected %>%
      dplyr::select(Coppicing_Plot, Stratum, ET_ID, num_unique_records) %>%
      distinct(),
    by = c("Coppicing_Plot","Stratum","ET_ID")
  )

#-----------------------------------------------------#
# Average Shannon across ET points per plot × stratum #
#-----------------------------------------------------#

shannon_index_normalised <- shannon_per_point %>%
  group_by(Coppicing_Plot, Stratum) %>%
  summarise(
    Shannon_Index_strata = mean(Shannon_ET, na.rm = TRUE),
    .groups = "drop"
  )

#----------------------------------------#
# Final modelling dataset (species-wide) #
#----------------------------------------#

model_data_species_wide <- sp_data_corrected %>%
  left_join(shannon_index_normalised, by = c("Coppicing_Plot", "Stratum")) %>%
  dplyr::select(
    Coppicing_Plot,
    Stratum,
    Species,
    Proportion_of_Species,
    Shannon_Index_strata,
    Mean_Dormice,
    num_unique_records,
    
  )

# Add Coppicing_Year to wide data
model_data_species_wide <- model_data_species_wide %>%
  left_join(
    final_data %>% dplyr::select(Coppicing_Plot, Coppicing_year) %>% distinct(),
    by = "Coppicing_Plot"
  )

# Save the dataframe
write.csv(model_data_species_wide,"Data/sp_data_shan.csv",row.names = F)

#------------------------------------------------------#
# Build modelling dataset (one row per plot × stratum) #
#------------------------------------------------------#

model_data_strata <- sp_data_corrected %>%
  left_join(shannon_index_normalised, by = c("Coppicing_Plot", "Stratum")) %>%
  group_by(Coppicing_Plot, Stratum) %>%
  summarise(
    Mean_Dormice = first(Mean_Dormice),
    Shannon_Index_strata = first(Shannon_Index_strata),
    num_unique_records = first(num_unique_records),
    Coppicing_year = first(Coppicing_year),
    .groups = "drop"
  )

head(model_data_strata)

write.csv(model_data_strata,"Data/model_data_strata.csv",row.names = F)

#/////////////////////#
#### Q2: DM Strata ####
#/////////////////////#

#//////////////////#
##### Emergent #####
#//////////////////#

#-------#
# Model #
#-------#

# Filter for Emergent stratum 
emergent <- model_data_strata %>%
  filter(Stratum == "Emergent")

mod_em<-glm(Mean_Dormice ~ Shannon_Index_strata + num_unique_records, data=emergent, family="poisson")
summary(mod_em)

#-------------------#
# Check assumptions #
#-------------------#

# Check dispersion
disp<- deviance(mod_em) / df.residual(mod_em)
disp # 1.5 Overdispersed

#-------#
# Model #
#-------#

mod_em_nb <- MASS::glm.nb(
  Mean_Dormice ~ Shannon_Index_strata + num_unique_records,
  data = emergent
)
summary(mod_em_nb) #use nb

# Call:
#   MASS::glm.nb(formula = Mean_Dormice ~ Shannon_Index_strata + 
#                  num_unique_records, data = emergent, init.theta = 4.063540198, 
#                link = log)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)          -0.05487    0.37698  -0.146 0.884278    
# Shannon_Index_strata -1.26018    0.63377  -1.988 0.046769 *  
#   num_unique_records    0.04135    0.01143   3.618 0.000297 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Negative Binomial(4.0635) family taken to be 1)
# 
# Null deviance: 43.880  on 26  degrees of freedom
# Residual deviance: 27.833  on 24  degrees of freedom
# AIC: 90.562
# 
# Number of Fisher Scoring iterations: 1
# 
# 
# Theta:  4.06 
# Std. Err.:  3.85 
# 
# 2 x log-likelihood:  -82.562 

#-------------------#
# Check assumptions #
#-------------------#

simres <- simulateResiduals(mod_em_nb)
plot(simres)

# Check dispersion
disp<- deviance(mod_em_nb) / df.residual(mod_em_nb)
disp
testZeroInflation(simres)
testUniformity(simres)
testOutliers(simres)

plot(simres, emergent$Shannon_Index_strata)
plot(simres, emergent$num_unique_records)

plot(cooks.distance(mod_em_nb), type = "h")
abline(h = 4 / nrow(emergent), col = 2, lty = 2)

car::vif(mod_em_nb)
cor(emergent[, c("Shannon_Index_strata", "num_unique_records")], use = "complete.obs")
# Assumptions met

#///////////////#
##### Upper #####
#///////////////#

#-------#
# Model #
#-------#

Upper <- model_data_strata %>%
  filter(Stratum == "Upper_Canopy")
mod_up<-glm(Mean_Dormice ~ Shannon_Index_strata+ num_unique_records, data=Upper, family="poisson")
summary(mod_up)

# Call:
#   glm(formula = Mean_Dormice ~ Shannon_Index_strata + num_unique_records, 
#       family = "poisson", data = Upper)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)          -0.189120   0.344997  -0.548  0.58357    
# Shannon_Index_strata -1.140566   0.394755  -2.889  0.00386 ** 
#   num_unique_records    0.034102   0.006928   4.922 8.56e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 89.241  on 44  degrees of freedom
# Residual deviance: 49.012  on 42  degrees of freedom
# AIC: Inf

#-------------------#
# Check assumptions #
#-------------------#

# Check dispersion
disp<- deviance(mod_up) / df.residual(mod_up)
disp # Not Overdispersed

simres <- simulateResiduals(mod_up)
plot(simres)
testZeroInflation(simres)
testUniformity(simres)
testOutliers(simres)

plot(simres, Upper$Shannon_Index_strata)
plot(simres, Upper$num_unique_records)

plot(cooks.distance(mod_up), type = "h")
abline(h = 4 / nrow(Upper), col = 2, lty = 2)

car::vif(mod_up)

# Assumptions met

#///////////////#
##### Lower #####
#///////////////#

#-------#
# Model #
#-------#

Lower <- model_data_strata %>%
  filter(Stratum == "Lower_Canopy")
mod_low<-glm(Mean_Dormice ~ Shannon_Index_strata+ num_unique_records, data=Lower, family="poisson")
summary(mod_low)

#-------------------#
# Check assumptions #
#-------------------#

# Check dispersion
disp<- deviance(mod_low) / df.residual(mod_low)
disp # Overdispersed

#-------#
# Model #
#-------#

#nb
mod_low_nb<-glm.nb(Mean_Dormice ~ Shannon_Index_strata + num_unique_records, data = Lower)
summary(mod_low_nb)

# Check dispersion

# Call:
#   glm.nb(formula = Mean_Dormice ~ Shannon_Index_strata + num_unique_records, 
#          data = Lower, init.theta = 1.05575622, link = log)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)          -1.54306    0.50426  -3.060  0.00221 ** 
#   Shannon_Index_strata  0.31857    0.54074   0.589  0.55577    
# num_unique_records    0.05872    0.01323   4.438 9.07e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Negative Binomial(1.0558) family taken to be 1)
# 
# Null deviance: 47.576  on 42  degrees of freedom
# Residual deviance: 32.889  on 40  degrees of freedom
# AIC: 102.67
# 
# Number of Fisher Scoring iterations: 1
# 
# 
# Theta:  1.056 
# Std. Err.:  0.512 
# Warning while fitting theta: alternation limit reached 
# 
# 2 x log-likelihood:  -94.672 

#-------------------#
# Check assumptions #
#-------------------#

simres <- simulateResiduals(mod_low_nb)
plot(simres)

# Check dispersion
disp<- deviance(mod_low_nb) / df.residual(mod_low_nb)
disp 
testZeroInflation(simres)
testUniformity(simres)
testOutliers(simres)

plot(cooks.distance(mod_low_nb), type = "h")
abline(h = 4 / nrow(Lower), col = 2, lty = 2)

plot(simres, Lower$Shannon_Index_strata)
plot(simres, Lower$num_unique_records)

car::vif(mod_low_nb)

# Assumptions met

#///////////////#
##### Shrub #####
#///////////////#

#-------#
# Model #
#-------#

Shrub <- model_data_strata %>%
  filter(Stratum == "Shrub")
mod_shrub<-glm(Mean_Dormice ~ Shannon_Index_strata+ num_unique_records, data=Shrub, family="poisson")
summary(mod_shrub)

# Call:
#   glm(formula = Mean_Dormice ~ Shannon_Index_strata + num_unique_records, 
#       family = "poisson", data = Shrub)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)          -1.590660   0.509888  -3.120  0.00181 ** 
#   Shannon_Index_strata  0.004549   0.365535   0.012  0.99007    
# num_unique_records    0.053014   0.009067   5.847    5e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 59.341  on 32  degrees of freedom
# Residual deviance: 23.610  on 30  degrees of freedom
# AIC: Inf

#-------------------#
# Check assumptions #
#-------------------#

simres <- simulateResiduals(mod_shrub)
plot(simres)

# Check dispersion
disp<- deviance(mod_shrub) / df.residual(mod_shrub)
disp 
testZeroInflation(simres)
testUniformity(simres)
testOutliers(simres)

plot(cooks.distance(mod_shrub), type = "h")
abline(h = 4 / nrow(Shrub), col = 2, lty = 2)

plot(simres, Shrub$Shannon_Index_strata)
plot(simres, Shrub$num_unique_records)

car::vif(mod_shrub)

#/////////////////////#
##### Understorey #####
#/////////////////////#

#-------#
# Model #
#-------#

Understorey <- model_data_strata %>%
  filter(Stratum == "Understorey")
mod_un<-glm(Mean_Dormice ~ Shannon_Index_strata+ num_unique_records, data=Understorey, family="poisson")
summary(mod_un)

# Call:
#   glm(formula = Mean_Dormice ~ Shannon_Index_strata + num_unique_records, 
#       family = "poisson", data = Understorey)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)          -0.779157   0.428729  -1.817   0.0692 .  
# Shannon_Index_strata -0.136844   0.246970  -0.554   0.5795    
# num_unique_records    0.040963   0.006346   6.455 1.08e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 111.995  on 55  degrees of freedom
# Residual deviance:  70.328  on 53  degrees of freedom
# AIC: Inf
# 
# Number of Fisher Scoring iterations: 6

#-------------------#
# Check assumptions #
#-------------------#

simres <- simulateResiduals(mod_un)
plot(simres)

dispersion <- deviance(mod_un) / df.residual(mod_un)
dispersion

testZeroInflation(simres)
testUniformity(simres)
testOutliers(simres)

plot(cooks.distance(mod_un), type = "h")
abline(h = 4 / nrow(Understorey), col = 2, lty = 2)

plot(simres, Understorey$Shannon_Index_strata)
plot(simres, Understorey$num_unique_records)

car::vif(mod_un)

#--------------------#
###### Plotting ######
#--------------------#

#---------#
# Palette #
#---------#

# Define consistent color palette (blues, purples, pinks)
all_species <- sort(unique(model_data_species_wide$Species))
num_species <- length(all_species)

# Use cohesive palette
cool_palette <- c(
  "#D34A68", "#E6C6DC", "pink", "#F5B856", 
  "hotpink", "#F7E7B6", "#843873", "#7A0016", 
  "#4C032A", "#D66F90", "#f48fb1", "#EAC2CA",
  "#C2185B", "#9C1C54", "#6E1E3A","#F26A6A",
  "#FF8A65","#F4A261","#F6C177")

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

#------------------------#
###### Plot Species ######
#------------------------#

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

#--------------------------------------------#
###### Plot Shannon Index across Strata ######
#--------------------------------------------#

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
         col=alpha("#D66F90",0.9),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(shan_strata$Shannon_Index_strata[shan_strata$Stratum=="Shrub"],at=2,add=T,
         col=alpha("hotpink",0.6),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(shan_strata$Shannon_Index_strata[shan_strata$Stratum=="Lower_Canopy"],at=3,add=T,
         col=alpha("#EAC2CA",0.9),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(shan_strata$Shannon_Index_strata[shan_strata$Stratum=="Upper_Canopy"],at=4,add=T,
         col=alpha("#9C1C54",0.9),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
wvioplot(shan_strata$Shannon_Index_strata[shan_strata$Stratum=="Emergent"],at=5,add=T,
         col=alpha("#D59AAE",0.9),border="grey20",wex=0.6,drawRect=T,lwd=1.5,adjust=1)
abline(h=par("usr")[4],col="white",lwd=4);abline(h=par("usr")[3],col="white",lwd=4)
abline(v=par("usr")[1],col="white",lwd=4);abline(v=par("usr")[2],col="white",lwd=4)
axis(2,at=seq(0,2.5,0.5),cex.axis=1,las=1,col="grey20",col.axis="grey20",lwd=1.5,line=-1)
axis(1, at = seq(1, 5, 1), labels = FALSE, col = "grey20", col.axis = "grey20", lwd = 1.5, line = 0.5)
mtext(side = 1, text = c("Under-\nstorey","\nShrub","Lower\ncanopy","Upper\ncanopy","\nEmergent"),
      at = seq(1, 5, 1), padj = 1.5, cex = 0.9, col = "grey20")

dev.off()

#-------------------------------------------#
###### Plot DM abundance across Strata ######
#-------------------------------------------#

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
custom_plot_colors <- c("#D66F90","hotpink","#EAC2CA","#9C1C54","#D59AAE")

#--------------#
# Store models #
#--------------#

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

#//////////////////////////////////////////#
#### Q3 DM abundance across key species ####
#//////////////////////////////////////////#

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

#-----------------#
# Fit Poisson GLM #
#-----------------#

# Already checked if interaction improved fit, it did not

model_data_filtered$Species <- relevel(
  factor(model_data_filtered$Species),
  ref = "Bramble"
)

mod <- glm(
  Mean_Dormice ~ Proportion_of_Species + Species + num_unique_records,
  family = poisson(link="log"),
  data = model_data_filtered
)

summary(mod)

# Call:
#   glm(formula = Mean_Dormice ~ Proportion_of_Species + Species + 
#         num_unique_records, family = poisson(link = "log"), data = model_data_filtered)
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)           -1.100280   0.205316  -5.359 8.37e-08 ***
#   Proportion_of_Species -0.245016   0.370946  -0.661    0.509    
# SpeciesHazel           0.131316   0.214038   0.614    0.540    
# SpeciesHoneysuckle     0.124146   0.278325   0.446    0.656    
# SpeciesOak             0.279744   0.231530   1.208    0.227    
# SpeciesSycamore        0.371633   0.268529   1.384    0.166    
# num_unique_records     0.042080   0.003616  11.639  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 377.56  on 194  degrees of freedom
# Residual deviance: 238.56  on 188  degrees of freedom
# AIC: Inf
# 
# Number of Fisher Scoring iterations: 6

#-------------------#
# Check assumptions #
#-------------------#

simres <- simulateResiduals(mod)
plot(simres)

# Check dispersion
disp<- deviance(mod) / df.residual(mod)
disp 
testZeroInflation(simres)
testUniformity(simres)
testOutliers(simres)

car::vif(mod)

plot(cooks.distance(mod), type = "h")
abline(h = 4 / nrow(model_data_filtered), col = 2, lty = 2) # Seems as though there could be an outlier

cd <- cooks.distance(mod)
which(cd > 4 / nrow(model_data_filtered))
sort(cd, decreasing = TRUE)[1:10]

# Sensitivity check:
mod_sens_1 <- glm(
  Mean_Dormice ~ Proportion_of_Species + Species + num_unique_records,
  family = poisson(link = "log"),
  data = model_data_filtered[-122, ]
)

summary(mod)
summary(mod_sens_1) # The most influential single point does not appear to be driving your main result.

#--------------------------#
# Prepare data for plotting
#--------------------------#

plot_data <- model_data_filtered %>%
  filter(!is.na(Mean_Dormice), !is.na(Proportion_of_Species))

# Custom colour palette
custom_colors <- c(
  "#E6C6DC",
  "#D66F90",
  "#f48fb1",
  "#6E1E3A",
  "#e85a95"
)

#--------------------#
###### Plotting ######
#--------------------#

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

#--------------#
# Sample sizes #
#--------------#

plot_data %>%
  group_by(Species) %>%
  summarise(n = n())

