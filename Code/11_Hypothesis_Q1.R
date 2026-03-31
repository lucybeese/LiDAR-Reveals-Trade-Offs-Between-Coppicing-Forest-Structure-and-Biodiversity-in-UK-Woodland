##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

### Bontuchel Biodiversity Analysis ###

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

##%%%%%%%%%%%%%%%%%%%%%%%%%##

### Hypothesis testing Q1 ###

##%%%%%%%%%%%%%%%%%%%%%%%%%##

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
library(patchwork)
library(TeachingDemos)

#/////////////////#
#### LOAD DATA ####
#/////////////////#

# Load data
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

#------------------#
# Diagnostics loop #
#------------------#

par(mfrow = c(2, 2))

for (rv in names(gam_models)) {
  
  model <- gam_models[[rv]]
  
  cat("\n\n---------------------------------\n")
  cat("Diagnostics for:", rv, "\n")
  cat("---------------------------------\n")
  
  # Standard GAM diagnostic output
  gam.check(model)
  
  # Residual plots / smooth plot
  plot(model, residuals = TRUE, pch = 16, cex = 0.7, main = rv)
  
  # Shapiro only for Gaussian models
  if (family(model)$family == "gaussian") {
    shapiro_p <- shapiro.test(residuals(model))$p.value
    cat("Shapiro-Wilk p-value:", shapiro_p, "\n")
  } else {
    cat("Shapiro-Wilk not used for", family(model)$family, "model.\n")
  }
}

#-----------#
# VCI model #
#-----------#

# VCI exceeds our shapiro-wilks threshold of .001, try Gamma
summary(sp_data_unique$VCI)
min(sp_data_unique$VCI, na.rm = TRUE)

gam_vci_gamma <- gam(
  VCI ~ s(recovery_period, k = 4),
  data = sp_data_unique,
  family = Gamma(link = "log"),
  method = "REML"
)

# Check assumptions
summary(gam_vci_gamma)
gam.check(gam_vci_gamma)

simres <- DHARMa::simulateResiduals(gam_vci_gamma)
plot(simres)
DHARMa::testDispersion(simres)
DHARMa::testUniformity(simres)
DHARMa::testOutliers(simres)

plot(simres, sp_data_unique$recovery_period)

# Gamma did not improve model adequacy, keep gaussian

#------------------------#
# Optional model summary #
#------------------------#

model_summary <- do.call(
  rbind,
  lapply(names(gam_models), function(rv) {
    
    m <- gam_models[[rv]]
    s <- summary(m)
    
    data.frame(
      Response = rv,
      Family = family(m)$family,
      Link = family(m)$link,
      EDF = round(s$s.table[1, "edf"], 3),
      Smooth_p = signif(s$s.table[1, "p-value"], 3),
      Dev_Explained = round(s$dev.expl * 100, 1),
      AIC = round(AIC(m), 2)
    )
  })
)

print(model_summary)

#-------------------#
##### Plotting ######
#-------------------#

# Named colours for each response variable
panel_colors <- c(
  Max_Canopy_Height_2020 = "#D34A68",
  Gap_Proportion_2020    = "hotpink",
  LAI                    = "#7A0016",
  Intensity              = "#D66F90",
  VCI                    = "#f48fb1",
  Canopy_Volume_2020     = "#C2185B"
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





