# ==========================================================================
#  Analyze Native Range (native distribution range)
#  Analyze scale.WorldCuP.n.tdwg3 (cultivation range outside china)
#  Analyze scale.planting.China.tdwg3 (cultivation range within China)
# ==========================================================================
# Clear memory
cat("\014")
rm(list = ls())

# Load libraries
library(broom)
library(dplyr)
library(ggeffects)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(purrr)
library(readr)
library(tidyverse)
library(stringr)
library(performance)
library(pscl)
library(broom.helpers)
library(tibble)
library(lmtest)
library(VGAM)
# Set working directory (please modify according to your actual path)
setwd("YOUR_PATH_HERE")
getwd()

# Create the target directory if it does not exist
processed_dir <- "./2.results"

if (!dir.exists(processed_dir)) {
  dir.create(processed_dir, recursive = TRUE)
  cat("The folder has been created：", processed_dir, "\n")
}

# Load data
native.flora <- read_csv("./data/20250530.native_plant.matched.cultivated_plant_DO11.csv", show_col_types = FALSE)

# Filter and preprocess data
native.flora01 <- native.flora %>%
  select(taxon_name,
         nat.extent,
         planting.China.tdwg3,
         WorldCuP.n.tdwg3,
         native.global.tdwg3,
         life.form.integrated)

native.flora02 <- native.flora01 %>%
  filter(!is.na(life.form.integrated)) %>% 
  mutate(life.form.integrated = factor(life.form.integrated, levels = c("annual herb", "perennial herb", "woody")))

str(native.flora02)

# Set 'woody' as reference level
native.flora02$life.form.integrated <- relevel(native.flora02$life.form.integrated, ref = "woody")

# Standardize predictor variables
native.flora02 <- native.flora02 %>%
  mutate(
    scale.WorldCuP.n.tdwg3 = as.numeric(scale(log(WorldCuP.n.tdwg3 + 1))),
    scale.native.global.tdwg3 = as.numeric(scale(log(native.global.tdwg3 + 1))),
    scale.planting.China.tdwg3 = as.numeric(scale(log(planting.China.tdwg3 + 1)))
  )

range(native.flora02$WorldCuP.n.tdwg3)
range(native.flora02$native.global.tdwg3)
range(native.flora02$planting.China.tdwg3)

# ==========================================================================
# hurdle models with more than one 2-way interactions
# ==========================================================================
model <- hurdle(nat.extent ~ scale.native.global.tdwg3 +
                  scale.WorldCuP.n.tdwg3 +
                  scale.planting.China.tdwg3 +
                  life.form.integrated +
                  scale.native.global.tdwg3:life.form.integrated +
                  scale.WorldCuP.n.tdwg3:life.form.integrated +
                  scale.planting.China.tdwg3:life.form.integrated |
                  scale.native.global.tdwg3 +
                  scale.WorldCuP.n.tdwg3 +
                  scale.planting.China.tdwg3 +
                  life.form.integrated +
                  scale.native.global.tdwg3:life.form.integrated +
                  scale.WorldCuP.n.tdwg3:life.form.integrated +
                  scale.planting.China.tdwg3:life.form.integrated,
                dist = "negbin", data = native.flora02)

summary_model <- summary(model)
summary_model

# Extract coefficients from the count part
result_count <- as.data.frame(summary_model$coefficients$count) %>%
  tibble::rownames_to_column(var = "term") %>%
  mutate(part = "count")

# Extract coefficients from the zero part
result_zero <- as.data.frame(summary_model$coefficients$zero) %>%
  tibble::rownames_to_column(var = "term") %>%
  mutate(part = "zero")

# Combine both parts
result_all <- bind_rows(result_count, result_zero)

write_csv(result_all, "./resutls/Table.S1.csv")

# ==========================================================================
# Figure.1
# The binomial model (zero part)
# The zero truncated negative binomial model
# ==========================================================================
# Model formula
formula <- as.formula(nat.extent ~
                        scale.native.global.tdwg3 + scale.WorldCuP.n.tdwg3 + scale.planting.China.tdwg3 + life.form.integrated +
                        scale.native.global.tdwg3:life.form.integrated +
                        scale.WorldCuP.n.tdwg3:life.form.integrated +
                        scale.planting.China.tdwg3:life.form.integrated)
# The binomial model (zero part)
zero_glm <- glm(formula, data = native.flora02 %>% mutate(nat.extent = as.numeric(nat.extent > 0)), family = binomial)
summary(zero_glm)

# The zero truncated negative binomial model
data_nonzero <- subset(native.flora02, nat.extent > 0)
count_vglm <- vglm(formula, data = data_nonzero, family = posnegbinomial)
summary(count_vglm)

# Extract life.form levels
life_levels <- levels(native.flora02$life.form.integrated)

# Create sequences for continuous predictors
native_seq <- seq(min(native.flora02$scale.native.global.tdwg3, na.rm = TRUE),
                  max(native.flora02$scale.native.global.tdwg3, na.rm = TRUE),
                  length.out = 100)
worldcup_seq <- seq(min(native.flora02$scale.WorldCuP.n.tdwg3, na.rm = TRUE),
                    max(native.flora02$scale.WorldCuP.n.tdwg3, na.rm = TRUE),
                    length.out = 100)
china_seq <- seq(min(native.flora02$scale.planting.China.tdwg3, na.rm = TRUE),
                 max(native.flora02$scale.planting.China.tdwg3, na.rm = TRUE),
                 length.out = 100)

# ----------------------------
# 1. Define global plot theme and color palette
# ----------------------------
theme_standard <- theme_classic(base_size = 14, base_family = "serif") +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 1.5),
    axis.ticks = element_line(colour = "black", linewidth = 1.5),
    axis.text = element_text(colour = "black", size = 14),
    axis.title = element_text(colour = "black", size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.position = "top"
  )

# Define color scheme
color_palette <- c("woody" = "#F8AF66", "perennial herb" = "#1F78B4", "annual herb" = "#E31A1C")

# ----------------------------
# 2. Define a general function for prediction and visualization
# ----------------------------

predict_hurdle_effect <- function(model, var_name, seq_values, data, life_levels) {
  
  ## Create prediction dataset: each life form gets the full sequence of var_name
  pred_data <- expand.grid(
    life.form.integrated = life_levels,
    x_seq = seq_values
  )
  
  # Add the three scaled variables (use mean for non-focal predictors)
  pred_data$scale.native.global.tdwg3 <- mean(data$scale.native.global.tdwg3, na.rm = TRUE)
  pred_data$scale.WorldCuP.n.tdwg3 <- mean(data$scale.WorldCuP.n.tdwg3, na.rm = TRUE)
  pred_data$scale.planting.China.tdwg3 <- mean(data$scale.planting.China.tdwg3, na.rm = TRUE)
  
  # Replace focal variable with the sequence
  pred_data[[var_name]] <- pred_data$x_seq
  pred_data$x_seq <- NULL
  
  #predict
  pred_data$response <- predict(model, newdata = pred_data, type = "response")
  
  
  pred_data <- pred_data %>%
    mutate(life.form.integrated = factor(life.form.integrated,
                                         levels = c("annual herb", "perennial herb", "woody")))
  
  return(pred_data)
}

# ----------------------------
# 3. Generate predictions for each focal variable
# ----------------------------
pred_native   <- predict_hurdle_effect(model, "scale.native.global.tdwg3", native_seq, native.flora02, life_levels)
pred_worldcup <- predict_hurdle_effect(model, "scale.WorldCuP.n.tdwg3", worldcup_seq, native.flora02, life_levels)
pred_china    <- predict_hurdle_effect(model, "scale.planting.China.tdwg3", china_seq, native.flora02, life_levels)

# 4. Create observed data set
# ----------------------------
observed_data <- native.flora02 %>%
  select(scale.native.global.tdwg3, scale.WorldCuP.n.tdwg3, scale.planting.China.tdwg3,
         life.form.integrated, nat.extent) %>%
  mutate(
    observed_count = nat.extent,
    observed_binary = as.numeric(nat.extent > 0),
    life.form.integrated = factor(life.form.integrated,
                                  levels = c("annual herb", "perennial herb", "woody"))
  )

# ----------------------------
# 5. Plotting function
# ----------------------------

plot_hurdle_parts <- function(pred_data, observed_data, x_var, x_label) {
  
  # response part: expected naturalization extent
  p_response <- ggplot(pred_data,
                       aes(x = .data[[x_var]], y = response, color = life.form.integrated)) +
    geom_point(data = observed_data, aes(x = .data[[x_var]], y = observed_count), alpha = 0.4, size = 4) +
    geom_line(size = 1.5) +
    scale_color_manual(values = color_palette, name = "Life form:") +
    labs(x = x_label, y = "Naturalization\n(no. of regions)") +
    scale_x_continuous(limits = c(-1.2, 4.5), breaks = seq(-1, 4, by = 1)) +
    scale_y_continuous(limits = c(-5, 200), breaks = seq(0, 200, by = 50)) +
    theme_standard +
    theme(axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          axis.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20))
  
  list(response = p_response)
}

# ----------------------------
# Generate three sets of plots
# ----------------------------
plots_native   <- plot_hurdle_parts(pred_native, observed_data, "scale.native.global.tdwg3", "Log(native range size + 1)(scaled)")
plots_worldcup <- plot_hurdle_parts(pred_worldcup, observed_data, "scale.WorldCuP.n.tdwg3", "Log(cultivation outside China + 1) (scaled)")
plots_china    <- plot_hurdle_parts(pred_china, observed_data,"scale.planting.China.tdwg3", "Log(cultivation within China + 1) (scaled)")

# ----------------------------
# Display or save plots
# ----------------------------

plot_response_native <- plots_native$response
plot_response_china <- plots_china$response
plot_response_worldcup <- plots_worldcup$response

#  Add inset plots to the response panels
plot_add_parts01 <- function(pred_data, x_var) {
  
  zoom_plot <- ggplot(pred_data, aes(x = .data[[x_var]], y = response,  color = life.form.integrated)) +
    geom_line(linewidth = 1.5) +
    scale_color_manual(values = color_palette) +
    scale_x_continuous(limits = c(-1, 4.5), breaks = seq(-1, 4, by = 1)) +
    scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2)) +
    labs(x = "", y = "", tag = "") +
    theme_standard +
    theme(legend.position = "none",
          axis.text = element_text(size = 14),
          title = element_text(size = 14),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5))
  
}

# Generate three inset plots
plots_zoom_native01   <- plot_add_parts01 (pred_native,"scale.native.global.tdwg3")
plots_zoom_worldcup01 <- plot_add_parts01 (pred_worldcup, "scale.WorldCuP.n.tdwg3")
plots_zoom_china01    <- plot_add_parts01 (pred_china, "scale.planting.China.tdwg3")+ scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, by = 0.1))

# Convert inset plots to grob objects
plots_zoom_native_response <- ggplotGrob(plots_zoom_native01)
plots_zoom_worldcup_response <- ggplotGrob(plots_zoom_worldcup01)
plots_zoom_china_response <- ggplotGrob(plots_zoom_china01)


# Add grobs directly onto the main plots
plot_response_native02 <- plot_response_native +
  annotation_custom(grob = plots_zoom_native_response,
                    xmin = -1, xmax = 1,  
                    ymin = 120, ymax = 220)

plot_response_china02 <- plot_response_china +
  annotation_custom(grob = plots_zoom_china_response,
                    xmin = -1, xmax = 1,  
                    ymin = 120, ymax = 220)

plot_response_worldcup02 <- plot_response_worldcup +
  annotation_custom(grob = plots_zoom_worldcup_response,
                    xmin = -1, xmax = 1,  
                    ymin = 120, ymax = 220)


### Visualizing the binomial model
# ----------------------------
# 1. Function to visualize the binomial model
# ----------------------------

predict_binomial_effect <- function(model, var_name, seq_values, data, life_levels) {
  
  pred_data <- expand.grid(
    life.form.integrated = life_levels,
    x_seq = seq_values
  )
  
  # Keep other predictors constant at their mean values
  pred_data$scale.native.global.tdwg3 <- mean(data$scale.native.global.tdwg3, na.rm = TRUE)
  pred_data$scale.WorldCuP.n.tdwg3 <- mean(data$scale.WorldCuP.n.tdwg3, na.rm = TRUE)
  pred_data$scale.planting.China.tdwg3 <- mean(data$scale.planting.China.tdwg3, na.rm = TRUE)
  
  # Replace the focal variable with its sequence
  pred_data[[var_name]] <- pred_data$x_seq
  pred_data$x_seq <- NULL
  
  # Predict probability of naturalization
  pred_data$prob_naturalized <- predict(model, newdata = pred_data, type = "response")
  
  # Reorder life form factor levels
  pred_data <- pred_data %>%
    mutate(life.form.integrated = factor(life.form.integrated,
                                         levels = c("annual herb", "perennial herb", "woody")))
  
  return(pred_data)
}

# Generate predictions for each focal variable
pred_binom_native   <- predict_binomial_effect(zero_glm, "scale.native.global.tdwg3", native_seq, native.flora02, life_levels)
pred_binom_worldcup <- predict_binomial_effect(zero_glm, "scale.WorldCuP.n.tdwg3", worldcup_seq, native.flora02, life_levels)
pred_binom_china    <- predict_binomial_effect(zero_glm, "scale.planting.China.tdwg3", china_seq, native.flora02, life_levels)

# ----------------------------
# 2. Plotting function for the binomial model
# ----------------------------

plot_binomial_parts <- function(pred_data, observed_data, x_var, x_label) {
  
  p_binom <- ggplot(pred_data,
                    aes(x = .data[[x_var]], y = prob_naturalized, color = life.form.integrated)) +
    geom_point(data = observed_data,
               aes(x = .data[[x_var]], y = observed_binary),
               alpha = 0.4, size = 4,
               position = position_jitter(height = 0.05, width = 0, seed = 123)) +
    geom_line(size = 1.5) +
    scale_color_manual(values = color_palette, name = "Life form:") +
    labs(x = x_label, y = "Natur. incidence\n(no, yes)") +
    scale_x_continuous(limits = c(-1.2, 4.5), breaks = seq(-1, 4, by = 1)) +
    scale_y_continuous(limits = c(-0.05, 1.05), breaks = seq(0, 1), labels = c("0", "1")) +
    theme_standard +
    theme(axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          axis.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20))
  
  return(p_binom)
}
# Generate individual plots
plot_binom_native   <- plot_binomial_parts(pred_binom_native, observed_data, "scale.native.global.tdwg3", "Log(native range size + 1) (scaled)")
plot_binom_worldcup <- plot_binomial_parts(pred_binom_worldcup, observed_data, "scale.WorldCuP.n.tdwg3", "Log(cultivation outside China + 1) (scaled)")
plot_binom_china    <- plot_binomial_parts(pred_binom_china, observed_data, "scale.planting.China.tdwg3", "Log(cultivation within China + 1) (scaled)")


### Visualizing the zero-truncated negative binomial model
# ----------------------------
# 1. Zero-truncated negative binomial model visualization function
# ----------------------------
predict_truncnegbin_effect <- function(model, var_name, seq_values, data, life_levels) {
  
  pred_data <- expand.grid(
    life.form.integrated = life_levels,
    x_seq = seq_values
  )
  
  # Keep other predictors constant at their mean values
  pred_data$scale.native.global.tdwg3 <- mean(data$scale.native.global.tdwg3, na.rm = TRUE)
  pred_data$scale.WorldCuP.n.tdwg3 <- mean(data$scale.WorldCuP.n.tdwg3, na.rm = TRUE)
  pred_data$scale.planting.China.tdwg3 <- mean(data$scale.planting.China.tdwg3, na.rm = TRUE)
  
  pred_data[[var_name]] <- pred_data$x_seq
  pred_data$x_seq <- NULL
  
  # Predict naturalization extent (for already-naturalized species)
  pred_data$expected_count <- predict(model, newdata = pred_data, type = "response")
  
  # Reorder factor levels for life forms
  pred_data <- pred_data %>%
    mutate(life.form.integrated = factor(life.form.integrated,
                                         levels = c("annual herb", "perennial herb", "woody")))
  
  return(pred_data)
}

# Generate predictions for each focal variable
pred_trunc_native   <- predict_truncnegbin_effect(count_vglm, "scale.native.global.tdwg3", native_seq, data_nonzero, life_levels)
pred_trunc_worldcup <- predict_truncnegbin_effect(count_vglm, "scale.WorldCuP.n.tdwg3", worldcup_seq, data_nonzero, life_levels)
pred_trunc_china    <- predict_truncnegbin_effect(count_vglm, "scale.planting.China.tdwg3", china_seq, data_nonzero, life_levels)

# ----------------------------
# 2. Plotting function for the zero-truncated negative binomial model
# ----------------------------

plot_truncnegbin_parts <- function(pred_data, observed_data, x_var, x_label) {
  
  p_trunc <- ggplot(pred_data,
                    aes(x = .data[[x_var]], y = expected_count, color = life.form.integrated)) +
    geom_point(data = filter(observed_data, observed_count > 0),
               aes(x = .data[[x_var]], y = observed_count),
               alpha = 0.4, size = 4) +
    geom_line(size = 1.5) +
    scale_color_manual(values = color_palette, name = "Life form:") +
    labs(x = x_label, y = "Natur. extent\n(no. of regions)") +
    scale_x_continuous(limits = c(-1.2, 4.5), breaks = seq(-1, 4, by = 1)) +
    scale_y_continuous(limits = c(-5, 200), breaks = seq(0, 200, by = 50)) +
    theme_standard +
    theme(axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          axis.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20))
  
  return(p_trunc)
}
# Generate individual plots
plot_trunc_native   <- plot_truncnegbin_parts(pred_trunc_native, observed_data, "scale.native.global.tdwg3", "Log(native range size + 1) (scaled)")
plot_trunc_worldcup <- plot_truncnegbin_parts(pred_trunc_worldcup, observed_data, "scale.WorldCuP.n.tdwg3", "Log(cultivation outside China + 1) (scaled)")
plot_trunc_china    <- plot_truncnegbin_parts(pred_trunc_china, observed_data, "scale.planting.China.tdwg3", "Log(cultivation within China + 1) (scaled)")

# ----------------------------
# 3. Add inset (zoomed subplots)
# ----------------------------

zoom_binom_plot <- ggplot(pred_binom_china, aes(x = scale.planting.China.tdwg3, y = prob_naturalized, color = life.form.integrated)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = color_palette) +
  scale_x_continuous(limits = c(-1, 4.5), breaks = seq(-1, 4, by = 1)) +
  scale_y_continuous(limits = c(0, 0.1), breaks = c(0, 0.1), labels = c("0", "0.1")) +
  labs(x = "", y = "", tag = "") +
  theme_standard +
  theme(legend.position = "none",
        axis.text = element_text(size = 14),
        title = element_text(size = 14),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5))

# Convert the zoomed plot to a grob object
zoom_binom_plot02 <- ggplotGrob(zoom_binom_plot)

# Add the inset grob directly to the main plot
plot_binom_china02 <- plot_binom_china +
  annotation_custom(grob = zoom_binom_plot02,
                    xmin = -1.5, xmax = 0.6,  
                    ymin = 0.4, ymax = 0.95)

# Function to generate zoomed subplots for truncated models
plot_add_parts <- function(pred_data, observed_data, x_var) {
  
  zoom_plot <- ggplot(filter(observed_data, observed_count > 0), aes(x = .data[[x_var]], color = life.form.integrated)) +
    geom_line(data = pred_data, aes(x = .data[[x_var]],y = expected_count), linewidth = 1.5) +
    scale_color_manual(values = color_palette) +
    scale_x_continuous(limits = c(-1, 4.5), breaks = seq(-1, 4, by = 1)) +
    scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2)) +
    labs(x = "", y = "", tag = "") +
    theme_standard +
    theme(legend.position = "none",
          axis.text = element_text(size = 14),
          title = element_text(size = 14),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5))
  
}

# Generate three zoomed subplots
plots_zoom_native   <- plot_add_parts(pred_trunc_native, observed_data, "scale.native.global.tdwg3")+ scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, by = 4))
plots_zoom_worldcup <- plot_add_parts(pred_trunc_worldcup, observed_data, "scale.WorldCuP.n.tdwg3") + scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, by = 4))
plots_zoom_china    <- plot_add_parts(pred_trunc_china, observed_data,"scale.planting.China.tdwg3") + scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, by = 6))

# Convert zoomed plots to grob objects
plots_zoom_native02 <- ggplotGrob(plots_zoom_native)
plots_zoom_worldcup02 <- ggplotGrob(plots_zoom_worldcup)
plots_zoom_china02 <- ggplotGrob(plots_zoom_china)

# Add zoomed insets directly to the main plots
plot_trunc_native02 <- plot_trunc_native +
  annotation_custom(grob = plots_zoom_native02,
                    xmin = -1, xmax = 1,  
                    ymin = 120, ymax = 220)

plot_trunc_china02 <- plot_trunc_china +
  annotation_custom(grob = plots_zoom_china02,
                    xmin = -1, xmax = 1, 
                    ymin = 120, ymax = 220)

plot_trunc_worldcup02 <- plot_trunc_worldcup +
  annotation_custom(grob = plots_zoom_worldcup02,
                    xmin = -1, xmax = 1, 
                    ymin = 120, ymax = 220)

# ----------------------------
# 4. Combine binomial and truncated model plots
# ----------------------------

plots <-
  plot_binom_native +
  plot_trunc_native02 +
  plot_binom_china02 +
  plot_trunc_china02 +
  plot_binom_worldcup +
  plot_trunc_worldcup02 +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "top") &
  scale_color_manual(
    values = c("annual herb" = "#E31A1C",
               "perennial herb" = "#1F78B4",
               "woody" = "#F8AF66"),
    name = "Life forms:",
    labels = c("annual herb" = "Annual herb.",
               "perennial herb" = "Perennial herb.",
               "woody" = "Woody")
  )
plots
# Save the plot
ggexport(plots, filename = "./resutls/Figure.1.png",
         width = 4200,
         height = 5500,
         pointsize = 12,
         res = 300)


# ==========================================================================
# Figure.S2
# Predicting the isolated effects of life forms 
# ==========================================================================
life_levels <- levels(native.flora02$life.form.integrated)

# ----------------------------
# 1. Predicting the life form effect for the Hurdle model
# ----------------------------
predict_lifeform_effect <- function(model, data, life_levels) {
  
  pred_data <- data.frame(
    life.form.integrated = factor(life_levels, 
                                  levels = c("annual herb", "perennial herb", "woody"))
  )
  pred_data$scale.native.global.tdwg3 <- mean(data$scale.native.global.tdwg3, na.rm = TRUE)
  pred_data$scale.WorldCuP.n.tdwg3 <- mean(data$scale.WorldCuP.n.tdwg3, na.rm = TRUE)
  pred_data$scale.planting.China.tdwg3 <- mean(data$scale.planting.China.tdwg3, na.rm = TRUE)
  pred_data$response <- predict(model, newdata = pred_data, type = "response")
  return(pred_data)
}

# life form effect predictions
pred_lifeform <- predict_lifeform_effect(model, native.flora02, life_levels)
pred_lifeform

# ----------------------------
# 2. Visualisation of the life form effect (Hurdle model)
# ----------------------------
theme_standard <- theme_classic(base_size = 20, base_family = "serif") +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 1.5),
    axis.ticks = element_line(colour = "black", linewidth = 1.5),
    axis.text = element_text(colour = "black", size = 20),
    axis.title = element_text(colour = "black", size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.position = "top"
  )
## ----------------------------
plot_lifeform_hurdle <- ggplot(pred_lifeform, 
                               aes(x = life.form.integrated, 
                                   y = response, 
                                   fill = life.form.integrated)) +
  geom_col(width = 0.6) +
  scale_fill_manual(
    values = c("annual herb" = "#E31A1C",
               "perennial herb" = "#1F78B4",
               "woody" = "#F8AF66"),
    name = "Life forms:"
  ) +
  labs(x = "Life form", 
       y = "Naturalization success\n(no. of regions)") +
  scale_x_discrete(labels = c(
    "annual herb" = "Annual herb.",
    "perennial herb" = "Perennial herb.",
    "woody" = "Woody"
  )) +
  scale_y_continuous(limits = c(0, 0.12), breaks = seq(0, 0.12, by = 0.03)) +
  theme_standard +
  theme(
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 20),
    legend.position = "none"
  )

plot_lifeform_hurdle

# ----------------------------
# 3. Predict life form effects for the separated model
# ----------------------------

# 3.1 Function to visualize the binomial model
predict_lifeform_binomial <- function(model, data, life_levels) {
  pred_data <- data.frame(
    life.form.integrated = factor(life_levels, 
                                  levels = c("annual herb", "perennial herb", "woody"))
  )
  
  pred_data$scale.native.global.tdwg3 <- mean(data$scale.native.global.tdwg3, na.rm = TRUE)
  pred_data$scale.WorldCuP.n.tdwg3 <- mean(data$scale.WorldCuP.n.tdwg3, na.rm = TRUE)
  pred_data$scale.planting.China.tdwg3 <- mean(data$scale.planting.China.tdwg3, na.rm = TRUE)
  
  pred_data$prob_naturalized <- predict(model, newdata = pred_data, type = "response")
  
  return(pred_data)
}

# 3.2 Function to visualize the Zero-truncated negative binomial model
predict_lifeform_truncnb <- function(model, data, life_levels) {
  pred_data <- data.frame(
    life.form.integrated = factor(life_levels, 
                                  levels = c("annual herb", "perennial herb", "woody"))
  )
  
  pred_data$scale.native.global.tdwg3 <- mean(data$scale.native.global.tdwg3, na.rm = TRUE)
  pred_data$scale.WorldCuP.n.tdwg3 <- mean(data$scale.WorldCuP.n.tdwg3, na.rm = TRUE)
  pred_data$scale.planting.China.tdwg3 <- mean(data$scale.planting.China.tdwg3, na.rm = TRUE)
  
  pred_data$expected_count <- predict(model, newdata = pred_data, type = "response")
  
  return(pred_data)
}

# predictions
pred_lifeform_binom <- predict_lifeform_binomial(zero_glm, native.flora02, life_levels)
pred_lifeform_trunc <- predict_lifeform_truncnb(count_vglm, data_nonzero, life_levels)

# ----------------------------
# 4. Visualisation of the life form effect in separation models
# ----------------------------
# 4.1 Binomial
plot_lifeform_binom <- ggplot(pred_lifeform_binom, 
                              aes(x = life.form.integrated, 
                                  y = prob_naturalized, 
                                  fill = life.form.integrated)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = color_palette) +
  labs(x = "Life form", 
       y = "Naturalization incidence\n(no, yes)") +
  scale_x_discrete(labels = c(
    "annual herb" = "Annual herb.",
    "perennial herb" = "Perennial herb.",
    "woody" = "Woody"
  )) +
  scale_y_continuous(limits = c(0, 0.04), breaks = seq(0, 0.04, 0.01)) +
  theme_standard +
  theme(
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.position = "none"
  )

# 4.2 Truncated NB 
plot_lifeform_trunc <- ggplot(pred_lifeform_trunc, 
                              aes(x = life.form.integrated, 
                                  y = expected_count, 
                                  fill = life.form.integrated)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = color_palette) +
  labs(x = "Life form", 
       y = "Naturalization extent\n(no. of regions)") +
  scale_x_discrete(labels = c(
    "annual herb" = "Annual herb.",
    "perennial herb" = "Perennial herb.",
    "woody" = "Woody"
  )) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 5)) +
  theme_standard +
  theme(
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.position = "none"
  )

# ----------------------------
# 5. combine plots
# ----------------------------

plots_lifeform_combined <- (
  plot_lifeform_binom + plot_lifeform_trunc
) +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 20, face = "bold"))

plots_lifeform_combined


plots_lifeform_combined02 <- (
  plot_lifeform_hurdle +
  plot_lifeform_binom + 
  plot_lifeform_trunc 
) +
  plot_layout(ncol = 1) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 20, face = "bold"))

plots_lifeform_combined02

# Save the plot
ggexport(plots_lifeform_combined02, filename = "./results/Figure.S2.png",
         width = 2500,
         height = 5500,
         pointsize = 12,
         res = 300)


#===============================================================================
# Table.1
# The likelihood ratio tests(Table.1)
# ===============================================================================
# ==========================================================================
# 1. Full model formula
# ==========================================================================
formula_full <- as.formula(
  "nat.extent ~ scale.native.global.tdwg3 +
   scale.WorldCuP.n.tdwg3 +
   scale.planting.China.tdwg3 +
   life.form.integrated +
   scale.native.global.tdwg3:life.form.integrated +
   scale.WorldCuP.n.tdwg3:life.form.integrated +
   scale.planting.China.tdwg3:life.form.integrated |
   scale.native.global.tdwg3 +
   scale.WorldCuP.n.tdwg3 +
   scale.planting.China.tdwg3 +
   life.form.integrated +
   scale.native.global.tdwg3:life.form.integrated +
   scale.WorldCuP.n.tdwg3:life.form.integrated +
   scale.planting.China.tdwg3:life.form.integrated"
)

model_full <- hurdle(formula_full, dist = "negbin", data = native.flora02)

# ==========================================================================
# 2. Reduced models for Zero part
# ==========================================================================
# Remove interaction terms from ZERO part
formula_zero_no_interaction <- as.formula(
  "nat.extent ~ scale.native.global.tdwg3 +
   scale.WorldCuP.n.tdwg3 +
   scale.planting.China.tdwg3 +
   life.form.integrated +
   scale.native.global.tdwg3:life.form.integrated +
   scale.WorldCuP.n.tdwg3:life.form.integrated +
   scale.planting.China.tdwg3:life.form.integrated |
   scale.native.global.tdwg3 +
   scale.WorldCuP.n.tdwg3 +
   scale.planting.China.tdwg3 +
   life.form.integrated"
)
model_zero_no_interaction <- hurdle(formula_zero_no_interaction, dist="negbin", data=native.flora02)
# Zero: remove native range
formula_zero_no_nat <- as.formula(
  "nat.extent ~ scale.native.global.tdwg3 +
   scale.WorldCuP.n.tdwg3 +
   scale.planting.China.tdwg3 +
   life.form.integrated +
   scale.native.global.tdwg3:life.form.integrated +
   scale.WorldCuP.n.tdwg3:life.form.integrated +
   scale.planting.China.tdwg3:life.form.integrated |
   scale.WorldCuP.n.tdwg3 +
   scale.planting.China.tdwg3 +
   life.form.integrated"
)
model_zero_no_nat <- hurdle(formula_zero_no_nat, dist="negbin", data=native.flora02)

# Zero: remove WorldCuP
formula_zero_no_worldcup <- as.formula(
  "nat.extent ~ scale.native.global.tdwg3 +
   scale.WorldCuP.n.tdwg3 +
   scale.planting.China.tdwg3 +
   life.form.integrated +
   scale.native.global.tdwg3:life.form.integrated +
   scale.WorldCuP.n.tdwg3:life.form.integrated +
   scale.planting.China.tdwg3:life.form.integrated |
   scale.native.global.tdwg3 +
   scale.planting.China.tdwg3 +
   life.form.integrated"
)
model_zero_no_worldcup <- hurdle(formula_zero_no_worldcup, dist="negbin", data=native.flora02)

# Zero: remove China planting
formula_zero_no_china <- as.formula(
  "nat.extent ~ scale.native.global.tdwg3 +
   scale.WorldCuP.n.tdwg3 +
   scale.planting.China.tdwg3 +
   life.form.integrated +
   scale.native.global.tdwg3:life.form.integrated +
   scale.WorldCuP.n.tdwg3:life.form.integrated +
   scale.planting.China.tdwg3:life.form.integrated |
   scale.native.global.tdwg3 +
   scale.WorldCuP.n.tdwg3 +
   life.form.integrated"
)
model_zero_no_china <- hurdle(formula_zero_no_china, dist="negbin", data=native.flora02)

# Remove each zero part interaction
formula_zero_no_int1 <- as.formula(
  "nat.extent ~ scale.native.global.tdwg3 +
   scale.WorldCuP.n.tdwg3 +
   scale.planting.China.tdwg3 +
   life.form.integrated +
   scale.native.global.tdwg3:life.form.integrated +
   scale.WorldCuP.n.tdwg3:life.form.integrated +
   scale.planting.China.tdwg3:life.form.integrated |
   scale.native.global.tdwg3 +
   scale.WorldCuP.n.tdwg3 +
   scale.planting.China.tdwg3 +
   life.form.integrated +

   scale.WorldCuP.n.tdwg3:life.form.integrated +
   scale.planting.China.tdwg3:life.form.integrated"
)
formula_zero_no_int2 <- as.formula(
  "nat.extent ~ scale.native.global.tdwg3 +
   scale.WorldCuP.n.tdwg3 +
   scale.planting.China.tdwg3 +
   life.form.integrated +
   scale.native.global.tdwg3:life.form.integrated +
   scale.WorldCuP.n.tdwg3:life.form.integrated +
   scale.planting.China.tdwg3:life.form.integrated |
   scale.native.global.tdwg3 +
   scale.WorldCuP.n.tdwg3 +
   scale.planting.China.tdwg3 +
   life.form.integrated +
   scale.native.global.tdwg3:life.form.integrated +

   scale.planting.China.tdwg3:life.form.integrated"
)
formula_zero_no_int3 <- as.formula(
  "nat.extent ~ scale.native.global.tdwg3 +
   scale.WorldCuP.n.tdwg3 +
   scale.planting.China.tdwg3 +
   life.form.integrated +
   scale.native.global.tdwg3:life.form.integrated +
   scale.WorldCuP.n.tdwg3:life.form.integrated +
   scale.planting.China.tdwg3:life.form.integrated |
   scale.native.global.tdwg3 +
   scale.WorldCuP.n.tdwg3 +
   scale.planting.China.tdwg3 +
   life.form.integrated +
   scale.native.global.tdwg3:life.form.integrated +
   scale.WorldCuP.n.tdwg3:life.form.integrated"
)

model_zero_no_int1 <- hurdle(formula_zero_no_int1, dist="negbin", data=native.flora02)
model_zero_no_int2 <- hurdle(formula_zero_no_int2, dist="negbin", data=native.flora02)
model_zero_no_int3 <- hurdle(formula_zero_no_int3, dist="negbin", data=native.flora02)

# Remove life.form.integrated from ZERO part
formula_zero_no_life <- as.formula(
  "nat.extent ~ scale.native.global.tdwg3 +
   scale.WorldCuP.n.tdwg3 +
   scale.planting.China.tdwg3 +
   life.form.integrated +
   scale.native.global.tdwg3:life.form.integrated +
   scale.WorldCuP.n.tdwg3:life.form.integrated +
   scale.planting.China.tdwg3:life.form.integrated |
   scale.native.global.tdwg3 +
   scale.WorldCuP.n.tdwg3 +
   scale.planting.China.tdwg3"
)
model_zero_no_life <- hurdle(formula_zero_no_life, dist="negbin", data=native.flora02)

# ==========================================================================
# 3. Reduced models for Count part
# ==========================================================================
# Remove interaction terms from COUNT part
formula_count_no_interaction <- as.formula(
  "nat.extent ~ scale.native.global.tdwg3 +
   scale.WorldCuP.n.tdwg3 +
   scale.planting.China.tdwg3 +
   life.form.integrated |
   scale.native.global.tdwg3 +
   scale.WorldCuP.n.tdwg3 +
   scale.planting.China.tdwg3 +
   life.form.integrated +
   scale.native.global.tdwg3:life.form.integrated +
   scale.WorldCuP.n.tdwg3:life.form.integrated +
   scale.planting.China.tdwg3:life.form.integrated"
)
model_count_no_interaction <- hurdle(formula_count_no_interaction, dist="negbin", data=native.flora02)

# Count: remove native range
formula_count_no_nat <- as.formula(
  "nat.extent ~ 
   scale.WorldCuP.n.tdwg3 +
   scale.planting.China.tdwg3 +
   life.form.integrated |
   scale.native.global.tdwg3 +
   scale.WorldCuP.n.tdwg3 +
   scale.planting.China.tdwg3 +
   life.form.integrated +
   scale.native.global.tdwg3:life.form.integrated +
   scale.WorldCuP.n.tdwg3:life.form.integrated +
   scale.planting.China.tdwg3:life.form.integrated"
)
model_count_no_nat <- hurdle(formula_count_no_nat, dist="negbin", data=native.flora02)

# Count: remove WorldCuP
formula_count_no_worldcup <- as.formula(
  "nat.extent ~ scale.native.global.tdwg3 +

   scale.planting.China.tdwg3 +
   life.form.integrated |
   scale.native.global.tdwg3 +
   scale.WorldCuP.n.tdwg3 +
   scale.planting.China.tdwg3 +
   life.form.integrated +
   scale.native.global.tdwg3:life.form.integrated +
   scale.WorldCuP.n.tdwg3:life.form.integrated +
   scale.planting.China.tdwg3:life.form.integrated"
)
model_count_no_worldcup <- hurdle(formula_count_no_worldcup, dist="negbin", data=native.flora02)

# Count: remove China planting
formula_count_no_china <- as.formula(
  "nat.extent ~ scale.native.global.tdwg3 +
   scale.WorldCuP.n.tdwg3 +

   life.form.integrated |
   scale.native.global.tdwg3 +
   scale.WorldCuP.n.tdwg3 +
   scale.planting.China.tdwg3 +
   life.form.integrated +
   scale.native.global.tdwg3:life.form.integrated +
   scale.WorldCuP.n.tdwg3:life.form.integrated +
   scale.planting.China.tdwg3:life.form.integrated"
)
model_count_no_china <- hurdle(formula_count_no_china, dist="negbin", data=native.flora02)

# Remove each count interaction
formula_count_no_int1 <- as.formula(
  "nat.extent ~ scale.native.global.tdwg3 +
   scale.WorldCuP.n.tdwg3 +
   scale.planting.China.tdwg3 +
   life.form.integrated +
   scale.WorldCuP.n.tdwg3:life.form.integrated +
   scale.planting.China.tdwg3:life.form.integrated |
   scale.native.global.tdwg3 +
   scale.WorldCuP.n.tdwg3 +
   scale.planting.China.tdwg3 +
   life.form.integrated +
   scale.native.global.tdwg3:life.form.integrated +
   scale.WorldCuP.n.tdwg3:life.form.integrated +
   scale.planting.China.tdwg3:life.form.integrated"
)
formula_count_no_int2 <- as.formula(
  "nat.extent ~ scale.native.global.tdwg3 +
   scale.WorldCuP.n.tdwg3 +
   scale.planting.China.tdwg3 +
   life.form.integrated +
   scale.native.global.tdwg3:life.form.integrated +

   scale.planting.China.tdwg3:life.form.integrated |
   scale.native.global.tdwg3 +
   scale.WorldCuP.n.tdwg3 +
   scale.planting.China.tdwg3 +
   life.form.integrated +
   scale.native.global.tdwg3:life.form.integrated +
   scale.WorldCuP.n.tdwg3:life.form.integrated +
   scale.planting.China.tdwg3:life.form.integrated"
)
formula_count_no_int3 <- as.formula(
  "nat.extent ~ scale.native.global.tdwg3 +
   scale.WorldCuP.n.tdwg3 +
   scale.planting.China.tdwg3 +
   life.form.integrated +
   scale.native.global.tdwg3:life.form.integrated +
   scale.WorldCuP.n.tdwg3:life.form.integrated |
   scale.native.global.tdwg3 +
   scale.WorldCuP.n.tdwg3 +
   scale.planting.China.tdwg3 +
   life.form.integrated +
   scale.native.global.tdwg3:life.form.integrated +
   scale.WorldCuP.n.tdwg3:life.form.integrated +
   scale.planting.China.tdwg3:life.form.integrated"
)
model_count_no_int1 <- hurdle(formula_count_no_int1, dist="negbin", data=native.flora02)
model_count_no_int2 <- hurdle(formula_count_no_int2, dist="negbin", data=native.flora02)
model_count_no_int3 <- hurdle(formula_count_no_int3, dist="negbin", data=native.flora02)

# Remove life.form.integrated from count
formula_count_no_life <- as.formula(
  "nat.extent ~ scale.native.global.tdwg3 +
   scale.WorldCuP.n.tdwg3 +
   scale.planting.China.tdwg3|
   scale.native.global.tdwg3 +
   scale.WorldCuP.n.tdwg3 +
   scale.planting.China.tdwg3 +
   life.form.integrated +
   scale.native.global.tdwg3:life.form.integrated +
   scale.WorldCuP.n.tdwg3:life.form.integrated +
   scale.planting.China.tdwg3:life.form.integrated"
)
model_count_no_life <- hurdle(formula_count_no_life, dist="negbin", data=native.flora02)

# ==========================================================================
# 4. Compile LR test results
# ==========================================================================
results_list <- list(
  "Zero: NatRange" = tidy(lmtest::lrtest(model_zero_no_interaction, model_zero_no_nat)),
  "Zero: WorldCuP" = tidy(lmtest::lrtest(model_zero_no_interaction, model_zero_no_worldcup)),
  "Zero: ChinaPlant" = tidy(lmtest::lrtest(model_zero_no_interaction, model_zero_no_china)),
  
  "Zero: NatRange × Life"  = tidy(lmtest::lrtest(model_full, model_zero_no_int1)),
  "Zero: WorldCuP × Life"  = tidy(lmtest::lrtest(model_full, model_zero_no_int2)),
  "Zero: ChinaPlant × Life" = tidy(lmtest::lrtest(model_full, model_zero_no_int3)),
  
  "Zero: Life Form" = tidy(lmtest::lrtest(model_zero_no_interaction, model_zero_no_life)),
  
  "Count: NatRange" = tidy(lmtest::lrtest(model_count_no_interaction, model_count_no_nat)),
  "Count: WorldCuP" = tidy(lmtest::lrtest(model_count_no_interaction, model_count_no_worldcup)),
  "Count: ChinaPlant" = tidy(lmtest::lrtest(model_count_no_interaction, model_count_no_china)),
  
  "Count: NatRange × Life"  = tidy(lmtest::lrtest(model_full, model_count_no_int1)),
  "Count: WorldCuP × Life"  = tidy(lmtest::lrtest(model_full, model_count_no_int2)),
  "Count: ChinaPlant × Life" = tidy(lmtest::lrtest(model_full, model_count_no_int3)),
  
  "Count: Life Form" = tidy(lmtest::lrtest(model_count_no_interaction, model_count_no_life))
)

native_hurdle_LRT <- bind_rows(results_list,.id = "Test")
print(native_hurdle_LRT)

write_csv(native_hurdle_LRT, "./resultS/Table.1.csv")



