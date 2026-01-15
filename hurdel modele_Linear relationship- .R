# -*- coding: utf-8 -*-
# @Author: dbc
# @Date: 2024-04-01 22:00:48
# @Last Modified by: dbc
# @Last Modified time: 2025-08-06 11:12:22
# @Description: Hurdle model——linear relationship

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

# Set working directory (please modify according to your actual path)
setwd("E:/中国本地植物归化/国内外栽培、原产地与归化强度")
#setwd("D:/我的坚果云/王颖/数据分析")
getwd()

# Create the target directory if it does not exist
processed_dir <- "./2.results"

if (!dir.exists(processed_dir)) {
  dir.create(processed_dir, recursive = TRUE)
  cat("The folder has been created：", processed_dir, "\n")
}

# Load data
native.flora <- read_csv("E:/中国本地植物归化/国内外栽培、原产地与归化强度/20250530.native_plant.matched.cultivated_plant_DO11.csv", show_col_types = FALSE)
#  native.flora <- read_csv("D:/我的坚果云/王颖/目前使用代码_数据20250603/主要数据/20250530.native_plant.matched.cultivated_plant_DO11.csv", show_col_types = FALSE)

# Filter and preprocess data
native.flora01 <- native.flora %>%
  select(taxon_name,
         nat.extent,
         planting.China.tdwg3,
         WorldCuP.n.tdwg3,
         native.global.tdwg3,
         life.form.integrated)

native.flora02 <- native.flora01 %>%
  filter(!is.na(life.form.integrated)) %>% # 去除没有生活型的数据
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
# Define global plot theme and color palette
# ==========================================================================
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


# Fit the hurdle model using Negative Binomial distribution
model <- hurdle( nat.extent ~ scale.native.global.tdwg3 + scale.WorldCuP.n.tdwg3 + scale.planting.China.tdwg3 + life.form.integrated + scale.native.global.tdwg3:life.form.integrated + scale.WorldCuP.n.tdwg3:life.form.integrated + scale.planting.China.tdwg3:life.form.integrated | scale.native.global.tdwg3 + scale.WorldCuP.n.tdwg3 + scale.planting.China.tdwg3 + life.form.integrated + scale.native.global.tdwg3:life.form.integrated + scale.WorldCuP.n.tdwg3:life.form.integrated + scale.planting.China.tdwg3:life.form.integrated , dist = "negbin", data = native.flora02)   
summary_model <- summary(model)

##=================================================================================
# Analyze scale.native.global.tdwg3 (Native range size)
##=================================================================================
#Generate a sequence of the predictor variable for smooth plotting
native_range_seq <- seq(
  min(native.flora02$scale.native.global.tdwg3),
  max(native.flora02$scale.native.global.tdwg3),
  length.out = 100
)


# Create a data frame for predictions
pred_data_native <- expand.grid(
  scale.native.global.tdwg3 = native_range_seq,
  life.form.integrated = levels(native.flora02$life.form.integrated)
)

# Hold other continuous variables at their mean values
pred_data_native$scale.WorldCuP.n.tdwg3   <- mean(native.flora02$scale.WorldCuP.n.tdwg3, na.rm = TRUE)
pred_data_native$scale.planting.China.tdwg3 <- mean(native.flora02$scale.planting.China.tdwg3, na.rm = TRUE)



# Extract design matrices for both model components (Zero and Count)
X_zero <- model.matrix(delete.response(terms(model, component = "zero")), 
                       data = pred_data_native)
X_count <- model.matrix(delete.response(terms(model, component = "count")), 
                        data = pred_data_native)

coef_zero <- coef(model, "zero")
coef_count <- coef(model, "count")

# Calculate linear predictors (link scale)
pred_data_native$linear_zero <- as.vector(X_zero %*% coef_zero)
pred_data_native$linear_count <- as.vector(X_count %*% coef_count)

# Sort data frame for consistent plotting lines
pred_data_native <- pred_data_native %>%
  arrange(life.form.integrated, scale.native.global.tdwg3)

# Visualization
# Plot A: Bernoulli/Zero component (Linear Predictors)

plot_native_linear_zero <- ggplot(pred_data_native, 
                                  aes(x = scale.native.global.tdwg3, 
                                      y = linear_zero, 
                                      color = life.form.integrated)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = color_palette, name = "Life forms:") +
  labs(x = "Log(native range size + 1)(scaled)",
       y = "Linear predictor \n(zero component)") +
  scale_x_continuous(limits = c(-1.2, 4.5), breaks = seq(-1, 4, by = 1)) +
  scale_y_continuous(limits = c(-6.5, 4), breaks = seq(-6, 4, by = 2)) +
  theme_standard +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)) +
  theme(legend.position = "none")


# Plot B: Count component (Linear Predictors)
plot_native_linear_count <- ggplot(pred_data_native, 
                                   aes(x = scale.native.global.tdwg3, 
                                       y = linear_count, 
                                       color = life.form.integrated)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = color_palette, name = "Life forms:") +
  labs(x = "Log(native range size + 1)(scaled)",
       y = "Linear predictor \n(count component)") +
  scale_x_continuous(limits = c(-1.2, 4.5), breaks = seq(-1, 4, by = 1)) +
  scale_y_continuous(limits = c(-2.5, 3), breaks = seq(-2, 3, by = 1)) +
  theme_standard +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)) +
  theme(legend.position = "none")

# Plot C: Combined layout with a single shared legend
combined_plot_lineart_native <- plot_native_linear_zero + plot_native_linear_count +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

combined_plot_lineart_native





##=================================================================================
# Analyze scale.WorldCuP.n.tdwg3 (cultivation range outside China)
##=================================================================================
# Generate a sequence of the predictor variable for smooth plotting
worldcup_seq <- seq(
  min(native.flora02$scale.WorldCuP.n.tdwg3),
  max(native.flora02$scale.WorldCuP.n.tdwg3),
  length.out = 100
)

# Create a data frame for predictions
pred_data_worldcup <- expand.grid(
  scale.WorldCuP.n.tdwg3 = worldcup_seq,
  life.form.integrated = levels(native.flora02$life.form.integrated)
)

# Hold other continuous variables at their mean values
pred_data_worldcup$scale.native.global.tdwg3   <- mean(native.flora02$scale.native.global.tdwg3, na.rm = TRUE)
pred_data_worldcup$scale.planting.China.tdwg3 <- mean(native.flora02$scale.planting.China.tdwg3, na.rm = TRUE)


# Extract design matrices for both model components (Zero and Count)
X_zero <- model.matrix(delete.response(terms(model, component = "zero")), 
                       data = pred_data_worldcup)
X_count <- model.matrix(delete.response(terms(model, component = "count")), 
                        data = pred_data_worldcup)

coef_zero <- coef(model, "zero")
coef_count <- coef(model, "count")

# Calculate linear predictors (link scale)
pred_data_worldcup$linear_zero <- as.vector(X_zero %*% coef_zero)
pred_data_worldcup$linear_count <- as.vector(X_count %*% coef_count)

# Sort data frame for consistent plotting lines
pred_data_worldcup <- pred_data_worldcup %>%
  arrange(life.form.integrated, scale.WorldCuP.n.tdwg3)

# Visualization
# Plot A: Bernoulli/Zero component (Linear Predictors)
plot_worldcup_linear_zero <- ggplot(pred_data_worldcup, 
                                    aes(x = scale.WorldCuP.n.tdwg3, 
                                        y = linear_zero, 
                                        color = life.form.integrated)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = color_palette, name = "Life forms:") +
  labs(x = "Log(cultivation outside China + 1) (scaled)",
       y = "Linear predictor \n(zero component)") +
  scale_x_continuous(limits = c(-1.2, 4.5), breaks = seq(-1, 4, by = 1)) +
  scale_y_continuous(limits = c(-6.5, 4), breaks = seq(-6, 4, by = 2)) +
  theme_standard +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)) +
  theme(legend.position = "none")

# Plot B: Count component (Linear Predictors)
plot_worldcup_linear_count <- ggplot(pred_data_worldcup, 
                                     aes(x = scale.WorldCuP.n.tdwg3, 
                                         y = linear_count, 
                                         color = life.form.integrated)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = color_palette, name = "Life forms:") +
  labs(x = "Log(cultivation outside China + 1) (scaled)",
       y = "Linear predictor \n(count component)") +
  theme_standard +
  scale_x_continuous(limits = c(-1.2, 4.5), breaks = seq(-1, 4, by = 1)) +
  scale_y_continuous(limits = c(-4, 6), breaks = seq(-4, 6, by = 2)) +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)) +
  theme(legend.position = "none")

# Plot C: Combined layout with a single shared legend
combined_plot_linea_worldcup <- plot_worldcup_linear_zero + plot_worldcup_linear_count +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

combined_plot_linea_worldcup



##=================================================================================
# Analyze scale.planting.China.tdwg3 (cultivation range within China)
##=================================================================================
# Generate a sequence of the predictor variable for smooth plotting
china_seq <- seq(
  min(native.flora02$scale.planting.China.tdwg3),
  max(native.flora02$scale.planting.China.tdwg3),
  length.out = 100
)

# Create a data frame for predictions
pred_data_china <- expand.grid(
  scale.planting.China.tdwg3 = china_seq,
  life.form.integrated = levels(native.flora02$life.form.integrated)
)


# Hold other continuous variables at their mean values
pred_data_china$scale.native.global.tdwg3   <- mean(native.flora02$scale.native.global.tdwg3, na.rm = TRUE)
pred_data_china$scale.WorldCuP.n.tdwg3 <- mean(native.flora02$scale.WorldCuP.n.tdwg3, na.rm = TRUE)



# Extract design matrices for both model components (Zero and Count)
X_zero <- model.matrix(delete.response(terms(model, component = "zero")), 
                       data = pred_data_china)
X_count <- model.matrix(delete.response(terms(model, component = "count")), 
                        data = pred_data_china)

coef_zero <- coef(model, "zero")
coef_count <- coef(model, "count")

# Calculate linear predictors (link scale)
pred_data_china$linear_zero <-  as.vector(X_zero %*% coef_zero)
pred_data_china$linear_count <-  as.vector(X_count %*% coef_count)

# Sort data frame for consistent plotting lines
pred_data_china <- pred_data_china %>%
  arrange(life.form.integrated, scale.planting.China.tdwg3)

# Visualization
# Plot A: Bernoulli/Zero component (Linear Predictors)
plot_china_linear_zero <- ggplot(pred_data_china, 
                                 aes(x = scale.planting.China.tdwg3, 
                                     y = linear_zero, 
                                     color = life.form.integrated)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = color_palette, name = "Life forms:") +
  labs(x = "Log(cultivation within China) (scaled)",
       y = "Linear predictor \n(zero component)") +
  scale_x_continuous(limits = c(-1.2, 4.5), breaks = seq(-1, 4, by = 1)) +
  scale_y_continuous(limits = c(-6, -1), breaks = seq(-6, -1, by = 1)) +
  theme_standard +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)) +
  theme(legend.position = "none")

# Plot B: Count component (Linear Predictors)
plot_china_linear_count <- ggplot(pred_data_china, 
                                  aes(x = scale.planting.China.tdwg3, 
                                      y = linear_count, 
                                      color = life.form.integrated)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = color_palette, name = "Life forms:") +
  labs(x = "Log(cultivation within China) (scaled)",
       y = "Linear predictor \n(count component)") +
  scale_x_continuous(limits = c(-1.2, 4.5), breaks = seq(-1, 4, by = 1)) +
  scale_y_continuous(limits = c(-3, 2), breaks = seq(-3, 2, by = 1)) +
  theme_standard +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)) +
  theme(legend.position = "none")

# Plot C: Combined layout with a single shared legend
combined_plot_linear_china <- plot_china_linear_zero + plot_china_linear_count +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

combined_plot_linear_china



# ==========================================================================
# combined all plots
# ==========================================================================
# Combine plots
plots <-
  plot_native_linear_zero +
  plot_native_linear_count +
  plot_china_linear_zero +
  plot_china_linear_count +
  plot_worldcup_linear_zero +
  plot_worldcup_linear_count+
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
               "woody" = "Woody"),
    breaks = c("annual herb", "perennial herb", "woody") )
plots


# Save the plot
ggexport(plots, filename = "./result1030/Combine_hurdle_linearized.png",
         width = 3800,
         height = 5500,
         pointsize = 12,
         res = 300)







