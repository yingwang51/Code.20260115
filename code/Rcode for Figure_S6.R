# -*- coding: utf-8 -*-
# Clear workspace memory
cat("\014")
rm(list = ls())

# loading packages
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
library(AER)

# Set working directory
setwd("YOUR_PATH_HERE")
getwd()

# Load native species dataset
native.species00 <- read_csv("./data/20250530.native_plant.matched.cultivated_plant_DO11.csv", show_col_types = FALSE)
str(native.species00)

# Select relevant variables for analysis
native.species <- native.species00 %>%
  select(taxon_name,
         nat.extent,
         planting.China.tdwg3,
         WorldCuP.n.tdwg3,
         native.global.tdwg3,
         life.form.integrated)

# ==========================================================================
# check data range
range(native.species$WorldCuP.n.tdwg3, na.rm = TRUE)
range(native.species$native.global.tdwg3, na.rm = TRUE)
range(native.species$planting.China.tdwg3, na.rm = TRUE)
# scale the data after log-transformation
# native.species$scale.Botanic_garden.n <- scale(log(native.species$Botanic_garden.n + 1))
# native.species$scale.Province.n <- scale(log(native.species$Province.n + 1))

native.species01 <- native.species %>%
  mutate(scale.native.global.tdwg3 = scale(log(native.species$native.global.tdwg3 + 1)),
         scale.WorldCuP.n.tdwg3 = scale(log(native.species$WorldCuP.n.tdwg3 + 1)),
         scale.planting.China.tdwg3 = scale(log(native.species$planting.China.tdwg3 + 1)))

str(native.species01)

# check the data structure
glimpse(native.species01)

# View(native.species)
range(native.species01$scale.native.global.tdwg3, na.rm = TRUE)
range(native.species01$scale.WorldCuP.n.tdwg3, na.rm = TRUE)
range(native.species01$scale.planting.China.tdwg3, na.rm = TRUE)


# Correlation Analysis
# ==========================================================================
library(Hmisc)
library(GGally)

# Prepare data for correlation plot
model_data <- native.species01 %>%
  select(scale.native.global.tdwg3,
         scale.WorldCuP.n.tdwg3,
         scale.planting.China.tdwg3) %>%
  mutate_all(~as.numeric(.))

colSums(is.na(model_data))
# How to customize lines in ggpairs [GGally]
# https://stackoverflow.com/questions/30858337/how-to-customize-lines-in-ggpairs-ggally
lowerFn <- function(data, mapping, method = "lm", ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_point(colour = "grey", alpha = 0.5, size = 3) +
    geom_smooth(method = method, color = "#8B1E3F", ...)
  p
}

# Create the pairwise correlation plot matrix
# https://allisonhorst.github.io/palmerpenguins/articles/intro.html?trk=public_post_comment-text
pair.plot <- ggpairs(model_data,
                     columns = c("scale.native.global.tdwg3", "scale.WorldCuP.n.tdwg3", "scale.planting.China.tdwg3"),
                     upper = list(continuous = wrap("cor",
                                                    size = 4,
                                                    color = "black")),
                     lower = list(continuous = wrap(lowerFn, method = "lm")),
                     diag = list(continuous = wrap("barDiag",
                                                   bins = 30,
                                                   fill = "#3F688C",     # NCS S 4040-R90B
                                                   color = "#1B3B6F")),
                     columnLabels = c("Native range size",
                                      "Cultivation outside China",
                                      "Cultivation within China"))
# Refine plot theme settings
pair.plot01 <- pair.plot +
  theme_bw() +
  theme(panel.background = element_rect(fill = NA),
        # axis.line = element_line(size = 1, linetype = "solid"),
        # axis.ticks = element_line(colour = "black", linetype = "solid"),
        axis.text = element_text(family = "serif", colour = "black", size = 12),
        strip.text = element_text(family = "serif", colour = "black", size = 12),
        axis.title = element_text(family = "serif", colour = "black", size = 12)) 

pair.plot01

# Export correlation plot to file
library(ggpubr)
ggexport(pair.plot01, filename = "./resutls/Figure_S6.png",
         width = 2000,
         height = 2000,
         pointsize = 12,
         res = 300)







