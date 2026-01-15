# clean 
cat("\014")
rm(list = ls())

# loading packages
library(broom)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(readr)
library(sf)
library(tidyverse)
library(scatterpie)
library(RColorBrewer) 

# Set the working directory
setwd("YOUR_PATH_HERE")
getwd()

# ==========================================================================
# Data Processing
# ==========================================================================
# Import Data
native.flora <- read_csv("./data/20250530.native_plant.matched.cultivated_plant_DO11.csv", locale = locale(encoding = "UTF-8"), show_col_types = FALSE)

# Filter data
# 2024-06-18
# Remove information on Glycine max nothosubsp. gracilis
# accepted_plant_name_id != "2995887"
native.flora <- native.flora %>% filter(accepted_plant_name_id != "2995887")

str(native.flora)
nrow(native.flora)

native.flora01 <- native.flora %>%
select(taxon_name, taxon_authors, accepted_plant_name_id, ABT:ZIM, cultivation.lin2019,life.form.integrated) %>%
distinct_all()

nrow(native.flora01)
# write_csv(native.flora01, "native.flora01.csv")

# Filter data suitable for analysis, as conflicting data exists; use conflicts to annotate.
native.flora02 <- native.flora01 %>%
  filter(cultivation.lin2019 %in% c("yes", "no")) %>%
  mutate(cultivation.lin2019 = factor(cultivation.lin2019, levels = c("no", "yes")))%>%
  filter(!is.na(life.form.integrated))
str(native.flora02)
levels(native.flora02$cultivation.lin2019)
# skimr::skim(native.flora02)
#--------------------------------------------------------------

native.flora02a <- native.flora02 %>%
pivot_longer(cols = ABT:ZIM, names_to = "tdwg_lvl03", values_to = "value") %>%
filter(value == 1)

native.flora02b <- native.flora02a %>%
group_by(tdwg_lvl03, cultivation.lin2019,life.form.integrated) %>%
summarise(n = n())

native.flora02c <- native.flora02b %>%
pivot_wider(names_from = cultivation.lin2019, values_from = n) %>%
mutate(not_culti = no, culti = yes)
native.flora02c

# Create a new logarithmic variable to avoid the log(0) issue.
native.flora03 <- native.flora02c %>%
  mutate_if(is.numeric, ~replace_na(., 0)) %>%
  mutate(total.nat = not_culti + culti,
         proportions = culti/total.nat) %>%
  mutate(culti_log = log10(culti + 1),
         not_culti_log = log10(not_culti + 1),
         total.nat_log = log10(total.nat + 1))

#Filter by life forms
native.flora03_an <- native.flora03 %>%
  filter(life.form.integrated=="annual herb")

native.flora03_pe <- native.flora03 %>%
  filter(life.form.integrated=="perennial herb")

native.flora03_wo <- native.flora03 %>%
  filter(life.form.integrated=="woody")

#Define functions: compute summary and global minimum/maximum
get_log_summary <- function(data) {
  summary_vals <- data %>%
    summarise(across(
      c(culti_log, not_culti_log, total.nat_log),
      list(min = ~min(.x, na.rm = TRUE), max = ~max(.x, na.rm = TRUE))
    ))
  
  global_min <- summary_vals %>% select(ends_with("_min")) %>% unlist() %>% min()
  global_max <- summary_vals %>% select(ends_with("_max")) %>% unlist() %>% max()
  
  return(list(summary = summary_vals, global_min = global_min, global_max = global_max))
}
res_an <- get_log_summary(native.flora03_an)
res_pe <- get_log_summary(native.flora03_pe)
res_wo <- get_log_summary(native.flora03_wo)

# ==========================================================================
# plots
# ==========================================================================
# Define file paths
tdwg_lvl03_path <- "./data/level3/level3.shp"
china_map_path <- "./data/China_map/bou2_4p.shp"
southsea_path <- "./data/China_map/south_sea.shp"

# Load maps using sf
tdwg_lvl03_map <- read_sf(tdwg_lvl03_path)
china_map <- read_sf(china_map_path)
southsea <- read_sf(southsea_path)

tdwg_lvl03_map
china_map
southsea

# 
st_crs(tdwg_lvl03_map) <- 4326
st_crs(china_map) <- 4326
st_crs(southsea) <- 4326

# Convert to data frame for ggplot2
tdwg_lvl03_map_df <- st_as_sf(tdwg_lvl03_map)
china_map_df <- st_as_sf(china_map)
southsea_df <- st_as_sf(southsea)

#naturalization information
tdwg_lvl03_map_df02a <- tdwg_lvl03_map_df %>%
  left_join(native.flora03_an, by = c("LEVEL3_COD"="tdwg_lvl03"))

tdwg_lvl03_map_df02b <- tdwg_lvl03_map_df %>%
  left_join(native.flora03_pe, by = c("LEVEL3_COD"="tdwg_lvl03"))

tdwg_lvl03_map_df02c <- tdwg_lvl03_map_df %>%
  left_join(native.flora03_wo, by = c("LEVEL3_COD"="tdwg_lvl03"))


# colour
color_palette <- colorRampPalette(brewer.pal(9, "YlOrRd"))(50)
color_palette <- alpha(color_palette, 0.8)  

##annual herb
# Global Distribution Map of Native Chinese Plants##annual herb
# ==========================================================================

# Create ggplot map of tdwg level-1 regions
tdwg_map_plot01A <- ggplot(data = tdwg_lvl03_map_df02a) +
  # geom_sf(data = tdwg_lvl01_map_df, fill = info.lvl01$map_color, color = "white", size = 0.2, alpha = 0.9) +
  geom_sf(aes(fill = total.nat_log)) +
  scale_fill_gradientn(
    colors = color_palette,
    na.value = "white",
    limits = c(res_an$global_min, res_an$global_max),
    breaks = c(res_an$global_min, res_an$global_max),
    labels = c(paste0("Min: ", round(10^res_an$global_min - 1, 0)),
               paste0("Max: ", round(10^res_an$global_max - 1, 0)))
  ) +
  geom_sf(data = china_map_df, fill = "gray90", color = "gray90", size = 0.3) +
  geom_sf(data = southsea_df, color = "gray90", size = 0.3) +
  theme_void() +
  guides(
    fill = guide_colorbar(
      title = NULL,
      barwidth = 0.5,   
      barheight = 5,    
      title.position = "top",
      label.position = "right"
    )
  ) +
  theme(
    legend.position = c(0.11, 0.3), 
    legend.direction = "vertical",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )
tdwg_map_plot01A

##perennial herb
# Global Distribution Map of Native Chinese Plants##perennial herb
# ==========================================================================

# Create ggplot map of tdwg level-1 regions
tdwg_map_plot01B <- ggplot(data = tdwg_lvl03_map_df02b) +
  # geom_sf(data = tdwg_lvl01_map_df, fill = info.lvl01$map_color, color = "white", size = 0.2, alpha = 0.9) +
  geom_sf(aes(fill = total.nat_log)) +
  scale_fill_gradientn(
    colors = color_palette,
    na.value = "white",
    limits = c(res_pe$global_min, res_pe$global_max),
    breaks = c(res_pe$global_min, res_pe$global_max),
    labels = c(paste0("Min: ", round(10^res_pe$global_min - 1, 0)),
               paste0("Max: ", round(10^res_pe$global_max - 1, 0)))
  ) +
  geom_sf(data = china_map_df, fill = "gray90", color = "gray90", size = 0.3) +
  geom_sf(data = southsea_df, color = "gray90", size = 0.3) +
  theme_void() +
  guides(
    fill = guide_colorbar(
      title = NULL,
      barwidth = 0.5,   
      barheight = 5,    
      title.position = "top",
      label.position = "right"
    )
  ) +
  theme(
    legend.position = c(0.11, 0.3), 
    legend.direction = "vertical",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )
tdwg_map_plot01B

##woody
# Global Distribution Map of Native Chinese Plants##woody
# ==========================================================================

# Create ggplot map of tdwg level-1 regions
tdwg_map_plot01C <- ggplot(data = tdwg_lvl03_map_df02c) +
  # geom_sf(data = tdwg_lvl01_map_df, fill = info.lvl01$map_color, color = "white", size = 0.2, alpha = 0.9) +
  geom_sf(aes(fill = total.nat_log)) +
  scale_fill_gradientn(
    colors = color_palette,
    na.value = "white",
    limits = c(res_wo$global_min, res_wo$global_max),
    breaks = c(res_wo$global_min, res_wo$global_max),
    labels = c(paste0("Min: ", round(10^res_wo$global_min - 1, 0)),
               paste0("Max: ", round(10^res_wo$global_max - 1, 0)))
  ) +
  geom_sf(data = china_map_df, fill = "gray90", color = "gray90", size = 0.3) +
  geom_sf(data = southsea_df, color = "gray90", size = 0.3) +
  theme_void() +
  guides(
    fill = guide_colorbar(
      title = NULL,
      barwidth = 0.5,  
      barheight = 5,    
      title.position = "top",
      label.position = "right"
    )
  ) +
  theme(
    legend.position = c(0.11, 0.3),  
    legend.direction = "vertical",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )
tdwg_map_plot01C

# life forms
# ==========================================================================
library(cowplot)
a.img_path <- "annual.png"
b.img_path <- "perennial.png"
c.img_path <- "woody.png"

tdwg_map_plot11A <- ggdraw(tdwg_map_plot01A) +
  draw_image(a.img_path,
             x = 0.17, y = 0.18,
             width = 0.35, height = 0.25, 
             hjust = 0.5, vjust = 0)
tdwg_map_plot11A

tdwg_map_plot11B <- ggdraw(tdwg_map_plot01B) +
  draw_image(b.img_path,
             x = 0.17, y = 0.18,  
             width = 0.35, height = 0.25, 
             hjust = 0.5, vjust = 0)

tdwg_map_plot11C <- ggdraw(tdwg_map_plot01C) +
  draw_image(c.img_path,
             x = 0.18, y = 0.18, 
             width = 0.25, height = 0.25, 
             hjust = 0.5, vjust = 0)


custom_labels <- c(
  "(A) Annual herbs", 
  "(B) Perennial herbs", 
  "(C) Woody plants"
)
combined_maps <- tdwg_map_plot11A +
                 tdwg_map_plot11B +  
                 tdwg_map_plot11C + 
                 plot_layout(ncol = 1) + 
                 plot_annotation(tag_levels = list(custom_labels))

combined_maps
ggexport(combined_maps, filename = "./results/figure.S1.png",
         width = 2500,
         height = 3600,
         pointsize = 12,
         res = 300)





