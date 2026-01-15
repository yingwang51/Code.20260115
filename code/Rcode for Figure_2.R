# --------------------------------------------------------------------------#
# --------------------------------------------------------------------------#
# Clear workspace memory
cat("\014")
rm(list = ls())


library(readr)
library(stringr)
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(purrr)

setwd("YOUR_PATH_HERE")
getwd()

# ==========================================================================
# function
# ==========================================================================
sem_models <- function(data, model_name, response_var) {

  data <- data %>%
  filter(model == model_name) %>%
  mutate(color = factor(color))

  color <- levels(data$color)

  data.g <- as_tbl_graph(data) %>%
    activate(nodes) %>%
    mutate(name = case_when(
      name == "culChinatdwg3scaled" ~ "cul.China\n(tdwg3)",
      name == "culforeigntdwg3scaled" ~ "cul.foreign\n(tdwg3)",
      name == "nativeglobaltdwg3scaled" ~ "native.global\n(tdwg3)",
      name == "natextentscaled" ~ "Naturalization extent\n(tdwg3)"
    ))

  nodes <- data.g %>%
    activate(nodes) %>%
    as_tibble() %>%
    mutate(x = c(0, 1, -1, 0),
           y = c(0.5, 0, 0, -1))

sem.plot <- ggraph(data.g, layout = nodes, x = x, y = y) +
    geom_edge_link(
      aes(
        label = label,
        color = factor(color),
        width = width,
        start_cap = circle(20, 'mm'),
        end_cap = circle(20, 'mm')
        ),
       angle_calc = "along",
       label_dodge = unit(8, "mm"),
       label_size = 4,
       arrow = arrow(length = unit(3, "mm"), type = "closed")) +
    geom_node_point(size = 0, color = "lightblue") +
    geom_node_text(aes(label = name),
                   size = 5,
                   color = 'black',
                   hjust = c(0.5, 1, 0, 0.5)) +
    scale_edge_color_manual(values = color) +
    theme_graph() +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          legend.position = "none"
    ) +
    labs(title = model_name)

return(sem.plot)

}


# ==========================================================================
# Loading naturalization data
# ==========================================================================
model_nat <- read.csv("./resutls/model_nat.tidy02.v202501.csv")
model_nat

#
model_nat <- model_nat %>%
filter(effect == "fixed" & term != "(Intercept)") %>%
  mutate(term = str_replace_all(term, "[.]", "")) %>%
mutate(estimate = round(estimate, 3),
       Post.Prob = round(Post.Prob, 4)) %>%
mutate(Star = case_when(Post.Prob > 0.999 ~ "***",
                        Post.Prob > 0.99 ~ "**",
                        Post.Prob > 0.95 ~ "*",
                        Post.Prob < 0.95 ~ "",
                        )) %>%
mutate(color = case_when(estimate > 0 & Post.Prob > 0.95 ~ "red",
                         estimate < 0 & Post.Prob > 0.95 ~ "blue",
                         Post.Prob < 0.95 ~ "grey",
                        ))

#
model_nat01 <- model_nat %>%
transmute(from = term,
          to = response,
          model = model,
          estimate = estimate,
          Post.Prob = Post.Prob,
          width = abs(estimate) * 5,
          color = color,
          label = paste(estimate, Star)
          )



model_names <- c("model_nat.all_sp","model_nat.annual_herb","model_nat.perennial_herb","model_nat.woody")
response_vars <- c("model_nat.all_sp","model_nat.annual_herb","model_nat.perennial_herb","model_nat.woody")


sem_models.nat <- purrr::map2(model_names, response_vars, ~sem_models(model_nat01, .x, .y))
sem_models.nat <- setNames(sem_models.nat, model_names)

sem_models.nat.all_sp <- sem_models.nat$model_nat.all_sp
sem_models.nat.annual_herb <- sem_models.nat$model_nat.annual_herb
sem_models.nat.perennial_herb <- sem_models.nat$model_nat.perennial_herb
sem_models.nat.woody <- sem_models.nat$model_nat.woody


# ==========================================================================
# Aggregate data
# ==========================================================================
library(patchwork)
library(ggpubr)

#all_sp
sem_models.nat.all_sp
ggexport(sem_models.nat.all_sp, filename = "./resutls/sem_models.nat.all_sp.png",
         width = 4000,
         height = 3000,
         pointsize = 12,
         res = 300)

#life forms_sp
sem_models01 <-
  sem_models.nat.annual_herb +
  sem_models.nat.perennial_herb +
  sem_models.nat.woody +
  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "A")

sem_models01
ggexport(sem_models01, filename = "./resutls/sem_models01.png",
         width = 8000,
         height = 3000,
         pointsize = 12,
         res = 300)




















