# ==========================================================================
# 1. Environment Setup and Data Preparation
# ==========================================================================
# Clear workspace memory and console
cat("\014")
rm(list = ls())

# loading packages
library(readr)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggeffects)
library(ggthemes)
library(ggpubr)
library(patchwork)
library(interactions)
library(glue)
library(piecewiseSEM)

# Set working directory
setwd("YOUR_PATH_HERE")
getwd()

# Load native species dataset
native.species00 <- read_csv("./data/20250530.native_plant.matched.cultivated_plant_DO11.csv", show_col_types = FALSE)
str(native.species00)

# Standardize and log-transform variables for the whole dataset
native.species <- native.species00 %>%
  mutate(cul.China.tdwg3.scaled = scale(log(planting.China.tdwg3 + 1)),
         cul.foreign.tdwg3.scaled = scale(log(WorldCuP.n.tdwg3 + 1)),
         native.global.tdwg3.scaled = scale(log(native.global.tdwg3)),
         nat.extent.scaled = scale(log(nat.extent + 1)))

# Check range of original variables
native.species00 %>%
  select(planting.China.tdwg3,
         WorldCuP.n.tdwg3,
         native.global.tdwg3,
         nat.extent) %>%
  purrr::map(range)

# --------------------------------------------------------------------------
# Data Filtering and Within-Group Scaling
# --------------------------------------------------------------------------

# Filter: Annual Herbs
native.species01 <- native.species %>%
  filter(life.form.integrated == "annual herb")

native.species01 <- native.species01 %>%
  mutate(cul.China.tdwg3.scaled = scale(log(planting.China.tdwg3 + 1)),
         cul.foreign.tdwg3.scaled = scale(log(WorldCuP.n.tdwg3 + 1)),
         native.global.tdwg3.scaled = scale(log(native.global.tdwg3)),
         nat.extent.scaled = scale(log(nat.extent + 1)))

# Filter: Perennial Herbs
native.species02 <- native.species %>%
  filter(life.form.integrated == "perennial herb")

native.species02 <- native.species02 %>%
  mutate(cul.China.tdwg3.scaled = scale(log(planting.China.tdwg3 + 1)),
         cul.foreign.tdwg3.scaled = scale(log(WorldCuP.n.tdwg3 + 1)),
         native.global.tdwg3.scaled = scale(log(native.global.tdwg3)),
         nat.extent.scaled = scale(log(nat.extent + 1)))

# Filter: Woody Plants
native.species03 <- native.species %>%
  filter(life.form.integrated == "woody")

native.species03 <- native.species03 %>%
  mutate(cul.China.tdwg3.scaled = scale(log(planting.China.tdwg3 + 1)),
         cul.foreign.tdwg3.scaled = scale(log(WorldCuP.n.tdwg3 + 1)),
         native.global.tdwg3.scaled = scale(log(native.global.tdwg3)),
         nat.extent.scaled = scale(log(nat.extent + 1)))
# ==========================================================================
# 2. Bayesian Structural Equation Modeling (BSEM) via brms
# ==========================================================================
library(brms)

# Detect CPU cores for parallel processing
total.cores <- parallel::detectCores(logical = FALSE)
total.cores_logical <- parallel::detectCores(logical = TRUE)

#Define model formulas for the path analysis
formula_native.cul.China <- brms::bf(cul.China.tdwg3.scaled ~ native.global.tdwg3.scaled, family = gaussian())
formula_native.cul.elsewhere <- brms::bf(cul.foreign.tdwg3.scaled ~ cul.China.tdwg3.scaled + native.global.tdwg3.scaled, family = gaussian())
formula_naturalized.global <- brms::bf(nat.extent.scaled ~ cul.foreign.tdwg3.scaled + cul.China.tdwg3.scaled + native.global.tdwg3.scaled, family = gaussian())

# Fit models for different plant groups
# Model 00: All species
model_naturalized.global00 <- brms::brm(formula_native.cul.China + formula_native.cul.elsewhere + formula_naturalized.global + set_rescor(FALSE),
                                        data = native.species,
                                        chains = 4,
                                        cores = total.cores/4,
                                        iter = 2000,
                                        seed = 2024)


# Model 01: Annual herbs
model_naturalized.global01 <- brms::brm(formula_native.cul.China + formula_native.cul.elsewhere + formula_naturalized.global + set_rescor(FALSE),
                                        data = native.species01,
                                        chains = 4,
                                        cores = total.cores/4,
                                        iter = 2000,
                                        seed = 2024)

# Model 02: Perennial herbs
model_naturalized.global02 <- brms::brm(formula_native.cul.China + formula_native.cul.elsewhere + formula_naturalized.global + set_rescor(FALSE),
                                        data = native.species02,
                                        chains = 4,
                                        cores = total.cores/4,
                                        iter = 2000,
                                        seed = 2024)

# Model 03: Woody plants
model_naturalized.global03 <- brms::brm(formula_native.cul.China + formula_native.cul.elsewhere + formula_naturalized.global + set_rescor(FALSE),
                                        data = native.species03,
                                        chains = 4,
                                        cores = total.cores/4,
                                        iter = 2000,
                                        seed = 2024)

# Display model summaries
summary(model_naturalized.global00)
summary(model_naturalized.global01)
summary(model_naturalized.global02)
summary(model_naturalized.global03)

# Save model objects
save(model_naturalized.global00,
     model_naturalized.global01,
     model_naturalized.global02,
     model_naturalized.global03,
     file = "./2025.2.results/model_naturalized.globals.v202501.Rdata")


# Calculate Bayesian R-squared
library(dplyr)
r2_01 <- bayes_R2(model_naturalized.global01)
r2_02 <- bayes_R2(model_naturalized.global02)
r2_03 <- bayes_R2(model_naturalized.global03)

# Consolidate R-squared values for comparison
bind_rows(
  r2_01 %>% as.data.frame() %>% mutate(model = "Model 01"),
  r2_02 %>% as.data.frame() %>% mutate(model = "Model 02"),
  r2_03 %>% as.data.frame() %>% mutate(model = "Model 03"),
  .id = "response"
)


# ==========================================================================
# 3. Posterior Path Analysis
# ==========================================================================

# Define path parameter names: ensure these match the output of parnames(fit)
b_China_native   <- "b_culChinatdwg3scaled_native.global.tdwg3.scaled"
b_Foreign_native <- "b_culforeigntdwg3scaled_native.global.tdwg3.scaled"
b_Foreign_China  <- "b_culforeigntdwg3scaled_cul.China.tdwg3.scaled"
b_Nat_native     <- "b_natextentscaled_native.global.tdwg3.scaled"
b_Nat_China      <- "b_natextentscaled_cul.China.tdwg3.scaled"
b_Nat_Foreign    <- "b_natextentscaled_cul.foreign.tdwg3.scaled"

fits <- list(
  Annual = model_naturalized.global01,
  Perennial = model_naturalized.global02,
  Woody = model_naturalized.global03
)

# Function to calculate direct and indirect effects from posterior draws
get_path_draws <- function(fit, group){
  draws <- as_draws_df(fit)
  
  effects <- draws %>%
    transmute(
      Direct = .data[[b_Nat_native]],
      Indirect_via_China = .data[[b_China_native]] * .data[[b_Nat_China]],
      Indirect_via_Foreign = .data[[b_Foreign_native]] * .data[[b_Nat_Foreign]],
      Indirect_via_Both = .data[[b_China_native]] * .data[[b_Foreign_China]] * .data[[b_Nat_Foreign]]
    ) %>%
    pivot_longer(everything(), names_to = "Path", values_to = "Value") %>%
    mutate(Group = group)
  
  effects
}

# Consolidate all posterior draws into long format
all_long <- bind_rows(lapply(names(fits), \(g) get_path_draws(fits[[g]], g)))

# --------------------------------------------------------------------------
# Visualization: Ridge Plots for Posterior Estimates
# --------------------------------------------------------------------------
library(ggplot2)
library(ggridges)
library(scales)


# 1) Standardize Group levels and labels
all_long$Group <- factor(
  all_long$Group,
  levels = c("Annual", "Perennial", "Woody"),
  labels = c("Annual herb.", "Perennial herb.", "Woody")
)

# 2) Define Path labels
all_long$Path <- factor(
  all_long$Path,
  levels = c(
    "Direct",
    "Indirect_via_China",
    "Indirect_via_Foreign",
    "Indirect_via_Both"
  ),
  labels = c(
    "NRS → NatSuc",
    "NRS → CWC → NatSuc",
    "NRS → COC → NatSuc",
    "NRS → CWC → COC → NatSuc"
  )
)# Generate facet labels (e.g., (A) NRS -> NatSuc)
path_levels <- levels(factor(all_long$Path))
path_labs <- paste0("(", LETTERS[seq_along(path_levels)], ") ", path_levels)
names(path_labs) <- path_levels

# Define standardized theme
theme_standard <- theme_classic(base_size = 16, base_family = "serif") +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
    axis.ticks = element_line(colour = "black", linewidth = 1.2),
    axis.text  = element_text(colour = "black", size = 16),
    axis.title = element_text(colour = "black", size = 18),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 16),
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 18, hjust = 0 )
  )

# Create the ridge plot
p <- ggplot(
  all_long,
  aes(x = Value, y = Group, fill = Group)
) +
  geom_density_ridges(
    alpha = 0.7,
    scale = 1.1,
    rel_min_height = 0.01,
    color = "grey30",
    size = 0.3
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey40") +
  facet_wrap(~ Path, scales = "free_x", ncol = 1, labeller = labeller(Path = path_labs)) +
  scale_fill_manual(
    name = NULL,
    values = c(
      "Annual herb."    = "#E31A1C",
      "Perennial herb." = "#1F78B4",
      "Woody"          = "#F8AF66"
    ),
    breaks = c("Annual herb.", "Perennial herb.", "Woody")  # 控制图例顺序
  ) +
  labs(x = "Posterior effect size", y = NULL) +
  theme_standard

p

# Export the plot
ggexport(p, filename = "./2.results/Posterior estimatel.png",
         width = 1800,
         height = 4200,
         pointsize = 12,
         res = 300)


# ==========================================================================
# 4. Tidying Model Outputs and Hypothesis Testing
# ==========================================================================
library(broom.mixed)
library(purrr)
library(tidyr)
#
model_naturalized.global.tidy <- model_naturalized.global00 %>%
  tidy() %>%
  mutate(
    dir = ifelse(estimate > 0, " > 0", " < 0"),
    hyp = paste0(response, "_", term, dir),
    Hypothesis = paste0("(", response, "_", term, ")", dir)
  )

model_naturalized.global.tidy01 <- model_naturalized.global.tidy %>%
  filter(effect == "fixed" & term != "(Intercept)") %>%
  mutate(hyp.test = purrr::map(hyp, ~ hypothesis(model_naturalized.global00, .x, seed = 2024))) %>%
  mutate(hypothesis = purrr::map(hyp.test, "hypothesis"))

model_naturalized.global.hyp <- model_naturalized.global.tidy01 %>%
  pull(hypothesis) %>%
  purrr::map_dfr(~.x)

model_naturalized.global.tidy02 <- model_naturalized.global.tidy %>%
  left_join(model_naturalized.global.hyp, by = "Hypothesis") %>%
  mutate(model = "model_nat.all_sp")

#
a.model_naturalized.global.tidy <- model_naturalized.global01 %>%
  tidy() %>%
  mutate(
    dir = ifelse(estimate > 0, " > 0", " < 0"),
    hyp = paste0(response, "_", term, dir),
    Hypothesis = paste0("(", response, "_", term, ")", dir)
  )

a.model_naturalized.global.tidy01 <- a.model_naturalized.global.tidy %>%
  filter(effect == "fixed" & term != "(Intercept)") %>%
  mutate(hyp.test = map(hyp, ~ hypothesis(model_naturalized.global01, .x, seed = 2024))) %>%
  mutate(hypothesis = map(hyp.test, "hypothesis"))

a.model_naturalized.global.hyp <- a.model_naturalized.global.tidy01 %>%
  pull(hypothesis) %>%
  map_dfr(~.x)

a.model_naturalized.global.tidy02 <- a.model_naturalized.global.tidy %>%
left_join(a.model_naturalized.global.hyp, by = "Hypothesis") %>%
mutate(model = "model_nat.annual_herb")

#
b.model_naturalized.global.tidy <- model_naturalized.global02 %>%
  tidy() %>%
  mutate(
    dir = ifelse(estimate > 0, " > 0", " < 0"),
    hyp = paste0(response, "_", term, dir),
    Hypothesis = paste0("(", response, "_", term, ")", dir)
  )

b.model_naturalized.global.tidy01 <- b.model_naturalized.global.tidy %>%
  filter(effect == "fixed" & term != "(Intercept)") %>%
  mutate(hyp.test = map(hyp, ~ hypothesis(model_naturalized.global02, .x, seed = 2024))) %>%
  mutate(hypothesis = map(hyp.test, "hypothesis"))

b.model_naturalized.global.hyp <- b.model_naturalized.global.tidy01 %>%
  pull(hypothesis) %>%
  map_dfr(~.x)

b.model_naturalized.global.tidy02 <- b.model_naturalized.global.tidy %>%
  left_join(b.model_naturalized.global.hyp, by = "Hypothesis") %>%
  mutate(model = "model_nat.perennial_herb")

#
c.model_naturalized.global.tidy <- model_naturalized.global03 %>%
  tidy() %>%
  mutate(
    dir = ifelse(estimate > 0, " > 0", " < 0"),
    hyp = paste0(response, "_", term, dir),
    Hypothesis = paste0("(", response, "_", term, ")", dir)
  )

c.model_naturalized.global.tidy01 <- c.model_naturalized.global.tidy %>%
  filter(effect == "fixed" & term != "(Intercept)") %>%
  mutate(hyp.test = map(hyp, ~ hypothesis(model_naturalized.global03, .x, seed = 2024))) %>%
  mutate(hypothesis = map(hyp.test, "hypothesis"))

c.model_naturalized.global.hyp <- c.model_naturalized.global.tidy01 %>%
  pull(hypothesis) %>%
  map_dfr(~.x)

c.model_naturalized.global.tidy02 <- c.model_naturalized.global.tidy %>%
  left_join(c.model_naturalized.global.hyp, by = "Hypothesis") %>%
  mutate(model = "model_nat.woody")


model_nat.tidy02 <- rbind(model_naturalized.global.tidy02,
                          a.model_naturalized.global.tidy02,
                          b.model_naturalized.global.tidy02,
                          c.model_naturalized.global.tidy02)

write_csv(model_nat.tidy02, "./2025.2.results/model_nat.tidy02.v202501.csv")



                        

