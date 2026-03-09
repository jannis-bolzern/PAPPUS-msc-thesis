# ============================================================
# PAPPUS MSc Thesis — Bee-only network ranking (minimal)
# Output: outputs/overview_site_period_key_metrics.csv
# Network unit: Site (UGS) × Period (A/B/C)
# ============================================================

library(tidyverse)
library(readxl)
library(janitor)
library(stringr)
library(here)
library(scales)

# ----------------------------
# Paths
# ----------------------------
path_flowers      <- here("data", "raw", "pappus_data_2023", "df_flowers.xlsx")
path_pollinators  <- here("data", "raw", "pappus_data_2023", "df_pollinator.xlsx")
path_ugs          <- here("data", "raw", "ugs_info", "UGS_list.xlsx")
path_floraltraits <- here("data", "raw", "trait_data", "floral_traits", "WP4_final_trait_matrix_imputed.xlsx")
path_beetraits    <- here("data", "raw", "trait_data", "pollinator_traits", "Bee_traits_all_CH.xlsx")

# Optional later:
path_botgarden_list <- here("data", "raw", "botanical_garden_plant_list.csv")

dir.create(here("outputs"), recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Helpers
# ----------------------------
clean_taxon <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("[\u00A0]", " ") %>%
    str_squish()
}

safe_divide <- function(n, d) ifelse(is.na(d) | d == 0, NA_real_, n / d)

rescale01_robust <- function(x) {
  q <- suppressWarnings(quantile(x, probs = c(0.05, 0.95), na.rm = TRUE))
  if (any(!is.finite(q)) || q[1] == q[2]) return(rep(0.5, length(x)))
  y <- scales::rescale(x, to = c(0, 1), from = q)
  pmin(1, pmax(0, y))
}

# ----------------------------
# 1) UGS metadata
# ----------------------------
ugs <- read_excel(path_ugs, sheet = "Official") %>%
  clean_names() %>%
  transmute(
    ugs_id = as.integer(ugs_idnew),
    ugs_type = clean_taxon(ugs_type),
    site_name = clean_taxon(site),
    
    contact_name  = clean_taxon(name),
    contact_email = clean_taxon(email),
    contact_tel   = clean_taxon(tel),
    
    lat = suppressWarnings(as.numeric(lat)),
    lon = suppressWarnings(as.numeric(long))
  ) %>%
  distinct(ugs_id, .keep_all = TRUE)

# ----------------------------
# 2) Plants recorded (df_flowers): full plant list per site×period
# ----------------------------
flowers <- read_excel(path_flowers, col_types = "text") %>%
  clean_names() %>%
  transmute(
    plant = clean_taxon(scientific_name),
    id = clean_taxon(id),
    period = str_extract(id, "^[ABC]"),
    ugs_id = as.integer(str_extract(id, "(?<=_)\\d+$"))
  ) %>%
  filter(!is.na(plant), plant != "", !is.na(period), !is.na(ugs_id))

plants_recorded_by_net <- flowers %>%
  group_by(period, ugs_id) %>%
  summarise(plants_recorded = n_distinct(plant), .groups = "drop")

# ----------------------------
# 3) Bee-only interactions (df_pollinator)
# ----------------------------
bee_families <- c("Apidae", "Andrenidae", "Halictidae", "Megachilidae", "Colletidae", "Melittidae")

poll_raw <- read_excel(path_pollinators, col_types = "text") %>%
  clean_names()

poll_bee <- poll_raw %>%
  transmute(
    plant = clean_taxon(scientific_name),
    ugsid_new = clean_taxon(ugsid_new),
    period = str_extract(ugsid_new, "^[ABC]"),
    ugs_id = as.integer(str_extract(ugsid_new, "(?<=^[ABC])\\d+(?=P)")),
    amount = suppressWarnings(as.numeric(amount)),
    pollinator_group = clean_taxon(pollinator_group),
    family = clean_taxon(family),
    species_binomial = clean_taxon(species),
    morphospecies = clean_taxon(morphospecies)
  ) %>%
  mutate(amount = ifelse(is.na(amount), 0, amount)) %>%
  filter(!is.na(period), !is.na(ugs_id), !is.na(plant), plant != "") %>%
  mutate(
    # bee filter
    is_bee = (family %in% bee_families) | str_detect(str_to_lower(pollinator_group), "wildbee|bombus|apis"),
    tax_resolution = case_when(
      !is.na(species_binomial) & species_binomial != "" ~ "species",
      !is.na(morphospecies) & morphospecies != "" ~ "morphospecies",
      TRUE ~ "group"
    ),
    bee_taxon = case_when(
      tax_resolution == "species" ~ species_binomial,
      tax_resolution == "morphospecies" ~ paste0("bee_morph_", morphospecies),
      TRUE ~ "bee_group"
    )
  ) %>%
  filter(is_bee) %>%
  filter(!is.na(bee_taxon), bee_taxon != "")

# ----------------------------
# 4) Build weighted bee interaction edges per site×period
# ----------------------------
edges_site_period <- poll_bee %>%
  group_by(period, ugs_id, plant, bee_taxon) %>%
  summarise(w = sum(amount, na.rm = TRUE), .groups = "drop") %>%
  filter(w > 0)

# ----------------------------
# 5) Core network metrics (minimal, used for ranking)
# ----------------------------
network_metrics <- edges_site_period %>%
  group_by(period, ugs_id) %>%
  summarise(
    N_interactions = sum(w, na.rm = TRUE),
    P_plants = n_distinct(plant),
    A_bee_taxa = n_distinct(bee_taxon),
    L_links = n(),                         # unique plant×bee_taxon links
    singleton_rate = mean(w == 1),         # fraction of links that occur once
    .groups = "drop"
  ) %>%
  mutate(network_id = paste0(period, "_", ugs_id))

# Bee ID resolution (species-level share), bee-only now
bee_id_quality <- poll_bee %>%
  group_by(period, ugs_id) %>%
  summarise(
    bee_species_id_share = safe_divide(sum(amount[tax_resolution == "species"], na.rm = TRUE), sum(amount, na.rm = TRUE)),
    .groups = "drop"
  )

# ----------------------------
# 6) Plant traits: weighted coverage (feasibility for trait analyses)
# ----------------------------
floral_traits <- read_excel(path_floraltraits, sheet = "in") %>%
  clean_names() %>%
  transmute(
    plant = clean_taxon(species_wfo),
    # keep only what we need to define "has traits"
    flower_color = clean_taxon(flower_color),
    flower_symmetry = clean_taxon(flower_symmetry),
    pollination_syndrome = clean_taxon(pollination_syndrome),
    plant_height_m = suppressWarnings(as.numeric(plant_height_vegetative_m))
  ) %>%
  filter(!is.na(plant), plant != "") %>%
  distinct(plant, .keep_all = TRUE)

plant_weights <- edges_site_period %>%
  group_by(period, ugs_id, plant) %>%
  summarise(plant_w = sum(w), .groups = "drop") %>%
  left_join(floral_traits, by = "plant") %>%
  mutate(has_plant_traits = !is.na(flower_color) | !is.na(flower_symmetry) | !is.na(pollination_syndrome) | !is.na(plant_height_m))

plant_trait_cov <- plant_weights %>%
  group_by(period, ugs_id) %>%
  summarise(
    plant_trait_weighted_cov = safe_divide(sum(plant_w[has_plant_traits], na.rm = TRUE), sum(plant_w, na.rm = TRUE)),
    .groups = "drop"
  )

# ----------------------------
# 7) Bee traits: species-level Bee_traits_all_CH only (no individual traits)
# ----------------------------
bee_raw <- read_excel(path_beetraits, sheet = "Original+Allometrics", col_names = FALSE, col_types = "text") %>%
  clean_names()

header_row <- 2
bee_colnames <- bee_raw %>%
  slice(header_row) %>%
  unlist(use.names = FALSE) %>%
  as.character() %>%
  janitor::make_clean_names()

bee_traits <- bee_raw %>%
  slice((header_row + 1):n()) %>%
  setNames(bee_colnames) %>%
  transmute(
    taxon = clean_taxon(full_name),
    itd_mean_f = suppressWarnings(as.numeric(itd_mean_f)),
    tongue_length = suppressWarnings(as.numeric(tongue_length_tongue)),
    foraging_distance = suppressWarnings(as.numeric(foraging_distance_mfd)),
    lecty = clean_taxon(lecty),
    nesting_trait = clean_taxon(nesting_trait),
    sociality = clean_taxon(sociality)
  ) %>%
  filter(!is.na(taxon), taxon != "") %>%
  distinct(taxon, .keep_all = TRUE) %>%
  mutate(has_bee_traits = !is.na(itd_mean_f) | !is.na(tongue_length) | !is.na(foraging_distance) |
           !is.na(lecty) | !is.na(nesting_trait) | !is.na(sociality))

bee_weights <- edges_site_period %>%
  group_by(period, ugs_id, bee_taxon) %>%
  summarise(bee_w = sum(w), .groups = "drop") %>%
  left_join(bee_traits %>% select(taxon, has_bee_traits), by = c("bee_taxon" = "taxon")) %>%
  mutate(has_bee_traits = ifelse(is.na(has_bee_traits), FALSE, has_bee_traits))

bee_trait_cov <- bee_weights %>%
  group_by(period, ugs_id) %>%
  summarise(
    bee_trait_weighted_cov = safe_divide(sum(bee_w[has_bee_traits], na.rm = TRUE), sum(bee_w, na.rm = TRUE)),
    .groups = "drop"
  )

# ----------------------------
# 8) Optional BG overlap (kept minimal; not used unless file exists)
# ----------------------------
if (file.exists(path_botgarden_list)) {
  botgarden <- read_csv(path_botgarden_list, show_col_types = FALSE) %>%
    clean_names() %>%
    mutate(plant = clean_taxon(plant)) %>%
    filter(!is.na(plant), plant != "") %>%
    distinct(plant)
  
  bg_overlap <- plant_weights %>%
    select(period, ugs_id, plant, plant_w) %>%
    left_join(botgarden %>% mutate(in_botgarden = TRUE), by = "plant") %>%
    mutate(in_botgarden = ifelse(is.na(in_botgarden), FALSE, in_botgarden)) %>%
    group_by(period, ugs_id) %>%
    summarise(
      botgarden_overlap_weighted = safe_divide(sum(plant_w[in_botgarden], na.rm = TRUE), sum(plant_w, na.rm = TRUE)),
      .groups = "drop"
    )
} else {
  bg_overlap <- tibble(period = character(), ugs_id = integer(), botgarden_overlap_weighted = double())
}

# ----------------------------
# 9) Assemble + score (only the minimal factors)
# ----------------------------
scored <- network_metrics %>%
  left_join(plants_recorded_by_net, by = c("period", "ugs_id")) %>%
  left_join(bee_id_quality, by = c("period", "ugs_id")) %>%
  left_join(plant_trait_cov, by = c("period", "ugs_id")) %>%
  left_join(bee_trait_cov, by = c("period", "ugs_id")) %>%
  left_join(bg_overlap, by = c("period", "ugs_id")) %>%
  left_join(ugs, by = "ugs_id") %>%
  mutate(
    # 1) Richness / power
    s_N = rescale01_robust(N_interactions),
    s_P = rescale01_robust(P_plants),
    s_A = rescale01_robust(A_bee_taxa),
    s_L = rescale01_robust(L_links),
    score_richness = rowMeans(cbind(s_N, s_P, s_A, s_L), na.rm = TRUE),
    
    # 2) Reliability (sparsity + taxonomic resolution)
    s_singletons = rescale01_robust(1 - singleton_rate),
    s_species_id = rescale01_robust(bee_species_id_share),
    score_reliability = rowMeans(cbind(s_singletons, s_species_id), na.rm = TRUE),
    
    # 3) Trait feasibility
    s_plant_traits = rescale01_robust(plant_trait_weighted_cov),
    s_bee_traits   = rescale01_robust(bee_trait_weighted_cov),
    score_trait_feasibility = rowMeans(cbind(s_plant_traits, s_bee_traits), na.rm = TRUE),
    
    # Overall
    score_overall = 0.50 * score_richness +
      0.30 * score_reliability +
      0.20 * score_trait_feasibility,
    
    rank_overall = dense_rank(desc(score_overall))
  ) %>%
  arrange(rank_overall)

# ----------------------------
# 10) Final output: overview table only
# ----------------------------
overview_site_period <- scored %>%
  transmute(
    rank_overall,
    network_id, period, ugs_id,
    ugs_type, site_name,
    
    contact_name, contact_email, contact_tel,
    lat, lon,
    
    # plant context
    plants_recorded,
    P_plants,
    plants_interaction_share = round(safe_divide(P_plants, plants_recorded), 3),
    
    # bee network size
    N_interactions,
    A_bee_taxa,
    L_links,
    
    # reliability + traits
    singleton_rate = round(singleton_rate, 3),
    bee_species_id_share = round(bee_species_id_share, 3),
    plant_trait_weighted_cov = round(plant_trait_weighted_cov, 3),
    bee_trait_weighted_cov = round(bee_trait_weighted_cov, 3),
    
    # scores (drivers of ranking)
    score_richness = round(score_richness, 3),
    score_reliability = round(score_reliability, 3),
    score_trait_feasibility = round(score_trait_feasibility, 3),
    score_overall = round(score_overall, 3),
    
    # optional
    botgarden_overlap_weighted = round(botgarden_overlap_weighted, 3)
  )

write_csv(overview_site_period, here("outputs", "overview_site_period_key_metrics.csv"))

message("Saved: outputs/overview_site_period_key_metrics.csv")
