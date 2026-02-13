# ============================================================
# PAPPUS MSc Thesis — Network scoring + trait matching pipeline
# Author: Jannis Bolzern
# Purpose: Compute data-quality / feasibility metrics per network
# Network unit: Site (UGS) × Period (A/B/C); also pooled per site
# ============================================================

# ----------------------------
# 0) Packages
# ----------------------------
required_pkgs <- c(
  "tidyverse", "readxl", "janitor", "stringr", "here",
  "bipartite", "scales"
)

to_install <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

library(tidyverse)
library(readxl)
library(janitor)
library(stringr)
library(here)
library(bipartite)
library(scales)

# ----------------------------
# 1) Paths (put raw files in data/raw/)
# ----------------------------
path_flowers      <- here("data", "raw", "pappus_data_2023", "df_flowers.xlsx")
path_pollinators  <- here("data", "raw", "pappus_data_2023", "df_pollinator.xlsx")
path_ugs          <- here("data", "raw", "ugs_info", "UGS_list.xlsx")
path_floraltraits <- here("data", "raw", "trait_data", "floral_traits", "WP4_final_trait_matrix_imputed.xlsx")
path_beetraits    <- here("data", "raw", "trait_data", "pollinator_traits", "Bee_traits_all_CH.xlsx")
path_indivtraits  <- here("data", "raw", "trait_data", "pollinator_traits", "BetterGardens", "06_trait_data", "individual_traits.csv")

dir.create(here("data", "processed"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("outputs"), recursive = TRUE, showWarnings = FALSE)

# Optional: Botanical garden list (CSV with column "plant")
path_botgarden_list <- here("data", "raw", "botanical_garden_plant_list.csv")

# ----------------------------
# 2) Helpers
# ----------------------------
clean_taxon <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_replace_all("[\u00A0]", " ") %>%  # non-breaking spaces
    stringr::str_squish()
}

safe_divide <- function(n, d) {
  ifelse(is.na(d) | d == 0, NA_real_, n / d)
}

shannon_evenness <- function(w) {
  w <- w[w > 0 & !is.na(w)]
  if (length(w) <= 1) return(NA_real_)
  p <- w / sum(w)
  H <- -sum(p * log(p))
  H / log(length(w))
}

rescale01_robust <- function(x) {
  q <- suppressWarnings(quantile(x, probs = c(0.05, 0.95), na.rm = TRUE))
  if (any(!is.finite(q)) || q[1] == q[2]) return(rep(0.5, length(x)))
  y <- scales::rescale(x, to = c(0, 1), from = q)
  pmin(1, pmax(0, y))
}

check_unique_keys <- function(df, keys, name = deparse(substitute(df))) {
  n_all <- nrow(df)
  n_u <- df %>% distinct(across(all_of(keys))) %>% nrow()
  if (n_all != n_u) {
    stop(
      sprintf("'%s' is not unique by keys: %s (n=%s, unique=%s). Fix before joining.",
              name, paste(keys, collapse = ", "), n_all, n_u),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# Build bipartite matrix for a single network
edges_to_web <- function(df_edges) {
  as.matrix(xtabs(w ~ plant + pollinator, data = df_edges))
}

# Compute core network metrics for one network
compute_network_metrics <- function(df_edges, period_val, ugs_id_val) {
  P <- dplyr::n_distinct(df_edges$plant)
  A <- dplyr::n_distinct(df_edges$pollinator)
  N <- sum(df_edges$w, na.rm = TRUE)
  L <- nrow(df_edges)
  
  singleton_rate <- mean(df_edges$w == 1, na.rm = TRUE)
  connectance <- ifelse(P > 0 && A > 0, L / (P * A), NA_real_)
  link_evenness <- shannon_evenness(df_edges$w)
  
  web <- edges_to_web(df_edges)
  
  # IMPORTANT: H2fun returns 4 values; we only want H2
  H2 <- tryCatch(
    unname(as.numeric(bipartite::H2fun(web, H2_integer = FALSE)["H2"])),
    error = function(e) NA_real_
  )
  
  # Only compute nestedness when network is sufficiently large
  calc_nested <- (P >= 5) && (A >= 5) && (N >= 20) && (L >= 20)
  
  # networklevel returns a named numeric; take first element
  WNODF <- if (calc_nested) {
    tryCatch(
      unname(as.numeric(bipartite::networklevel(web, index = "weighted NODF")[1])),
      error = function(e) NA_real_
    )
  } else {
    NA_real_
  }
  
  tibble::tibble(
    period = period_val,
    ugs_id = ugs_id_val,
    network_id = paste0(period_val, "_", ugs_id_val),
    N_interactions = N,
    P_plants = P,
    A_pollinators = A,
    L_links = L,
    singleton_rate = singleton_rate,
    connectance = connectance,
    link_evenness = link_evenness,
    H2_specialization = H2,
    WNODF_nestedness = WNODF
  )
}

# ----------------------------
# 3) Read UGS metadata
# ----------------------------
ugs <- read_excel(path_ugs, sheet = "Official") %>%
  clean_names() %>%
  transmute(
    ugs_id = as.integer(ugs_idnew),
    ugs_type = clean_taxon(ugs_type),
    site_name = clean_taxon(site),
    lat = lat,
    lon = long
  )

check_unique_keys(ugs, c("ugs_id"), "ugs")

# ----------------------------
# 4) Read + clean plant observations (df_flowers)
# ----------------------------
flowers_raw <- read_excel(path_flowers, sheet = 1, col_types = "text") %>%
  clean_names()

flowers <- flowers_raw %>%
  transmute(
    plant = clean_taxon(scientific_name),
    plant_raw = clean_taxon(old_species),
    id = clean_taxon(id),
    period = str_extract(id, "^[ABC]"),
    ugs_id = as.integer(str_extract(id, "(?<=_)\\d+$")),
    availability = suppressWarnings(as.numeric(unnamed_2)),
    microhabitat = clean_taxon(unnamed_3),
    quadrat_id = clean_taxon(unnamed_4)
  ) %>%
  filter(!is.na(period), !is.na(ugs_id), !is.na(plant), plant != "")

plants_recorded_by_net <- flowers %>%
  group_by(period, ugs_id) %>%
  summarise(
    plants_recorded = n_distinct(plant),
    plant_records_rows = n(),
    availability_sum = sum(availability, na.rm = TRUE),
    availability_n_nonmissing = sum(!is.na(availability)),
    .groups = "drop"
  )

check_unique_keys(plants_recorded_by_net, c("period", "ugs_id"), "plants_recorded_by_net")

# ----------------------------
# 5) Read + clean pollinator interactions (df_pollinator)
# ----------------------------
poll_raw <- read_excel(path_pollinators, sheet = 1, col_types = "text") %>%
  clean_names()

bee_families <- c("Apidae", "Andrenidae", "Halictidae", "Megachilidae", "Colletidae", "Melittidae")

poll <- poll_raw %>%
  mutate(
    is_duplicate = tolower(clean_taxon(is_duplicate)) == "true",
    is_duplicate_across = tolower(clean_taxon(is_duplicate_across)) == "true"
  ) %>%
  filter(!is_duplicate, !is_duplicate_across) %>%
  transmute(
    plant = clean_taxon(scientific_name),
    ugsid_new = clean_taxon(ugsid_new),
    period = str_extract(ugsid_new, "^[ABC]"),
    ugs_id = as.integer(str_extract(ugsid_new, "(?<=^[ABC])\\d+(?=P)")),
    plant_obs_id = as.integer(str_extract(ugsid_new, "(?<=P)\\d+$")),
    amount = suppressWarnings(as.numeric(amount)),
    pollinator_group = clean_taxon(pollinator_group),
    order = clean_taxon(order),
    family = clean_taxon(family),
    species_binomial = clean_taxon(species),
    morphospecies = clean_taxon(morphospecies),
    tax_resolution = case_when(
      !is.na(species_binomial) & species_binomial != "" ~ "species",
      !is.na(morphospecies) & morphospecies != "" ~ "morphospecies",
      TRUE ~ "group"
    ),
    pollinator = case_when(
      tax_resolution == "species" ~ species_binomial,
      tax_resolution == "morphospecies" ~ paste0(pollinator_group, "_", morphospecies),
      TRUE ~ pollinator_group
    ),
    is_bee = (family %in% bee_families) | str_detect(str_to_lower(pollinator_group), "wildbee|bombus|apis")
  ) %>%
  filter(!is.na(period), !is.na(ugs_id), !is.na(plant), plant != "", !is.na(pollinator), pollinator != "") %>%
  mutate(amount = ifelse(is.na(amount), 0, amount))

# Aggregate to edges (plant x pollinator)
edges_site_period <- poll %>%
  group_by(period, ugs_id, plant, pollinator) %>%
  summarise(w = sum(amount, na.rm = TRUE), .groups = "drop") %>%
  filter(w > 0)

check_unique_keys(edges_site_period, c("period", "ugs_id", "plant", "pollinator"), "edges_site_period")

# Pollinator metadata (unique per pollinator label)
pollinator_meta <- poll %>%
  group_by(pollinator) %>%
  summarise(
    tax_resolution = case_when(
      any(tax_resolution == "species", na.rm = TRUE) ~ "species",
      any(tax_resolution == "morphospecies", na.rm = TRUE) ~ "morphospecies",
      TRUE ~ "group"
    ),
    is_bee = any(is_bee, na.rm = TRUE),
    order  = if (all(is.na(order)))  NA_character_ else names(sort(table(order),  decreasing = TRUE))[1],
    family = if (all(is.na(family))) NA_character_ else names(sort(table(family), decreasing = TRUE))[1],
    .groups = "drop"
  )

check_unique_keys(pollinator_meta, c("pollinator"), "pollinator_meta")

# ----------------------------
# 6) Compute network metrics (site × period)
# ----------------------------
net_tbl <- edges_site_period %>%
  group_by(period, ugs_id) %>%
  group_nest()

network_metrics_site_period <- purrr::pmap_dfr(
  list(net_tbl$data, net_tbl$period, net_tbl$ugs_id),
  ~ compute_network_metrics(..1, period_val = ..2, ugs_id_val = ..3)
)

check_unique_keys(network_metrics_site_period, c("period", "ugs_id"), "network_metrics_site_period")

# ID-resolution metrics
poll_resolution_site_period <- poll %>%
  group_by(period, ugs_id) %>%
  summarise(
    share_interactions_species_id = safe_divide(sum(amount[tax_resolution == "species"], na.rm = TRUE), sum(amount, na.rm = TRUE)),
    share_taxa_species_id = safe_divide(n_distinct(pollinator[tax_resolution == "species"]), n_distinct(pollinator)),
    share_bee_interactions = safe_divide(sum(amount[is_bee], na.rm = TRUE), sum(amount, na.rm = TRUE)),
    .groups = "drop"
  )

check_unique_keys(poll_resolution_site_period, c("period", "ugs_id"), "poll_resolution_site_period")

network_metrics_site_period <- network_metrics_site_period %>%
  left_join(plants_recorded_by_net, by = c("period", "ugs_id")) %>%
  left_join(poll_resolution_site_period, by = c("period", "ugs_id")) %>%
  left_join(ugs, by = "ugs_id")

check_unique_keys(network_metrics_site_period, c("period", "ugs_id"), "network_metrics_site_period_after_joins")

# ----------------------------
# 7) Read + prep floral traits
# ----------------------------
floral_traits_raw <- read_excel(path_floraltraits, sheet = "in") %>%
  clean_names()

floral_traits <- floral_traits_raw %>%
  transmute(
    plant = clean_taxon(species_wfo),
    is_alien = as.integer(is_alien),
    is_invasive = as.integer(is_invasive),
    is_commercial = as.integer(is_commercial),
    flower_color = clean_taxon(flower_color),
    flower_symmetry = clean_taxon(flower_symmetry),
    pollination_syndrome = clean_taxon(pollination_syndrome),
    flowering_start = suppressWarnings(as.numeric(flowering_start)),
    flowering_duration = suppressWarnings(as.numeric(flowering_duration)),
    plant_height_m = suppressWarnings(as.numeric(plant_height_vegetative_m))
  ) %>%
  filter(!is.na(plant), plant != "")

check_unique_keys(floral_traits %>% distinct(plant, .keep_all = TRUE), c("plant"), "floral_traits (expected unique plant)")

plant_weights <- edges_site_period %>%
  group_by(period, ugs_id, plant) %>%
  summarise(plant_w = sum(w, na.rm = TRUE), .groups = "drop")

plant_trait_summary <- plant_weights %>%
  left_join(floral_traits, by = "plant") %>%
  mutate(has_plant_traits = !is.na(is_alien) | !is.na(flower_color) | !is.na(plant_height_m)) %>%
  group_by(period, ugs_id) %>%
  summarise(
    plant_trait_node_cov = mean(has_plant_traits, na.rm = TRUE),
    plant_trait_weighted_cov = safe_divide(sum(plant_w[has_plant_traits], na.rm = TRUE), sum(plant_w, na.rm = TRUE)),
    alien_share_unweighted = mean(is_alien == 1, na.rm = TRUE),
    alien_share_weighted = safe_divide(sum(plant_w[is_alien == 1], na.rm = TRUE), sum(plant_w, na.rm = TRUE)),
    flower_color_richness = n_distinct(flower_color[!is.na(flower_color) & flower_color != ""]),
    symmetry_richness = n_distinct(flower_symmetry[!is.na(flower_symmetry) & flower_symmetry != ""]),
    pollination_syndrome_richness = n_distinct(pollination_syndrome[!is.na(pollination_syndrome) & pollination_syndrome != ""]),
    plant_height_range = ifelse(all(is.na(plant_height_m)), NA_real_, max(plant_height_m, na.rm = TRUE) - min(plant_height_m, na.rm = TRUE)),
    flowering_start_range = ifelse(all(is.na(flowering_start)), NA_real_, max(flowering_start, na.rm = TRUE) - min(flowering_start, na.rm = TRUE)),
    .groups = "drop"
  )

check_unique_keys(plant_trait_summary, c("period", "ugs_id"), "plant_trait_summary")

network_metrics_site_period <- network_metrics_site_period %>%
  left_join(plant_trait_summary, by = c("period", "ugs_id"))

check_unique_keys(network_metrics_site_period, c("period", "ugs_id"), "network_metrics_site_period_after_plant_traits")

# ----------------------------
# 8) Read + prep pollinator traits (bee backbone + individual enrichment)
# ----------------------------
bee_raw <- read_excel(
  path_beetraits,
  sheet = "Original+Allometrics",
  col_names = FALSE,
  col_types = "text"
) %>% clean_names()

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
    genus = clean_taxon(genus),
    family = clean_taxon(family),
    itd_mean_f = suppressWarnings(as.numeric(itd_mean_f)),
    tongue_length_tongue = suppressWarnings(as.numeric(tongue_length_tongue)),
    foraging_distance_mfd = suppressWarnings(as.numeric(foraging_distance_mfd)),
    lecty = clean_taxon(lecty),
    nesting_trait = clean_taxon(nesting_trait),
    sociality = clean_taxon(sociality)
  ) %>%
  filter(!is.na(taxon), taxon != "")

bee_traits_genus <- bee_traits %>%
  group_by(genus) %>%
  summarise(
    itd_mean_f_genus = ifelse(all(is.na(itd_mean_f)), NA_real_, mean(itd_mean_f, na.rm = TRUE)),
    tongue_length_tongue_genus = ifelse(all(is.na(tongue_length_tongue)), NA_real_, mean(tongue_length_tongue, na.rm = TRUE)),
    foraging_distance_mfd_genus = ifelse(all(is.na(foraging_distance_mfd)), NA_real_, mean(foraging_distance_mfd, na.rm = TRUE)),
    .groups = "drop"
  )

indiv <- read_csv(path_indivtraits, show_col_types = FALSE) %>%
  mutate(
    taxon = clean_taxon(taxon),
    intertegular_distance = na_if(intertegular_distance, 0),
    proboscis_length = na_if(proboscis_length, 0),
    prementum_length = na_if(prementum_length, 0),
    forewing_length = na_if(forewing_length, 0)
  )

indiv_species_means <- indiv %>%
  group_by(taxon) %>%
  summarise(
    itd_indiv = ifelse(all(is.na(intertegular_distance)), NA_real_, mean(intertegular_distance, na.rm = TRUE)),
    proboscis_indiv = ifelse(all(is.na(proboscis_length)), NA_real_, mean(proboscis_length, na.rm = TRUE)),
    prementum_indiv = ifelse(all(is.na(prementum_length)), NA_real_, mean(prementum_length, na.rm = TRUE)),
    forewing_indiv = ifelse(all(is.na(forewing_length)), NA_real_, mean(forewing_length, na.rm = TRUE)),
    n_individuals = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(taxon), taxon != "")

poll_traits <- full_join(bee_traits, indiv_species_means, by = "taxon") %>%
  mutate(
    itd_final = coalesce(itd_mean_f, itd_indiv),
    tongue_final = coalesce(tongue_length_tongue, proboscis_indiv),
    forage_final = foraging_distance_mfd
  )

poll_traits_u <- poll_traits %>%
  group_by(taxon) %>%
  summarise(
    genus = first(na.omit(genus)),
    family = first(na.omit(family)),
    lecty = first(na.omit(lecty)),
    nesting_trait = first(na.omit(nesting_trait)),
    sociality = first(na.omit(sociality)),
    itd_final = if (all(is.na(itd_final))) NA_real_ else mean(itd_final, na.rm = TRUE),
    tongue_final = if (all(is.na(tongue_final))) NA_real_ else mean(tongue_final, na.rm = TRUE),
    forage_final = if (all(is.na(forage_final))) NA_real_ else mean(forage_final, na.rm = TRUE),
    n_individuals = if (all(is.na(n_individuals))) NA_integer_ else max(n_individuals, na.rm = TRUE),
    .groups = "drop"
  )

check_unique_keys(poll_traits_u, c("taxon"), "poll_traits_u")

# ----------------------------
# 9) Trait coverage + variability per network (pollinators)
# ----------------------------
poll_weights <- edges_site_period %>%
  group_by(period, ugs_id, pollinator) %>%
  summarise(poll_w = sum(w, na.rm = TRUE), .groups = "drop") %>%
  left_join(pollinator_meta, by = "pollinator") %>%
  mutate(
    genus = if_else(str_detect(pollinator, " "), str_extract(pollinator, "^[^ ]+"), NA_character_)
  )

poll_trait_joined <- poll_weights %>%
  left_join(poll_traits_u %>% select(-genus), by = c("pollinator" = "taxon")) %>%
  left_join(bee_traits_genus, by = "genus") %>%
  mutate(
    itd_final2    = coalesce(itd_final, itd_mean_f_genus),
    tongue_final2 = coalesce(tongue_final, tongue_length_tongue_genus),
    forage_final2 = coalesce(forage_final, foraging_distance_mfd_genus),
    has_poll_traits = !is.na(itd_final2) | !is.na(tongue_final2) | !is.na(forage_final2)
  )

poll_trait_summary <- poll_trait_joined %>%
  group_by(period, ugs_id) %>%
  summarise(
    poll_trait_node_cov = mean(has_poll_traits, na.rm = TRUE),
    poll_trait_weighted_cov = safe_divide(sum(poll_w[has_poll_traits], na.rm = TRUE), sum(poll_w, na.rm = TRUE)),
    bee_trait_weighted_cov = {
      bee_w <- sum(poll_w[is_bee], na.rm = TRUE)
      bee_w_traits <- sum(poll_w[is_bee & has_poll_traits], na.rm = TRUE)
      safe_divide(bee_w_traits, bee_w)
    },
    bee_itd_range = {
      v <- itd_final2[is_bee & has_poll_traits]
      if (length(v[is.finite(v)]) == 0) NA_real_ else max(v, na.rm = TRUE) - min(v, na.rm = TRUE)
    },
    bee_tongue_range = {
      v <- tongue_final2[is_bee & has_poll_traits]
      if (length(v[is.finite(v)]) == 0) NA_real_ else max(v, na.rm = TRUE) - min(v, na.rm = TRUE)
    },
    lecty_richness = n_distinct(lecty[is_bee & !is.na(lecty) & lecty != ""]),
    nesting_richness = n_distinct(nesting_trait[is_bee & !is.na(nesting_trait) & nesting_trait != ""]),
    sociality_richness = n_distinct(sociality[is_bee & !is.na(sociality) & sociality != ""]),
    .groups = "drop"
  )

check_unique_keys(poll_trait_summary, c("period", "ugs_id"), "poll_trait_summary")

network_metrics_site_period <- network_metrics_site_period %>%
  left_join(poll_trait_summary, by = c("period", "ugs_id"))

check_unique_keys(network_metrics_site_period, c("period", "ugs_id"), "network_metrics_site_period_after_poll_traits")

# ----------------------------
# 10) Optional: Botanical Garden overlap (if list exists)
# ----------------------------
if (!is.null(path_botgarden_list) && file.exists(path_botgarden_list)) {
  botgarden <- read_csv(path_botgarden_list, show_col_types = FALSE) %>%
    clean_names() %>%
    mutate(plant = clean_taxon(plant)) %>%
    filter(!is.na(plant), plant != "") %>%
    distinct(plant)
  
  bg_overlap <- plant_weights %>%
    left_join(botgarden %>% mutate(in_botgarden = TRUE), by = "plant") %>%
    mutate(in_botgarden = ifelse(is.na(in_botgarden), FALSE, in_botgarden)) %>%
    group_by(period, ugs_id) %>%
    summarise(
      botgarden_overlap_node = mean(in_botgarden, na.rm = TRUE),
      botgarden_overlap_weighted = safe_divide(sum(plant_w[in_botgarden], na.rm = TRUE), sum(plant_w, na.rm = TRUE)),
      .groups = "drop"
    )
  
  check_unique_keys(bg_overlap, c("period", "ugs_id"), "bg_overlap")
  
  network_metrics_site_period <- network_metrics_site_period %>%
    left_join(bg_overlap, by = c("period", "ugs_id"))
} else {
  network_metrics_site_period <- network_metrics_site_period %>%
    mutate(
      botgarden_overlap_node = NA_real_,
      botgarden_overlap_weighted = NA_real_
    )
}

check_unique_keys(network_metrics_site_period, c("period", "ugs_id"), "network_metrics_site_period_final")

# ----------------------------
# 11) Scoring (site × period)
# ----------------------------
scored_site_period <- network_metrics_site_period %>%
  mutate(
    # --- Bucket 1
    s_N = rescale01_robust(N_interactions),
    s_P = rescale01_robust(P_plants),
    s_A = rescale01_robust(A_pollinators),
    s_L = rescale01_robust(L_links),
    s_singletons = rescale01_robust(1 - singleton_rate),
    s_species_id = rescale01_robust(share_interactions_species_id),
    
    score_richness = rowMeans(cbind(s_N, s_P, s_A, s_L), na.rm = TRUE),
    score_completeness = rowMeans(cbind(s_singletons, s_species_id), na.rm = TRUE),
    
    # --- Bucket 2
    s_connectance = rescale01_robust(connectance),
    s_evenness = rescale01_robust(link_evenness),
    s_H2 = rescale01_robust(H2_specialization),
    score_structure = rowMeans(cbind(s_connectance, s_evenness, s_H2), na.rm = TRUE),
    
    # --- Bucket 3a: trait COVERAGE (data completeness for traits)
    s_plant_traits_cov = rescale01_robust(plant_trait_weighted_cov),
    s_bee_traits_cov   = rescale01_robust(bee_trait_weighted_cov),
    score_trait_coverage = rowMeans(cbind(s_plant_traits_cov, s_bee_traits_cov), na.rm = TRUE),
    
    # --- Bucket 3b: trait BREADTH (representativeness proxy)
    s_color_rich = rescale01_robust(flower_color_richness),
    s_symm_rich  = rescale01_robust(symmetry_richness),
    s_syn_rich   = rescale01_robust(pollination_syndrome_richness),
    s_height_rng = rescale01_robust(plant_height_range),
    s_flow_rng   = rescale01_robust(flowering_start_range),
    
    score_plant_trait_breadth = rowMeans(
      cbind(s_color_rich, s_symm_rich, s_syn_rich, s_height_rng, s_flow_rng),
      na.rm = TRUE
    ),
    
    s_itd_rng    = rescale01_robust(bee_itd_range),
    s_tong_rng   = rescale01_robust(bee_tongue_range),
    s_lecty_rich = rescale01_robust(lecty_richness),
    s_nest_rich  = rescale01_robust(nesting_richness),
    s_soc_rich   = rescale01_robust(sociality_richness),
    
    score_bee_trait_breadth = rowMeans(
      cbind(s_itd_rng, s_tong_rng, s_lecty_rich, s_nest_rich, s_soc_rich),
      na.rm = TRUE
    ),
    
    score_trait_breadth = rowMeans(cbind(score_plant_trait_breadth, score_bee_trait_breadth), na.rm = TRUE),
    
    # --- Overall score (now includes representativeness proxy)
    score_overall_noBG = 0.40 * score_richness +
      0.20 * score_completeness +
      0.15 * score_structure +
      0.15 * score_trait_coverage +
      0.10 * score_trait_breadth
  ) %>%
  arrange(desc(score_overall_noBG))


check_unique_keys(scored_site_period, c("period", "ugs_id"), "scored_site_period")

# ----------------------------
# 12) Pooled networks per site (A+B+C combined)
# ----------------------------
edges_site <- poll %>%
  group_by(ugs_id, plant, pollinator) %>%
  summarise(w = sum(amount, na.rm = TRUE), .groups = "drop") %>%
  filter(w > 0)

check_unique_keys(edges_site, c("ugs_id", "plant", "pollinator"), "edges_site")

network_metrics_site <- edges_site %>%
  group_by(ugs_id) %>%
  group_nest() %>%
  mutate(
    metrics = map2(data, ugs_id, ~ compute_network_metrics(.x, period_val = "ABC", ugs_id_val = .y) %>%
                     select(-ugs_id))  # ugs_id already present outside
  ) %>%
  select(-data) %>%
  unnest(metrics) %>%
  left_join(ugs, by = "ugs_id")

check_unique_keys(network_metrics_site, c("ugs_id"), "network_metrics_site")

scored_site <- network_metrics_site %>%
  mutate(
    s_N = rescale01_robust(N_interactions),
    s_P = rescale01_robust(P_plants),
    s_A = rescale01_robust(A_pollinators),
    s_L = rescale01_robust(L_links),
    s_singletons = rescale01_robust(1 - singleton_rate),
    
    score_richness = rowMeans(cbind(s_N, s_P, s_A, s_L), na.rm = TRUE),
    score_completeness = s_singletons,
    
    s_connectance = rescale01_robust(connectance),
    s_evenness = rescale01_robust(link_evenness),
    s_H2 = rescale01_robust(H2_specialization),
    score_structure = rowMeans(cbind(s_connectance, s_evenness, s_H2), na.rm = TRUE),
    
    score_overall_site = 0.55 * score_richness +
      0.20 * score_completeness +
      0.25 * score_structure
  ) %>%
  arrange(desc(score_overall_site))

check_unique_keys(scored_site, c("ugs_id"), "scored_site")

# ----------------------------
# 13) Clean overview table for scored_site_period (shareable)
# ----------------------------
overview_site_period <- scored_site_period %>%
  distinct(network_id, .keep_all = TRUE) %>%  # safety
  mutate(
    plants_interaction_share = safe_divide(P_plants, plants_recorded),
    rank_overall = dense_rank(desc(score_overall_noBG))
  ) %>%
  select(
    rank_overall,
    network_id, period, ugs_id, ugs_type, site_name,
    plants_recorded, P_plants, plants_interaction_share,
    N_interactions, A_pollinators, L_links,
    singleton_rate, share_interactions_species_id, share_bee_interactions,
    plant_trait_weighted_cov, bee_trait_weighted_cov,
    score_richness, score_completeness, score_trait_coverage, score_trait_breadth,
    score_overall_noBG
  ) %>%
  mutate(
    across(c(singleton_rate, share_interactions_species_id, share_bee_interactions,
             plant_trait_weighted_cov, bee_trait_weighted_cov,
             plants_interaction_share,
             score_richness, score_completeness, score_trait_coverage, score_trait_breadth,
             score_overall_noBG), ~ round(.x, 3))
  ) %>%
  arrange(rank_overall)

write_csv(overview_site_period, here("outputs", "overview_site_period_key_metrics.csv"))


# ----------------------------
# 13) Site overview table (simple, shareable)
# ----------------------------

plants_recorded_site <- flowers %>%
  group_by(ugs_id) %>%
  summarise(
    plants_recorded_total = n_distinct(plant),
    availability_sum_total = sum(availability, na.rm = TRUE),
    .groups = "drop"
  )

poll_quality_site <- poll %>%
  group_by(ugs_id) %>%
  summarise(
    share_interactions_species_id = safe_divide(sum(amount[tax_resolution == "species"], na.rm = TRUE), sum(amount, na.rm = TRUE)),
    share_bee_interactions = safe_divide(sum(amount[is_bee], na.rm = TRUE), sum(amount, na.rm = TRUE)),
    .groups = "drop"
  )

site_overview <- scored_site %>%
  left_join(plants_recorded_site, by = "ugs_id") %>%
  left_join(poll_quality_site, by = "ugs_id") %>%
  transmute(
    ugs_id, ugs_type, site_name,
    plants_recorded_total,
    P_plants_interacting = P_plants,
    A_pollinators,
    N_interactions,
    L_links,
    singleton_rate,
    share_interactions_species_id,
    share_bee_interactions,
    score_overall_site
  ) %>%
  arrange(desc(score_overall_site))

check_unique_keys(site_overview, c("ugs_id"), "site_overview")

# ----------------------------
# 14) Plots (saved to outputs/)
# ----------------------------
p1 <- ggplot(site_overview, aes(P_plants_interacting, N_interactions)) +
  geom_point(aes(shape = ugs_type)) +
  labs(x = "Interacting plant richness (pooled)", y = "Total interactions (pooled)")

ggsave(here("outputs", "plot_sites_plants_vs_interactions.png"), p1, width = 7, height = 5, dpi = 300)

p2 <- ggplot(site_overview, aes(ugs_type, score_overall_site)) +
  geom_boxplot() +
  geom_jitter(width = 0.15, height = 0, alpha = 0.6) +
  labs(x = "UGS type", y = "Site score (pooled)")

ggsave(here("outputs", "plot_score_by_ugs_type.png"), p2, width = 7, height = 5, dpi = 300)

p3 <- scored_site_period %>%
  mutate(ugs_id_f = fct_reorder(factor(ugs_id), score_overall_noBG, .fun = max, na.rm = TRUE)) %>%
  ggplot(aes(period, ugs_id_f, fill = score_overall_noBG)) +
  geom_tile() +
  labs(x = "Period", y = "Site (ugs_id)", fill = "Score (no BG)")

ggsave(here("outputs", "plot_heatmap_site_period_scores.png"), p3, width = 7, height = 10, dpi = 300)

# ----------------------------
# 15) Save outputs
# ----------------------------
write_csv(network_metrics_site_period, here("outputs", "network_metrics_site_period_raw.csv"))
write_csv(scored_site_period, here("outputs", "network_metrics_site_period_scored.csv"))
write_csv(network_metrics_site, here("outputs", "network_metrics_site_pooled.csv"))
write_csv(scored_site, here("outputs", "network_metrics_site_scored.csv"))
write_csv(site_overview, here("outputs", "site_overview_key_metrics.csv"))

write_csv(edges_site_period, here("data", "processed", "edges_site_period.csv"))

# ----------------------------
# 16) Quick console summary
# ----------------------------
message("Saved:")
message(" - outputs/network_metrics_site_period_raw.csv")
message(" - outputs/network_metrics_site_period_scored.csv")
message(" - outputs/network_metrics_site_pooled.csv")
message(" - outputs/network_metrics_site_scored.csv")
message(" - outputs/site_overview_key_metrics.csv")
message(" - outputs/plot_sites_plants_vs_interactions.png")
message(" - outputs/plot_score_by_ugs_type.png")
message(" - outputs/plot_heatmap_site_period_scores.png")

message("\nTop 10 networks by overall score (site × period):")
print(scored_site_period %>%
        select(network_id, ugs_type, N_interactions, P_plants, A_pollinators, L_links,
               share_interactions_species_id, plant_trait_weighted_cov, bee_trait_weighted_cov,
               H2_specialization) %>%
        slice(1:10))

message("\nTop 10 sites by pooled score:")
print(site_overview %>%
        select(ugs_id, ugs_type, plants_recorded_total, P_plants_interacting,
               A_pollinators, N_interactions, score_overall_site) %>%
        slice(1:10))

# End.
