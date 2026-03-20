library(tidyverse)
library(readxl)
library(readr)
library(janitor)
library(stringr)
library(here)

# ----------------------------
# Paths
# ----------------------------
path_flowers <- here("data", "raw", "pappus_data_2023", "df_flowers.xlsx")
path_pollinators <- here("data", "raw", "pappus_data_2023", "df_pollinator.xlsx")
path_botgarden <- here("data", "raw", "Bot_Garden_Avaialble.csv")
path_traits <- here("data", "raw", "trait_data", "floral_traits", "20200811_Plant_traits_database_EXTENDED.csv")
path_manual_traits <- here("data", "raw", "trait_data", "floral_traits", "trait_candidates_for_5_missing_species.csv")

out_dir <- here("outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_site_period <- file.path(out_dir, "plant_site_period_table.csv")
out_ranked <- file.path(out_dir, "plant_feasible_ranked.csv")
out_grouped <- file.path(out_dir, "plant_candidates_grouped.csv")
out_selected <- file.path(out_dir, "plant_candidates_selected.csv")
out_selection_overview <- file.path(out_dir, "plant_selection_group_overview.csv")

# ----------------------------
# Thresholds: hard feasibility filter
# ----------------------------
# minimum thresholds
MIN_SITE_PERIODS_PRESENT <- 5
MIN_SITES_PRESENT <- 2
MIN_PERIODS_PRESENT <- 2
MIN_SITE_PERIODS_WITH_BEE <- 2
MIN_SITE_PERIODS_WITHOUT_BEE <- 2
MIN_TOTAL_BEE_VISITS <- 5

# preferred thresholds
PREF_SITE_PERIODS_PRESENT <- 10
PREF_SITES_PRESENT <- 5
PREF_PERIODS_PRESENT <- 3
PREF_SITE_PERIODS_WITH_BEE <- 5
PREF_SITE_PERIODS_WITHOUT_BEE <- 5
PREF_TOTAL_BEE_VISITS <- 10

# ----------------------------
# Helpers
# ----------------------------
clean_taxon <- function(x) {
  x %>%
    as.character() %>%
    na_if("NA") %>%
    str_replace_all("[\u00A0]", " ") %>%
    str_squish() %>%
    na_if("")
}

make_binomial <- function(genus, species) {
  out <- paste(clean_taxon(genus), clean_taxon(species))
  out <- str_squish(out)
  out[out %in% c("", "NA NA", "NA")] <- NA_character_
  out
}

lookup_named_count <- function(x, counts, default = 0) {
  out <- rep(default, length(x))
  if (length(counts) == 0) return(out)

  hits <- !is.na(x) & x %in% names(counts)
  out[hits] <- unname(counts[x[hits]])
  out
}

period_penalty <- function(flag, label, counts) {
  ifelse(flag, unname(counts[label]), max(counts) + 1)
}

collapse_unique <- function(x, sep = "; ") {
  x <- unique(x[!is.na(x) & x != ""])
  if (length(x) == 0) return(NA_character_)
  paste(sort(x), collapse = sep)
}

select_balanced_within_group <- function(df_group, n_target) {
  remaining <- df_group %>%
    arrange(feasibility_tier, feasibility_rank)

  chosen <- tibble()
  period_counts <- c(A = 0, B = 0, C = 0)
  family_counts <- c()

  while (nrow(chosen) < n_target && nrow(remaining) > 0) {
    scored <- remaining %>%
      mutate(
        tier_score = case_when(
          feasibility_tier == "preferred" ~ 0,
          feasibility_tier == "good" ~ 1,
          TRUE ~ 2
        ),
        family_use = lookup_named_count(family, family_counts, default = 0),
        period_score =
          period_penalty(period_a, "A", period_counts) +
          period_penalty(period_b, "B", period_counts) +
          period_penalty(period_c, "C", period_counts)
      ) %>%
      arrange(
        tier_score,
        family_use,
        period_score,
        feasibility_rank
      )

    pick <- scored %>% slice(1)
    chosen <- bind_rows(chosen, pick)

    if (isTRUE(pick$period_a[[1]])) period_counts["A"] <- period_counts["A"] + 1
    if (isTRUE(pick$period_b[[1]])) period_counts["B"] <- period_counts["B"] + 1
    if (isTRUE(pick$period_c[[1]])) period_counts["C"] <- period_counts["C"] + 1

    fam <- pick$family[[1]]
    if (!is.na(fam)) {
      if (fam %in% names(family_counts)) {
        family_counts[fam] <- family_counts[fam] + 1
      } else {
        family_counts[fam] <- 1
      }
    }

    remaining <- remaining %>%
      filter(plant != pick$plant[[1]])
  }

  chosen
}

# ----------------------------
# 1) Flowers -> occupied plant x site x period units
#    Presence = abundance > 0
# ----------------------------
flowers <- read_excel(path_flowers, col_types = "text") %>%
  clean_names() %>%
  transmute(
    plant = clean_taxon(scientific_name),
    id = clean_taxon(id),
    period = str_extract(id, "^[ABC]"),
    ugs_id = as.integer(str_extract(id, "(?<=_)\\d+$")),
    floral_abundance = suppressWarnings(as.numeric(unnamed_2))
  ) %>%
  filter(!is.na(plant), !is.na(period), !is.na(ugs_id))

occupied_units <- flowers %>%
  group_by(plant, ugs_id, period) %>%
  summarise(
    floral_abundance = ifelse(all(is.na(floral_abundance)), NA_real_, sum(floral_abundance, na.rm = TRUE)),
    plant_present = any(replace_na(floral_abundance, 0) > 0),
    .groups = "drop"
  ) %>%
  filter(plant_present)

period_summary <- flowers %>%
  filter(period %in% c("A", "B", "C")) %>%
  group_by(plant, period) %>%
  summarise(period_present = any(replace_na(floral_abundance, 0) > 0), .groups = "drop") %>%
  filter(period_present) %>%
  group_by(plant) %>%
  summarise(
    periods_present = paste(sort(unique(period)), collapse = "/"),
    period_a = any(period == "A"),
    period_b = any(period == "B"),
    period_c = any(period == "C"),
    n_periods_from_flowers = n_distinct(period),
    .groups = "drop"
  )

# ----------------------------
# 2) Pollinators -> bees only, then plant x site x period bee visits
# ----------------------------
bee_families <- c("Apidae", "Andrenidae", "Halictidae", "Megachilidae", "Colletidae", "Melittidae")

poll_bee <- read_excel(path_pollinators, col_types = "text") %>%
  clean_names() %>%
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
    bee_visits = suppressWarnings(as.numeric(amount)),
    pollinator_group = clean_taxon(pollinator_group),
    family = clean_taxon(family)
  ) %>%
  mutate(
    bee_visits = replace_na(bee_visits, 0),
    is_bee = family %in% bee_families |
      str_detect(str_to_lower(pollinator_group), "wildbee|bombus|apis")
  ) %>%
  filter(is_bee, !is.na(plant), !is.na(period), !is.na(ugs_id))

bee_units <- poll_bee %>%
  group_by(plant, ugs_id, period) %>%
  summarise(
    bee_visits = sum(bee_visits, na.rm = TRUE),
    .groups = "drop"
  )

# ----------------------------
# 3) Core analysis table: occupied units + bee visits
# ----------------------------
plant_site_period <- occupied_units %>%
  left_join(bee_units, by = c("plant", "ugs_id", "period")) %>%
  mutate(
    bee_visits = replace_na(bee_visits, 0),
    bee_present = bee_visits > 0
  )

# ----------------------------
# 4) Botanical Garden lookup
# ----------------------------
in_botgarden <- tibble(plant = character(), in_botgarden = logical())

botgarden_names <- read_csv(path_botgarden, show_col_types = FALSE) %>%
  clean_names() %>%
  bind_rows(
    transmute(., plant = make_binomial(genus_lookup, species_lookup)),
    transmute(., plant = make_binomial(synonym_genus_lookup, synonym_species_lookup)),
    transmute(., plant = str_squish(word(clean_taxon(taxon_name_without_author), 1, 2)))
  ) %>%
  filter(!is.na(plant)) %>%
  distinct(plant)

in_botgarden <- botgarden_names %>%
  mutate(in_botgarden = TRUE)

# ----------------------------
# 5) Plant-level feasibility summary
# ----------------------------
max_unit_share <- plant_site_period %>%
  group_by(plant) %>%
  summarise(
    max_site_period_visits = max(bee_visits, na.rm = TRUE),
    .groups = "drop"
  )

plant_feasibility <- plant_site_period %>%
  group_by(plant) %>%
  summarise(
    n_site_periods_present = n(),
    n_sites_present = n_distinct(ugs_id),
    n_periods_present = n_distinct(period),
    n_site_periods_with_bee = sum(bee_present),
    n_site_periods_without_bee = sum(!bee_present),
    total_bee_visits = sum(bee_visits, na.rm = TRUE),
    mean_floral_abundance = mean(floral_abundance, na.rm = TRUE),
    sd_floral_abundance = sd(floral_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(max_unit_share, by = "plant") %>%
  mutate(
    balance_score = pmin(n_site_periods_with_bee, n_site_periods_without_bee),
    max_site_period_share = ifelse(total_bee_visits > 0, max_site_period_visits / total_bee_visits, NA_real_),
    abundance_cv = ifelse(mean_floral_abundance > 0, sd_floral_abundance / mean_floral_abundance, NA_real_)
  ) %>%
  left_join(in_botgarden, by = "plant") %>%
  mutate(in_botgarden = replace_na(in_botgarden, FALSE))

# ----------------------------
# 6) Hard filter + ranked feasible pool
# ----------------------------
ranked_plants <- plant_feasibility %>%
  filter(
    n_site_periods_present >= MIN_SITE_PERIODS_PRESENT,
    n_sites_present >= MIN_SITES_PRESENT,
    n_periods_present >= MIN_PERIODS_PRESENT,
    n_site_periods_with_bee >= MIN_SITE_PERIODS_WITH_BEE,
    n_site_periods_without_bee >= MIN_SITE_PERIODS_WITHOUT_BEE,
    total_bee_visits >= MIN_TOTAL_BEE_VISITS
  ) %>%
  arrange(
    desc(balance_score),
    desc(n_sites_present),
    desc(n_periods_present),
    desc(n_site_periods_present),
    desc(total_bee_visits),
    max_site_period_share
  ) %>%
  mutate(
    feasibility_rank = row_number(),
    preferred_tier = (
      n_site_periods_present >= PREF_SITE_PERIODS_PRESENT &
        n_sites_present >= PREF_SITES_PRESENT &
        n_periods_present >= PREF_PERIODS_PRESENT &
        n_site_periods_with_bee >= PREF_SITE_PERIODS_WITH_BEE &
        n_site_periods_without_bee >= PREF_SITE_PERIODS_WITHOUT_BEE &
        total_bee_visits >= PREF_TOTAL_BEE_VISITS
    ),
    flag_low_balance = balance_score <= 2,
    flag_low_periods = n_periods_present == 2,
    flag_high_concentration = max_site_period_share > 0.7,
    flag_low_sites = n_sites_present == 2,
    fragility_score = flag_low_balance + flag_low_periods + flag_high_concentration + flag_low_sites,
    feasibility_tier = case_when(
      preferred_tier ~ "preferred",
      fragility_score >= 2 ~ "reserve",
      TRUE ~ "good"
    )
  ) %>%
  relocate(feasibility_tier, .after = feasibility_rank)

# ----------------------------
# 7) Trait database + manual supplement
# ----------------------------
trait_db <- read_delim(
  path_traits,
  delim = ";",
  locale = locale(encoding = "Latin1"),
  show_col_types = FALSE
) %>%
  clean_names() %>%
  transmute(
    plant = clean_taxon(taxon_name),
    family = clean_taxon(family),
    structural_blossom_class = clean_taxon(structural_blossom_class),
    tube_or_spur_present = suppressWarnings(as.numeric(tube_or_spur_present)),
    displaying_pollen = clean_taxon(displaying_pollen),
    pollinator_behaviour = clean_taxon(pollinator_behaviour),
    tube_spur_length_mm = suppressWarnings(as.numeric(tube_spur_length_mm)),
    corolla_or_substitute_length_mm = suppressWarnings(as.numeric(corolla_or_substitute_length_mm)),
    alighting = clean_taxon(alighting),
    trait_source_dataset = "extended_db"
  ) %>%
  filter(!is.na(plant))

manual_traits <- read_csv(
  path_manual_traits,
  locale = locale(encoding = "Latin1"),
  show_col_types = FALSE
) %>%
  clean_names() %>%
  transmute(
    plant = clean_taxon(taxon_name),
    family = clean_taxon(family),
    structural_blossom_class = clean_taxon(structural_blossom_class),
    tube_or_spur_present = suppressWarnings(as.numeric(tube_or_spur_present)),
    displaying_pollen = clean_taxon(displaying_pollen),
    pollinator_behaviour = clean_taxon(pollinator_behaviour),
    tube_spur_length_mm = suppressWarnings(as.numeric(tube_spur_length_mm)),
    corolla_or_substitute_length_mm = suppressWarnings(as.numeric(corolla_or_substitute_length_mm)),
    alighting = clean_taxon(alighting),
    trait_source_dataset = "manual_fill"
  ) %>%
  filter(!is.na(plant))

all_traits <- bind_rows(trait_db, manual_traits) %>%
  group_by(plant) %>%
  slice_tail(n = 1) %>%
  ungroup()

trait_lookup <- tribble(
  ~plant, ~trait_lookup_name, ~trait_match_type, ~trait_match_note,
  "Valeriana rubra", "Centranthus ruber", "synonym", "Trait row taken from accepted synonym / taxonomic equivalent.",
  "Achillea millefolium", "Achillea millefolium aggr.", "aggregate", "Trait row taken from aggregate taxon.",
  "Leucanthemum vulgare", "Leucanthemum vulgare aggr.", "aggregate", "Trait row taken from aggregate taxon.",
  "Taraxacum officinale", "Taraxacum officinale aggr.", "aggregate", "Trait row taken from aggregate taxon."
)

# ----------------------------
# 8) Join ranked pool with period and trait data
# ----------------------------
candidate_traits <- ranked_plants %>%
  left_join(period_summary, by = "plant") %>%
  left_join(trait_lookup, by = "plant") %>%
  mutate(
    trait_lookup_name = coalesce(trait_lookup_name, plant),
    trait_match_type = coalesce(trait_match_type, "exact"),
    trait_match_note = coalesce(trait_match_note, "Exact species-level trait match."),
    periods_present = replace_na(periods_present, ""),
    period_a = replace_na(period_a, FALSE),
    period_b = replace_na(period_b, FALSE),
    period_c = replace_na(period_c, FALSE),
    n_periods_from_flowers = replace_na(n_periods_from_flowers, 0L)
  ) %>%
  left_join(all_traits, by = c("trait_lookup_name" = "plant"))

# ----------------------------
# 9) Primary grouping
#    Primary groups:
#    - open_access
#    - restrictive_access
#
#    Restrictive subtypes:
#    - concealed_nectar
#    - restricted_pollen_handling
#    - both_restrictive
# ----------------------------
ranked_grouped <- candidate_traits %>%
  mutate(
    structural_blossom_class = str_to_lower(structural_blossom_class),
    displaying_pollen = str_to_lower(displaying_pollen),
    pollinator_behaviour = str_to_lower(pollinator_behaviour),
    alighting = str_to_lower(alighting),

    flag_restricted_pollen_handling = (
      displaying_pollen == "hidden" |
        pollinator_behaviour %in% c("alighting_forcing", "buzz_pollination") |
        structural_blossom_class %in% c("flag", "closed_blossom")
    ),

    flag_concealed_nectar = (
      structural_blossom_class %in% c("tube", "gullet") |
        (tube_or_spur_present == 1 & !is.na(tube_spur_length_mm) & tube_spur_length_mm >= 5)
    ),

    flag_open_override = (
      structural_blossom_class %in% c("dish_bowl", "brush") &
        displaying_pollen == "open"
    ),

    restrictive_subtype = case_when(
      flag_open_override & !flag_restricted_pollen_handling ~ NA_character_,
      flag_restricted_pollen_handling & flag_concealed_nectar ~ "both_restrictive",
      flag_restricted_pollen_handling ~ "restricted_pollen_handling",
      flag_concealed_nectar ~ "concealed_nectar",
      TRUE ~ NA_character_
    ),

    primary_group = if_else(is.na(restrictive_subtype), "open_access", "restrictive_access"),

    decision_rule = case_when(
      restrictive_subtype == "both_restrictive" ~ "restrictive_both",
      restrictive_subtype == "restricted_pollen_handling" ~ "restrictive_pollen_handling",
      restrictive_subtype == "concealed_nectar" ~ "restrictive_concealed_nectar",
      flag_open_override ~ "open_override_blossom_plus_open_pollen",
      displaying_pollen == "open" ~ "open_pollen",
      structural_blossom_class %in% c("dish_bowl", "brush", "bell_trumpet", "stielteller") ~ "open_blossom_class",
      TRUE ~ "default_open"
    ),

    primary_group = factor(primary_group, levels = c("open_access", "restrictive_access")),
    restrictive_subtype = factor(
      restrictive_subtype,
      levels = c("concealed_nectar", "restricted_pollen_handling", "both_restrictive")
    ),
    feasibility_tier = factor(feasibility_tier, levels = c("preferred", "good", "reserve"))
  ) %>%
  arrange(primary_group, restrictive_subtype, feasibility_tier, feasibility_rank)

plant_candidates_grouped <- ranked_grouped %>%
  select(
    feasibility_rank,
    feasibility_tier,
    preferred_tier,
    plant,
    family,
    primary_group,
    restrictive_subtype,
    decision_rule,
    periods_present,
    period_a,
    period_b,
    period_c,
    n_periods_from_flowers,
    n_site_periods_present,
    n_site_periods_with_bee,
    n_site_periods_without_bee,
    n_sites_present,
    n_periods_present,
    total_bee_visits,
    balance_score,
    max_site_period_share,
    in_botgarden,
    structural_blossom_class,
    displaying_pollen,
    tube_or_spur_present,
    tube_spur_length_mm,
    pollinator_behaviour,
    corolla_or_substitute_length_mm,
    alighting,
    trait_lookup_name,
    trait_match_type,
    trait_source_dataset,
    trait_match_note
  )

# ----------------------------
# 10) Concrete species selection
#     - eligible main pool = preferred + good only
#     - main pool size derived as the largest 50/50 balance
#       across the two primary groups
#     - within each group, selection prioritises tier, family spread,
#       period spread, and feasibility rank
# ----------------------------
main_eligible <- plant_candidates_grouped %>%
  mutate(
    family = na_if(family, ""),
    period_a = replace_na(as.logical(period_a), FALSE),
    period_b = replace_na(as.logical(period_b), FALSE),
    period_c = replace_na(as.logical(period_c), FALSE)
  ) %>%
  filter(feasibility_tier %in% c("preferred", "good"))

eligible_counts <- main_eligible %>%
  count(primary_group, name = "n_eligible") %>%
  arrange(primary_group)

derived_n_main_per_group <- min(eligible_counts$n_eligible)
derived_n_main_total <- derived_n_main_per_group * 2

main_selected <- main_eligible %>%
  group_by(primary_group) %>%
  group_modify(~ select_balanced_within_group(.x, derived_n_main_per_group)) %>%
  ungroup() %>%
  select(plant)

plant_candidates_selected <- plant_candidates_grouped %>%
  left_join(
    main_selected %>% mutate(selection_pool = "main"),
    by = "plant"
  ) %>%
  mutate(selection_pool = replace_na(selection_pool, "backup")) %>%
  arrange(selection_pool, primary_group, feasibility_tier, feasibility_rank)

plant_selection_group_overview <- plant_candidates_selected %>%
  group_by(primary_group) %>%
  summarise(
    n_species_total = n(),
    n_species_main = sum(selection_pool == "main"),
    n_species_backup = sum(selection_pool == "backup"),
    n_preferred = sum(feasibility_tier == "preferred", na.rm = TRUE),
    n_good = sum(feasibility_tier == "good", na.rm = TRUE),
    n_reserve = sum(feasibility_tier == "reserve", na.rm = TRUE),
    n_period_a = sum(period_a, na.rm = TRUE),
    n_period_b = sum(period_b, na.rm = TRUE),
    n_period_c = sum(period_c, na.rm = TRUE),
    n_main_period_a = sum(selection_pool == "main" & period_a, na.rm = TRUE),
    n_main_period_b = sum(selection_pool == "main" & period_b, na.rm = TRUE),
    n_main_period_c = sum(selection_pool == "main" & period_c, na.rm = TRUE),
    families_present = collapse_unique(family),
    .groups = "drop"
  ) %>%
  arrange(primary_group)

# ----------------------------
# 11) Save outputs
# ----------------------------
write_csv(plant_site_period, out_site_period)
write_csv(ranked_plants, out_ranked)
write_csv(plant_candidates_grouped, out_grouped)
write_csv(plant_candidates_selected, out_selected)
write_csv(plant_selection_group_overview, out_selection_overview)
