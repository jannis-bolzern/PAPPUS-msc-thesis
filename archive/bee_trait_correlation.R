library(tidyverse)
library(readxl)
library(janitor)
library(stringr)

# -------------------------
# Paths (adjust if needed)
# -------------------------
path_bee_traits <- "data/raw/trait_data/pollinator_traits/Bee_traits_all_CH.xlsx"
path_indiv      <- "data/raw/trait_data/pollinator_traits/BetterGardens/06_trait_data/individual_traits.csv"

# -------------------------
# Helpers
# -------------------------
clean_taxon <- function(x) {
  x %>% as.character() %>% str_replace_all("[\u00A0]", " ") %>% str_squish()
}

# -------------------------
# 1) Load species-level bee traits (Bee_traits_all_CH)
# -------------------------
bee_raw <- read_excel(
  path_bee_traits,
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
    itd_mean_f = suppressWarnings(as.numeric(itd_mean_f)),
    tongue_length_tongue = suppressWarnings(as.numeric(tongue_length_tongue))
  ) %>%
  filter(!is.na(taxon), taxon != "") %>%
  distinct(taxon, .keep_all = TRUE)

# -------------------------
# 2) Load individual traits and aggregate to species means
# -------------------------
indiv <- read_csv(path_indiv, show_col_types = FALSE) %>%
  mutate(
    taxon = clean_taxon(taxon),
    intertegular_distance = na_if(intertegular_distance, 0),
    proboscis_length = na_if(proboscis_length, 0),
    prementum_length = na_if(prementum_length, 0)
  ) %>%
  filter(!is.na(taxon), taxon != "")

indiv_means <- indiv %>%
  group_by(taxon) %>%
  summarise(
    itd_indiv = mean(intertegular_distance, na.rm = TRUE),
    proboscis_indiv = mean(proboscis_length, na.rm = TRUE),
    prementum_indiv = mean(prementum_length, na.rm = TRUE),
    n_individuals = n(),
    .groups = "drop"
  ) %>%
  mutate(
    itd_indiv = ifelse(is.nan(itd_indiv), NA_real_, itd_indiv),
    proboscis_indiv = ifelse(is.nan(proboscis_indiv), NA_real_, proboscis_indiv),
    prementum_indiv = ifelse(is.nan(prementum_indiv), NA_real_, prementum_indiv)
  )

# -------------------------
# 3) Join and compare
# -------------------------
comp <- indiv_means %>%
  inner_join(bee_traits, by = "taxon")

message("Matched species: ", nrow(comp))

# Trait pairs to compare:
# - ITD: species-level itd_mean_f vs individual mean intertegular_distance
# - Tongue/proboscis proxy: species-level tongue_length_tongue vs individual proboscis / prementum (if present)

# Correlation helper
corr_pair <- function(df, x, y) {
  d <- df %>% select(all_of(c(x, y, "n_individuals"))) %>% drop_na()
  if (nrow(d) < 5) return(tibble(x = x, y = y, n = nrow(d), pearson = NA, spearman = NA))
  tibble(
    x = x, y = y,
    n = nrow(d),
    pearson = cor(d[[x]], d[[y]], method = "pearson"),
    spearman = cor(d[[x]], d[[y]], method = "spearman")
  )
}

corr_summary <- bind_rows(
  corr_pair(comp, "itd_mean_f", "itd_indiv"),
  corr_pair(comp, "tongue_length_tongue", "proboscis_indiv"),
  corr_pair(comp, "tongue_length_tongue", "prementum_indiv")
)

print(corr_summary)

# -------------------------
# 4) Quick visual check
# -------------------------
# ITD plot (size by number of individuals)
p_itd <- comp %>%
  ggplot(aes(itd_mean_f, itd_indiv)) +
  geom_point(aes(size = n_individuals), alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Species-level ITD (Bee_traits_all_CH)", y = "Mean ITD (individual_traits)")

ggsave("trait_compare_ITD.png", p_itd, width = 6, height = 4, dpi = 300)

# Tongue/proboscis plot (if available)
p_tongue <- comp %>%
  drop_na(tongue_length_tongue, proboscis_indiv) %>%
  ggplot(aes(tongue_length_tongue, proboscis_indiv)) +
  geom_point(aes(size = n_individuals), alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Species-level tongue length (Bee_traits_all_CH)",
       y = "Mean proboscis length (individual_traits)")

ggsave("trait_compare_TONGUE_vs_PROBOSCIS.png", p_tongue, width = 6, height = 4, dpi = 300)

message("Saved plots: trait_compare_ITD.png, trait_compare_TONGUE_vs_PROBOSCIS.png")
