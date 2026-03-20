library(tidyverse)
library(readxl)
library(janitor)
library(stringr)

# paths (adjust if needed)
path_overview <- "outputs/overview_site_period_key_metrics.csv"
path_flowers  <- "data/raw/pappus_data_2023/df_flowers.xlsx"

TOP_N <- 20

# 1) Read overview table and pick top N networks
overview <- read_csv(path_overview, show_col_types = FALSE)

if ("rank_overall" %in% names(overview)) {
  top_keys <- overview %>%
    arrange(rank_overall) %>%
    slice(1:TOP_N) %>%
    transmute(
      period = as.character(period),
      ugs_id = as.integer(ugs_id)
    ) %>%
    distinct()
} else if ("score_overall" %in% names(overview)) {
  top_keys <- overview %>%
    arrange(desc(score_overall)) %>%
    slice(1:TOP_N) %>%
    transmute(
      period = as.character(period),
      ugs_id = as.integer(ugs_id)
    ) %>%
    distinct()
} else {
  stop("overview table must contain either 'rank_overall' or 'score_overall'.")
}

# 2) Read raw df_flowers and extract period + ugs_id from id (e.g., 'B_33')
flowers_raw <- read_excel(path_flowers, col_types = "text") %>% clean_names()

flowers_clean <- flowers_raw %>%
  transmute(
    plant  = str_squish(as.character(scientific_name)),
    id     = str_squish(as.character(id)),
    period = str_extract(id, "^[ABC]"),
    ugs_id = as.integer(str_extract(id, "(?<=_)\\d+$"))
  ) %>%
  filter(!is.na(plant), plant != "", !is.na(period), !is.na(ugs_id))

# 3) Plants recorded in the top N networks (full plant list per site×period)
plants_topN <- flowers_clean %>%
  semi_join(top_keys, by = c("period", "ugs_id")) %>%
  distinct(plant) %>%
  arrange(plant)

# 4) Total unique species count
n_unique_plants_topN <- nrow(plants_topN)
cat("Top", TOP_N, "networks contain", n_unique_plants_topN, "unique recorded plant taxa.\n")

# Optional: save the plant list
out_file <- paste0("outputs/plants_in_top", TOP_N, "_networks_from_df_flowers.csv")
write_csv(plants_topN, out_file)
cat("Saved plant list to:", out_file, "\n")

