library(tidyverse)

setwd("//wsl.localhost/Ubuntu-20.04/home/jonr/Prosjekter/NIPH_NGS_vis/data/SC2_samples")
df <- read_csv("2023-11-09_variants_per_month_Torstein.csv")

# Ønsker denne strukturen
# PANGOLIN_NOM,MONTH,ANTALL,PROSENT,ENDRING,SPVFLAGG
df %>% 
  mutate(MONTH = case_when(
    CATEGORY == "PREV_MONTH_ANTALL" ~ "September",
    CATEGORY == "PREV_MONTH_PERC"   ~ "September",
    CATEGORY == "CURR_MONTH_ANTALL" ~ "Oktober",
    CATEGORY == "CURR_MONTH_PERC"   ~ "Oktober",
  )) %>%
  pivot_wider(names_from = CATEGORY, values_from = ANTALL) %>% 
  # Må samle alle "Antall" i samme kolonne ANTALL
  pivot_longer(!c(PANGOLIN_NOM, PREV_MONTH_PERC, CURR_MONTH_PERC, PERC_CHANGE, SPVFLAGG, MONTH), names_to = "ANTALL", values_to = "OBS")
