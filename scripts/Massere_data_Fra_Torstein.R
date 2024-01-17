library(tidyverse)


# Topp 10 -----------------------------------------------------------------



Topp10November_og_oktober <- readxl::read_excel("C:/Users/jonr/Downloads/Topp10November og oktober.xlsx")

okt <- Topp10November_og_oktober %>% 
  tail(n=10)
nov <- Topp10November_og_oktober %>% 
  head(n=10)

okt <- okt %>% 
  select(PANGOLIN_NOM, "Antall" = AntallPrøver, "Prosent" = `Andel prøver(%)`) %>% 
  add_column("MONTH" = "Oktober")

nov <- nov %>% 
  select(PANGOLIN_NOM, "Antall" = AntallPrøver, "Prosent" = `Andel prøver(%)`) %>% 
  add_column("MONTH" = "November")

# Need the variants from the last month only

okt <- okt %>% 
  filter(PANGOLIN_NOM %in% nov$PANGOLIN_NOM)

df <- bind_rows(okt, nov) %>% 
  # Fill in the missing months in the previous month
  complete(PANGOLIN_NOM, MONTH) %>% 
  mutate(Prosent = as.numeric(Prosent),
         Antall = as.numeric(Antall)) %>% 
  mutate(Prosent = round(Prosent * 100, digits = 2))
  

df <- df %>% 
  rename("VARIANT" = PANGOLIN_NOM) %>% 
  pivot_longer(!c(VARIANT, MONTH), names_to = "MALTALL", values_to = "OBS") %>% 
  select(MONTH, MALTALL, VARIANT, OBS) %>% 
  arrange(VARIANT, MONTH, MALTALL) %>% 
  # Replace NA with zero
  mutate(OBS = replace_na(OBS, 0)) %>% 
  # Legg til prikking
  mutate(SPVFLAGG = 0)

write_csv(df, "C:/Users/jonr/Downloads/2023-12-06_variants_per_month.csv")

# MONTH,MALTALL,VARIANT,OBS,SPVFLAGG
# Oktober,Antall,JG3,22,0
# Oktober,Prosent,JG3,17.5,0
# September,Antall,JG3,7,0
# September,Prosent,JG3,4.2,0
# Oktober,Antall,HK3,13,0
# Oktober,Prosent,HK3,10.3,0
# September,Antall,HK3,19,0
# September,Prosent,HK3,11.5,0
# Oktober,Antall,EG.5.1.1,11,0
# Oktober,Prosent,EG.5.1.1,8.7,0
# September,Antall,EG.5.1.1,19,0
# September,Prosent,EG.5.1.1,11.5,0
# Oktober,Antall,EG.5.1.3,7,0
# Oktober,Prosent,EG.5.1.3,5.6,0
# September,Antall,EG.5.1.3,10,0
# September,Prosent,EG.5.1.3,6.1,0
# Oktober,Antall,EG.5.1,6,0
# Oktober,Prosent,EG.5.1,4.8,0
# September,Antall,EG.5.1,17,0
# September,Prosent,EG.5.1,10.3,0
# Oktober,Antall,FL.1.5.1,6,0
# Oktober,Prosent,FL.1.5.1,4.8,0
# September,Antall,FL.1.5.1,3,0
# September,Prosent,FL.1.5.1,1.8,0
# Oktober,Antall,HV.1,6,0
# Oktober,Prosent,HV.1,4.8,0
# September,Antall,HV.1,6,0
# September,Prosent,HV.1,3.6,0
# Oktober,Antall,JD.1.1,6,0
# Oktober,Prosent,JD.1.1,4.8,0
# September,Antall,JD.1.1,1,0
# September,Prosent,JD.1.1,0.6,0
# Oktober,Antall,XBB.1.16.11,5,0
# Oktober,Prosent,XBB.1.16.11,4.0,0
# September,Antall,XBB.1.16.11,3,0
# September,Prosent,XBB.1.16.11,1.8,0
# Oktober,Antall,GK.1,3,0
# Oktober,Prosent,GK.1,2.4,0
# September,Antall,GK.1,1,0
# September,Prosent,GK.1,0.6,0


# VOI ---------------------------------------------------------------------


# Read lineage descriptions from GitHub
pango <- read_delim(file = "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/lineage_notes.txt")

# 2023.11.08: Include BA.2.86 abbreviations
pango_str <- pango %>% 
  # Get the BA.2.75's
  filter(str_detect(Description, "B.1.1.529.2.75") | str_detect(Lineage, "^BA.2.75")) %>% 
  # Remove any withdrawn lineages
  filter(str_detect(Lineage, "\\*", negate = TRUE)) %>% 
  # Pull all the aliases into a character vector
  pull(Lineage)

pango_str <- pango %>% 
  # Get the BA.2.86's
  filter(str_detect(Description, "B.1.1.529.2.86") | str_detect(Lineage, "^BA.2.86")) %>% 
  # Remove any withdrawn lineages
  filter(str_detect(Lineage, "\\*", negate = TRUE)) %>% 
  # Pull all the aliases into a character vector
  pull(Lineage)

BN %>% 
  filter(PANGOLIN_NOM %in% pango_str)
