library(tidyverse)
library(lubridate)

# Read the data
df <- read_csv("data/HCV_NGS.csv")


HCV_genotypes_per_year <- 
df %>% 
  # Create date column
  mutate("Date" = as.Date(SAMPLED_DATE)) %>% 
  # Create year and month column 
  mutate("Year"  = lubridate::year(Date),
         "Month" = lubridate::month(Date)) %>% 
  unite("YearMonth", c(Year, Month), sep = "-", remove = FALSE) %>% 
  # remove completely duplicated rows
  distinct() %>% 
  # Remove samples with no sampling date
  filter(!is.na(Date)) %>%
  # Keep only rows with "Resultat" in NAME
  filter(NAME == "Resultat") %>% 
  # Clean up genotype names
  mutate(ENTRY = case_when(
    str_detect(ENTRY, "IKKE")        ~ "IKKE",
    str_detect(ENTRY, "Ikke")        ~ "IKKE",
    str_detect(ENTRY, "\\.")         ~ "IKKE",
    str_detect(ENTRY, "kansellert")  ~ "IKKE",
    str_detect(ENTRY, "A")           ~ "a",
    str_detect(ENTRY, "B")           ~ "b",
    str_detect(ENTRY, "C")           ~ "c",
    str_detect(ENTRY, "D")           ~ "d",
    str_detect(ENTRY, "H")           ~ "H",
    str_detect(ENTRY, "N")           ~ "n",
    str_detect(ENTRY, "2k1b")        ~ "2k1b",
    .default = ENTRY
  )) %>% 
  # Count the number of each genotype per year
  count(Year, ENTRY, name = "Genotypes_per_year") %>% 
  # Group by year to calculate total samples per year
  group_by(Year) %>%
  mutate(Total_per_year = sum(Genotypes_per_year), # Total samples per year
         Percentage = round(100 * Genotypes_per_year / Total_per_year, digits = 2)) %>% # Percentage of each genotype per year
  ungroup() %>%

  # Reorder and rename columns
  select(YEAR = Year,
         GENOTYPE = ENTRY,
         ANTALL = Genotypes_per_year,
         PERCENT = Percentage) %>% 
  add_column("FLAGG" = "0") %>% 
  # Complete the series for all genotypes per year
  complete(YEAR, GENOTYPE) %>% 
  replace_na(list(ANTALL = 0, PERCENT = 0, FLAGG = "0")) 
  
write_csv(HCV_genotypes_per_year, "HCV_genotypes_per_year.csv")

# Set up data structures

# Create the data structure that will be filled with values. Using strings only to get factors later in crossing()
month_numbers <- as.character(1:12)
month_names   <- lubridate::month(1:12, label = T, abbr = T, locale = "Norwegian_Norway") 
years         <- c("2023")

# Create variables to use in filtering
current_year <- year(Sys.Date())
#current_month <- month(Sys.Date())
# Hardcode current month.
current_month <- 9
current_month_label <- lubridate::month(current_month, label = T, abbr = T, locale = "Norwegian_Norway")
current_yearmonth <- paste(current_year, current_month_label, sep = "-")

if (current_month == 1) {
  previous_month <- 12
  previous_month_label <- lubridate::month(previous_month, label = T, abbr = T, locale = "Norwegian_Norway")
  previous_year <- current_year-1
  previous_yearmonth <- paste(previous_year, previous_month_label, sep = "-")
} else {
  previous_month <- current_month - 1
  previous_month_label <- lubridate::month(previous_month, label = T, abbr = T, locale = "Norwegian_Norway")
  previous_yearmonth <- paste(current_year, previous_month_label, sep = "-")
}


# Create all possible combinations of month and year 
final_data <- crossing(month_names, years) %>% 
  rename(
    "MONTH"  = "month_names",
    "YEAR"   = "years"
  ) %>% 
  # Create YEARMONTH column for combining later
  unite("YEARMONTH", c("YEAR", "MONTH"), sep = "-", remove = T) %>% 
  mutate(YEARMONTH = as.factor(YEARMONTH))


#### OLD 

library(tidyverse)
library(lubridate)


# Get data from BN --------------------------------------------------------


# Source the file from Andreas Rohringer to get the latest data and variants
#source("N:/Virologi/Influensa/ARoh/SARS_COV2/Scripts/BN_COVID-SQLquery.R")

# Establish connection to the BN_Covid19 database on SQL Server
con <- odbc::dbConnect(odbc::odbc(),
                       Driver = "SQL Server",
                       Server = "sql-bn-covid19",
                       Database = "BN_Covid19")

# Query the ENTRYFLD table to extract relevant data fields
entryfld <- tbl(con, "ENTRYFLD") %>%
  filter(FIELDID %in% c(272, 360)) %>%
  pivot_wider(names_from = FIELDID, values_from = CONTENT) %>%
  mutate(CONTENT = paste0(`272`, `360`, sep = ";")) %>%
  select(KEY, CONTENT)  %>%
  collect()

# Query the ENTRYFLD table to extract NSP5 mutation data
NSP5mut <- tbl(con, "ENTRYFLD") %>%
  filter(FIELDID %in% "388") %>%
  pivot_wider(names_from = FIELDID, values_from = CONTENT) %>%
  rename(NSP5mut = `388`) %>%
  collect()

ORF1A <- tbl(con, "ENTRYFLD") %>%
  filter(FIELDID %in% "270") %>%
  pivot_wider(names_from = FIELDID, values_from = CONTENT) %>%
  rename(ORF1A = `270`) %>%
  collect()


# Query the ENTRYTABLE table to extract relevant sequencing data
entrytable <- tbl(con, "ENTRYTABLE") %>%
  select(
    KEY,
    PROVE_TATT,
    P,
    FYLKENAVN,
    ST,
    MATERIALE,
    PROSENTDEKNING_GENOM,
    DEKNING_NANOPORE,
    SEKV_OPPSETT_NANOPORE,
    SEKV_OPPSETT_SWIFT7,
    SEQUENCEID_NANO29,
    SEQUENCEID_SWIFT,
    COVERAGE_BREADTH_SWIFT,
    GISAID_PLATFORM,
    GENOTYPE_SVART_I_LABWARE,
    COVERAGE_BREATH_EKSTERNE,
    SAMPLE_CATEGORY,
    INNSENDER,
    COVERAGE_DEPTH_SWIFT,
    COVARAGE_DEPTH_NANO,
    S,
    PANGOLIN_NOM,
    FULL_PANGO_LINEAGE,
    AR,
    N,
    K,
    ORF1A,
    ORF1B,
    ORF3A,
    E,
    ORF6,
    ORF7A,
    ORF8,
    ORF9B,
    ORF14,
    ORF10,
    M
  ) %>%
  rename(
    "Dekning_Artic" = PROSENTDEKNING_GENOM,
    "Dekning_Swift" = COVERAGE_BREADTH_SWIFT,
    "Dekning_Eksterne" = COVERAGE_BREATH_EKSTERNE,
    "Dekning_Nano" = DEKNING_NANOPORE,
    "Project" = P,
    "HSP" = ST
  ) %>%
  collect()

# Join the ENTRYFLD and NSP5mut tables with the ENTRYTABLE table to create a single dataset
# Assuming you have a data frame named 'replacement_data' that contains the replacement values
# Replace 'replacement_data' with the actual name of your replacement data frame

allvariants <- entrytable %>%
  left_join(entryfld, by = "KEY") %>%
  mutate_all(list(~ na_if(., ""))) %>%
  left_join(NSP5mut, by = "KEY") %>%
  left_join(ORF1A, by = "KEY") %>%
  mutate(ORF1A = coalesce(ORF1A.x, ORF1A.y)) %>%
  select(-ends_with(".x"), -ends_with(".y"))%>%  # Remove temporary columns
  
  # Now 'ORF1A' column in 'allvariants' contains values from 'replacement_data' if empty
  
  mutate("S_merged" = coalesce(S, CONTENT), "flag" = if_else(is.na(S), 1, 0)) %>%
  #NB! Here are changes from Andreas' code:
  #mutate(week = week(as.Date(PROVE_TATT)), year = lubridate::year(as.Date(PROVE_TATT)), 
  #       wy = tsibble::yearweek(as.Date(PROVE_TATT)),
  #       my = tsibble::yearmonth(as.Date(PROVE_TATT))) %>% 
  
  filter(PROVE_TATT != "")

# Convert all column values to UTF-8 encoding to deal with Norwegian letters
allvariants[] <- lapply(allvariants, function(x)
  iconv(x, from = "ISO-8859-1", to = "UTF-8"))
allvariants <- as.data.frame(allvariants)

# Filter out low quality samples
allvariants_v <- allvariants %>%
  filter((Dekning_Swift >= 70) |
           (Dekning_Artic >= 70) |
           (Dekning_Nano >= 70)  |
           (Dekning_Eksterne >= 70)
  ) %>%
  filter (PANGOLIN_NOM != "#BESTILT#") %>%
  filter (PANGOLIN_NOM != "Inkonklusiv") %>%
  filter (PANGOLIN_NOM != "inkonklusiv") %>%
  filter (PANGOLIN_NOM != "Se kommentar") %>%
  filter (PANGOLIN_NOM != "Seekom") %>%
  filter (PANGOLIN_NOM != "") %>%
  filter (PANGOLIN_NOM != "Failed") %>%
  filter (PANGOLIN_NOM != "failed") %>%
  #filter (PANGOLIN_NOM != "Unassigned")%>%
  filter (PANGOLIN_NOM != "NA")

# Clean up
rm(entrytable)
rm(entryfld)
rm(NSP5mut)
rm(ORF1A)

# Close connection 
closeAllConnections()


######################   Collapsing Pangolins  ##################################
# https://www.who.int/activities/tracking-SARS-CoV-2-variants go here to update the list 
# Add a Collapsed pango column to the dataset based on the long Pangolin lineage 

allvariants_v <- allvariants_v %>%
  mutate(Collapsed_pango = case_when(
    data.table::like(FULL_PANGO_LINEAGE, "B.1.1.529.2.86") ~ "BA.2.86.X",
    data.table::like(FULL_PANGO_LINEAGE, "XBB.1.9.2.5.1") ~ "EG.5.1.X",
    data.table::like(FULL_PANGO_LINEAGE, "BJ.1BM.1.1.1].1.16") ~ "XBB.1.16.X",
    data.table::like(FULL_PANGO_LINEAGE, "BJ.1BM.1.1.1].2.3") ~ "XBB.2.3.X",
    data.table::like(FULL_PANGO_LINEAGE, "XBB.2.3") ~ "XBB.2.3.X",
    data.table::like(FULL_PANGO_LINEAGE, "XBB.1.9.2.5") ~ "EG.5.X",
    data.table::like(FULL_PANGO_LINEAGE, "XBB.1.16") ~ "XBB.1.16.X",
    data.table::like(FULL_PANGO_LINEAGE, "XBB.1.9.1") ~ "XBB.1.9.1.X",
    data.table::like(FULL_PANGO_LINEAGE, "XBB.1.9.2") ~ "XBB.1.9.2.X",
    data.table::like(FULL_PANGO_LINEAGE, "XBB.1.5") ~ "XBB.1.5.X",
    data.table::like(FULL_PANGO_LINEAGE, "BJ.1BM.1.1.1].1.9.1") ~ "XBB.1.9.1.X",
    data.table::like(FULL_PANGO_LINEAGE, "BJ.1BM.1.1.1].1.9.2") ~ "XBB.1.9.2.X",
    data.table::like(FULL_PANGO_LINEAGE, "B.1.1.529.2.75.3.4.1.1.1.1") ~ "CH.1.1.X",
    data.table::like(FULL_PANGO_LINEAGE, "BJ.1BM.1.1.1].1.5") ~ "XBB.1.5.X",
    data.table::like(FULL_PANGO_LINEAGE, "BJ.1BM.1.1.1].1.9") ~ "XBB.1.9.X",
    data.table::like(FULL_PANGO_LINEAGE, "XBB.1.9") ~ "XBB.1.9.X",
    data.table::like(FULL_PANGO_LINEAGE, "XBB.") ~ "XBB.X",
    data.table::like(FULL_PANGO_LINEAGE, "BJ.1BM.1.1.1") ~ "XBB.X",
    grepl("^\\[", FULL_PANGO_LINEAGE) ~ "Andre recombinanter",
    data.table::like(FULL_PANGO_LINEAGE, "X") ~ "Andre recombinanter",
    data.table::like(FULL_PANGO_LINEAGE, "B.1.1.529.5.3.1.1.1.1.1") ~ "BQ.1.X",
    data.table::like(FULL_PANGO_LINEAGE, "B.1.1.529.2.75") ~ "BA.2.75.X",
    data.table::like(FULL_PANGO_LINEAGE, "B.1.1.529.5") ~ "BA.5.X",
    data.table::like(FULL_PANGO_LINEAGE, "B.1.1.529.4") ~ "BA.4.X",
    data.table::like(FULL_PANGO_LINEAGE, "B.1.1.529.3") ~ "BA.3.X",
    data.table::like(FULL_PANGO_LINEAGE, "B.1.1.529.2") ~ "BA.2.X",
    data.table::like(FULL_PANGO_LINEAGE, "B.1.1.529.1") ~ "BA.1.X",
    data.table::like(FULL_PANGO_LINEAGE, "B.1.1.529") ~ "B.1.1.529",
    data.table::like(FULL_PANGO_LINEAGE, "B.1.617") ~ "Delta",
    data.table::like(FULL_PANGO_LINEAGE, "B.1.1.7") ~ "Alpha",
    data.table::like(FULL_PANGO_LINEAGE, "B.1.351") ~ "Beta",
    data.table::like(FULL_PANGO_LINEAGE, "B.1.1.28.1") ~ "Gamma",
    data.table::like(FULL_PANGO_LINEAGE, "B.1.427") | data.table::like(FULL_PANGO_LINEAGE, "B.1.429") ~ "Epsilon",
    data.table::like(FULL_PANGO_LINEAGE, "B.1.1.28.2") ~ "Zeta",
    data.table::like(FULL_PANGO_LINEAGE, "B.1.525") ~ "Eta",
    data.table::like(FULL_PANGO_LINEAGE, "B.1.1.28.3") ~ "Theta",
    data.table::like(FULL_PANGO_LINEAGE, "B.1.526") ~ "Iota",
    data.table::like(FULL_PANGO_LINEAGE, "B.1.617.1") ~ "Kappa",
    data.table::like(FULL_PANGO_LINEAGE, "B.1.1.1.37.1") ~ "Lambda",
    data.table::like(FULL_PANGO_LINEAGE, "B.1.621") ~ "Mu",
    TRUE ~ "Andre SARS CoV 2"
  ))


# Load the CSV file containing VirusVariant and included.sub.lineages
variant_mappings_url <- "https://www.ecdc.europa.eu/sites/default/files/documents/PathogenVariant_public_mappings.csv"
variant_mappings <- read.csv(variant_mappings_url)

# Split the included.sub.lineages into separate variants
variant_mappings$included.sub.lineages <- strsplit(variant_mappings$included.sub.lineages, "\\|")

# Ensure that included.sub.lineages are character vectors
variant_mappings$included.sub.lineages <- lapply(variant_mappings$included.sub.lineages, as.character)

# Modify the find_matched_variant function to return a single value (if any) rather than a vector
find_matched_variant <- function(PANGOLIN_NOM) {
  matched_variants <- variant_mappings$VirusVariant[sapply(variant_mappings$included.sub.lineages, function(variants) PANGOLIN_NOM %in% variants)]
  if (length(matched_variants) > 0) {
    return(matched_variants[1])  # Return the first matched variant
  } else {
    return(NA)
  }
}

# Match included.sub.lineages and create a new column "Tessy" in the dataframe
allvariants_v$Tessy <- sapply(allvariants_v$PANGOLIN_NOM, find_matched_variant)


allvariants_v <- allvariants_v %>%
  mutate(
    Tessy = case_when(
      data.table::like(FULL_PANGO_LINEAGE, "XBF") ~ "BA.2.75",
      data.table::like(FULL_PANGO_LINEAGE, "XBK") ~ "BA.2.75",
      grepl("XBB.1.5-like\\+F456L", Tessy) & grepl("L455F", S_merged) ~ "XBB.1.5-like+L455F+F456L",
      grepl("BA.2.75", Tessy) & grepl("DV.7.1", PANGOLIN_NOM) ~ "DV.7.1",
      grepl("BA.2", Tessy) & grepl("BA.2.86", PANGOLIN_NOM) ~ "BA.2.86",
      grepl("XBB.1.5-like", Tessy) & grepl("[BJ.1BM.1.1.1].1.16", FULL_PANGO_LINEAGE) ~ "XBB.1.16",
      grepl("XBB.1.5-like", Tessy) & grepl("XBB.1.16", FULL_PANGO_LINEAGE) ~ "XBB.1.16",
      TRUE ~ Tessy  # Keep the original value for other cases
    )
  )


# Define custom colors for each Collapsed_pango name
custom_colors <- c(
  "BA.2.86.X" = "#83e4da",
  "XBB.1.16.X" = "#CD6090",
  "XBB.2.3.X" = "#EE799F",
  "EG.5.X" = "#EEA9B8",
  "EG.5.1.X" = "#ff5733",
  "XBB.1.9.X" = "#481B6DFF",
  "XBB.1.9.1.X" = "#9932CC",
  "XBB.1.9.2.X" = "#DDA0DD",
  "XBB.1.5.X" = "#EE2C2C",
  "CH.1.1.X" = "#2F6B8EFF",
  "XBB.X" = "#8B1A1A",
  "BQ.1.X" = "#31B57BFF",
  "BA.2.75.X" = "#FFF68F",
  "BA.5.X" = "#556B2F",
  "BA.4.X" = "#8BD646FF",
  "BA.3.X" = "#C0FF3E",
  "BA.2.X" = "#CDC673",
  "BA.1.X" = "#00EE00",
  "B.1.1.529" = "#26AD81FF",
  "Omicron"= "#00ffc4",
  "Delta" = "#FF8C00",
  "Alpha" = "#97FFFF",
  "Beta" = "#009ACD",
  "Gamma" = "#104E8B",
  "Epsilon" = "#CAE1FF",
  "Andre SARS CoV 2" = "darkslategray",
  "Mu" = "#668B8B",
  "Lambda" = "#836FFF",
  "Kappa" = "#27408B",
  "Zeta" =  "#5D00D9",
  "Theta" = "#0000ff",
  "Kappa"= "#3399ff",
  "Andre recombinanter" = "black"
)

# Top 10 variants ---------------------------------------------------------


# Set up data structures

# Create the data structure that will be filled with values. Using strings only to get factors later in crossing()
month_numbers <- as.character(1:12)
month_names   <- lubridate::month(1:12, label = T, abbr = T, locale = "Norwegian_Norway") 
years         <- c("2023")

# Create variables to use in filtering
current_year <- year(Sys.Date())
#current_month <- month(Sys.Date())
# Hardcode current month.
current_month <- 9
current_month_label <- lubridate::month(current_month, label = T, abbr = T, locale = "Norwegian_Norway")
current_yearmonth <- paste(current_year, current_month_label, sep = "-")

if (current_month == 1) {
  previous_month <- 12
  previous_month_label <- lubridate::month(previous_month, label = T, abbr = T, locale = "Norwegian_Norway")
  previous_year <- current_year-1
  previous_yearmonth <- paste(previous_year, previous_month_label, sep = "-")
} else {
  previous_month <- current_month - 1
  previous_month_label <- lubridate::month(previous_month, label = T, abbr = T, locale = "Norwegian_Norway")
  previous_yearmonth <- paste(current_year, previous_month_label, sep = "-")
}


# Create all possible combinations of month and year 
final_data <- crossing(month_names, years) %>% 
  rename(
    "MONTH"  = "month_names",
    "YEAR"   = "years"
  ) %>% 
  # Create YEARMONTH column for combining later
  unite("YEARMONTH", c("YEAR", "MONTH"), sep = "-", remove = T) %>% 
  mutate(YEARMONTH = as.factor(YEARMONTH))

allvariants_v <- tibble(allvariants_v)

# Filter the data
tmp <- allvariants_v %>% 
  # Fjerne evt positiv controll
  filter(str_detect(KEY, "pos", negate = TRUE)) %>%
  # Fix date format
  mutate("PROVE_TATT" = ymd(PROVE_TATT)) %>% 
  # Keep samples from current year only
  filter(PROVE_TATT >= paste0(current_year, "-01-01")) %>% 
  # select relevant columns
  select(PROVE_TATT, PANGOLIN_NOM, FULL_PANGO_LINEAGE, Collapsed_pango)

# Create year and month column
tmp <- tmp %>% 
  mutate(YEAR = year(PROVE_TATT),
         MONTH = lubridate::month(PROVE_TATT, label = T, abbr = T, locale = "Norwegian_Norway")) %>% 
  unite("YEARMONTH", c(YEAR, MONTH), sep = "-", remove = FALSE) %>% 
  mutate(YEARMONTH = factor(YEARMONTH))


# Calculate variant frequency

# Uncomment if using the Collapsed pangos
#tmp <- tmp %>% select(-PANGOLIN_NOM) %>% rename("PANGOLIN_NOM" = Collapsed_pango)

freqs <- tmp %>% 
  # Count number of lineages per month per year
  count(PANGOLIN_NOM, YEARMONTH) %>% 
  # Add total number of sequences per month per year
  group_by(YEARMONTH) %>% 
  mutate("Total_seqs_pr_month" = sum(n)) %>% 
  ungroup() %>% 
  # Calculate percentages
  mutate(Percent = round(100 * n / Total_seqs_pr_month, 2)) %>% 
  arrange(desc(Percent))

# Pull out current and previous month
data <- freqs %>% 
  filter(YEARMONTH == current_yearmonth | YEARMONTH == previous_yearmonth) %>% 
  select(n, PANGOLIN_NOM, YEARMONTH, Percent) %>% 
  pivot_wider(names_from = YEARMONTH, values_from = c(n, Percent)) %>% 
  rename("CURR_MONTH_ANTALL" = paste0("n_", current_yearmonth),
         "PREV_MONTH_ANTALL" = paste0("n_", previous_yearmonth),
         "CURR_MONTH_PERC" = paste0("Percent_", current_yearmonth),
         "PREV_MONTH_PERC" = paste0("Percent_", previous_yearmonth)) %>% 
  arrange(desc(CURR_MONTH_PERC)) %>% # Sorterer etter høyest andel inneværende måned. Altså topp 10 inneværende måned 
  # Calculate percent change
  mutate(PERC_CHANGE = format((CURR_MONTH_PERC - PREV_MONTH_PERC) / PREV_MONTH_PERC * 100, digits = 2)) %>% 
  head(n=10) 
# Add n to percentage column in brackets
#mutate("CURR_MONTH_ANTALL" = paste0('(', CURR_MONTH_ANTALL, ')'),
#       "PREV_MONTH_ANTALL" = paste0('(', PREV_MONTH_ANTALL, ')')) %>% 
#unite("CURR_MONTH_PERC", c(CURR_MONTH_PERC, CURR_MONTH_ANTALL), sep = " ") %>%
#unite("PREV_MONTH_PERC", c(PREV_MONTH_PERC, PREV_MONTH_ANTALL), sep = " ")


#%>% 
# clean up column names. 
#  rename(!!current_yearmonth := CURR_MONTH_PERC,
#         !!previous_yearmonth := PREV_MONTH_PERC)

# Jeg må ha tabellen i long format
final_data <- data %>% 
  mutate(PERC_CHANGE = as.numeric(PERC_CHANGE)) %>% 
  pivot_longer(!PANGOLIN_NOM, names_to = "CATEGORY", values_to = "ANTALL") %>%  
  # Add flags for FHI Statistikk
  add_column("SPVFLAGG" = 0) %>%  # 0 er default og betyr at verdien finnes i tabellen
  mutate(SPVFLAGG = case_when(
    is.na(ANTALL) & str_detect(CATEGORY, "PERC") ~ 2, # Sett SPVFLAGG til 2 hvis NA for prosenter. Lar seg ikke beregne.
    is.na(ANTALL) & str_detect(CATEGORY, "ANTALL") ~ 1, # Sett SPVFLAGG til 1 hvis NA for antall. Altså at vi ikke har samlet inn data for denne perioden
    #ANTALL < 5 ~ 3, # Vi setter alle verdier lavere enn 5 til 0
    .default = 0
  )) %>% 
  # Endre NA til null
  mutate(ANTALL = replace_na(ANTALL, 0)) 

# Write as csv
write_csv(final_data, 
          file = paste0("data/SC2_samples/", Sys.Date(), "_variants_per_month.csv"))


# Jeg bør etterhvert generere en kategorifil automatisk? Det er vel bare variantene som vil endres?
# Bør også endre visningen av nåværende og forrige måned til den faktiske måneden
# Prepare to make the categories
variants <- final_data %>% distinct(PANGOLIN_NOM) %>% pull(PANGOLIN_NOM)
variants_title <- variants
v_levelName <- "Variant"

categories <- final_data %>% distinct(CATEGORY) %>% pull(CATEGORY)
categories_title <- c("Nåværende måned antall", "Forrige måned antall", "Nåværende måned %", "Forrige måned %", "% endring")
c_levelName <- "Kategori"

measure_type <- "ANTALL" 
measure_type_title <- "Antall"

flags <- 0:4
flags_title <- c("", "..", ".", ":", "-")
flags_text <- c("Verdi finnes i tabell", "Manglende data", "Lar seg ikke beregne", "Anonymisert", "Data ikke tilgjengelig")

# Tester å lage listen med Pangolins
variants_list <- vector("list", length = length(variants))
# Populere listen
for (i in length(variants_list)) {
  variants_list[[i]] <- "TEST"
}

# Hardcoded category list
category_file_list <- list(
  "Dimensions" = list(
    list("Value" = "PANGOLIN_NOM", 
         "Title" = "Variant",
         "Categories" = list(
           list(
             "Value" = "EG.5.1.1",
             "Title" = "EG.5.1.1",
             "LevelName" = "Variant"
           ),
           list(
             "Value" = "EG.5.1.3",
             "Title" = "EG.5.1.3",
             "LevelName" = "Variant"
           ),
           list(
             "Value" = "XBB.1.16",
             "Title" = "XBB.1.16",
             "LevelName" = "Variant"
           ),
           list(
             "Value" = "XBB.1.16.12",
             "Title" = "XBB.1.16.12",
             "LevelName" = "Variant"
           ),
           list(
             "Value" = "EG.6.1",
             "Title" = "EG.6.1",
             "LevelName" = "Variant"
           ),
           list(
             "Value" = "FE.1.1",
             "Title" = "FE.1.1",
             "LevelName" = "Variant"
           ),
           list(
             "Value" = "XBB.1.16.6",
             "Title" = "XBB.1.16.6",
             "LevelName" = "Variant"
           ),
           list(
             "Value" = "GE.1",
             "Title" = "GE.1",
             "LevelName" = "Variant"
           ),
           list(
             "Value" = "XBB.1.16.15",
             "Title" = "XBB.1.16.15",
             "LevelName" = "Variant"
           ),
           list(
             "Value" = "XBB.1.42.1",
             "Title" = "XBB.1.42.1",
             "LevelName" = "Variant"
           )
         )),
    list(
      "Value" = "CATEGORY",
      "Title" = "Kategori",
      "Categories" = list(
        list(
          "Value" = "CURR_MONTH_ANTALL",
          "Title" = "Nåværende måned antall",
          "LevelName" = "Kategori"
        ),
        list(
          "Value" = "PREV_MONTH_ANTALL",
          "Title" = "Forrige måned antall",
          "LevelName" = "Kategori"
        ),
        list(
          "Value" = "CURR_MONTH_PERC",
          "Title" = "Nåværende måned %",
          "LevelName" = "Kategori"
        ),
        list(
          "Value" = "PREV_MONTH_PERC",
          "Title" = "Forrige måned %",
          "LevelName" = "Kategori"
        ),
        list(
          "Value" = "PERC_CHANGE",
          "Title" = "% endring",
          "LevelName" = "Kategori"
        )
      )
    )),
  "MeasureTypes" = list(
    list(
      "Value" = "ANTALL",
      "Title" = "Antall"
    )
  ),
  "Flags" = list(
    list(
      "Value" = 0,
      "Title" = "",
      "Text" = "Verdi finnes i tabell"
    ),
    list(
      "Value" = 1,
      "Title" = "..",
      "Text" = "Manglende data"
    ),
    list(
      "Value" = 2,
      "Title" = ".",
      "Text" = "Lar seg ikke beregne"
    ),
    list(
      "Value" = 3,
      "Title" = ":",
      "Text" = "Anonymisert"
    ),
    list(
      "Value" = 4,
      "Title" = "-",
      "Text" = "Data ikke tilgjengelig"
    )
  )
)


jsonlite::write_json(category_file_list, "/home/jonr/Prosjekter/NIPH_NGS_vis/data/script_test.json")


# Variants of Interest ----------------------------------------------------


# Change current month if you need to
#current_month_label <- lubridate::month(current_month, label = T, abbr = T, locale = "Norwegian_Norway")
#current_yearmonth <- paste(current_year, current_month_label, sep = "-")

# Need to include last 4 months
if (current_month == 1) {
  previous_month             <- 12
  previous_month_label       <- lubridate::month(previous_month, label = T, abbr = T, locale = "Norwegian_Norway")
  previous_year              <- current_year-1
  previous_yearmonth         <- paste(previous_year, previous_month_label, sep = "-")
  previous_two_month         <- 11
  previous_two_month_label   <- lubridate::month(previous_two_month, label = T, abbr = T, locale = "Norwegian_Norway")
  previous_two_yearmonth     <- paste(previous_year, previous_two_month_label, sep = "-")
  previous_three_month       <- 10
  previous_three_month_label <- lubridate::month(previous_three_month, label = T, abbr = T, locale = "Norwegian_Norway")
  previous_three_yearmonth   <- paste(previous_year, previous_three_month_label, sep = "-")
} else if (current_month == 2) {
  previous_month             <- 1
  previous_month_label       <- lubridate::month(previous_month, label = T, abbr = T, locale = "Norwegian_Norway")
  previous_yearmonth         <- paste(current_year, previous_month_label, sep = "-")
  previous_year              <- current_year-1
  previous_two_month         <- 12
  previous_two_month_label   <- lubridate::month(previous_two_month, label = T, abbr = T, locale = "Norwegian_Norway")
  previous_two_yearmonth     <- paste(previous_year, previous_two_month_label, sep = "-")
  previous_three_month       <- 11
  previous_three_month_label <- lubridate::month(previous_three_month, label = T, abbr = T, locale = "Norwegian_Norway")
  previous_three_yearmonth   <- paste(previous_year, previous_three_month_label, sep = "-")
} else if (current_month == 3) {
  previous_month             <- 2
  previous_month_label       <- lubridate::month(previous_month, label = T, abbr = T, locale = "Norwegian_Norway")
  previous_yearmonth         <- paste(current_year, previous_month_label, sep = "-")
  previous_two_month         <- 1
  previous_two_month_label   <- lubridate::month(previous_two_month, label = T, abbr = T, locale = "Norwegian_Norway")
  previous_two_yearmonth     <- paste(current_year, previous_two_month_label, sep = "-")
  previous_year              <- current_year-1
  previous_three_month       <- 12
  previous_three_month_label <- lubridate::month(previous_three_month, label = T, abbr = T, locale = "Norwegian_Norway")
  previous_three_yearmonth   <- paste(previous_year, previous_three_month_label, sep = "-")
} else {
  previous_month             <- current_month - 1
  previous_month_label       <- lubridate::month(previous_month, label = T, abbr = T, locale = "Norwegian_Norway")
  previous_yearmonth         <- paste(current_year, previous_month_label, sep = "-")
  previous_two_month         <- current_month - 2
  previous_two_month_label   <- lubridate::month(previous_two_month, label = T, abbr = T, locale = "Norwegian_Norway")
  previous_two_yearmonth     <- paste(current_year, previous_two_month_label, sep = "-")
  previous_three_month       <- current_month - 3
  previous_three_month_label <- lubridate::month(previous_three_month, label = T, abbr = T, locale = "Norwegian_Norway")
  previous_three_yearmonth   <- paste(current_year, previous_three_month_label, sep = "-")
}

# Filter the data

# For the VOI I will take out the lineages called "BA.2.75", "XBB.1.5-like+F456L" and "XBB.1.5-like+L455F+F456L" from
# the Tessy column and count these

ba2.75 <- allvariants_v %>% 
  filter(Tessy == "BA.2.75") %>% 
  # Fix date format
  mutate("PROVE_TATT" = ymd(PROVE_TATT)) %>% 
  # select relevant columns
  select(PROVE_TATT, PANGOLIN_NOM, FULL_PANGO_LINEAGE, Collapsed_pango, Tessy) 

xbb1.5 <- allvariants_v %>% 
  filter(Tessy == "XBB.1.5-like+F456L") %>% 
  # Fix date format
  mutate("PROVE_TATT" = ymd(PROVE_TATT)) %>% 
  # select relevant columns
  select(PROVE_TATT, PANGOLIN_NOM, FULL_PANGO_LINEAGE, Collapsed_pango, Tessy)

xbb1.5_2 <- allvariants_v %>% 
  filter(Tessy == "XBB.1.5-like+L455F+F456L") %>% 
  # Fix date format
  mutate("PROVE_TATT" = ymd(PROVE_TATT)) %>% 
  # select relevant columns
  select(PROVE_TATT, PANGOLIN_NOM, FULL_PANGO_LINEAGE, Collapsed_pango, Tessy)

# combine data
tmp <- bind_rows(
  ba2.75,
  xbb1.5,
  xbb1.5_2
)

# Get total observations for each VOI
ba2.75_tot   <- tmp %>% count(Tessy, name = "ANTALL") %>% filter(Tessy == "BA.2.75")  %>% add_column("CATEGORY" = "TOTAL_OBS")
xbb1.5_tot   <- tmp %>% count(Tessy, name = "ANTALL") %>% filter(Tessy == "XBB.1.5-like+F456L")  %>% add_column("CATEGORY" = "TOTAL_OBS")
xbb1.5_2_tot   <- tmp %>% count(Tessy, name = "ANTALL") %>% filter(Tessy == "XBB.1.5-like+L455F+F456L")  %>% add_column("CATEGORY" = "TOTAL_OBS")

# Count per month of interest
pr_month <- tmp %>% 
  # Create year and month column
  mutate(YEAR = lubridate::year(PROVE_TATT),
         MONTH = lubridate::month(PROVE_TATT, label = T, abbr = T, locale = "Norwegian_Norway")) %>% 
  unite("YEARMONTH", c(YEAR, MONTH), sep = "-", remove = FALSE) %>% 
  mutate(YEARMONTH = factor(YEARMONTH)) %>% 
  # Filter out months of interest
  filter(YEARMONTH == current_yearmonth | YEARMONTH == previous_yearmonth | YEARMONTH == previous_two_yearmonth | YEARMONTH == previous_three_yearmonth) %>% 
  # Count VOI's per month
  count(YEARMONTH, Tessy, name = "ANTALL") %>% 
  # Rename YEARMONTH to CATEGORY for joining with total observations
  rename("CATEGORY" = YEARMONTH)

# Join the data
final_data <- bind_rows(
  ba2.75_tot,
  xbb1.5_tot,
  xbb1.5_2_tot,
  pr_month
) %>%  
  # Add flags for FHI Statistikk
  add_column("SPVFLAGG" = 0) %>%  # 0 er default og betyr at verdien finnes i tabellen
  # Complete the series to fill in any missing months for any variant
  complete(CATEGORY, Tessy) %>% 
  # NA in ANTALL means that we have no observations of that variant for that month
  mutate(SPVFLAGG = case_when(
    is.na(ANTALL) ~ 1,
    .default = SPVFLAGG
  )) %>% 
  # Change NA to 0
  mutate(ANTALL = replace_na(ANTALL, 0))

# Write as csv
write_csv(final_data, 
          file = paste0("data/SC2_samples/", Sys.Date(), "_variants_of_interest.csv"))



# Variants under monitoring -----------------------------------------------


# Filter the data

xbb1.16 <- allvariants_v %>% 
  filter(str_detect(Tessy, "XBB.1.6")) %>% 
  # Fix date format
  mutate("PROVE_TATT" = ymd(PROVE_TATT)) %>% 
  # select relevant columns
  select(PROVE_TATT, PANGOLIN_NOM, FULL_PANGO_LINEAGE, Collapsed_pango, Tessy) 

xbb1.5 <- allvariants_v %>% 
  filter(Tessy == "XBB.1.5-like+F456L") %>% 
  # Fix date format
  mutate("PROVE_TATT" = ymd(PROVE_TATT)) %>% 
  # select relevant columns
  select(PROVE_TATT, PANGOLIN_NOM, FULL_PANGO_LINEAGE, Collapsed_pango, Tessy)

xbb1.5_2 <- allvariants_v %>% 
  filter(Tessy == "XBB.1.5-like+L455F+F456L") %>% 
  # Fix date format
  mutate("PROVE_TATT" = ymd(PROVE_TATT)) %>% 
  # select relevant columns
  select(PROVE_TATT, PANGOLIN_NOM, FULL_PANGO_LINEAGE, Collapsed_pango, Tessy)

# combine data
tmp <- bind_rows(
  ba2.75,
  xbb1.5,
  xbb1.5_2
)

# Get total observations for each VOI
ba2.75_tot   <- tmp %>% count(Tessy, name = "ANTALL") %>% filter(Tessy == "BA.2.75")  %>% add_column("CATEGORY" = "TOTAL_OBS")
xbb1.5_tot   <- tmp %>% count(Tessy, name = "ANTALL") %>% filter(Tessy == "XBB.1.5-like+F456L")  %>% add_column("CATEGORY" = "TOTAL_OBS")
xbb1.5_2_tot   <- tmp %>% count(Tessy, name = "ANTALL") %>% filter(Tessy == "XBB.1.5-like+L455F+F456L")  %>% add_column("CATEGORY" = "TOTAL_OBS")

# Count per month of interest
pr_month <- tmp %>% 
  # Create year and month column
  mutate(YEAR = lubridate::year(PROVE_TATT),
         MONTH = lubridate::month(PROVE_TATT, label = T, abbr = T, locale = "Norwegian_Norway")) %>% 
  unite("YEARMONTH", c(YEAR, MONTH), sep = "-", remove = FALSE) %>% 
  mutate(YEARMONTH = factor(YEARMONTH)) %>% 
  # Filter out months of interest
  filter(YEARMONTH == current_yearmonth | YEARMONTH == previous_yearmonth | YEARMONTH == previous_two_yearmonth | YEARMONTH == previous_three_yearmonth) %>% 
  # Count VOI's per month
  count(YEARMONTH, Tessy, name = "ANTALL") %>% 
  # Rename YEARMONTH to CATEGORY for joining with total observations
  rename("CATEGORY" = YEARMONTH)

# Join the data
final_data <- bind_rows(
  ba2.75_tot,
  xbb1.5_tot,
  xbb1.5_2_tot,
  pr_month
) %>%  
  # Add flags for FHI Statistikk
  add_column("SPVFLAGG" = 0) %>%  # 0 er default og betyr at verdien finnes i tabellen
  # Complete the series to fill in any missing months for any variant
  complete(CATEGORY, Tessy) %>% 
  # NA in ANTALL means that we have no observations of that variant for that month
  mutate(SPVFLAGG = case_when(
    is.na(ANTALL) ~ 1,
    .default = SPVFLAGG
  )) %>% 
  # Change NA to 0
  mutate(ANTALL = replace_na(ANTALL, 0))

# Write as csv
write_csv(final_data, 
          file = paste0("data/SC2_samples/", Sys.Date(), "_variants_under_monitoring.csv"))

# DELETE ------------------------------------------------------------------




# MASSERT UTSLUSNINGSFIL LAGET I SIKKER SONE ------------------------------

# Les filen
df <- read_tsv("data/SC2_samples/2023.08.15-SARS-CoV-2_samples_for_FHI_Statistikk.tsv.zip")

# Denne filen inneholder alle SARS-CoV-2-prøver som har "Test_status" "A". 
# ToDo: Er alle disse relevante å ta med? Noen som skal ut?
# TODO: bør sjekke om det er noen fylkenavn som ikke er med i listen over
# TODO: Trenger vi dele inn i aldersgruppe? Blir fort for små tall...

# Hva er det jeg trenger å ta med videre?
# Kanskje enklere å trekke ut info separat og binde sammen senere enn å pivotere først?

# I MSIS i dag på FHI Statistikk kan man filtrere på Aldersgruppe, Kjønn, Fylke, og Måned (i tillegg til smittestd Norge eller utlandet).
# Alt dette bør derfor være regnet ut og tilgjengelig i samme csv-fil.

# Først regne ut antall prøver per måned per fylke 
# Denne informasjonen står i alle radene 

# Date formats have changed at some point. Harmonize this and merge before further processing
dates_1 <- df %>% 
  # Pull out dates om day.month.year format
  filter(str_detect(Sample_date, "\\.")) %>% 
  # Convert to year-month-day
  mutate(Sample_date = lubridate::dmy_hms(Sample_date))

dates_2 <- df %>% 
  # Pull out dates om day.month.year format
  filter(str_detect(Sample_date, "-")) %>% 
  # Convert to year-month-day
  mutate(Sample_date = lubridate::ymd(Sample_date)) 

dates <- rbind(dates_1, dates_2)

df_mod <- dates %>% 
  select("PROVE_TATT" = Sample_date, 
         "FYLKE" = Fylke, 
         "FODT_AAR" = Born_year, 
         "FODT_MND" = Born_month, 
         "KJONN" = Gender) %>% 
  # Remove identical rows
  distinct() %>% 
  # Rename NA til "Ukjent fylke" and convert to factor
  mutate(FYLKE = factor(replace_na(FYLKE, "Ukjent fylke"))) %>% 
  # Drop samples with no date 
  filter(!is.na(PROVE_TATT)) %>% 
  # Add year and month columns
  mutate(MONTH = month(PROVE_TATT),
         month_name = month(PROVE_TATT, label = T, abbr = T, locale = "Norwegian_Norway"),
         YEAR = year(PROVE_TATT)) %>% 
  unite("YEARMONTH", c(YEAR, month_name), sep = "-", remove = FALSE) %>% 
  mutate(YEARMONTH = factor(YEARMONTH)) %>% 
  # Remove samples with wrongly formatted dates
  filter(!is.na(PROVE_TATT))

# Calculate samples per month per year per Fylke and total for the country
fylke_month <- df_mod %>% 
  #count(MONTH, YEAR, Fylke) %>% 
  count(YEARMONTH, FYLKE) %>% 
  rename("ANTALL" = "n")

# Calculate samples per month per year for all Fylker
total_month <- df_mod %>% 
  count(YEARMONTH) %>% 
  rename("ANTALL" = "n") %>% 
  add_column("FYLKE" = "Totalt") %>% 
  select(YEARMONTH, FYLKE, ANTALL)

# Calculate per age group per fylke per month per year

# Combine data
data <- bind_rows(fylke_month, total_month) 

# Join with final data structure to get all possible combinations of month, year and fylke

final_data <- left_join(final_data, data, by = join_by("YEARMONTH", "FYLKE")) %>% 
  # Create new  columns
  separate(YEARMONTH, into = c("YEAR", "MONTH"), sep = "-", remove = T) %>% 
  mutate(YEAR = as.factor(YEAR),
         MONTH = as.factor(MONTH),
         FYLKE = as.factor(FYLKE)) %>% 
  # Add flags for FHI Statistikk
  add_column("SPVFLAGG" = 0) %>% # 0 er default og betyr at verdien finnes i tabellen
  mutate(SPVFLAGG = case_when(
    is.na(ANTALL) ~ 1, # Sett SPVFLAGG til 1 hvis NA. Altså at vi ikke har samlet inn data for denne perioden
    ANTALL < 5 ~ 3, # Vi setter alle verdier lavere enn 5 til 0
    .default = 0
  )) %>% 
  # Replace NA with 0
  mutate(ANTALL = replace_na(ANTALL, 0)) %>% 
  # Change all numbers < 5 to zero
  mutate(ANTALL = case_when(
    ANTALL < 5 ~ 0,
    .default = ANTALL
  )) 


# Write as csv
write_csv(final_data, 
          file = paste0("data/SC2_samples/", Sys.Date(), "_seqs_per_month_per_fylke.csv"))



# Pangolin lineages -------------------------------------------------------

# Tanker:
# På hvilket pango-nivå skal vi legge oss?
# Skal vi stratifisere etter fylke? Blir mye prikking kanskje?
# Hvordan skal kategori-filen utformes? Blir det samme som fylke bare for pango?


# UTSLUSNINGSFILENE FRA SIKKER SONE ---------------------------------------


# Les alle utslusningsfilene fra LabWare
files <- list.files(path = "N:/Virologi/Influensa/2223/LabwareUttrekk/", 
                    pattern = "csv$",
                    full.names = TRUE)

# Trekk ut SARS-CoV-2 prøver som er svar ut (godkjent)

df <- tibble()
for (i in 1:length(files)) {
  
  # Clear old objects
  if (exists("tmp")) {
    rm(tmp)
  }
  
  # Read data files
  x <- read_delim(files[i], delim = ";", col_names = TRUE, locale=locale(encoding="latin1")) %>% 
    filter(`AGENS Agens` == "SARS-CoV-2") %>%
    filter(Test_status == "A")
  
  # Continue only if relevant columns exist
  if (all(c("Key", "Prøvedato2", "Fylke", "FylkeNr", "Kjønn", "Innsender", "Landsdel", "Stammenavn") %in% colnames(x))) {
    tmp <- x %>%
      # Select relevant columns
      select("Key", "Prøvedato2", "Fylke", "FylkeNr", "Kjønn", "Innsender", "Landsdel", "Stammenavn") %>%
      # Convert to character for combining later
      mutate(Key = as.character(Key),
             FylkeNr = as.character(FylkeNr),
             Innsender = as.character(Innsender))
  }
  
  # Combine files
  if (exists("tmp")) {
    df <- bind_rows(df, tmp)
  }
}

# Remove duplicate entries (identical rows)
df <- df %>% distinct()

# Drop samples and create columns
df_mod <- df %>% 
  # Drop samples with no date or Fylke
  filter(!is.na(Prøvedato2)) %>% 
  filter(Fylke != "0") %>% 
  # Fix date format
  rename("PROVE_TATT" = Prøvedato2) %>% 
  rename("FYLKE" = "Fylke") %>% 
  # Convert to date
  mutate(PROVE_TATT = dmy(PROVE_TATT)) %>% 
  # Add year and month columns
  mutate(MONTH = month(PROVE_TATT),
         month_name = month(PROVE_TATT, label = T, abbr = T, locale = "Norwegian_Norway"),
         YEAR = year(PROVE_TATT)) %>% 
  unite("YEARMONTH", c(YEAR, month_name), sep = "-", remove = FALSE) %>% 
  mutate(YEARMONTH = factor(YEARMONTH)) %>% 
  # Remove samples with wrongly formatted dates
  filter(!is.na(PROVE_TATT))


# Calculate samples per month per year per Fylke and total for the country
fylke_month <- df_mod %>% 
  #count(MONTH, YEAR, Fylke) %>% 
  count(YEARMONTH, FYLKE) %>% 
  rename("ANTALL" = "n")

# Calculate samples per month per year for all Fylker
total_month <- df_mod %>% 
  count(YEARMONTH) %>% 
  rename("ANTALL" = "n") %>% 
  add_column("FYLKE" = "Totalt")

# Combine data
data <- bind_rows(fylke_month, total_month) 

# TODO: How and where to order the YEARMONTH column correctly? 
# Should be independent of the data
https://stackoverflow.com/questions/15103562/sort-year-month-column-by-year-and-month

# Join with final data structure to get all possible combinations of month, year and fylke
final_data <- left_join(final_data, data) %>% 
  rename("ANTALL" = "n") %>% 
  # Create new factor columns
  mutate(
    "MONTH"  = factor(MONTH),
    "YEAR"   = factor(YEAR),
    "FYLKE"  = factor(FYLKE)
  )


# Add flags for FHI Statistikk
final_data <- final_data %>% 
  add_column("SPVFLAGG" = 0)

# Replace NA with 0
final_data <- final_data %>% 
  mutate(ANTALL = replace_na(ANTALL, 0)) 

# Write as csv
write_csv(final_data, file = "data/SC2_samples/seqs_per_month.csv")


# Create the category file

# First creating a list with the different levels necessary

# Need three top levels: 
## Dimensions
### [[1]]
#### Value
#### Title
#### Categories

### [[2]]
### [[3]]
## MeasureTypes
## Flags

# Read example file
example <- read_json("data/SC2_samples/seqs_per_month.category.json")

my_list <- list()

list(levels(final_data$YEAR))

{
  "Dimensions": [
    {
      "Value": "YEAR",
      "Title": "År",
      "Categories": [
        {
          "Value": "2022",
          "Title": "2022",
          "LevelName": "År"
        },
        {
          "Value": "2023",
          "Title": "2023",
          "LevelName": "År"
        }
      ]
    },
    {
      "Value": "Fylke",
      "Title": "Fylke",
      "Categories": [
        {
          "Value": "Agder",
          "Title": "Agder",
          "LevelName": "Fylke"
        },
        {
          "Value": "Innlandet",
          "Title": "Innlandet",
          "LevelName": "Fylke"
        },
        {
          "Value": "Møre og Romsdal",
          "Title": "Møre og Romsdal",
          "LevelName": "Fylke"
        },
        {
          "Value": "Nordland",
          "Title": "Nordland",
          "LevelName": "Fylke"
        },
        {
          "Value": "Oslo",
          "Title": "Oslo",
          "LevelName": "Fylke"
        },
        {
          "Value": "Rogaland",
          "Title": "Rogaland",
          "LevelName": "Fylke"
        },
        {
          "Value": "Troms og Finnmark",
          "Title": "Troms og Finnmark",
          "LevelName":"Fylke"
        },
        {
          "Value": "Trøndelag",
          "Title": "Trøndelag",
          "LevelName":"Fylke"
        },
        {
          "Value": "Vestfold og Telemark",
          "Title": "Vestfold og Telemark",
          "LevelName":"Fylke"
        },
        {
          "Value": "Vestland",
          "Title": "Vestland",
          "LevelName":"Fylke"
        },
        {
          "Value": "Viken",
          "Title": "Viken",
          "LevelName":"Fylke"
        }
      ]
    }
  ],
  "MeasureTypes": [
    {
      "Value": "antall",
      "Title": "Antall"
    }
  ],
  "Flags": [
    {
      "Value": 0,
      "Title": "",
      "Text": ""
    },
    {
      "Value": 1,
      "Title": "..",
      "Text": "Manglende data"
    },
    {
      "Value": 2,
      "Title": ".",
      "Text": "Lar seg ikke beregne"
    },
    {
      "Value": 3,
      "Title": ":",
      "Text": "Anonymisert"
    },
    {
      "Value": 4,
      "Title": "-",
      "Text": "Data ikke tilgjengelig"
    }
  ]
}
