library(tidyverse)
library(lubridate)


# Set up data structures --------------------------------------------------

# Create the data structure that will be filled with values. Using strings only to get factors later in crossing()
month_numbers <- as.character(1:12)
month_names   <- lubridate::month(1:12, label = T, abbr = T, locale = "Norwegian_Norway") 
years         <- c("2023")

# Create variables to use in filtering
current_year <- year(Sys.Date())
current_month <- month(Sys.Date())
current_month_label <- month(current_month, label = T, abbr = T, locale = "Norwegian_Norway")
current_yearmonth <- paste(current_year, current_month_label, sep = "-")
  
if (current_month == 12) {
  previous_month <- 1
} else {
  previous_month <- current_month - 1
}
previous_month_label <- month(previous_month, label = T, abbr = T, locale = "Norwegian_Norway")
previous_yearmonth <- paste(current_year, previous_month_label, sep = "-")

# Create all possible combinations of month and year 
final_data <- crossing(month_names, years) %>% 
  rename(
    "MONTH"  = "month_names",
    "YEAR"   = "years"
  ) %>% 
  # Create YEARMONTH column for combining later
  unite("YEARMONTH", c("YEAR", "MONTH"), sep = "-", remove = T) %>% 
  mutate(YEARMONTH = as.factor(YEARMONTH))


# Get data from BN. Needs to be changed from Andreas ----------------------

con <- odbc::dbConnect(odbc::odbc(),
                 Driver = "SQL Server",
                 Server = "sql-bn-covid19",
                 Database = "BN_Covid19")

# Query the ENTRYTABLE table to extract relevant sequencing data
allvariants <- tbl(con, "ENTRYTABLE") %>%
  select(
    KEY,
    PROVE_TATT,
    PROSENTDEKNING_GENOM,
    DEKNING_NANOPORE,
    COVERAGE_BREADTH_SWIFT,
    COVERAGE_BREATH_EKSTERNE,
    PANGOLIN_NOM,
    FULL_PANGO_LINEAGE,
    
  ) %>%
  rename(
    "Dekning_Artic" = PROSENTDEKNING_GENOM,
    "Dekning_Swift" = COVERAGE_BREADTH_SWIFT,
    "Dekning_Eksterne" = COVERAGE_BREATH_EKSTERNE,
    "Dekning_Nano" = DEKNING_NANOPORE,
  ) %>%
  collect()

# Convert all column values to UTF-8 encoding
allvariants[] <- lapply(allvariants, function(x)
  iconv(x, from = "ISO-8859-1", to = "UTF-8"))

# Filter out low quality samples
allvariants_v <- allvariants  %>%
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
  filter (PANGOLIN_NOM != "Unassigned")%>%
  filter (PANGOLIN_NOM != "NA") %>%
  filter(PROVE_TATT != "")

######################   Collapsing Pangolins  ##################################
# https://www.who.int/activities/tracking-SARS-CoV-2-variants go here to update the list 
# Add a Collapsed pango column to the dataset based on the long Pangolin lineage 

# https://raw.githubusercontent.com/cov-lineages/pango-designation/master/lineage_notes.txt

allvariants_v <- allvariants_v %>%
  mutate(Collapsed_pango = case_when(
    FULL_PANGO_LINEAGE == "B.1.1.529.2.86" ~ "BA.2.86.X",
    FULL_PANGO_LINEAGE == "XBB.1.9.2.5.1" ~ "EG.5.1.X",
    FULL_PANGO_LINEAGE == "BJ.1BM.1.1.1].1.16" ~ "XBB.1.16.X",
    FULL_PANGO_LINEAGE == "BJ.1BM.1.1.1].2.3" ~ "XBB.2.3.X",
    FULL_PANGO_LINEAGE == "XBB.2.3" ~ "XBB.2.3.X",
    FULL_PANGO_LINEAGE == "XBB.1.9.2.5" ~ "EG.5.X",
    FULL_PANGO_LINEAGE == "XBB.1.16" ~ "XBB.1.16.X",
    FULL_PANGO_LINEAGE == "XBB.1.9.1" ~ "XBB.1.9.1.X",
    FULL_PANGO_LINEAGE == "XBB.1.9.2" ~ "XBB.1.9.2.X",
    FULL_PANGO_LINEAGE == "XBB.1.5" ~ "XBB.1.5.X",
    FULL_PANGO_LINEAGE == "BJ.1BM.1.1.1].1.9.1" ~ "XBB.1.9.1.X",
    FULL_PANGO_LINEAGE == "BJ.1BM.1.1.1].1.9.2" ~ "XBB.1.9.2.X",
    FULL_PANGO_LINEAGE == "B.1.1.529.2.75.3.4.1.1.1.1" ~ "CH.1.1.X",
    FULL_PANGO_LINEAGE == "BJ.1BM.1.1.1].1.5" ~ "XBB.1.5.X",
    FULL_PANGO_LINEAGE == "BJ.1BM.1.1.1].1.9" ~ "XBB.1.9.X",
    FULL_PANGO_LINEAGE == "XBB.1.9" ~ "XBB.1.9.X",
    FULL_PANGO_LINEAGE == "XBB." ~ "XBB.X",
    FULL_PANGO_LINEAGE == "BJ.1BM.1.1.1" ~ "XBB.X",
    grepl("^\\[", FULL_PANGO_LINEAGE) ~ "Andre rekombinanter",
    FULL_PANGO_LINEAGE == "X" ~ "Andre rekombinanter",
    FULL_PANGO_LINEAGE == "B.1.1.529.5.3.1.1.1.1.1" ~ "BQ.1.X",
    FULL_PANGO_LINEAGE == "B.1.1.529.2.75" ~ "BA.2.75.X",
    FULL_PANGO_LINEAGE == "B.1.1.529.5" ~ "BA.5.X",
    FULL_PANGO_LINEAGE == "B.1.1.529.4" ~ "BA.4.X",
    FULL_PANGO_LINEAGE == "B.1.1.529.3" ~ "BA.3.X",
    FULL_PANGO_LINEAGE == "B.1.1.529.2" ~ "BA.2.X",
    FULL_PANGO_LINEAGE == "B.1.1.529.1" ~ "BA.1.X",
    FULL_PANGO_LINEAGE == "B.1.1.529" ~ "B.1.1.529",
    FULL_PANGO_LINEAGE == "B.1.617" ~ "Delta",
    FULL_PANGO_LINEAGE == "B.1.1.7" ~ "Alpha",
    FULL_PANGO_LINEAGE == "B.1.351" ~ "Beta",
    FULL_PANGO_LINEAGE == "B.1.1.28.1" ~ "Gamma",
    FULL_PANGO_LINEAGE == "B.1.427" | FULL_PANGO_LINEAGE == "B.1.429" ~ "Epsilon",
    FULL_PANGO_LINEAGE == "B.1.1.28.2" ~ "Zeta",
    FULL_PANGO_LINEAGE == "B.1.525" ~ "Eta",
    FULL_PANGO_LINEAGE == "B.1.1.28.3" ~ "Theta",
    FULL_PANGO_LINEAGE == "B.1.526" ~ "Iota",
    FULL_PANGO_LINEAGE == "B.1.617.1" ~ "Kappa",
    FULL_PANGO_LINEAGE == "B.1.1.1.37.1" ~ "Lambda",
    FULL_PANGO_LINEAGE == "B.1.621" ~ "Mu",
    TRUE ~ "Andre SARS CoV 2"
  ))

# Clean up
rm(allvariants)

# Close connection 
closeAllConnections()

# Convert comma to dot in the coverage
allvariants_v <- allvariants_v %>% 
  # Replace a few double commas
  mutate(Dekning_Nano = str_replace(Dekning_Nano, ",,", ",")) %>% 
  # Change comma to decimal for the coverage
  mutate(Dekning_Artic = str_replace(Dekning_Artic, ",", "."),
         Dekning_Swift = str_replace(Dekning_Swift, ",", "."),
         Dekning_Nano = str_replace(Dekning_Nano, ",", "."),
         Dekning_Eksterne = str_replace(Dekning_Eksterne, ",", "."))

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
         MONTH = month(PROVE_TATT, label = T, abbr = T, locale = "Norwegian_Norway")) %>% 
  unite("YEARMONTH", c(YEAR, MONTH), sep = "-", remove = FALSE) %>% 
  mutate(YEARMONTH = factor(YEARMONTH))
  

# Calculate variant frequency ---------------------------------------------

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
  arrange(desc(CURR_MONTH_PERC)) %>% 
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

# Create a list
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
