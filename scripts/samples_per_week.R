library(tidyverse)
library(lubridate)
library(jsonlite)

#### THOUGHTS ####
# Maybe don't always include all months and all years etc.
# Only include what we have data for 
# Then it will be important to generate the category file from the actual data

# Maybe better to just operate with year-months together

# Create the data structure that will be filled with values. Using strings only to get factors later in crossing()
month_numbers <- as.character(1:12)
fylke_names   <- c("Agder", "Innlandet", "Møre og Romsdal", "Nordland", "Oslo", "Rogaland", "Troms og Finnmark", "Trøndelag", "Vestfold og Telemark", "Vestland", "Viken", "Ukjent fylke", "Totalt")
years         <- c("2022", "2023")

# Create all possible combinations of fylke number, fylke name and year
final_data <- crossing(month_numbers, years, fylke_names) %>% 
  # Create new factor columns
  rename(
    "MONTH"  = "month_numbers",
    "YEAR"   = "years",
    "FYLKE"  = "fylke_names"
  )


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
