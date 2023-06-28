library(tidyverse)
library(lubridate)
library(jsonlite)

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

# Drop samples with no sample date
df <- df %>% 
  filter(!is.na(Prøvedato2)) %>% 
  # Fix date format
  rename("PROVE_TATT" = Prøvedato2) %>% 
  # Convert to date
  mutate(PROVE_TATT = dmy(PROVE_TATT)) %>% 
  # Add year and month columns
  mutate(MONTH = month(PROVE_TATT),
         YEAR = year(PROVE_TATT))


# Calculate samples per month per year per Fylke and total for the country
fylke_month <- df %>% 
  group_by(MONTH, YEAR, Fylke) %>% 
  count() %>% 
  # Remove empty Fylkenavn
  filter(str_detect(Fylke, "")) %>% 
  # Some Fylke are named zero
  filter(Fylke != "0") %>% 
  # Drop empty year info
  filter(!is.na(YEAR))

# Calculate samples per month per year for all Fylker
total_month <- df %>% 
  group_by(MONTH, YEAR) %>% 
  count() %>% 
  add_column("Fylke" = "Totalt") %>% 
  # Drop empty year info
  filter(!is.na(YEAR))

# Combine data
data <- bind_rows(fylke_month, total_month) %>% 
  # Complete series
  ungroup() %>% 
  complete(MONTH, YEAR, Fylke)

# Add flags for FHI Statistikk
data <- data %>% 
  add_column("SPVFLAGG" = 0)

# Replace NA with 0
data <- data %>% 
  mutate(n = replace_na(n, 0)) 

# Write as csv
write_csv(data, file = "data/SC2_samples/seqs_per_month.csv")


# Create the category file


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
