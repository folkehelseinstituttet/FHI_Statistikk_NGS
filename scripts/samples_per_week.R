library(tidyverse)
library(lubridate)
library(jsonlite)

#### THOUGHTS ####
# Maybe don't always include all months and all years etc.
# Only include what we have data for 
# Then it will be important to generate the category file from the actual data

# Create the data structure that will be filled with values. Using strings only to get factors later in crossing()
month_numbers <- as.character(1:12)
fylke_names   <- c("Agder", "Innlandet", "Møre og Romsdal", "Nordland", "Oslo", "Rogaland", "Troms og Finnmark", "Trøndelag", "Vestfold og Telemark", "Vestland", "Viken", "Totalt")
years         <- c("2022", "2023")

# Create all possible combinations of fylke number, fylke name and year
final_data <- crossing(month_numbers, years, fylke_names) %>% 
  # Create new factor columns
  rename(
    "MONTH"  = "month_numbers",
    "YEAR"   = "years",
    "FYLKE"  = "fylke_names"
  )

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
  ungroup() %>% 
  # Remove empty Fylkenavn
  filter(str_detect(Fylke, "")) %>% 
  # Some Fylke are named zero
  filter(Fylke != "0") %>% 
  # Drop empty year info
  filter(!is.na(YEAR)) %>% 
  rename("FYLKE" = "Fylke") %>% 
  mutate(MONTH = as.character(MONTH),
         YEAR = as.character(YEAR))

# Calculate samples per month per year for all Fylker
total_month <- df %>% 
  group_by(MONTH, YEAR) %>% 
  count() %>% 
  ungroup() %>% 
  add_column("Fylke" = "Totalt") %>% 
  # Drop empty year info
  filter(!is.na(YEAR)) %>% 
  rename("FYLKE" = "Fylke") %>% 
  mutate(MONTH = as.character(MONTH),
         YEAR = as.character(YEAR))

# Combine data
data <- bind_rows(fylke_month, total_month) 

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
