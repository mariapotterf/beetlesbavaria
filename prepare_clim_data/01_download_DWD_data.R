# Get climate data from the DWD: Deutsche weather institute:

library(rvest)
library(stringr)

page_link <- "https://opendata.dwd.de/climate_environment/CDC/grids_germany/monthly/air_temperature_mean/"
page_link_precip <- 'https://opendata.dwd.de/climate_environment/CDC/grids_germany/monthly/precipitation/'

month_name <- c("01_Jan",
                "02_Feb",
                "03_Mar",
                "04_Apr",
                "05_May",
                "06_Jun",
                "07_Jul",
                "08_Aug",
                "09_Sep",
                "10_Oct",
                "11_Nov",
                "12_Dec")
                
# you can set month_name as "05_May" to get the data from 05_May
download_data <- function(x, ...) {
  print(x)
  
  month_name = x
  print(month_name)
  
  # getting the html page for 01_Jan folder
  page <- read_html(paste0(page_link_precip, month_name, "/"))

  # getting the link text
  link_text <- page %>%
    html_elements("a") %>%
    html_text()

  # creating links
  links <- paste0(page_link_precip, month_name, "/", link_text)[-1]

  # extracting the numbers for filename
  filenames <- stringr::str_extract(pattern = "\\d+", string = link_text[-1])

  # creating a directory
  dir.create(month_name)

  # setting the option for maximizing time limits for downloading
  options(timeout = max(600, getOption("timeout")))

  # downloading the file
  for (i in seq_along(links)) {
    download.file(links[i], paste0(month_name, "/", filenames[i], "asc.gz"))
  }
  
}


lapply(month_name, download_data)
