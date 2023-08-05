# R HElp


# filetr only locations that occurs oevr all years:

# Example: keep only records that are present over all categories
df <- data.frame(loc = c('a', "b", "c",
                         'a', "b",
                         "b", "c",
                         "b", "c"
),
year = c(1,1,1,
         2,2,
         3,3,
         4,4))

# keep only rows that have the same amount of observations
loc.unique <- unique(df$loc)
year.unique <- unique(df$year)

df %>% 
  # group_by(years, loc) %>% 
  group_by(loc) %>%
  filter(n_distinct(year) == n_distinct(df$year))
