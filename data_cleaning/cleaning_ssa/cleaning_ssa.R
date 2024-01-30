### Script to read in data from the SSA downloaded from the Bureau of 
### Meteorology (http://www.bom.gov.au/australia/stormarchive/) and clean it 
### into a form acceptable for Python.
###
### Last modified: 2023-01-31
### Author: Isabelle Greco

# get args from command line
# 1 = library for r installs if need be
# 2 = location of original data file
# 3 = location of output data file
args = commandArgs(trailingOnly = TRUE)

# assert correct number of args
if (length(args) != 3) {
  stop(paste("Must give library for R package installs, data location, and",
             "output location."), call. = FALSE)
}

# load libraries + install if need be
for (package in c("tidyverse", "janitor", "lubridate")) {
    if (!require(package, character.only = TRUE)) {
	# if not installed, install in the directory given by first argument
        install.packages(package, lib = args[1], 
                         repos = "https://cran.csiro.au", dependencies = TRUE)
        library(package, character.only = TRUE)
    }
}


# read in data - location is second argument
ssa <- read_csv(args[2], col_types = "iTddccdcc") %>%
  clean_names()

# row 961 is causing the extra column - merge these two columns after cleaning 
# data for extraneous symbols/punctuation
ssa_clean_comment <- ssa %>% 
  # firstly getting rid of %D
  mutate(comments_tidy1 = gsub('%D', "", comments)) %>%
  mutate(x9_tidy1 = gsub('%D', "", x9)) %>% 
  # then only keeping a-z, 0-9, and some punctuation
  mutate(comments_tidy2 = gsub('[^-<>\\.\\:\\(\\)[:alnum:] ]', "", 
                               comments_tidy1)) %>%
  mutate(x9_tidy2 = gsub('[^-<>\\.\\:\\(\\)[:alnum:] ]', "", x9_tidy1)) %>%
  # joining two colums
  mutate(comments_clean = case_when(
    is.na(x9_tidy2) ~ comments_tidy2,
    TRUE ~ paste(comments_tidy2, x9_tidy2, sep = " ")
    )) %>%
  # empty strings to na
  mutate(comments_clean = case_when(comments_clean == "" ~ as.character(NA),
                                    TRUE ~ comments_clean)) %>%
  # drop unecessary intermediary columns
  select(-c(comments, comments_tidy1, comments_tidy2, x9, x9_tidy1, x9_tidy2))

# parses ID 2086 effectively but also adds the extra row of all NA
ssa_valid_id <- ssa_clean_comment %>%
  filter(!is.na(hail_id))

# all towns to lower case for consistency
ssa_lower_town <- ssa_valid_id %>%
  mutate(nearest_town = tolower(nearest_town))

# updating impossible locations
ssa_better_locations <- ssa_lower_town %>%
  # use nearest town
  rows_update(tibble(hail_id = 4035, longitude = 153.459), by = "hail_id") %>% 
  rows_update(tibble(hail_id = 4046, latitude = -32.245, longitude = 148.619),
              by = "hail_id") %>% 
  # use description to approximate
  # taking point between metro area and reported snow 
  rows_update(tibble(hail_id = 4076, latitude = -33.035, longitude = 117.005), 
              by = "hail_id") %>% 
  # unclear where exactly Dunrock is but using Hyden
  rows_update(tibble(hail_id = 4077, latitude = -32.446, longitude = 119.103), 
              by = "hail_id") %>% 
  # state was also wrong from description
  rows_update(tibble(hail_id = 3789, latitude = -33.912, longitude = 138.883, 
                     state = "SA"), by = "hail_id") %>%
  # the following had a location of lat = 0, lon = 0
  # taking half way between Kempsey and Crescent head
  rows_update(tibble(hail_id = 2871, latitude = -31.150, longitude = 152.910),
              by = "hail_id") %>%
  # halfway between moree and collarenebri as per comments
  rows_update(tibble(hail_id = 2881, latitude = -29.463, longitude = 149.192), 
              by = "hail_id") %>%
  # yarralumla, canberra as per nearest town - no ACT included in the data 
  rows_update(tibble(hail_id = 2912, latitude = -35.305, longitude = 149.101),
              by = "hail_id")


# quality control from manual examination: all 00:00:00 UTC cases unless stated
# otherwise
ssa_qc_study_area <- ssa_better_locations %>%
  # comment + ABC suggests late afternoon for related reports
  # abc.net.au/news/2010-12-17/storms-knock-out-power-to-thousands/2378530
  rows_update(tibble(hail_id = 3230, 
                     date_time = as_datetime("2010-12-17 09:00:00")), 
	      by = "hail_id") %>%
  rows_update(tibble(hail_id = 3231, 
                     date_time = as_datetime("2010-12-17 09:00:00")),
	      by = "hail_id") %>%
  # reports in Ramsay, Tawoomba, Laidley all very close (<60km)
  # Laidley report has the most accurate time (2013-11-13 06:06:00) but ABC says
  # after 4pm local (06:00:00 UTC)
  # abc.net.au/news/2013-11-14/qlds-darling-downs-cleans-up-after-severe-storms/5090704
  rows_update(tibble(hail_id = 4096, 
                     date_time = as_datetime("2013-11-13 06:00:00")),
	      by = "hail_id") %>%
  rows_update(tibble(hail_id = 4097, 
                     date_time = as_datetime("2013-11-13 06:00:00")),
	      by = "hail_id") %>%
  rows_update(tibble(hail_id = 4098, 
                     date_time = as_datetime("2013-11-13 06:00:00")),
	      by = "hail_id") %>%
  # duplication of the Laidley report at same time of previous day
  filter(hail_id != 4137) %>%
  # midnight UTC recorded but ABC article and radar suggests afternoon
  # abc.net.au/news/2011-02-21/severe-storms-sweep-southern-queensland/1951762
  rows_update(tibble(hail_id = 3241, 
                     date_time = as_datetime("2011-02-21 03:00:00")),
	      by = "hail_id") %>%
  # one of many afternoon reports in the area - this one had 00:00:00 UTC and 
  # the rest were 0310 to 0815 UTC (1310 to 1815 UTC), we take the centre
  # abc.net.au/news/2012-11-18/severe-storms-sweep-through-se-queensland/4378430?nw=0&r=HtmlFragment
  rows_update(tibble(hail_id = 3892, 
                     date_time = as_datetime("2012-11-18 05:45:00")),
	      by = "hail_id") %>%
  # coordinates exact same as second report from Coffs Harbour but six hours 
  # earlier. likely to be an error, as ABC implies damage throughout the 
  # afternoon and warning put out at 1230 AEST
  # abc.net.au/news/2014-12-09/violent-storm-murwillumbah/5955544
  # ewn.com.au/alerts/2014-12-09-034200-82267-432.weather
  rows_update(tibble(hail_id = 4024, 
                     date_time = as_datetime("2014-12-09 06:00:00")),
	      by = "hail_id") %>%
  # no media attention, but radar suggests more likely to be afternoon
  rows_update(tibble(hail_id = 3927, 
                     date_time = as_datetime("2012-02-13 03:00:00")),
	      by = "hail_id") %>%
  # given radar and time precision likely a mix up between UTC time and local
  # time
  rows_update(tibble(hail_id = 3882, 
                     date_time = as_datetime("2012-02-13 04:14:00")),
	      by = "hail_id") %>%
  # photographic evidence of evening lightning
  # couriermail.com.au/news/queensland/warwick/warwick-storms-2013/image-gallery/85ed3426757b0ad24beeb590ce20c6c1
  # radar suggests likely to be later afternoon 
  rows_update(tibble(hail_id = 3928, 
                     date_time = as_datetime("2013-01-19 03:00:00")),
	      by = "hail_id") %>%
  # given timing of neighbouring report, likely to have been a mix up again
  # between local and UTC time
  rows_update(tibble(hail_id = 3883, 
                     date_time = as_datetime("2012-10-14 10:02:00")),
	      by = "hail_id")

# just need to save to csv - third argument
write_csv(ssa_qc_study_area, args[3])
