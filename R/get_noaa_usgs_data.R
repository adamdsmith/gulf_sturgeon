if (!requireNamespace("remotes", quietly = T)) install.packages("remotes")
if (!requireNamespace("pacman", quietly = T)) install.packages("pacman")
if (!requireNamespace("nrsmisc", quietly = T)) remotes::install_github("adamdsmith/nrsmisc")
pacman::p_load(readr, nrsmisc, dplyr, purrr, lubridate, dataRetrieval)

# Water temperature from a couple of NOAA tide/weather stations
noaa_stns <- c(8761305, 8761927) # Shell Beach & New Canal Station, LA
noaa <- lapply(noaa_stns, function(stn) {
  tmp <- get_coop(start = "2017-06-01", end = Sys.Date(),
                  station_name = stn, product = "water_temperature",
                  verbose = FALSE)
})
names(noaa) <- noaa_stns
noaa <- bind_rows(noaa, .id = "site_no") %>%
  select(site_no, dateTime = t, tempC = v) %>%
  filter(!is.na(tempC))
attr(noaa$dateTime, "tzone") <- "America/Chicago"
saveRDS(noaa, file = "./Output/noaa_temps.rds")

# USGS gages of interest (note that two have only gage height data)
usgs_stns <- c("02479310", "02489500", "02490500", "02492000", "02492111", 
               "02492600", "02492620", "02492700", "07375000", "07375280", "07375500", 
               "07375800", "07375960", "07376000", "07376500", "300722089150100", 
               "301001089442600", "301104089253400", "301141089320300", "301429089145600")

# Parameters of interest
# water temp (C), daily discharge (ft3/s), salinity (ppt), mean water velocity (ft/s)
water_parms <- c("00010", "00060", "00480", "72255")
water_parms_nm <- c("tempC", "discharge_ft3s", "salinity_ppt", "velocity_fts")

# See what unit/daily data are available for each station
avail <- whatNWISdata(siteNumber = usgs_stns,
                      parameterCd = water_parms, service = c("uv", "dv")) %>%
  select(site_no, station_nm, data_type_cd:stat_cd, begin_date:count_nu) %>%
  # Keep unit values or daily means only (not max/min, etc).
  # Drop data that ends before the study period
  filter(is.na(stat_cd) | stat_cd == "00003",
         end_date > as.Date("2017-06-01")) %>%
  arrange(site_no, parm_cd, data_type_cd) %>% as.data.frame()

dailys <- filter(avail, data_type_cd == "dv") %>%
  select(site_no, parm_cd)
daily_data <- lapply(unique(dailys$parm_cd), function(cd) {
  cd_nm <- water_parms_nm[which(water_parms == cd)]
  stns <- filter(dailys, parm_cd == cd) %>% 
    pull(site_no) %>% unique()
  tmp <- readNWISdv(siteNumbers = stns,
                    parameterCd = cd,
                    startDate = "2017-06-01")
  names(tmp) <- gsub(paste0("X_", cd, "_00003"), paste0("dly_", cd_nm), names(tmp))
  tmp
})
daily_data <- reduce(daily_data, full_join, by = c("agency_cd", "site_no", "Date"))
saveRDS(daily_data, file = "./Output/usgs_daily.rds")

units <- filter(avail, data_type_cd == "uv", parm_cd != "") %>%
  select(site_no, parm_cd)
unit_data <- lapply(unique(units$parm_cd), function(cd) {
  cd_nm <- water_parms_nm[which(water_parms == cd)]
  stns <- filter(units, parm_cd == cd) %>% 
    pull(site_no) %>% unique()
  tmp <- readNWISuv(siteNumbers = stns,
                    parameterCd = cd,
                    startDate = "2017-06-01",
                    tz = "America/Chicago")
  names(tmp) <- gsub(paste0("X_", cd, "_00000"), cd_nm, names(tmp))
  tmp
})
unit_data <- reduce(unit_data, full_join, by = c("agency_cd", "site_no", "dateTime", "tz_cd"))
saveRDS(unit_data, file = "./Output/usgs_unit.rds")

# CRMS station data
CRMS <- list.files("./Data/CRMS", pattern = "*.zip$", full.names = TRUE)
crms <- lapply(CRMS, function(i) {
  tmp <- read_csv(i, locale = locale(tz = "America/Chicago", encoding = "latin1")) %>%
    select(site_no = `Station ID`, date = `Date (mm/dd/yyyy)`, time = `Time (hh:mm:ss)`,
           tempC = `Adjusted Water Temperature (°C)`,
           salinity_ppt = `Adjusted Salinity (ppt)`) %>%
    mutate(dateTime = mdy_hms(paste(date, time), tz = "Etc/GMT-6")) %>%
    select(-date, -time)
  attr(tmp$dateTime, "tzone") <- "America/Chicago"
  tmp
})
crms <- bind_rows(crms) %>%
  filter(dateTime >= as.Date("2017-06-01", tz = "America/Chicago"))
saveRDS(crms, file = "./Output/crms_data.rds")
