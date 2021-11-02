library(daymetr)
library(lubridate)
require(dplyr)

mac_daymet_list <- download_daymet(site = "Maricopa Agricultural Center",
                                   lat = 33.068941,
                                   lon =  -111.972244,
                                   start = 2000,
                                   end = 2020, internal = TRUE)


# rename variables, create a date column
mac_daymet <- mac_daymet_list$data %>% 
  transmute(date = ymd(paste0(year, '01-01'))+ days(yday) -1, 
            precip = prcp..mm.day., 
            tmax = tmax..deg.c.,
            tmin = tmin..deg.c.,
            tmean = (tmax + tmin) / 2,
            trange = tmax - tmin,
            srad = srad..W.m.2.,
            vpd = vp..Pa.)%>%
  select(date, precip, tmean, srad, vpd)

head(mac_daymet)

tmean.ts <- ts(mac_daymet$tmean, 
               start = c(2000, 1), 
               end = c(2020, 365), 
               deltat = 1/365)

mac_ts <- ts(mac_daymet, 
             start = c(2000, 1), 
             end = c(2020, 365), 
             deltat = 1/365)

plot(tmean.ts, ylab = "Daily mean T", xlab = "Year")

lag.plot(tmean.ts, set.lags = c(1, 10, 100, 180, 360))

tmean_mo <- mac_daymet %>% 
  mutate(year = year(date), month = month(date)) %>% 
  group_by(year, month) %>% 
  summarise(tmean = mean(tmean), .groups = 'keep') %>% 
  ungroup() %>% 
  select(tmean)

tmean.mo.ts <- ts(tmean_mo, start = c(2000, 1), end = c(2020, 12), frequency = 12)

lag.plot(tmean.mo.ts, lags = 12)

plot(acf(tmean.mo.ts))
