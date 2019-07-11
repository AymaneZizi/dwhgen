library(nasapower)
library(rio)
library(jsonlite)
library(dplyr)
library(tidyr)
library(naniar)
library(lubridate)
# install_formats() instalar formatos adicionales de rio
library(purrr)
library(anytime)

read_ideam <- function(filename, type){
  
  switch(type, "rjson" = {
    
    x = rjson::fromJSON(file = filename, simplify = TRUE)
    x = tibble::as_tibble(x) %>%
      mutate_at(.vars = c('Start Date', 'End Date', 'init_date', 'last_date'), .funs = anydate) 
    
  }, "rio" = {
    
    x = rio::import(filename, setclass =  "tibble")
    
  }
  )
  return(x)
}

write_ideam <- function(x, filename, type = "rio"){
  
  switch(type, "rjson" = {
    
    x = x %>%
      dplyr::mutate_if(.predicate = is.Date, .funs = as.character) %>%
      rjson::toJSON()
    
    write_lines(x, path = filename)
    rm(x)
    gc()
    gc(reset = T)
    
    write_json(x, path = filename)
    
  }, "rio" = {
    
    rio::export(x, filename)
    gc(reset = T)
  }
  )
}

x <- read_ideam("fill_prec_IDEAM_1981_2010", type = "rio")
y <- read_ideam("IDEAM_qc_and_chirps.json", type = "rio")
coords <- dplyr::select(y, Code, lat, lon, qc_climate)

z <- left_join(x, coords, by = "Code")
coords_precip <- dplyr::select(z, lat, lon)

library(future)

get_nasapower <- function(lat, lon, from, to){
  
  # daily_region_ag <- get_power(community = "AG",
  #                              lonlat = c(-76.1, 5.88),
  #                              pars = c("PRECTOT" ,
  #                                       "ALLSKY_SFC_SW_DWN", 
  #                                       "RH2M",
  #                                       "T2M_MAX",
  #                                       "T2M_MIN",
  #                                       "WS2M"),
  #                              dates = c("1981-01-01", "2010-12-31"),
  #                              temporal_average = "DAILY")
  
  daily_region_ag <- get_power(community = "AG",
                               lonlat = c(lon, lat),
                               pars = c("PRECTOT" ,
                                        "ALLSKY_SFC_SW_DWN", 
                                        "RH2M",
                                        "T2M_MAX",
                                        "T2M_MIN",
                                        "WS2M"),
                               dates = c(from, to),
                               temporal_average = "DAILY")
  
}
library(furrr)
options(future.globals.maxSize= 891289600)
plan(future::cluster, workers = 8)


distribute_load <- function(x, n) {
  assertthat::assert_that(assertthat::is.count(x),
                          assertthat::is.count(n),
                          isTRUE(x > 0),
                          isTRUE(n > 0))
  if (n == 1) {
    i <- list(seq_len(x))
  } else if (x <= n) {
    i <- as.list(seq_len(x))
  } else {
    j <- as.integer(floor(seq(1, n + 1 , length.out = x + 1)))
    i <- list()
    for (k in seq_len(n)) {
      i[[k]] <- which(j == k)
    }
  }
  
  i
}

l = distribute_load(x = dim(z)[1], n = 32)

w <- purrr::map(.x = l, .f = function(x, y){
  
  y %>% 
    filter(row_number() %in% x )
}, y = x)


plan(sequential)
plan(list(tweak(multisession, workers = 2), tweak(multisession, workers = 4)))

tic("nasapower")
jeferson <- furrr::future_map(.x = w, .f = function(x){
  
  get_information <- x %>% 
    mutate(nasapower = furrr::future_map2(.x = Latitude, .y = Longitude, .f = get_nasapower,
                                          from = "1985-01-01", "2019-02-28"))
  
  ## agregar un export csv
  
})
toc()

jeferson %>% 
  bind_rows() -> jeferson

jeferson <- dplyr::select(jeferson, Code, data, lat, lon, nasapower)
# w <- z %>% 
#   # filter(row_number() <= 32) %>% 
#   mutate(nasapower = furrr::future_map2(.x = lat, .y = lon, .f = get_nasapower,
#                                         from = "1981-01-01", "2010-12-31"))

write_ideam(jeferson, filename = "fill_prec_IDEAM_1981_2010_nasapower.json", type = "rio")

### descargar datos Nicaragua
library(readr)
x <-  read_csv("data/datos_arroz_final taller_2.csv") %>% 
  filter(Proyecto == "AVT")



l = distribute_load(x = dim(x)[1], n = 6)

w <- purrr::map(.x = l, .f = function(x,y){
  
  y %>% 
    filter(row_number() %in% x )
}, y = z)


jeferson %>% 
  bind_rows() -> jeferson

write_ideam(jeferson, filename = "nasapower_nicaragua.json", type = "rio")



install.packages("nasapower")
library(nasapower)
daily_region_ag <- get_power(community = "AG", # agroclimatology
                             lonlat = c(-76.1, 5.88),
                             pars = c("PRECTOT", ## precipitation
                                      "ALLSKY_SFC_SW_DWN", #srad
                                      "RH2M", ## relative humidity
                                      "T2M_MAX", ## Maximum temp
                                      "T2M_MIN", ## minimun temp
                                      "WS2M"),   ## wind
                             dates = c("2010-01-01", "2010-12-31"),
                             temporal_average = "DAILY")
library(readr) ## to load and save information
write_csv(daily_region_ag, path = "D:/nasapower.csv")



library(ggplot2)

ggplot() +
  geom_line(data = daily_region_ag, aes(x= YYYYMMDD,
                                        y = PRECTOT)) +
  theme_bw()

sse_i <- get_power(
    community = "SSE",
    lonlat = c(112.5, -55.5, 115.5, -50.5),
    dates = c("1984", "1985"),
    temporal_average = "INTERANNUAL",
    pars = c("T2M_MAX"))



library(raster)

daily_region_ag <- get_power(community = "AG",
                             lonlat = c(150.5, -28.5 , 153.5, -25.5),
                             pars = c("T2M_MAX"),
                             dates = c("1985-01-01", "1985-01-02"),
                             temporal_average = "DAILY")

#> Loading required package: sp
# Use split to create a list of data frames split by YYYYMMDD
daily_region_ag <- split(daily_region_ag, daily_region_ag$YYYYMMDD)

# Remove date information from data frame, list names will carry YYYYMMDD
daily_region_ag <-
  lapply(daily_region_ag, function(x)
    x[(!names(x) %in% c("YEAR", "MM", "DD", "DOY", "YYYYMMDD"))])

# Create a list of raster bricks from each YYYYMMDD data frame
raster_list <- lapply(daily_region_ag, rasterFromXYZ,
                      crs = "+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

stack_names <- paste0(names(raster_list), rep(c("T2M_maximum"), 2))

raster_stack <- stack(unlist(raster_list))
names(raster_stack) <- stack_names
plot(raster_stack)
writeRaster(raster_stack, 
            filename=paste0("D:/", names(raster_stack)),
            bylayer=TRUE,format="GTiff")












