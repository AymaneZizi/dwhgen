## time series clustering on rainfall

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
filter_date <- function(x){
  
  filter(x, year >= 1981, year <= 2010) %>%
    mutate(Date = lubridate::ymd(Date))
}

climate_var <- function(x, var = "prec_qc"){
  
  x %>%
    dplyr::select(!!var)
}
library(rio)
library(jsonlite)
# install.packages("dtwclust")
library(dplyr)
library(dtwclust)
# install.packages("nngeo")
library(nngeo)
library(tidyr)
# https://github.com/topics/time-series-clustering?l=r
# https://www.r-bloggers.com/tsrepr-use-case-clustering-time-series-representations-in-r/
# https://github.com/PetoLau/TSrepr
# https://github.com/lnferreira/time_series_clustering_via_community_detection/tree/gh-pages/R
# https://hal.archives-ouvertes.fr/tel-01394280v2/document
# https://pdfs.semanticscholar.org/e181/49ae37583332f06e8288726a786077934ed9.pdf
# http://sisifospage.tech/2017-05-15-time-series-clustering-pulsi.html
# https://rpubs.com/Hailstone/346625
# https://cran.r-project.org/web/packages/HiClimR/HiClimR.pdf
# https://stackoverflow.com/questions/53440909/polygon-from-cluster-of-lat-long-points-in-r
# https://gis.stackexchange.com/questions/249762/calculating-distances-between-two-geometry-columns-using-r
# https://cran.r-project.org/web/packages/nngeo/vignettes/intro.pdf
# https://www.youtube.com/watch?v=Tf2CSDDTezI
# http://rstudio-pubs-static.s3.amazonaws.com/398402_abe1a0343a4e4e03977de8f3791e96bb.html
# https://rpubs.com/janoskaz/10351
IDEAM_qc_chirps <- read_ideam(filename  = 'IDEAM_qc_and_chirps.json',
                              type = "rio")



IDEAM_qc_chirps <- IDEAM_qc_chirps %>%
  mutate(qc_climate = purrr::map(.x = qc_climate, .f = filter_date))

IDEAM_qc_chirps <- IDEAM_qc_chirps %>%
  sf::st_as_sf(coords = c("lon", "lat")) %>%
  st_set_crs( "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>%
  st_transform("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

near_stations <- IDEAM_qc_chirps %>% 
  dplyr::select(qc_climate)

# nn %>%
#   group_by(Code)
# near_stations
library(here)
install.packages("HiClimR")
colombia <- read_sf(dsn = here('data', 'municipios_wgs84.shp')) %>%
  st_transform("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>%
  group_by(NOM_DEPART) %>%
  summarize()

nn =st_nn(IDEAM_qc_chirps, near_stations, progress = FALSE,  k = 100, maxdist = 5000)
l =st_connect(IDEAM_qc_chirps, near_stations, ids = nn, progress = FALSE)

nn =st_join(IDEAM_qc_chirps, near_stations, join = st_nn, k = 100, maxdist = 5000)



ggplot() +
  geom_sf(data = IDEAM_qc_chirps, aes(color = "red"))+
  # geom_sf(data = near_stations, aes(color = "black"))+
  geom_sf(data =l, fill = NA, color = gray(.01) ,size=0.125)+
  geom_sf(data = colombia, fill = NA, color = gray(.01) ,size=0.125) +
  theme_bw()

ggsave('mapa.pdf', width = 15, height = 10, dpi = 300)
plot(st_geometry(IDEAM_qc_chirps[1:5, ]), col = "darkgrey")
plot(st_geometry(l), add = TRUE)
plot(st_geometry(IDEAM_qc_chirps[1:5, ]), col = "red", add = TRUE)
text(st_coordinates(IDEAM_qc_chirps[1:5, ])[, 1],st_coordinates(IDEAM_qc_chirps[1:5, ])[, 2],1:3, col = "red", pos = 4)
text(st_coordinates(IDEAM_qc_chirps[1:5, ])[, 1],st_coordinates(towns[1:5, ])[, 2],1:5, pos = 4)


x <- IDEAM_qc_chirps %>%
  # filter(percent_na <= 25)
  geometry_to_lonlat %>%
  mutate(rainfall = purrr::map(.x = qc_climate, 
                               .f = climate_var, var = "prec_qc")) %>% 
  # filter(row_number() <= 10) %>%
  dplyr::select(Code, rainfall) %>%
  # group_split(Code, keep = F)
  unnest() %>%
  group_split(Code, keep = F) %>% 
  purrr::map(.x = ., .f = function(x){pull(x) %>% na.omit})
  # group_by(Code) %>%
  # mutate(id = row_number()) %>% 
  # spread(Code, prec_qc)

series <- reinterpolate(x, new.length = max(lengths(x)))

series <- series
labels <- IDEAM_qc_chirps %>%
  # filter(percent_na <= 25)
  geometry_to_lonlat %>%
  mutate(rainfall = purrr::map(.x = qc_climate, 
                               .f = climate_var, var = "prec_qc")) %>% 
  # filter(row_number() <= 10) %>%
  dplyr::select(Code) %>%
  # group_split(Code, keep = F)
  unnest() %>%
  group_keys(Code) %>%
  pull

pc.l2 <- tsclust(series, k = 20L, 
                 centroid = "pam",seed = 3247, 
                 trace = TRUE,
                 distance = "dtw_basic")
pc <- tsclust(x, type = "partitional", k = 5, 
              distance = "dtw_basic", centroid = "pam", 
              seed = 3247L, trace = TRUE)
pc@cluster 

x <- IDEAM_qc_chirps %>%
  # filter(percent_na <= 25)
  geometry_to_lonlat %>%
  mutate(rainfall = purrr::map(.x = qc_climate, 
                               .f = climate_var, var = "prec_qc")) %>% 
  filter(row_number() %in% 3:6) %>%
  dplyr::select(Code, qc_climate) %>%
  unnest %>%
  arrange(Date)
  group_by(Code, year, month) %>%
  summarise(prec_acum = sum(prec_qc, na.rm = T)) %>%
  arrange(month)
  
ggplot()+
  geom_line(data = x, aes(x = Date, y = prec_qc, color= as.factor(Code)))

ggsave('lines.pdf', width = 15, height = 10, dpi = 300)

pc <- dtwclust(x, k = 4L,
                  distance = "dtw_basic", centroid = "pam",
                  seed = 3247, control = list(trace = TRUE,
                                              nrep = 10L))