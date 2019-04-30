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

geometry_to_lonlat <- function(x) {
  if (any(sf::st_geometry_type(x) != "POINT")) {
    stop("Selecting non-points isn't implemented.")
  }
  coord_df <- sf::st_transform(x, sf::st_crs("+proj=longlat +datum=WGS84")) %>%
    sf::st_coordinates() %>%
    dplyr::as_tibble() %>%
    dplyr::select(X, Y) %>%
    dplyr::rename(lon = X, lat = Y)
  out <- sf::st_set_geometry(x, NULL) %>%
    dplyr::bind_cols(coord_df)
  return(out)
}

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

library(rio)
library(jsonlite)
# install.packages("dtwclust")
library(dplyr)
library(dtwclust)
# install.packages("nngeo")
library(nngeo)
library(tidyr)
library(forcats)
library(lubridate)
library(here)
library(dbscan)
library(glue)

## 819 estaciones para precipitacion
IDEAM_qc_chirps <- read_ideam(filename  = 'IDEAM_qc_and_chirps.json',
                              type = "rio")

## se filtran las series para los años 1981-2010

IDEAM_qc_chirps <- IDEAM_qc_chirps %>%
  mutate(qc_climate = purrr::map(.x = qc_climate, .f = filter_date))

## convertir estaciones de IDEAM a un objeto espacial

IDEAM_qc_chirps <- IDEAM_qc_chirps %>%
  sf::st_as_sf(coords = c("lon", "lat")) %>%
  st_set_crs( "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>%
  st_transform("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

## Cada estacion climatica como una lista (solo precipitacion)
x <- IDEAM_qc_chirps %>%
  geometry_to_lonlat %>% 
  mutate(rainfall = purrr::map(.x = qc_climate, 
                               .f = climate_var, var = "prec_qc")) %>% 
  dplyr::select(Code, rainfall) %>%
  unnest() %>%
  group_split(Code, keep = F) %>% 
  purrr::map(.x = ., .f = function(x){
    pull(x) %>% na.omit
    }
    )

## interpolar para tener el mismo tamaño en todas las series de tiempo
series <- reinterpolate(x, new.length = max(lengths(x)))

series <- series
labels <- IDEAM_qc_chirps %>%
  mutate(rainfall = purrr::map(.x = qc_climate, 
                               .f = climate_var, var = "prec_qc")) %>% 
  dplyr::select(Code) %>%
  unnest() %>%
  group_keys(Code) #%>%
  # pull()

## realizar los cluster para 20 medoides (averiguar bien que son medoides)
pc.l2 <- tsclust(series, k = 20L, 
                 centroid = "pam",seed = 3247, 
                 trace = TRUE,
                 distance = "dtw_basic")

id_cluster <- pc.l2@cluster %>%
  tibble::enframe(value = "cluster_time") %>%
  dplyr::select(cluster_time) %>%
  mutate(cluster_time = as_factor(cluster_time))

length_cluster <- pc.l2@clusinfo %>%
  as_tibble() %>%
  dplyr::select(size) %>%
  mutate(cluster = as_factor(row_number()))

time_series_cl <- bind_cols(labels, id_cluster) %>%
  left_join(length_cluster, by = c("cluster_time" = "cluster"))

write_csv(time_series_cl, "cluster_dtw_prec.csv")

## cluster para las zonas espaciales dentro de cada cluster determinado por el paso anterior

## extraer las coordenadas para realizar el cluster


cluster_sf <- IDEAM_qc_chirps %>%
  left_join(length_cluster, by = "cluster") %>%
  filter(size > 10 )

coords <- cluster_sf %>%
  geometry_to_lonlat %>%
  dplyr::select(lat,lon, Code, cluster)


make_stdbscan <- function(x, vars = c("lat", "lon")){
  
  # x <- coords %>%
  #   group_split(cluster) %>% 
  #   purrr::pluck(1)
  res <- x %>%
    dplyr::select(!!vars) %>%
    dbscan::dbscan(eps = .33, minPts = 10)
  
  clust <- res$cluster %>%
    tibble::enframe(value = "cluster_sp")  %>%
    dplyr::select(cluster_sp)
  
  x <- x %>%
    bind_cols(clust) %>%
    filter(cluster_sp >0) %>%
    mutate(cluster_inside = glue::glue("time_{cluster}_sp_{cluster_sp}")) %>%
    mutate(cluster_inside = as_factor(cluster_inside))
  
  return(x)
}

data_spatial <- coords %>%
  group_split(cluster) %>% 
  purrr::map(.x = ., .f = make_stdbscan) %>%
  bind_rows() %>%
  dplyr::select(Code, cluster_inside)

# data_spatial <- coords %>%
#   group_split(cluster) %>%
#   purrr::map(.x = ., .f = function(x){
#     res <- x %>%
#       dplyr::select(-cluster) %>%
#       dbscan(eps = .33, minPts = 8)
# 
#     clust <- res$cluster %>%
#       tibble::enframe(value = "cluster_sp")  %>%
#       dplyr::select(cluster_sp)
# 
#     x <- x %>%
#       bind_cols(clust) %>%
#       mutate(cluster_inside = glue::glue("cluster_{cluster}_{cluster_sp}"))
# 
#     return(x)
#   }) %>%
#   bind_rows() %>%
#   dplyr::select(cluster_inside)



cluster_sf <- cluster_sf %>%
  inner_join(data_spatial, by = "Code")

size_cluster <- cluster_sf %>%
  as_tibble() %>%
  count(cluster_inside) 

cluster_sf <- cluster_sf %>%
  left_join(size_cluster, by = "cluster_inside") %>%
  filter(n > 10)


write_ideam(cluster_sf %>% geometry_to_lonlat, filename = "IDEAM_clustering_time_spatial_1981-2010.json", type = "rio")

##################################################33


near_stations <- IDEAM_qc_chirps %>% 
  dplyr::select(qc_climate)

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





id_cluster <- pc.l2@cluster %>%
  tibble::enframe(value = "cluster") %>%
  dplyr::select(cluster)

length_cluster <- pc.l2@clusinfo %>%
  as_tibble() %>%
  dplyr::select(size) %>%
  mutate(cluster = as_factor(row_number()))

write_csv(id_cluster, "id_cluster.csv")
write_csv(length_cluster, "length_cluster.csv")
IDEAM_qc_chirps <- bind_cols(IDEAM_qc_chirps, id_cluster)

IDEAM_qc_chirps <- IDEAM_qc_chirps %>%
  # sf::st_as_sf(coords = c("LatitudeDD", "LongitudeDD")) %>%
  sf::st_as_sf(coords = c("lon", "lat")) %>%
  st_set_crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>%
  st_transform("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

IDEAM_qc_chirps <- IDEAM_qc_chirps %>%
  mutate(cluster = as_factor(cluster))

cluster_sf <- IDEAM_qc_chirps %>%
  left_join(length_cluster, by = "cluster") %>%
  # dplyr::select(size, cluster, id)  %>%
  filter(size > 10 )

coords <- cluster_sf %>%
  geometry_to_lonlat %>%
  dplyr::select(lat,lon, cluster)

res <- dbscan(coords, eps = .33, minPts = 10)
clust <- res$cluster %>%
  tibble::enframe(value = "cluster") %>%
  dplyr::select(cluster)

cluster_sf <- cluster_sf %>%
  bind_cols(data_spatial)

size_cluster <- cluster_sf %>%
  as_tibble() %>%
  count(cluster_inside) 

cluster_sf <- cluster_sf %>%
  left_join(size_cluster, by = "cluster_inside") %>%
  filter(n > 10)

clust_Dist <- tibble::enframe(res$Best.partition, value = "cluster")%>%
  dplyr::select(cluster)

clust <- clust_Dist
cluster_sf <- bind_cols(cluster_sf, clust) 
cluster_sf <- cluster_sf %>%
  # filter(cluster1>0) %>%
  mutate(cluster1 =as_factor(cluster1))

cluster <- cluster_sf %>%
  geometry_to_lonlat %>%
  dplyr::select(cluster, cluster1) %>%
  mutate_all(~as.numeric(.)) %>%
  mutate(cluster_s_t = cluster*cluster1) 

cluster_sf <- bind_cols(cluster_sf, cluster) 

cluster_sf <- cluster_sf %>%
  mutate(cluster_s_t = as_factor(cluster_s_t))
cluster_sf %>%
  as_tibble() %>%
  mutate(Elevation = as.numeric(Elevation)) %>%
  group_by(cluster) %>%
  summarise(avg_elev = mean(Elevation), sd_elev = sd(Elevation))

size_cluster <- cluster_sf %>%
  as_tibble() %>%
  count(cluster1) 

cluster_sf <- cluster_sf %>%
  left_join(size_cluster, by = "cluster1") %>%
  filter(n > 10)
  
# cluster_sf <- IDEAM_qc_chirps %>%
#   filter(cluster != "1")

library(here)
colombia <- read_sf(dsn = here('data', 'municipios_wgs84.shp')) %>%
  st_transform("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>%
  group_by(NOM_DEPART) %>%
  summarize()

y = ggplot()+
  geom_sf(data = cluster_sf, aes(color = cluster_inside, fill = cluster_inside)) +
  geom_sf(data = colombia, fill = NA, color = gray(.01) ,size=0.125) +
  theme_bw() +
  # ggthemes::theme_map() +
  # guides(color = FALSE) +
  scale_color_viridis_d() +
  scale_fill_viridis_d()
  # ggtitle(label = "28 Weather Stations with Quality Control (with, tmax, tmin, prec, sbright)", 
          # subtitle = "1981-2010")
ggsave('cluster_inside_2.pdf', width = 15, height =  10, dpi = 300) 
# pc <- tsclust(x, type = "partitional", k = 5, 
#               distance = "dtw_basic", centroid = "pam", 
#               seed = 3247L, trace = TRUE)


wheater_filter <- IDEAM_qc_chirps # %>%
  # filter(cluster > 1)

x = wheater_filter %>%
  filter(cluster == 2)


x <- x %>%
  dplyr::select(Code, qc_climate) %>%
  unnest() %>%
  group_by(Code) %>%
  mutate(id = row_number()) %>%
  spread(Code, prec_qc) %>%
  dplyr::select(-id, -Date, -month, -year, -day, -prec_chirps)

pc <- pca(x.na, nPcs=5, method="nlpca")

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