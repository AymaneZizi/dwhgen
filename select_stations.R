library(rio)
library(jsonlite)
library(dplyr)
library(tidyr)
library(naniar)
library(lubridate)
# install_formats() instalar formatos adicionales de rio
library(purrr)
library(broom)
library(NbClust)

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


IDEAM_qc_chirps <- read_ideam(filename  = 'IDEAM_qc_and_chirps.json',
                              type = "rio")

# IDEAM_qc_chirps <- IDEAM_qc_chirps %>%
  # dplyr::select(-chirps)

# write_ideam(IDEAM_qc_chirps, filename = "IDEAM_qc_and_chirps.json", type = "rio")

### Selecting a "good" Wheather stations
## seleccionar años completos

## years with good information
IDEAM_qc_chirps <- IDEAM_qc_chirps %>%
  mutate(qc_climate = purrr::map(.x = qc_climate, .f = filter_date))

years_good <- function(x){
  
  # x <- IDEAM_qc_chirps %>%
  #   filter(Category %in% c("AM", "CP")) %>%
  #   filter(row_number()==1) %>%
  #   dplyr::select(qc_climate) %>%
  #   unnest()

  x %>%
    drop_na() %>%
    group_by(year(Date)) %>%
    count() %>%
    ungroup() %>%
    top_n(2, wt = "n") %>%
    rename(year_good = `year(Date)`, n_good = n)
  
}

wheather_station <- IDEAM_qc_chirps %>%
  filter(Category %in% c("AM", "CP"))

omit_na_df <- function(x){
  
  x <- x %>%
    drop_na()
    
  
}

s_without_na <- function(x){
  
  x <- x %>%
    mutate(good_station = purrr::map(.x = qc_climate, .f = omit_na_df)) %>%
    dplyr::select(Code, good_station) %>%
    mutate(n_row = purrr::map_dbl(.x = good_station, .f = nrow)) %>%
    filter(n_row > 0)

}

wheather_station <- s_without_na(wheather_station)

## write the wheather stations

write_ideam(x = wheather_station, filename = "good_stations.json", type = "rio")

## methods to fill wheather stations
make_na <- function(x){
  
  
  make_na <- function(x){
    
     x %>%
      dplyr::select(-Date, -day, -month, -year) %>%
      dplyr::select(-prec_chirps) %>%
      prodNA(noNA = 0.20) %>%
      rename(prec_na = prec_qc) 
  }
  library(missForest)
  ## generar NA por estacion climatica
  x.na <- IDEAM_qc_chirps %>%
    filter(Departament == "Antioquia") %>%
    filter(percent_na <= 5) %>%
    filter(row_number()<=10) %>%
    # filter(year(init_date) == 1980) %>%
    dplyr::select(Code, qc_climate) %>%
    mutate(data_na = purrr::map(.x = qc_climate, .f = make_na)) %>%
    dplyr::select(Code, data_na) %>%
    # filter(row_number()==1) %>%
    unnest()
  
    
  x <- IDEAM_qc_chirps %>%
    filter(Departament == "Antioquia") %>%
    filter(percent_na <= 5) %>%
    filter(row_number()<=10) %>%
    # filter(year(init_date) == 1980) %>%
    dplyr::select(Code, qc_climate) %>%
    unnest %>%
    dplyr::select(Date, day, month, year, prec_qc)
  
  x <- bind_cols(x, x.na) %>%
    filter(year(Date)==1981)
z = parallelPCA(lcounts, min.rank=0, value="n")

  x.na <- x %>%
    dplyr::select(Code, prec_na) %>%
    group_by(Code) %>%
    mutate(id = row_number()) %>%
    spread(Code, prec_na) %>%
    dplyr::select(-id)
  
  pc <- pca(x.na, nPcs=5, method="nlpca")
  
  imputed <- completeObs(pc) %>%
    as_tibble() %>%
    gather(Code, prec_imp)
  
  y <- IDEAM_qc_chirps %>%
    filter(Departament == "Antioquia") %>%
    filter(percent_na <= 5) %>%
    filter(row_number()<=10) %>%
    dplyr::select(Code, qc_climate) %>%
    unnest %>%
    filter(year(Date)==1981) %>%
    mutate(Code = as.character(Code))
  
  z = bind_cols(y, imputed) %>%
    mutate(diff = abs(prec_qc - prec_imp)) %>% 
    filter(Code == "11070020") %>%
    dplyr::select(Date, prec_qc, prec_imp) %>%
    # dplyr::select(prec_qc, prec_imp)
    bind_cols(bias_correction) %>%
    gather(key, value, -Date) %>%
    filter(key != "prec_imp")
  
  
    
  library(qmap)
  qm.fit <- fitQmapRQUANT(z %>% pull(prec_qc),z %>% pull(prec_imp),
                          qstep=0.1, nboot = 10, wet.day = TRUE)
  
  qm.fit$par$fitq
  doQmapRQUANT(modprecip[,2],qm.fit,type="linear")
  bias_correction <- doQmapRQUANT(z %>% pull(prec_imp),qm.fit,type="linear") %>%
    tibble(bias_correction = .)
  
  qm.b <- doQmapRQUANT(z %>% pull(prec_imp),qm.fit,type="linear2")  %>%
    tibble(bias_correction = .)
  
  qm1.fit <- fitQmap(z %>% pull(prec_qc),z %>% pull(prec_imp),
                     method="PTF",
                     transfun="expasympt",
                     cost="RSS",wett.day=TRUE)
  qm1 <- doQmap(z %>% pull(prec_imp),qm1.fit)
  
  
  
  qm2.fit <- fitQmap(sqrt(z %>% pull(prec_qc)),sqrt(abs(z %>% pull(prec_imp))),
                     method="DIST",qstep=0.001,
                     transfun="berngamma")
  
  qm2 <- doQmap(sqrt(abs(z %>% pull(prec_imp))),qm2.fit)^2
  
  grafico <- ggplot()+
    geom_boxplot(data = z, aes(x = Code, y = diff)) +
    theme_bw()
  grafico <- ggplot()+
    geom_line(data = z, aes(x = Date, y = value, color = key)) +
    theme_bw()
  
  z <- bind_cols(y, imputed) %>%
    mutate(diff = abs(prec_qc - prec_imp)) %>% 
    filter(Code == "11070020") %>%
    dplyr::select(Date, prec_qc, prec_imp, prec_chirps) %>%
    # dplyr::select(prec_qc, prec_imp)
    bind_cols(bias_correction) %>%
    gather(key, value, -Date) %>%
    filter(key != "prec_imp")
    
  z <- z %>%
    drop_na()
  cor(z$prec_qc, z$bias_correction, method = "kendall")
  ggsave('imp_prec.pdf', width = 15, height = 10, dpi = 300)
  # inner_join(imputed, y, by = "Code")
  library(pcaMethods)
  z = parallelPCA(x, min.rank=0, value="n")
  BiocManager::install("scran")

  pc <- pca(x, nPcs=15, method="nlpca",  maxSteps = 500)
  imputed <- completeObs(pc)
  
  BiocManager::install("scran")
  library(scran)
  imputed <- imputed %>%
    as_tibble() %>%
    mutate(id = row_number()) %>%
    gather(Code, prec_imp, -id) %>%
    dplyr::select(-id) %>%
    nest(-Code)
    
  antioquia <- IDEAM_qc_chirps %>%
    filter(Departament == "Antioquia") %>%
    filter(percent_na <= 5) %>%
    dplyr::select(Code, qc_climate) %>%
    mutate(Code = as.character(Code))
  
  x= left_join(antioquia, imputed, by = "Code")
  # x <- wheather_station %>%
  #   dplyr::select(good_station) %>%
  #   filter(row_number()==1) %>%
  #   unnest

  x.na <- x %>%
    dplyr::select(-Date, -day, -month, -year) %>%
    prodNA(noNA = 0.25) 
  
  prec <- x.na %>%
    dplyr::select(prec_qc) %>%
    drop_na()
    
  
  kclust <- kmeans(prec, centers = 6) 
  prec_clust <- augment(kclust, prec)
  tidy(kclust)
  glance(kclust)
  prec_clust %>%
    group_by(.cluster) %>%
    summarise(avg = mean(prec_qc))
  library(pcaMethods)
  
  data(metaboliteData)
  mD <- metaboliteData
  pc <- pca(x, nPcs=3, method="nlpca", maxSteps = 400)
  imputed <- completeObs(pc)
  
  x <- x %>%
    dplyr::select(-Date)
  nlpca(x, nPcs=3)
  prec_clust %>%
    filter(.cluster ==6 )
  x %>%
    filter(prec_qc >= 40)
  nc <- NbClust(prec, min.nc=2, max.nc=9, method="kmeans")
  NbClust(prec,distance = "euclidean", 
          min.nc=2, max.nc=10, method = "kmeans", 
          index = "all", alphaBeale = 0.1)
  library(cluster)
  gap_stat <- clusGap(prec, FUN = kmeans, nstart = 25,
                      K.max = 10, B = 50)
  print(gap_stat, method = "firstmax")
  gap_stat <- clusGap(prec, FUN = kmeans, K.max = 10, spaceH0 = "original")
  gap_stat %>%
    as_tibble
  k <- maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"], method="firstmax")
  kmean_withinss <- function(k) {
    cluster <- kmeans(prec, k)
    return (cluster$tot.withinss)
  }
  
  wss <- sapply(2:9, kmean_withinss)
  elbow <-data.frame(2:9, wss)
}

unlink("/home/jeison/R/x86_64-pc-linux-gnu-library/3.4/00LOCK-RcppEigen", recursive = TRUE) 
