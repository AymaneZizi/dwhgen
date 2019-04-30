### implementacion de NLPCA (componentes principales no lineales)

library(pcaMethods)


cluster_sf <- cluster_sf %>% 
  geometry_to_lonlat

library(pcaMethods)

QM <- function(x){
  obs <- pull(x, prec_qc_qm)
  mod <- pull(x, prec_imp_qm)
  qm_fit <- fitQmapRQUANT(obs, mod,
                          qstep = 0.01, nboot = 1, wet.day = 0.2)
  
  qm_corrected <- doQmapRQUANT(mod, qm_fit,type="tricub") %>%
    tibble(prec_bias_correction = .) %>% 
    mutate(prec_bias_correction = if_else(prec_bias_correction <= 0.2, 0, prec_bias_correction))
  
  return(qm_corrected)
}

make_nlpca <- function(x){
  
  # x <- cluster_sf %>%
    # group_split(cluster_inside) %>%
    # purrr::pluck(1) # %>%
  
  
    # # filter(row_number()==1) %>%
    # dplyr::select(Code, qc_climate) %>%
    # unnest()
  x <- x %>% 
  dplyr::select(Code, qc_climate) %>%
    unnest()
  
  dates <- x %>%
    dplyr::select(Code, Date, year, month, day)
  
  prec_obs <- dplyr::select(x, prec_qc)
  
  x <- x %>%
    group_by(Code) %>%
    mutate(id = row_number()) %>%
    dplyr::select(Code, id, prec_qc) %>% 
    spread(Code, prec_qc) %>%
    dplyr::select(-id)
  
  # n_Pcs <- dim()
  pc <- pca(x, nPcs = 3, method="nlpca",  maxSteps = 50)
  

  imputed <- completeObs(pc) %>% 
    as_tibble() %>%
    gather(Code, prec_imp) %>%
    bind_cols(prec_obs) %>%
    mutate(prec_imp_qm = if_else(prec_imp <= 0.2, 0, prec_imp),
           prec_qc_qm = if_else(prec_qc <= 0.2, 0, prec_qc)) %>%
    mutate(prec_imp_qm = if_else(prec_imp_qm == 0, runif(n(), 0, 0.2), prec_imp_qm),
           prec_qc_qm = if_else(prec_qc_qm == 0, runif(n(), 0, 0.2), prec_qc_qm)) %>%
    nest(-Code) %>%
    mutate(QM_processing = purrr::map(.x = data, .f = QM)) 
  
  
  dates <- dates %>%
    nest(-Code) # %>%
    # filter(row_number()==1)
  
  imputed <- imputed %>%
    # filter(row_number()==1) %>%
    # unnest() %>%
    mutate(Code = as.integer(Code)) %>% 
    left_join(x = ., y = dates, by = "Code")  %>%
    unnest() %>%
    nest(-Code)
  
  return(imputed)
    # dplyr::select(Date, prec_qc, prec_bias_correction) %>% 
    # gather(key, value, -Date)
    

  
  # grafico <- ggplot() +
  #   geom_line(data = prueba, aes(x = Date, y = value, color = key)) +
  #   theme_bw()
  #   
  # ggsave('comp_obs_mod.pdf', width = 15, height = 10, dpi = 300)

  # library(qmap)
  # 
  # qm.fit <- fitQmapRQUANT(prec_qc_qm, prec_imp_qm,
  #                         qstep = 0.99, nboot = 10, wet.day = 0.2)
  
}
cluster_sf <- bind_rows(cluster_sf)
x <- cluster_sf %>% 
  # filter(row_number() %in% 1:5) %>% 
  group_split(cluster_inside)
# cluster_sf <- cluster_sf %>%
#   group_split(cluster_inside)

library(furrr)
library(future)
x <- pluck(x, 1)
options(future.globals.maxSize= 891289600)
plan(future::cluster, workers = 15)

x <-  furrr::future_map(.x = x, make_nlpca)

for(i in 1:)
2 * prod(dim(x))
x <- cluster_sf %>%
  group_split(cluster_inside) %>%
  purrr::pluck(1) %>%
  # filter(row_number()==1) %>%
  dplyr::select(Code, qc_climate) %>%
  unnest()

x.na <- x %>%
  dplyr::select(Code, prec_na) %>%
  group_by(Code) %>%
  mutate(id = row_number()) %>%
  spread(Code, prec_na) %>%
  dplyr::select(-id)

pc <- pca(x.na, nPcs=5, method="nlpca")



  
  mutate(imp_prec = furrr::future_map(.x = ))

simPopN <- data.frame(slope = 0.00237, 
                      intercept=115.767,
                      sigma = 2:6) %>%
  crossing(n=10:100) 
  
