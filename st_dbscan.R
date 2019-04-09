########################################################################
# ST-DBSCAN : An algorithm for clustering spatial-temporal data        #
# (Birant and Kut, 2006)                                 #   
# Application on a trajectory                           #
########################################################################


########################################################################
# INPUTS :                                                             #
# traj = traj gps (x, y and time)                                      #
# eps = distance minimum for longitude and latitude                    #
# eps2 =  distance minimum for date                                    #
# minpts = number of points to consider a cluster                      #
########################################################################

# data = spatio-temporal data 
# x = data longitude 
# y = data latitude 
# time = data timestamps 
# 
# eps = distance minimum for longitude and latitude 
# eps2 = temporal window 
# minpts = number of points to consider a cluster 
# cldensity = TRUE or FALSE to display the number of points reachables for every point within a cluster 
# 
# OUTPUTS :

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

traj <- IDEAM_qc_chirps %>%
  # filter(row_number() <= 5) %>%
  geometry_to_lonlat %>%
  dplyr::select(lat,lon) #%>%
  # unnest() %>%
  # dplyr::select(-year,-month,-year)

x <- dplyr::select(traj, lon) %>%
  pull

y <- dplyr::select(traj, lat) %>%
  pull

time <- seq(ymd('1981-01-01'),ymd('2010-12-31'),by='day')
time <- strptime(time, "%Y-%m-%d")
  # dplyr::select(traj, Date) %>%
  # mutate(Date = anytime::anydate(Date)) 
  
stdbscan = function (traj, 
                     x, 
                     y, 
                     time, 
                     eps = 0.3, 
                     eps2 = 121190.58, 
                     minpts = 10, 
                     cldensity = TRUE) { 
  
  countmode = 1:length(x)
  seeds = TRUE
  
  data_spatial <- as.matrix(dist(cbind(y, x)))
  
  data_spatial <- coords %>%
    group_split(cluster) %>%
    purrr::map(.x = ., .f = function(x){
      res <- x %>%
        dplyr::select(-cluster) %>%
        dbscan(eps = .33, minPts = 10)
      
      clust <- res$cluster %>%
        tibble::enframe(value = "cluster_sp")  %>%
        dplyr::select(cluster_sp)
      
      x <- x %>%
        bind_cols(clust) %>%
        mutate(cluster_inside = glue::glue("cluster_{cluster}_{cluster_sp}"))
        
        return(x)
    }) %>%
    bind_rows() %>%
    dplyr::select(cluster_inside)
  
  data_temporal<- as.matrix(pc.l2@distmat)
  dd <- fuse(data_spatial, data_temporal, weights = c(0.5, 0.5))
  # hc <- hclust(dd, "ave")
  res <- NbClust(diss=dd, distance = NULL, min.nc=20, max.nc=50, 
               method = "ward.D2", index = "silhouette")
  
  n <- nrow(data_spatial)
  
  classn <- cv <- integer(n)
  isseed <- logical(n)
  cn <- integer(1)
  
  for (i in 1:n) {
    if (i %in% countmode)
      #cat("Processing point ", i, " of ", n, ".\n")
      unclass <- (1:n)[cv < 1]
    
    if (cv[i] == 0) {
      reachables <- intersect(unclass[data_spatial[i, unclass] <= eps],  unclass[data_temporal[i, unclass] <=eps2])
      if (length(reachables) + classn[i] < minpts)
        cv[i] <- (-1)                    
      else {
        cn <- cn + 1                   
        cv[i] <- cn
        isseed[i] <- TRUE
        reachables <- setdiff(reachables, i)
        unclass <- setdiff(unclass, i)       
        classn[reachables] <- classn[reachables] + 1
        while (length(reachables)) {
          cv[reachables] <- cn           
          ap <- reachables                           
          reachables <- integer()
          
          for (i2 in seq(along = ap)) {
            j <- ap[i2]
            
            jreachables <- intersect(unclass[data_spatial[j, unclass] <= eps], unclass[data_temporal[j, unclass] <= eps2])
            
            if (length(jreachables) + classn[j] >= minpts) {
              isseed[j] <- TRUE
              cv[jreachables[cv[jreachables] < 0]] <- cn
              reachables <- union(reachables, jreachables[cv[jreachables] == 0])
            }
            classn[jreachables] <- classn[jreachables] + 1
            unclass <- setdiff(unclass, j)
          }
        }
      }
    }
    if (!length(unclass))
      break
    
  }
  
  
  if (any(cv == (-1))) {
    cv[cv == (-1)] <- 0
  }
  out <- list(cluster = cv, eps = eps, minpts = minpts, density = classn)
  rm(classn)
  if (seeds && cn > 0) {
    out$isseed <- isseed
  }
  class(out) <- "stdbscan"
  return(out)
}