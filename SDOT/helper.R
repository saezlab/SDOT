###################################################################################################################
########################### SDOT: Spatial Cell Type Deconvolution by Optimal Transport ############################
## Preliminary work. Under review by the International Conferenceon Machine Learning (ICML). Do not distribute. ###
###################################################################################################################

library(igraph)

find_spatial_neighborhood_radious <- function(coordinates, neighbors = 8, return_dists = FALSE)
{
  N <- nrow(coordinates)
  
  if(N > 5000 & !return_dists)
  {
    center <- as.data.frame(t(colMeans(coordinates)))
    center_dists <- fields::rdist(center, coordinates)
    
    centered <- order(center_dists[1, ])[1:min(5000, N)]
    dists <- as.matrix(dist(coordinates[centered,]))
  }else
  {
    dists <- as.matrix(dist(coordinates))
  }
  
  if(neighbors > 1)
  {
    min_dists <- matrix(NA, nrow = nrow(dists), ncol = neighbors)
    for(r in 1:nrow(dists))
      min_dists[r, ] <- sort(dists[r, which(dists[r,]>0)])[1:neighbors]
  }else
  {
    diag(dists) <- NA
    min_dists <- apply(dists, 1, min, na.rm = TRUE)
  }
  
  radius <- as.numeric(quantile(min_dists, 0.99)) # 99% percentile
  
  if(return_dists)
  {
    return(list(radius = radius, dists = dists))
  }else
  {
    return(radius)  
  }
}

find_spatial_distance_matrix <- function(features, coordinates, radius = "auto")
{
  # features is assumed to be normalized so that l2-norm of each row is 1
  
  N <- nrow(features)
  ST_DS <- NULL
  if(radius == "auto")
  {
    radius_info <- find_spatial_neighborhood_radious(coordinates, neighbors = 8, return_dists = TRUE)
    if(!is.null(radius_info$dists))
      ST_DS <- radius_info$dists
    radius <- radius_info$radius
  }
  
  if(is.null(ST_DS))
    ST_DS <- as.matrix(dist(coordinates))
  
  ST_DS <- 0.5*(ST_DS > radius)
  ST_DX <- 1 - features %*% t(features)
  ST_DX[which(ST_DX <0)] <- 0
  ST_DX <- sqrt(ST_DX)
  
  ST_D <- ST_DS + 0.5*(ST_DX)
  
  return(list(distance_matrix = ST_D, spatial_radius = radius))
}

load_sample_data <- function(data_directory = "Data")
{
  sc_file <- sprintf("%s/sample_scRNAseq.rds", data_directory)
  if(!file.exists(sc_file))
  {
    stop("scRNAseq file not found")
  }
  
  st_file <- sprintf("%s/sample_multicell_spatial_data.rds", data_directory)
  if(!file.exists(st_file))
  {
    stop("Spatial data file not found")
  }
  
  return(list(sc_data = readRDS(sc_file), st_data = readRDS(st_file)))
}