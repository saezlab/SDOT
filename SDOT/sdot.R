###################################################################################################################
########################### SDOT: Spatial Cell Type Deconvolution by Optimal Transport ############################
## Preliminary work. Under review by the International Conferenceon Machine Learning (ICML). Do not distribute. ###
###################################################################################################################

setClass("SDOT", slots = list(st_X = "matrix", st_coordinates = "matrix", st_D = "matrix", st_radius = "numeric",
                            sc_X = "matrix", sc_ratios = "numeric", sc_D = "matrix"))


sdot.setup <- function(st_X, st_coordinates, celltype_signatures, celltype_ratios = NULL, st_radius = "auto")
{
  l2 <- function(x) norm(x, type="2")
  
  normalize <- function(x)
  {
    x <- x / matrix(apply(x, 1, l2), byrow = FALSE,
                    nrow = nrow(x), ncol = ncol(x))
    x[which(is.nan(x))] <- 0
    
    return(x)
  }
  
  if(!is.null(celltype_ratios))
    celltype_signatures <- celltype_signatures[which(celltype_ratios > 0), ]
  
  message("Computing the distance matrices")
  SC_X <- normalize(celltype_signatures)
  SC_D <- 1 - SC_X %*% t(SC_X)
  SC_D[which(SC_D <0)] <- 0
  SC_D <- sqrt(SC_D)
  
  ST_D_info <- find_spatial_distance_matrix(features = normalize(st_X), coordinates = st_coordinates, radius = st_radius)
  
  common_genes <- intersect(colnames(st_X), colnames(celltype_signatures))
  
  new("SDOT", st_X = normalize(st_X[, common_genes]), st_coordinates = st_coordinates, 
      st_D = ST_D_info$distance_matrix, st_radius = ST_D_info$spatial_radius,
      sc_X = normalize(celltype_signatures[, common_genes]), sc_ratios = celltype_ratios, sc_D = SC_D)
}

sdot.solve <- function(object,
                       iterations = -1, gap_threshold = 0.01, initial_solution = c("ratios"),
                       lambda_GW = 0.1, lambda_A = 1, lambda_R = 1, 
                       verbose = TRUE, result_saver = NULL, save_iterations = -1)
{
  G <- ncol(object@st_X)
  if(G != ncol(object@sc_X))
    stop("Invalid arguments")
  
  C <- nrow(object@sc_X)
  S <- nrow(object@st_X)
  
  SC_ratios <- object@sc_ratios
  if(is.null(SC_ratios))
  {
    SC_ratios <- rep(1/C, C)
    names(SC_ratios) <- rownames(object@sc_X)
    
    lambda_A <- 0 #cannot penalize abundance if ratios are not available
  }else
  {
    SC_ratios <- SC_ratios/sum(SC_ratios)
  }
  
  r_SC <- SC_ratios * S # expected abundance of cell types
  r_ST <- rep(1, S)
  
  if(iterations <= 0)
    iterations <- 1000
  
  iteration_start <- Sys.time()
  
  if(is(initial_solution, "matrix"))
  {
    Yt <- initial_solution
  }else if(initial_solution == "ratios")
  {
    Yt <- outer(SC_ratios, r_ST)
  }
  
  l2 <- function(x) norm(x, type="2")
  
  linear_gw <- NULL
  if(lambda_GW > 0)
  {
    SC_D2 <- object@sc_D^2  
    ST_D2 <- (object@st_D^2) %*% r_ST # rowSums(object@st_D^2)
    
    linear_gw <- matrix(ST_D2, nrow = C, ncol = S, byrow = TRUE)
  }
  
  f <- Inf #best value (upper bound)
  LB <- -Inf #best lower bound
  Y <- NULL #best solution
  
  iteration_cols <- c("ft", "UB", "LB", "Gap", "Time", "d_ST", "d_SC", "d_GW", "d_A", "d_R", "step_size")
  iteration_info <- matrix(NA, nrow = iterations, ncol = length(iteration_cols))
  colnames(iteration_info) <- iteration_cols
  
  iteration <- 1
  converged <- FALSE
  while(!converged)
  {
    rho_t <- rowSums(Yt)
    rho_t[which(rho_t < 1e-5)] <- 1e-5
    
    delta_ratio <- -lambda_R/rho_t
    
    d_A <- 0
    if(lambda_A > 0)
    {
      d_A <- l2(rho_t-r_SC)  
      
      if(d_A > 1e-5)
      {
        delta_ratio <- lambda_A*(rho_t-r_SC)/d_A - lambda_R/rho_t
      }
    }
    
    Delta <- matrix(delta_ratio, nrow = C, ncol = S, byrow = F)
    
    dcosine_ST <- 0
    dcosine_SC <- 0
    st_Xt <- t(Yt) %*% object@sc_X
    ST_De <- matrix(NA, nrow = S, ncol = G)
    for (i in 1:S) {
      p <- sum(object@st_X[i, ]*st_Xt[i, ])
      n <- l2(st_Xt[i, ])
      d <- 1-p/n
      if(d < 1e-5)
        d <- 1e-5
      d <- sqrt(d)
      
      ST_De[i, ] <- -0.5*r_ST[i]/d*(object@st_X[i, ]/n - st_Xt[i, ]/n^3*p)
      
      dcosine_ST <- dcosine_ST + d*r_ST[i]
    }
    
    sc_Xt <- Yt %*% object@st_X
    SC_De <- matrix(NA, nrow = C, ncol = G)
    for (k in 1:C) {
      p <- sum(object@sc_X[k, ]*sc_Xt[k, ])
      n <- l2(sc_Xt[k, ])
      if(n <= 1e-5)
      {
        d <- 1
        n <- 1e-5
      }else
      {
        d <- 1-p/n
        if(d < 1e-5)
          d <- 1e-5
        d <- sqrt(d)
      }
      
      SC_De[k, ] <- -0.5*r_SC[k]/d*(object@sc_X[k, ]/n - sc_Xt[k, ]/n^3*p)
      
      dcosine_SC <- dcosine_SC + d*r_SC[k]
    }
    
    Delta <- Delta + object@sc_X %*% t(ST_De) + SC_De %*% t(object@st_X)
    
    d_GW <- 0
    if(lambda_GW > 0)
    {
      Z <- object@sc_D %*% Yt %*% object@st_D
      
      SC_D2t <- matrix(SC_D2 %*% rho_t, nrow = C, ncol = S, byrow = FALSE)
      d_GW <- sum(Yt * (SC_D2t + linear_gw - 2*Z))
      if(d_GW > 1e-5)
      {
        d_GW <- sqrt(d_GW)
        Delta <- Delta + (0.5*lambda_GW/d_GW) *(2*SC_D2t + linear_gw - 4*Z) 
      }
    }
    
    d_R <- 0
    if(lambda_R > 0)
      d_R <- -sum(log(rho_t))
    
    ft <- dcosine_SC + dcosine_ST + 
      lambda_GW * d_GW + lambda_A * d_A + lambda_R * d_R
    
    if(ft < f)
    {
      f <- ft
      Y <- Yt
    }
    
    # atom:
    Yt_h <- matrix(0, nrow = C, ncol = S)
    for (i in 1:S)
    {
      Yt_h[which.min(Delta[, i]), i] <- r_ST[i]
    }
    
    Y_diff <- Yt-Yt_h
    
    step <- 2/(iteration + 1)
    
    #gap:
    gap <- sum(Delta * Y_diff)
    if(LB < ft - gap)
      LB <- ft - gap
    
    gap <- f - LB
    
    if(abs(f) > 1e-10)
      gap <- abs(gap / f)
    
    time <- as.numeric(Sys.time() - iteration_start, units = "secs")
    if(verbose)
    {
      if(iteration %% 10 == 1)
        message("Iteration: obj-val, UB, LB, d_ST, d_SC, d_GW, d_A, d_R, gap, step, time")
      message(sprintf("%d: %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g",
                      iteration, ft, f, LB, 
                      dcosine_ST, dcosine_SC, d_GW, d_A, d_R, 
                      gap, step, time))
    }
    
    
    iteration_info[iteration, ] <- c(ft, f, LB, gap, time, dcosine_ST, dcosine_SC, d_GW, d_A, d_R, step)
    
    iteration_start <- Sys.time()
    
    if(save_iterations > 0 & iteration %% save_iterations == 0 & !is.null(result_saver))
    {
      message("... Saving mapping ...")
      result_saver(Y)
    }
    
    if(gap <= gap_threshold)
      converged <- TRUE
    
    if(step <= 1e-5)
      converged <- TRUE
    
    if(iterations > 0 & iteration >= iterations)
      converged <- TRUE
    
    if(converged)
      break
    
    Yt <- Yt - step * Y_diff
    iteration <- iteration + 1
  }
  
  return(list(Y = Y, f = f, history = as.data.frame(iteration_info[1:iteration, ])))
}

