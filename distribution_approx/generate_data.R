library(tidyverse)
library(furrr)

source("simulation_functions.R")
seed <- 12345

# Number of simulations
B <- 1000
mc.cores <- 2

parameter_list <- expand_grid(p = c(500, 1000, 1500, 2000), 
                              rho = c(0, 0.2))

plan(multisession, workers = mc.cores)

# We don't need to save the output of future_map
tmp <- future_map(seq_len(nrow(parameter_list)), function(ind) {
  p <- parameter_list$p[ind]
  rho <- parameter_list$rho[ind]
  
  dist <- generate_largestRootDist(p, B, rho)
  
  saveRDS(dist, paste0("cache/largRootDist_",
                       p, "_", rho, ".rds"))
  
  return(NULL)
})
