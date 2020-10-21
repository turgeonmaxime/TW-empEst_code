library(tidyverse)

set.seed(12345)

source("simulation_functions.R")

# Number of simulations
B <- 1000

# Simulation parameters are encoded in the file names
list_params <- list.files("cache", full.names = TRUE) %>%
  tibble::tibble("filename" = .) %>% 
  filter(str_detect(filename, "largRoot"),
         !str_detect(filename, "shrink")) %>% # Exclude shrinkage estimates
  separate(filename, c("dir", "p", "rho"), sep = "[_]", remove = FALSE) %>% 
  mutate(rho = str_remove(rho, ".rds"),
         p = as.integer(p),
         rho = as.numeric(rho))

# Analyse the generated data
results <- purrr::map_df(seq_len(nrow(list_params)), 
                         function(index) {
  params <- list_params[index,]
  
  dist <- readRDS(params$filename)
  ls_param25 <- fit_heuristic(sample(dist, size = 25, replace = FALSE))
  ls_param50 <- fit_heuristic(sample(dist, size = 50, replace = FALSE))
  ls_param75 <- fit_heuristic(sample(dist, size = 75, replace = FALSE))
  ls_param100 <- fit_heuristic(sample(dist, size = 100, replace = FALSE))
  
  result_orig <- data.frame(X = dist) %>%
    mutate(Y = cume_dist(X), Type = "True CDF") %>%
    arrange(X)

  result_heuristic25 <- data.frame(X = seq(min(dist), max(dist), length.out = B),
                                   Y = ptw_ls(log(seq(min(dist), max(dist), length.out = B)),
                                                   ls_param25[1], ls_param25[2])) %>%
    mutate(Type = "EE.25")
  result_heuristic50 <- data.frame(X = seq(min(dist), max(dist), length.out = B),
                                   Y = ptw_ls(log(seq(min(dist), max(dist), length.out = B)),
                                                   ls_param50[1], ls_param50[2])) %>%
    mutate(Type = "EE.50")
  result_heuristic75 <- data.frame(X = seq(min(dist), max(dist), length.out = B),
                                   Y = ptw_ls(log(seq(min(dist), max(dist), length.out = B)),
                                                   ls_param75[1], ls_param75[2])) %>%
    mutate(Type = "EE.75")
  result_heuristic100 <- data.frame(X = seq(min(dist), max(dist), length.out = B),
                                    Y = ptw_ls(log(seq(min(dist), max(dist), length.out = B)),
                                                    ls_param100[1], ls_param100[2])) %>%
    mutate(Type = "EE.100")
  
  bind_rows(result_orig,
            result_heuristic25,
            result_heuristic50,
            result_heuristic75,
            result_heuristic100) %>%
    mutate("p" = params$p, 
           "rho" = params$rho)
  
})

plot_results <- results %>%
  mutate(Type = factor(Type, levels = c("True CDF", 
                                        "EE.25",
                                        "EE.50",
                                        "EE.75",
                                        "EE.100"))) %>% 
  ggplot(aes(X, Y, colour = Type)) + geom_line() +
  facet_grid(rho ~ p, labeller = pryr::partial(label_both, sep = " = "), 
             scales = "free") + 
  xlab("Largest root") + ylab("CDF") +
  scale_colour_manual(values = c("black", gg_color_hue(4))) + 
  theme(legend.position = "top", axis.text.x = element_text(angle = 90, 
                                                            hjust = 1))

ggsave(plot_results, file = "figures/approx_results.pdf", 
       width = 7,  height = 3.5)
