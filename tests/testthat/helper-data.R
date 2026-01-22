library(testthat)

set.seed(123)
sim_data_two_sample <- get_data(n = 200, seed = 123)
dat_2s <- sim_data_two_sample$data
covars_2s <- paste0("X", 0:9)

sim_data_single_sample <- get_data(n = 200, seed = 456)
dat_1s <- sim_data_single_sample$data
dat_1s$S <- 1
covars_1s <- paste0("X", 0:9)

# Tiny dataset: 50 rows to be safe for GLM/Ridge
# Use "Y" not "Yobs" - that's what get_data() returns
dat_tiny <- dat_2s[1:50, c("Y", "Tr", "S", paste0("X", 0:4))]
