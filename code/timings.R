all_mod <- c(list.files("PLR-UKBB", full.names = TRUE),
             list.files("PLR-UKBB-binary", full.names = TRUE))

library(doParallel)
registerDoParallel(cl <- makeCluster(15))

res <- foreach(res_file = all_mod) %dopar% {

  library(bigstatsr)
  mod <- readRDS(res_file)
  summ <- summary(mod)[2:4]

  library(foreach)
  summ$nb_var <- foreach(chosen_mod = mod) %:%
    foreach(mod_sub = chosen_mod, .combine = 'c') %do% {
      sum(mod_sub$beta != 0)
    }
  summ$time <- foreach(chosen_mod = mod) %:%
    foreach(mod_sub = chosen_mod, .combine = 'c') %do% {
      mod_sub$time
    }

  summ
}
stopCluster(cl)

res2 <- do.call(dplyr::bind_rows, res)
res3 <- tidyr::unnest(res2, cols = c(nb_var, time))

library(ggplot2)
ggplot(res3, aes(nb_var, time / 3600)) +
  bigstatsr::theme_bigstatsr() +
  # geom_point(aes(x, y), alpha = 0, data = data.frame(x = 0, y = 0)) +
  geom_point(alpha = 0.2) +
  scale_y_sqrt(breaks = c(0:6 * 10, 1, 5), minor_breaks = NULL) +
  geom_smooth(method = "lm", formula = "y ~ x + 0") +
  labs(x = "Number of non-zero variables", y = "Time (hours, sqrt-scale)")
# ggsave("figures/timings.png", width = 9, height = 6)


## LDpred2
length(files <- list.files("assess-pred-ldpred2", "\\.rds$", full.names = TRUE))
library(furrr)
(NCORES <- availableCores() - 1L)
plan("multisession", workers = NCORES)
all_res0 <- future_map_dfr(files, readRDS)

library(ggplot2)
theme_set(bigstatsr::theme_bigstatsr())
qplot(node, time / 3600, geom = "boxplot", data = all_res0) +
  coord_flip() +
  geom_hline(yintercept = 12, color = "red") +
  scale_y_continuous(breaks = 0:12, limits = c(0, 12)) +
  labs(x = "Compute node", y = "Time (in hours)") +
  theme(plot.margin = margin(r = 1, 0.5, 0.5, 0.5, unit = "lines"))
# ggsave("figures/timings-ldpred2.pdf", width = 7, height = 9)
