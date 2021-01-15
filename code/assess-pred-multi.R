library(bigsnpr)
ukbb <- snp_attach("data/UKBB_HM3.rds")
G <- ukbb$genotypes
POP <- c("United Kingdom", "Caribbean")
ind_test <- which(ukbb$fam$set == "test" & ukbb$fam$group %in% POP)
fam_test <- ukbb$fam[ind_test, ]
ind_csv <- match(fam_test$eid, readRDS("data/csv_eid.rds"))
covar <- dplyr::select(fam_test, -eid, -group, -set)
pop <- factor(fam_test$group, levels = POP)


NCORES <- 15
library(future.batchtools)
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES, mem = "125g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

bigassertr::assert_dir("assess-pred-multi")

all_mod <- c(list.files("PLR-UKBB-multi",        full.names = TRUE),
             list.files("PLR-UKBB-binary-multi", full.names = TRUE))
(all_mod_todo <- all_mod[
  !file.exists(file.path("assess-pred-multi", basename(all_mod)))])

furrr::future_walk(all_mod_todo, function(res_file) {

  res_file2 <- file.path("assess-pred-multi", basename(res_file))
  if (!file.exists(res_file2)) {

    dir <- `if`(dirname(res_file) == "PLR-UKBB-multi", "data/ukbb-quant-pheno",
                "data/ukbb-binary-pheno")
    pheno <- readRDS(file.path(dir, basename(res_file)))
    y <- pheno[ind_csv]
    mod <- readRDS(res_file)
    summ <- summary(mod)
    K <- length(summ$beta[[1]]) - ncol(G)

    summ$pcor <- lapply(seq_along(mod), function(k) {
      print(k)
      pred <- if (attr(mod, "family") == "gaussian") {
        predict(mod[k], G, ind.row = ind_test, ncores = NCORES,
                covar.row = matrix(0, length(ind_test), K))
      } else {
        predict(mod[k], G, ind.row = ind_test, ncores = NCORES, proba = FALSE,
                covar.row = matrix(0, length(ind_test), K))
      }

      lapply(split(seq_along(pop), pop), function(ind) {
        if (sum(y[ind] != 0, na.rm = TRUE) == 0) return(rep(NA_real_, 3))
        bigstatsr::pcor(pred[ind], y[ind], covar[ind, ])
      })
    })

    library(dplyr)
    summ %>%
      select(-alpha, -intercept, -beta, -message, -all_conv) %>%
      tidyr::unnest_wider("pcor") %>%
      saveRDS(res_file2)
  }

  NULL
})


length(files <- list.files("assess-pred-multi", full.names = TRUE))
library(dplyr)
res_multi <- purrr::map_dfr(files, function(file) {
  readRDS(file) %>%
    mutate(pheno = sub("\\.rds$", "", basename(file))) %>%
    arrange(validation_loss) %>%
    select(-validation_loss) %>%
    slice(1)
}) %>%
  mutate(train = "multi")

res <- file.path("assess-pred", basename(files)) %>%
  purrr::map_dfr(function(file) {
    readRDS(file) %>%
      mutate(pheno = sub("\\.rds$", "", basename(file))) %>%
      arrange(validation_loss) %>%
      select(-validation_loss) %>%
      slice(1)
  }) %>%
  mutate(train = "UK") %>%
  select(names(res_multi)) %>%
  bind_rows(res_multi)

res

library(ggplot2)

POP <- c("United Kingdom", "Caribbean")

res %>%
  tidyr::pivot_longer(all_of(POP), names_to = "pop", values_to = "pcor") %>%
  tidyr::unnest_wider("pcor", names_sep = "_") %>%
  mutate(pop = factor(pop, POP), train = factor(train, c("UK", "multi"))) %>%
  ggplot(aes(pop, pcor_1, fill = train)) +
  bigstatsr::theme_bigstatsr() +
  geom_col(position = position_dodge(), alpha = 0.5, color = "black") +
  geom_errorbar(aes(ymin = pcor_2, ymax = pcor_3),
                position = position_dodge(width = 0.9), width = 0.3) +
  facet_wrap(~ pheno, scales = "free_y", ncol = 3) +
  theme(panel.spacing = unit(1, "lines")) +
  labs(x = "Population", y = "Partial correlation", fill = "Training")
# ggsave("figures/lasso_multi_pcor.pdf", width = 12, height = 8)
