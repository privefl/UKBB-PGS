library(bigsnpr)
ukbb <- snp_attach("data/UKBB_HM3.rds")
G <- ukbb$genotypes
ind_test <- which(ukbb$fam$set == "test")
fam_test <- ukbb$fam[ind_test, ]
ind_csv <- match(fam_test$eid, readRDS("data/csv_eid.rds"))
covar <- dplyr::select(fam_test, -eid, -group, -set)
# center <- bigreadr::fread2("UKBB/ukb41181.csv", select = "54-0.0")[ind_csv, 1]
# table(center, exclude = NULL)
# covar$center <- as.factor(center)
POP <- c("United Kingdom", "Poland", "Italy", "Iran",
         "India", "China", "Caribbean", "Nigeria")
pop <- factor(fam_test$group, levels = POP)


NCORES <- 15
library(future.batchtools)
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES, mem = "125g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

bigassertr::assert_dir("assess-pred")

all_mod <- c(list.files("PLR-UKBB",        full.names = TRUE),
             list.files("PLR-UKBB-binary", full.names = TRUE))
(all_mod_todo <- all_mod[!file.exists(file.path("assess-pred", basename(all_mod)))])

furrr::future_walk(all_mod_todo, function(res_file) {

  res_file2 <- file.path("assess-pred", basename(res_file))
  if (!file.exists(res_file2)) {

    dir <- `if`(dirname(res_file) == "PLR-UKBB", "data/ukbb-quant-pheno",
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


length(files <- grep("/(?!LTFH_).*$", perl = TRUE, value = TRUE,
                     list.files("assess-pred", full.names = TRUE)))
library(dplyr)
res <- purrr::map_dfr(files, function(file) {
  readRDS(file) %>%
    mutate(pheno = sub("\\.rds$", "", basename(file))) %>%
    arrange(validation_loss) %>%
    select(-validation_loss) %>%
    slice(1)
})
res

with(res, table(power_adaptive, power_scale))
#               power_scale
# power_adaptive  0 0.5  1
#            0   18  97 17
#            0.5  2  10 25
#            1.5 11  19 46

library(ggplot2)

POP <- c("United Kingdom", "Poland", "Italy", "Iran",
         "India", "China", "Caribbean", "Nigeria")

OUTLIERS <- c("less_tanned", "darker_hair", "darker_skin", "red_hair",
              "log_bilirubin", "log_lipoA")
stopifnot(all(OUTLIERS %in% res$pheno))

filter(res, pheno %in% OUTLIERS) %>%
  select(-(1:3)) %>%
  filter(pheno != "log_bilirubin") %>%
  tidyr::unnest()

slopes <- list("United Kingdom" = 1)

bigstatsr::plot_grid(plotlist = lapply(POP[-1], function(pop) {

  df <- tibble(pheno = res$pheno, x = res$`United Kingdom`, y = res[[pop]]) %>%
    filter(purrr::map_dbl(y, ~ .[3] - .[2]) < 0.2,
           !startsWith(pheno, "LTFH")) %>%
    tidyr::unnest_wider("x", names_sep = "_") %>%
    tidyr::unnest_wider("y", names_sep = "_") %>%
    na.omit() %>%
    mutate(., eps = abs(resid(lm(y_1 ~ x_1 + 0, data = .))))

  robust_slope <- deming::deming(y_1 ~ x_1 + 0,
                                 data = filter(df, !pheno %in% OUTLIERS, y_2 != y_3),
                                 xstd = x_3 - x_2,
                                 ystd = y_3 - y_2)$coef[2]

  cat(pop, ": ", round(100 * (slopes[[pop]] <<- unname(robust_slope^2)), 1), "%\n", sep = "")

  ggplot(df, aes(x_1, y_1, label = pheno)) +
    bigstatsr::theme_bigstatsr(0.85) +
    geom_abline(color = "red", linetype = 2) +
    geom_hline(yintercept = 0, color = "red", linetype = 3) +
    geom_vline(xintercept = 0, color = "red", linetype = 3) +
    geom_abline(slope = robust_slope, color = "blue") +
    geom_errorbar(aes(xmin = x_2, xmax = x_3), width = 0, color = "green") +
    geom_errorbar(aes(ymin = y_2, ymax = y_3), width = 0, color = "green") +
    geom_point() +
    ggrepel::geom_text_repel(data = slice_max(df, eps, n = 5),
                             min.segment.length = 0, seed = 42, box.padding = 0.5) +
    labs(x = "United Kingdom", y = pop) +
    xlim(-0.1, 0.65) + ylim(-0.2, 0.65) +
    ggtitle(paste0(round(100 * unname(robust_slope^2), 1), "%")) +
    coord_equal()
}), nrow = 2, scale = 0.95)
# Poland: 93.8%
# Italy: 85.6%
# Iran: 72.2%
# India: 64.7%
# China: 48.6%
# Caribbean: 25.2%
# Nigeria: 18%
# ggsave("figures/lasso-ancestry.pdf", width = 12, height = 7)

library(bigsnpr)
ukbb <- snp_attach("data/UKBB_HM3.rds")
PC <- dplyr::select(ukbb$fam, PC1:PC16)
centers <- bigutilsr::geometric_median(PC, by_grp = ukbb$fam$group)
dist_to_UK <- as.matrix(dist(centers))[, "United Kingdom"]

qplot(dist_to_UK[names(slopes)], unlist(slopes), size = I(2)) +
  theme_bigstatsr() +
  geom_smooth(method = "lm") +
  ggrepel::geom_label_repel(aes(label = names(slopes)), min.segment.length = 0) +
  # scale_y_log10() +
  scale_y_continuous(breaks = 0:5 / 5, minor_breaks = 0:10 / 10) +
  labs(x = "PC distance to UK", y = "Relative predictive performance with UK")
# ggsave("figures/ratio-dist.pdf", width = 8, height = 6)
