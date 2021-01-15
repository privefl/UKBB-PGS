library(bigsnpr)
fam <- snp_attach("data/UKBB_HM3.rds")$fam
G <- snp_attach("data/UKBB_geno.rds")$genotypes
ind_test <- which(fam$set == "test")
fam_test <- fam[ind_test, ]
ind_csv <- match(fam_test$eid, readRDS("data/csv_eid.rds"))
covar <- dplyr::select(fam_test, -eid, -group, -set)

POP <- c("United Kingdom", "Poland", "Italy", "Iran",
         "India", "China", "Caribbean", "Nigeria")
pop <- factor(fam_test$group, levels = POP)


NCORES <- 15
library(future.batchtools)
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES, mem = "125g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

bigassertr::assert_dir("assess-pred-geno")

all_mod <- list.files("PLR-UKBB-geno", full.names = TRUE)
(all_mod_todo <- all_mod[
  !file.exists(file.path("assess-pred-geno", basename(all_mod)))])

furrr::future_walk(all_mod_todo, function(res_file) {

  res_file2 <- file.path("assess-pred-geno", basename(res_file))
  if (!file.exists(res_file2)) {

    dir <- `if`(dirname(res_file) == "PLR-UKBB-geno",
                "data/ukbb-quant-pheno",
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


length(files <- list.files("assess-pred-geno", full.names = TRUE))
library(dplyr)
res <- purrr::map_dfr(files, function(file) {
  readRDS(file) %>%
    mutate(pheno = sub("\\.rds$", "", basename(file))) %>%
    arrange(validation_loss) %>%
    select(-validation_loss) %>%
    slice(1)
})
res

POP <- c("United Kingdom", "Poland", "Italy", "Iran",
         "India", "China", "Caribbean", "Nigeria")

OUTLIERS <- c("less_tanned", "darker_hair", "darker_skin", "log_bilirubin", "log_lipoA")
stopifnot(all(OUTLIERS %in% res$pheno))


library(ggplot2)

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
# Poland: 92.5%
# Italy: 86.5%
# Iran: 71.3%
# India: 56.5%
# China: 49.6%
# Caribbean: 25.5%
# Nigeria: 18.5%
# ggsave("figures/lasso-ancestry-geno.pdf", width = 12, height = 7)
