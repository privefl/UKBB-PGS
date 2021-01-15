library(bigsnpr)
ukbb <- snp_attach("data/UKBB_HM3.rds")
G <- ukbb$genotypes
is_train <- (ukbb$fam$set == "train")
covar <- as.matrix(dplyr::select(ukbb$fam, -eid, -group, -set))
ind_csv <- match(ukbb$fam$eid, readRDS("data/csv_eid.rds"))

pheno_files <- c(list.files("data/ukbb-quant-pheno", full.names = TRUE),
                 list.files("data/ukbb-binary-pheno", full.names = TRUE))

library(dplyr)
files <- tibble(
  pheno_file = pheno_files,
  res_file = file.path("GWAS", basename(pheno_file))
) %>%
  filter(!file.exists(res_file)) %>%
  print(n = Inf)

bigassertr::assert_dir("GWAS")

NCORES <- 15
library(future.batchtools)
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 1, mem = "125g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

furrr::future_pmap(files, function(pheno_file, res_file) {

  y <- readRDS(print(pheno_file))[ind_csv] + 0
  ind.train <- which(!is.na(y) & is_train)

  if (length(table(covar[ind.train, 1])) == 1) {
    COVAR <- covar[, -1]  # rm sex as covariate
  } else {
    COVAR <- covar
  }

  gwas <- big_univLinReg(G, y[ind.train], ind.train = ind.train,
                         covar.train = COVAR[ind.train, ],
                         ncores = NCORES)

  saveRDS(gwas, res_file)
})


library(ggplot2)
library(bigsnpr)
ukbb <- snp_attach("data/UKBB_HM3.rds")
CHR <- as.integer(ukbb$map$chromosome)
POS <- ukbb$map$physical.pos

gwas <- readRDS("GWAS/log_bilirubin.rds")
gwas <- readRDS("GWAS/less_tanned.rds")
gwas <- readRDS("GWAS/log_lipoA.rds")
gwas <- readRDS("GWAS/red_hair.rds")

snp_manhattan(gwas, CHR, POS, npoints = 100e3) +
  scale_y_log10() +
  geom_hline(yintercept = -log10(5e-8), color = "red") +
  ggtitle(NULL, NULL)
# ggsave("figures/bilirubin_manhattan.png", width = 13, height = 6)
# ggsave("figures/less_tanned_manhattan.png", width = 13, height = 5)
