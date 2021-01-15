library(bigsnpr)
ukbb_large <- snp_attach("data/UKBB_large.rds")
G <- ukbb_large$genotypes
CHR <- as.integer(ukbb_large$map$chromosome)
POS <- ukbb_large$map$physical.pos

fam <- snp_attach("data/UKBB_HM3.rds")$fam
is_train <- (fam$set == "train")
ind_csv <- match(fam$eid, readRDS("data/csv_eid.rds"))

pheno_files <- c(
  paste0("data/ukbb-quant-pheno/",
         c("log_BMI", "height", "log_bilirubin", "log_lipoA", "apoB"), ".rds"),
  paste0("data/ukbb-binary-pheno/",
         c("174.1", "185", "250.2", "401", "411.4"), ".rds"))
file.exists(pheno_files)

bigassertr::assert_dir("LD")

library(dplyr)
files <- expand.grid(
  pheno_file = pheno_files,
  chr = 1:22,
  stringsAsFactors = FALSE
) %>%
  mutate(pheno = sub("\\.rds$", "", basename(pheno_file)),
         res_file = paste0("LD/", pheno, "_chr", chr, ".rds")) %>%
  filter(!file.exists(res_file)) %>%
  print()

NCORES <- 15
library(future.batchtools)
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 1, mem = "125g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

furrr::future_pwalk(files, function(pheno_file, chr, pheno, res_file) {

  # pheno_file <- "data/ukbb-quant-pheno/log_BMI.rds"
  # chr <- 22
  # pheno <- "log_BMI"
  # res_file <- "LD/log_BMI_chr22.rds"

  gwas <- purrr::map_dfr(
    paste0("gwas-large/", pheno, "_part", 1:83, ".rds"), readRDS)
  lpval <- -bigstatsr:::predict.mhtest(gwas)
  thr <- sort(lpval, decreasing = TRUE)[1e6]
  ind_top <- which(lpval >= thr)

  y <- readRDS(pheno_file)[ind_csv] + 0
  ind.train <- which(!is.na(y) & is_train)

  if (chr == 22) {
    gwas2 <- gwas[ind_top, ]
    gwas3 <- dplyr::transmute(gwas2, beta_se = std.err, n_eff = length(ind.train))
    gwas3$beta <- bigsnpr::snp_thr_correct(gwas2$estim, gwas2$std.err, thr_lpS = thr)
    gwas3$id <- ind_top
    saveRDS(gwas3, paste0("gwas-large/", pheno, "_top.rds"))
  }

  ## indices in 'ind_top'
  ind.chr <- which(CHR[ind_top] == chr)
  ## indices in 'G'
  ind.chr2 <- ind_top[ind.chr]

  POS2 <- snp_asGeneticPos(CHR[ind.chr2], POS[ind.chr2], dir = "tmp-data")
  corr_chr <- snp_cor(G, ind.row = ind.train, ind.col = ind.chr2,
                      size = 3 / 1000, infos.pos = POS2, ncores = NCORES)

  saveRDS(corr_chr, res_file)
})

#### LDpred2-auto ####

bigassertr::assert_dir("ldpred2-large")

files <- tibble(
  pheno_file = pheno_files,
  res_file = file.path("ldpred2-large", basename(pheno_file))
) %>%
  filter(!file.exists(res_file)) %>%
  print()

NCORES <- 30
library(future.batchtools)
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 2, mem = "250g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

furrr::future_pmap(files, function(pheno_file, res_file) {

  # pheno_file <- "data/ukbb-binary-pheno/174.1.rds"

  pheno <- sub("\\.rds$", "", basename(pheno_file))
  df_beta <- readRDS(paste0("gwas-large/", pheno, "_top.rds"))

  print(Sys.getenv("SLURMD_NODENAME"))

  ## Assemble correlations
  tmp <- tempfile(tmpdir = "tmp-data")
  for (chr in 1:22) {
    cat(chr, ".. ", sep = "")
    corr_chr <- readRDS(paste0("LD/", pheno, "_chr", chr, ".rds"))
    if (chr == 1) {
      corr <- bigsparser::as_SFBM(corr_chr, tmp)
    } else {
      corr$add_columns(corr_chr, nrow(corr))
    }
  }

  # LDpred2-auto
  multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = 0.2,
                                 vec_p_init = seq_log(1e-4, 0.9, 30),
                                 sparse = TRUE, ncores = NCORES)

  # cleanup
  rm(corr); gc(); file.remove(paste0(tmp, ".sbk"))

  saveRDS(multi_auto, res_file)
})
