library(bigsnpr)
ukbb_large <- snp_attach("data/UKBB_large.rds")
G <- ukbb_large$genotypes

fam <- snp_attach("data/UKBB_HM3.rds")$fam
is_train <- (fam$set == "train")
covar <- as.matrix(dplyr::select(fam, -eid, -group, -set))
ind_csv <- match(fam$eid, readRDS("data/csv_eid.rds"))

pheno_files <- c(
  paste0("data/ukbb-quant-pheno/",
         c("log_BMI", "height", "log_bilirubin", "log_lipoA", "apoB"), ".rds"),
  paste0("data/ukbb-binary-pheno/",
         c("174.1", "185", "250.2", "401", "411.4"), ".rds"))
file.exists(pheno_files)

intervals <- bigparallelr::split_len(ncol(G), block_len = 100e3)

NCORES <- 15
library(future.batchtools)
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 1, mem = "125g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

bigassertr::assert_dir("gwas-large")

furrr::future_walk(rows_along(intervals), function(ic) {

  ind <- seq(intervals[ic, "lower"], intervals[ic, "upper"])

  lapply(pheno_files, function(pheno_file) {

    pheno <- sub("\\.rds$", "", basename(pheno_file))
    res_file <- paste0("gwas-large/", pheno, "_part", ic, ".rds")

    if (!file.exists(res_file)) {

      y <- readRDS(pheno_file)[ind_csv] + 0
      ind.train <- which(!is.na(y) & is_train)

      if (length(table(covar[ind.train, 1])) == 1) {
        COVAR <- covar[, -1]  # rm sex as covariate
      } else {
        COVAR <- covar
      }

      gwas <- big_univLinReg(G, y[ind.train], ind.train = ind.train,
                             covar.train = COVAR[ind.train, ],
                             ind.col = ind, ncores = NCORES)

      saveRDS(gwas, res_file)
    }
  })
})
