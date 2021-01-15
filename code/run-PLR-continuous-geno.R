library(bigsnpr)
G <- snp_attach("data/UKBB_geno.rds")$genotypes
fam <- snp_attach("data/UKBB_HM3.rds")$fam
is_train <- (fam$set == "train")
covar <- as.matrix(dplyr::select(fam, -eid, -group, -set))
ind_csv <- match(fam$eid, readRDS("data/csv_eid.rds"))

SLOW_ONES <- paste0(c(
  "log_HDL", "log_IGF1", "log_height", "log_BMI", "log_cystatinC", "log_erythrocyte",
  "log_platelet", "log_leukocyte", "log_haemoglobin", "log_haematocrit",
  "log_log_gammaGT",  "MCH", "log_platelet_crit", "apoA", "log_urate",
  "protein", "MCV", "log_ALP", "log_creatinine",
  "log_potassium_urine", "log_HLR_reticulocyte", "log_HLR_reticulocyte_perc",
  "sphered_cell_volume", "log_triglycerides", "log_HbA1c", "log_CRP",
  "log_platelet_volume", "log_platelet_width", "log_lymphocyte", "log_monocyte",
  "log_reticulocyte", "log_reticulocyte_perc", "log_basophil_perc",
  "log_eosinophil_perc", "neutrophil_perc", "logit_monocyte_perc",
  "log_neutrophil", "log_impedance", "log_water_mass", "log_fat_free_mass",
  "log_fat_mass", "logit_fat_perc", "less_happy_with_health", "log_FVC",
  "haematocrit_perc", "reticulocyte_volume", "height", "LTFH_296.2",
  "erythrocyte", "fat_perc", "log_waist_circ", "sitting_height"
), ".rds")

library(dplyr)
files <- tibble(
  basename = list.files("data/ukbb-quant-pheno"),
  pheno_file = file.path("data/ukbb-quant-pheno", basename),
  res_file = file.path("PLR-UKBB-geno", basename)
) %>%
  filter(!file.exists(res_file)) %>%
  filter(!basename %in% SLOW_ONES) %>%
  print(n = Inf)

bigassertr::assert_dir("PLR-UKBB-geno")
bigassertr::assert_dir("log")

# dput(sub("\\.rds$", "", files$basename))


library(future.batchtools)
## For fast ones
NCORES <- 30
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 2, mem = "250g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))
## For slow ones
NCORES <- 60
plan(batchtools_slurm(resources = list(
  t = "5-00:00", c = NCORES + 4, mem = "500g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))
## For queued ones
NCORES <- 23
plan(batchtools_slurm(resources = list(
  t = "7-00:00", c = NCORES + 1, mem = "250g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

furrr::future_pmap(files[-1], function(pheno_file, res_file) {

  y <- readRDS(print(pheno_file))[ind_csv]
  ind.train <- which(!is.na(y) & is_train)

  if (length(table(covar[ind.train, 1])) == 1) {
    COVAR <- covar[, -1]  # rm sex as covariate
  } else {
    COVAR <- covar
  }

  ## Penalized reg
  mod <- big_spLinReg(G, y[ind.train], ind.train = ind.train, K = 10,
                      covar.train = COVAR[ind.train, ],
                      pf.covar = rep(0, ncol(COVAR)),
                      power_scale = c(0, 0.5, 1),
                      power_adaptive = c(0, 0.5, 1.5),
                      lambda.min.ratio = 1e-6, nlam.min = 30,
                      n.abort = 2, dfmax = 200e3, ncores = NCORES)

  saveRDS(mod, res_file)
})
