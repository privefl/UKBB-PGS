library(bigsnpr)
ukbb <- snp_attach("data/UKBB_HM3.rds")
G <- ukbb$genotypes
# replace some UK people in the training by other ancestries
table(set <- ukbb$fam$set)
ind <- which(ukbb$fam$group %in%
               c("Poland", "Italy", "Iran", "India", "China", "Nigeria"))
set.seed(1); set[sample(which(set == "train"), length(ind))] <- "skip"
set[ind] <- "train"
table(set)
is_train <- (set == "train")
covar <- as.matrix(dplyr::select(ukbb$fam, -eid, -group, -set))
ind_csv <- match(ukbb$fam$eid, readRDS("data/csv_eid.rds"))

focus <- c("years_of_edu", "log_lipoA", "log_BMI", "darker_skin", "log_bilirubin")

library(dplyr)
files <- tibble(
  basename = list.files("data/ukbb-quant-pheno"),
  pheno_file = file.path("data/ukbb-quant-pheno", basename),
  res_file = file.path("PLR-UKBB-multi", basename)
) %>%
  filter(!file.exists(res_file)) %>%
  filter(basename %in% paste0(focus, ".rds")) %>%
  print(n = Inf)

bigassertr::assert_dir("PLR-UKBB-multi")
bigassertr::assert_dir("log")


library(future.batchtools)
NCORES <- 30
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 2, mem = "250g",
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
