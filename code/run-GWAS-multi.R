library(bigsnpr)
ukbb <- snp_attach("data/UKBB_HM3.rds")
G <- ukbb$genotypes
# replace some UK people in the training by other ancestries
table(set <- ukbb$fam$set)
pop <- ukbb$fam$group
ind <- which(pop %in% c("Poland", "Italy", "Iran", "India", "China", "Nigeria"))
set.seed(1); set[sample(which(set == "train"), length(ind))] <- "skip"
set[ind] <- "train"
table(set)

pop[set != "train" | pop == "Iran"] <- NA
pop2 <- factor(pop, levels = c("United Kingdom", "Poland", "Italy",
                               "India", "China", "Nigeria"),
               labels = c("eur", "eur", "eur", "sas", "eas", "afr"))
table(pop2)
#    eur    sas    eas    afr
# 377859   6331   1810   3924

covar <- as.matrix(dplyr::select(ukbb$fam, -eid, -group, -set))
ind_csv <- match(ukbb$fam$eid, readRDS("data/csv_eid.rds"))


pheno_files <- c(
  file.path("data/ukbb-quant-pheno",
            c("years_of_edu", "log_lipoA", "darker_skin", "log_bilirubin")),
  file.path("data/ukbb-binary-pheno", c("174.1", "185", "250.2", "401", "411.4")))

library(dplyr)
files <- tidyr::expand_grid(
  pheno_file = paste0(pheno_files, ".rds"),
  pop = c("eur", "sas", "eas", "afr")
) %>%
  mutate(res_file = paste0("GWAS-multi/", pop, "_", basename(pheno_file))) %>%
  filter(!file.exists(res_file)) %>%
  print(n = Inf)

bigassertr::assert_dir("GWAS-multi")

NCORES <- 15
library(future.batchtools)
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 1, mem = "100g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

furrr::future_pmap(files, function(pheno_file, pop, res_file) {

  y <- readRDS(print(pheno_file))[ind_csv] + 0
  ind.train <- which(!is.na(y) & pop2 == pop)

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
