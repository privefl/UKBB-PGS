library(bigsnpr)
ukbb <- snp_attach("data/UKBB_HM3_AJ.rds")
G <- ukbb$genotypes
fam_test <- ukbb$fam
ind_csv <- match(fam_test$eid, readRDS("data/csv_eid.rds"))
covar <- dplyr::select(fam_test, -eid, -pop, -country, -continent)
dim(covar)  # 1709 x 20

library(future.batchtools)
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = 3, mem = "25g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

bigassertr::assert_dir("assess-pred-AJ")

all_mod <- c(list.files("PLR-UKBB",        full.names = TRUE),
             list.files("PLR-UKBB-binary", full.names = TRUE))
(all_mod_todo <- all_mod[!file.exists(file.path("assess-pred-AJ", basename(all_mod)))])

furrr::future_walk(all_mod_todo, function(res_file) {

  res_file2 <- file.path("assess-pred-AJ", basename(res_file))
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
        predict(mod[k], G, covar.row = matrix(0, nrow(G), K))
      } else {
        predict(mod[k], G, covar.row = matrix(0, nrow(G), K), proba = FALSE)
      }

      bigstatsr::pcor(pred, y, covar)
    })

    library(dplyr)
    summ %>%
      select(-alpha, -intercept, -beta, -message, -all_conv) %>%
      saveRDS(res_file2)
  }

  NULL
})
