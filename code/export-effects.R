library(bigsnpr)
ukbb <- snp_attach("data/UKBB_HM3.rds")
G <- ukbb$genotypes
ind_test <- which(ukbb$fam$set == "test")
fam_test <- ukbb$fam[ind_test, ]
ind_csv <- match(fam_test$eid, readRDS("data/csv_eid.rds"))
covar <- dplyr::select(fam_test, -eid, -group, -set)
POP <- c("United Kingdom", "Poland", "Italy", "Iran",
         "India", "China", "Caribbean", "Nigeria")
pop <- factor(fam_test$group, levels = POP)
ind_pop <- split(seq_along(pop), pop)


NCORES <- 15
library(future.batchtools)
plan(batchtools_slurm(workers = 250, finalize = FALSE, resources = list(
  t = "12:00:00", c = NCORES, mem = "125g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

bigassertr::assert_dir("export-effects")

library(dplyr)
all_mod <- c(list.files("PLR-UKBB",        full.names = TRUE),
             list.files("PLR-UKBB-binary", full.names = TRUE)) %>%
  grep("/(?!LTFH_).*$", perl = TRUE, value = TRUE, .) %>%
  { .[!file.exists(file.path("export-effects", basename(.)))] } %>%
  print()

furrr::future_walk(all_mod, function(res_file) {

  res_file2 <- file.path("export-effects", basename(res_file))
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

      sapply(ind_pop, function(ind) {
        if (sum(y[ind] != 0, na.rm = TRUE) == 0) return(NA_real_)
        bigstatsr::pcor(pred[ind], y[ind], covar[ind, ])[1]
      })
    })

    library(dplyr)
    beta_PLR <- summ %>%
      mutate(stat = sapply(summ$pcor, function(.) {
        weighted.mean(., 1 / lengths(ind_pop), na.rm = TRUE)
      })) %>%
      arrange(desc(stat)) %>%
      slice(1) %>%
      mutate(beta = list(head(beta[[1]], ncol(G)))) %>%
      select(stat, beta, pcor) %>%
      tidyr::unnest_wider("pcor")

    mod_ldpred2 <- readRDS(file.path("assess-pred-ldpred2", basename(res_file)))
    if (mod_ldpred2$n_keep > 0) {
      pcor_ldpred2 <- sapply(ind_pop, function(ind) {
        if (sum(y[ind] != 0, na.rm = TRUE) == 0) return(NA_real_)
        bigstatsr::pcor(mod_ldpred2$pred_sp[[1]][ind], y[ind], covar[ind, ])[1]
      })
      stat_ldpred2 <- weighted.mean(pcor_ldpred2, 1 / lengths(ind_pop), na.rm = TRUE)
    } else {
      stat_ldpred2 <- -Inf
    }

    res <- if (stat_ldpred2 > beta_PLR$stat) {
      if (stat_ldpred2 < 0.01) NULL else {
        bind_cols(select(mod_ldpred2, beta = beta_sp), as.list(pcor_ldpred2))
      }
    } else {
      if (beta_PLR$stat < 0.01) NULL else beta_PLR[-1]
    }

    saveRDS(res, res_file2)
  }

  NULL
})


length(files <- list.files("export-effects", full.names = TRUE))
library(dplyr)
res <- purrr::map_dfr(files, function(file) {
  readRDS(file) %>%
    bind_cols(pheno = sub("\\.rds$", "", basename(file)), .)
})
res <- filter(res, !pheno %in% c("darker_hair", "darker_skin"))  # bad transforms, use 0 versions instead
res2 <- res[lengths(res$beta) > 0, ]

all_betas <- do.call("cbind", res2$beta)
colnames(all_betas) <- res2$pheno
map <- ukbb$map %>%
  transmute(rsid, chr = as.integer(chromosome), pos = physical.pos,
            a0 = allele1, a1 = allele2) # reversed in the BGEN file somehow

# Verification (sign)
sumstats <- bigreadr::fread2("../paper-ldpred2/tmp-data/sumstats_BRCA.txt",
                             na.strings = "NULL",
                             select = c("chr", "position_b37", "a0", "a1",
                                        "bcac_onco_icogs_gwas_beta"),
                             col.names = c("chr", "pos", "a0", "a1", "beta"))
info_snp <- snp_match(sumstats, bind_cols(map, beta2 = all_betas[, "174.1"]))
# 11,792,542 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 1,040,087 variants have been matched; 0 were flipped and 1,040,087 were reversed.
ind <- which(info_snp$beta2 != 0)
library(ggplot2)
qplot(beta, beta2, data = info_snp[ind, ]) + theme_bw(15)

sumstats <- bigreadr::fread2(
  "tmp-data/meta_v3_onco_euro_overall_ChrAll_1_release.txt",
  select = c("Chr", "position", "Allele1", "Allele2", "Effect"),
  col.names = c("chr", "pos", "a1", "a0", "beta")
) %>%
  mutate(a0 = toupper(a0), a1 = toupper(a1))

info_snp <- snp_match(sumstats, bind_cols(map, beta2 = all_betas[, "185"]))
# 20,370,946 variants to be matched.
# 0 ambiguous SNPs have been removed.
# Some duplicates were removed.
# 1,039,112 variants have been matched; 0 were flipped and 534,105 were reversed.
ind <- which(info_snp$beta2 != 0)
qplot(beta, beta2, data = info_snp[ind, ]) + theme_bw(15)


#### TO EXPORT ####

write(res$pheno, "pheno.txt", ncolumns = 1)

bigreadr::fwrite2(mutate_if(res2[-2], is.numeric, round, digits = 4), "pred-cor.csv")
readLines("pred-cor.csv", n = 3)

csv <- bigreadr::fwrite2(bind_cols(map, as.data.frame(all_betas)), "PGS-effects.csv")
substr(readLines(csv, n = 3), 0, 100)
R.utils::gzip(csv)


#### TO USE INTERNALLY IN IPSYCH ####

ipsych <- snp_attach("../data_ipsych/dosage_ipsych2015.rds")
map_ipsych <- ipsych$map %>%
  transmute(chr = CHR, pos = POS, a1, a0 = a2)

info_snp <- snp_match(bind_cols(map, beta = 1), map_ipsych)
# all there and same sign (or all reversed)

all_pgs <- runonce::save_run(
  big_prodMat(ipsych$genotypes, all_betas, ind.col = info_snp$`_NUM_ID_`, ncores = 15),
  "tmp-data/UKBB_PGS_IPSYCH.rds"
)
colnames(all_pgs) <- colnames(all_betas)
all_cor <- cor(all_pgs, ipsych$fam$affection)
all_cor[order(abs(all_cor), decreasing = TRUE)[1:12], ]
#          poorer_health      log_age_first_sex            neuroticism
#            -0.10608606             0.09673061            -0.07771761
#                    318     self_harm_thoughts                 income
#            -0.07034071            -0.06825544             0.06746521
#              geek_time         log_waist_circ                log_BMI
#            -0.06365020            -0.06347000            -0.06317889
# less_happy_with_health                    496           less_alcohol
#            -0.06013341            -0.05900470            -0.05786728

## -> NEED TO REVERSE INDEED
csv <- bigreadr::fwrite2(bind_cols(setNames(ipsych$fam[1:2], c("ID1", "ID2")),
                                            as.data.frame(-all_pgs)),
                         "UKBB-PGS-iPSYCH.csv")
substr(readLines(csv, n = 3), 0, 100)
R.utils::gzip(csv)
