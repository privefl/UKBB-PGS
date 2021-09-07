library(dplyr)
library(bigsnpr)
fam <- snp_attach("data/UKBB_HM3.rds")$fam
is_train <- (fam$set == "train")
ind_csv <- match(fam$eid, readRDS("data/csv_eid.rds"))

sample_size_binary <- purrr::map_dfr(
  list.files("data/ukbb-binary-pheno", full.names = TRUE), ~ {
    y <- readRDS(.)[ind_csv] + 0
    ind.train <- which(!is.na(y) & is_train)
    tibble(pheno = sub("\\.rds$", "", basename(.)),
           N_case = sum(y[ind.train]),
           N_control = sum(!y[ind.train]))
  })

sample_size_quant <- purrr::map_dfr(
  list.files("data/ukbb-quant-pheno", full.names = TRUE), ~ {
    y <- readRDS(.)[ind_csv] + 0
    ind.train <- which(!is.na(y) & is_train)
    tibble(pheno = sub("\\.rds$", "", basename(.)),
           N = length(ind.train))
  })

length(files <- list.files("assess-pred-ldpred2", "\\.rds$", full.names = TRUE))
all_res0 <- purrr::map_dfr(files, function(file) {
  readRDS(print(file)) %>%
    mutate(pheno = sub("\\.rds$", "", basename(file)))
})

# LDSC intercept
all_ldsc_int <- all_res0 %>%
  pull(ldsc) %>%
  purrr::map_dbl("int")

hist(all_ldsc_int, "FD")
all_res0[all_ldsc_int > 1.2, ]  # standing and sitting heights


# Heritability
all_res <- all_res0 %>%
  filter(!startsWith(pheno, "LTFH_")) %>%
  mutate(ldpred2_h2 = purrr::map_dbl(h2, 1),
         ldpred2_p  = purrr::map_dbl(p,  1)) %>%
  select(-time, -node, -h2_ldsc, -h2, -p,
         -postp, -beta, -beta_sp, -pred, -pred_sp) %>%
  left_join(sample_size_binary, by = "pheno") %>%
  left_join(sample_size_quant, by = "pheno") %>%
  relocate(c(N, N_case, N_control), .after = pheno) %>%
  tidyr::unnest_wider("ldsc", names_sep = "_") %>%
  rowwise() %>%
  mutate_at(c("ldsc_h2", "ldsc_h2_se", "ldpred2_h2"), ~ {
    K <- N_case / (N_case + N_control)
    `if`(is.na(K), ., . * coef_to_liab(K, K))
  }) %>%
  mutate_at(c(5:8, 10:11), signif, digits = 3) %>%
  rename(ldpred2_n_keep = n_keep) %>%
  ungroup()

writexl::write_xlsx(all_res, "phenotype-info.xlsx")
bigreadr::fwrite2(all_res, "phenotype-info.csv")
