library(bigsnpr)
ukbb <- snp_attach("data/UKBB_large.rds")
G <- ukbb$genotypes
fam <- snp_attach("data/UKBB_HM3.rds")$fam
ind_test <- which(fam$set == "test")

NCORES <- 3
library(future.batchtools)
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES, mem = "25g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

bigassertr::assert_dir("assess-pred-ldpred2-large")

all_mod <- list.files("ldpred2-large", full.names = TRUE)
(all_mod_todo <- all_mod[
  !file.exists(file.path("assess-pred-ldpred2-large", basename(all_mod)))])

furrr::future_walk(all_mod_todo, function(.) {

  # . <- all_mod_todo[11]
  res <- readRDS(.)
  print(pheno <- sub("\\.rds", "", basename(.)))
  # str(res$ldpred2[[1]])

  library(dplyr)

  # should probably remove this filtering on ACF
  all_acf <- sapply(res, function(auto) {
    acf(tail(auto$path_p_est, 500), lag.max = 50, plot = FALSE)$acf
  })
  keep <- apply(abs(all_acf), 2, min) < 0.05

  all_h2 <- sapply(res, function(auto) auto$h2_est)
  h2 <- median(all_h2[keep])
  keep <- keep & between(all_h2, 0.7 * h2, 1.4 * h2)
  all_p <- sapply(res, function(auto) auto$p_est)
  p <- median(all_p[keep])
  keep <- keep & between(all_p, 0.5 * p, 2 * p)

  res2 <- tibble(
    pheno   = pheno,
    n_keep  = sum(keep)
  )

  if (res2$n_keep > 0) {

    path_h2 <- sapply(res[keep], function(auto) tail(auto$path_h2_est, 500))
    path_p <- sapply(res[keep], function(auto) tail(auto$path_p_est, 500))

    postp <- rowMeans(
      do.call("cbind", lapply(res[keep], function(auto) auto$postp_est)))
    beta_sp <- rowMeans(
      do.call("cbind", lapply(res[keep], function(auto) auto$beta_est_sparse)))

    ind <- which(beta_sp != 0)
    ind_top <- readRDS(paste0("gwas-large/", pheno, "_top.rds"))$id
    pred_sp <- big_prodVec(G, beta_sp[ind], ind.row = ind_test,
                           ind.col = ind_top[ind], ncores = NCORES)

    res2 <- res2 %>%
      mutate(
        h2      = list(unname(quantile(path_h2, probs = c(0.5, 0.025, 0.975)))),
        p       = list(unname(quantile(path_p, probs = c(0.5, 0.025, 0.975)))),
        postp   = list(postp),
        beta_sp = list(beta_sp),
        pred_sp = list(pred_sp)
      )
  }

  saveRDS(res2, file.path("assess-pred-ldpred2-large", basename(.)))
})


library(dplyr)
length(files <- list.files("assess-pred-ldpred2-large", "\\.rds$", full.names = TRUE))

res <- purrr::map_dfr(files, function(file) {
  file2 <- sub("-large", "", file)
  res <- bind_rows(
    mutate(readRDS(file),  varset = "top1M"),
    mutate(readRDS(file2), varset = "HM3")
  ) %>%
    mutate(pheno = sub("\\.rds$", "", basename(file)))
  `if`(any(res$n_keep == 0), NULL, res)
})

res_pred <- res %>%
  mutate(nb_var_ldpred2 = sapply(beta_sp, function(beta) sum(beta != 0))) %>%
  select(pheno, n_keep, h2, p, pred_sp, nb_var_ldpred2, varset) %>%
  group_by(pheno) %>%
  filter(any(n_keep > 0))

library(bigsnpr)
fam <- snp_attach("data/UKBB_HM3.rds")$fam
fam_test <- fam[fam$set == "test", ]
covar <- dplyr::select(fam_test, -eid, -group, -set)
ind_csv_test <- match(fam_test$eid, readRDS("data/csv_eid.rds"))

pheno_files <- c(list.files("data/ukbb-quant-pheno", full.names = TRUE),
                 list.files("data/ukbb-binary-pheno", full.names = TRUE))
all_pheno <- lapply(pheno_files, function(.) readRDS(.)[ind_csv_test])
names(all_pheno) <- sub("\\.rds$", "", basename(pheno_files))

POP <- c("United Kingdom", "Poland", "Italy", "Iran",
         "India", "China", "Caribbean", "Nigeria")
pop <- factor(fam_test$group, levels = POP)

res_pred$pcor <- purrr::pmap(res_pred[c("pred_sp", "pheno")], function(pred_sp, pheno) {
  # pred <- res_pred$pred[[29]]
  # pheno <- res_pred$pheno[[29]]
  y <- all_pheno[[print(pheno)]] + 0
  lapply(split(seq_along(pop), pop), function(ind) {
    if (sum(y[ind] != 0, na.rm = TRUE) == 0) return(rep(NA_real_, 3))
    bigstatsr::pcor(pred_sp[ind], y[ind], covar[ind, ])
  })
})

res2 <- res_pred %>%
  select(-pred_sp) %>%
  tidyr::unnest_longer("pcor", indices_to = "pop") %>%
  mutate(pop = factor(pop, levels = POP)) %>%
  tidyr::unnest_wider("pcor", names_sep = "_") %>%
  print()

library(ggplot2)

ggplot(res2, aes(pop, pcor_1, fill = varset)) +
  bigstatsr::theme_bigstatsr(0.9) +
  geom_col(position = position_dodge(), alpha = 0.5, color = "black") +
  geom_errorbar(aes(ymin = pcor_2, ymax = pcor_3),
                position = position_dodge(width = 0.9), width = 0.3) +
  facet_wrap(~pheno, scales = "free_y", ncol = 3) +
  geom_hline(yintercept = 0) +
  labs(x = "Population", y = "Partial correlation", fill = "Set of variants") +
  theme(axis.text.x = element_text(angle = 40, hjust = 0.7, vjust = 0.9),
        legend.position = c(0.85, 0.12))
# ggsave("figures/ldpred2-large.pdf", width = 13, height = 9)

res_pred %>%
  select(pheno, varset, h2, p) %>%
  tidyr::unnest_wider("h2", names_sep = "_") %>%
  tidyr::unnest_wider("p", names_sep = "_") %>%
  print(n = Inf) %>%
  mutate_if(is.numeric, ~ as.character(signif(., 3))) %>%
  transmute(pheno, varset,
            h2 = glue::glue("{h2_1} [{h2_2}-{h2_3}]"),
            p = glue::glue("{p_1} [{p_2}-{p_3}]")) %>%
  xtable::xtable(align = "|l|c|c|c|c|") %>%
  print(include.rownames = FALSE)
#   pheno         varset   h2_1   h2_2   h2_3      p_1      p_2      p_3
# 1 174.1         top1M  0.0889 0.0860 0.0920 0.00760  0.00678  0.00841
# 2 174.1         HM3    0.0299 0.0264 0.0334 0.000881 0.000636 0.00117
# 3 185           top1M  0.113  0.109  0.116  0.00819  0.00743  0.00906
# 4 185           HM3    0.0381 0.0343 0.0423 0.000784 0.000588 0.00105
# 5 411.4         top1M  0.0641 0.0624 0.0659 0.0152   0.0138   0.0168
# 6 411.4         HM3    0.0401 0.0379 0.0422 0.00457  0.00397  0.00526
# 7 apoB          top1M  0.269  0.265  0.272  0.0533   0.0498   0.0568
# 8 apoB          HM3    0.163  0.160  0.166  0.00132  0.00119  0.00145
# 9 height        top1M  0.482  0.479  0.486  1.00     1.00     1.00
# 10 height        HM3    0.546  0.541  0.552  0.0226   0.0218   0.0235
# 11 log_bilirubin top1M  0.301  0.267  0.363  0.214    0.195    0.227
# 12 log_bilirubin HM3    0.361  0.357  0.365  0.000481 0.000423 0.000545
# 13 log_BMI       top1M  0.173  0.171  0.176  1.00     1.00     1.00
# 14 log_BMI       HM3    0.263  0.260  0.266  0.0426   0.0404   0.0446
# 15 log_lipoA     top1M  0.696  0.689  0.702  0.0116   0.0110   0.0122
# 16 log_lipoA     HM3    0.340  0.336  0.345  0.000229 0.000192 0.000268

# % latex table generated in R 3.6.1 by xtable 1.8-4 package
# \begin{table}[ht]
# \centering
# \begin{tabular}{|c|c|c|c|}
# \hline
# pheno & varset & h2 & p \\
# \hline
# 174.1 & top1M & 0.0889 [0.086-0.092] & 0.0076 [0.00678-0.00841] \\
# 174.1 & HM3 & 0.0299 [0.0264-0.0334] & 0.000881 [0.000636-0.00117] \\
# 185 & top1M & 0.113 [0.109-0.116] & 0.00819 [0.00743-0.00906] \\
# 185 & HM3 & 0.0381 [0.0343-0.0423] & 0.000784 [0.000588-0.00105] \\
# 411.4 & top1M & 0.0641 [0.0624-0.0659] & 0.0152 [0.0138-0.0168] \\
# 411.4 & HM3 & 0.0401 [0.0379-0.0422] & 0.00457 [0.00397-0.00526] \\
# apoB & top1M & 0.269 [0.265-0.272] & 0.0533 [0.0498-0.0568] \\
# apoB & HM3 & 0.163 [0.16-0.166] & 0.00132 [0.00119-0.00145] \\
# height & top1M & 0.482 [0.479-0.486] & 1 [1-1] \\
# height & HM3 & 0.546 [0.541-0.552] & 0.0226 [0.0218-0.0235] \\
# log\_bilirubin & top1M & 0.301 [0.267-0.363] & 0.214 [0.195-0.227] \\
# log\_bilirubin & HM3 & 0.361 [0.357-0.365] & 0.000481 [0.000423-0.000545] \\
# log\_BMI & top1M & 0.173 [0.171-0.176] & 1 [1-1] \\
# log\_BMI & HM3 & 0.263 [0.26-0.266] & 0.0426 [0.0404-0.0446] \\
# log\_lipoA & top1M & 0.696 [0.689-0.702] & 0.0116 [0.011-0.0122] \\
# log\_lipoA & HM3 & 0.34 [0.336-0.345] & 0.000229 [0.000192-0.000268] \\
# \hline
# \end{tabular}
# \end{table}
