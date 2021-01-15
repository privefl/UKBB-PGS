library(bigsnpr)
ukbb <- snp_attach("data/UKBB_HM3.rds")
G <- ukbb$genotypes
ind_test <- which(ukbb$fam$set == "test")

NCORES <- 3
library(future.batchtools)
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES, mem = "25g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

bigassertr::assert_dir("assess-pred-ldpred2")

all_mod <- list.files("ldpred2", full.names = TRUE)
(all_mod_todo <- all_mod[
  !file.exists(file.path("assess-pred-ldpred2", basename(all_mod)))])

furrr::future_walk(all_mod_todo, function(.) {

  # . <- all_mod_todo[11]
  res <- readRDS(.)
  print(pheno <- sub("\\.rds", "", basename(.)))
  # str(res$ldpred2[[1]])

  library(dplyr)

  all_acf <- sapply(res$ldpred2, function(auto) {
    acf(tail(auto$path_p_est, 500), lag.max = 50, plot = FALSE)$acf
  })
  keep <- apply(abs(all_acf), 2, min) < 0.05

  all_h2 <- sapply(res$ldpred2, function(auto) auto$h2_est)
  h2 <- median(all_h2[keep])
  keep <- keep & between(all_h2, 0.7 * h2, 1.4 * h2)
  all_p <- sapply(res$ldpred2, function(auto) auto$p_est)
  p <- median(all_p[keep])
  keep <- keep & between(all_p, 0.5 * p, 2 * p)

  res2 <- tibble(
    pheno   = pheno,
    time    = res$time,
    node    = res$node,
    ldsc    = list(res$ldsc),
    h2_ldsc = list(res$ldsc[["h2"]] + c(0, -1.96, 1.96) * res$ldsc[["h2_se"]]),
    n_keep  = sum(keep)
  )

  if (res2$n_keep > 0) {

    path_h2 <- sapply(res$ldpred2[keep], function(auto) tail(auto$path_h2_est, 500))
    path_p <- sapply(res$ldpred2[keep], function(auto) tail(auto$path_p_est, 500))

    postp <- rowMeans(
      do.call("cbind", lapply(res$ldpred2[keep], function(auto) auto$postp_est)))
    beta <- rowMeans(
      do.call("cbind", lapply(res$ldpred2[keep], function(auto) auto$beta_est)))
    beta_sp <- rowMeans(
      do.call("cbind", lapply(res$ldpred2[keep], function(auto) auto$beta_est_sparse)))

    ind <- which(beta_sp != 0)
    pred_sp <- big_prodVec(G, beta_sp[ind], ind.row = ind_test, ind.col = ind,
                           ncores = NCORES)
    pred <- big_prodVec(G, beta, ind.row = ind_test, ncores = NCORES)

    res2 <- res2 %>%
      mutate(
        h2      = list(unname(quantile(path_h2, probs = c(0.5, 0.025, 0.975)))),
        p       = list(unname(quantile(path_p, probs = c(0.5, 0.025, 0.975)))),
        postp   = list(postp),
        beta    = list(beta),
        beta_sp = list(beta_sp),
        pred    = list(pred),
        pred_sp = list(pred_sp)
      )
  }

  saveRDS(res2, file.path("assess-pred-ldpred2", basename(.)))
})


library(dplyr)
length(files <- list.files("assess-pred-ldpred2", "\\.rds$", full.names = TRUE))
library(furrr)
(NCORES <- availableCores() - 1L)
plan("multisession", workers = NCORES)
res <- future_map_dfr(files, function(file) {
  readRDS(print(file)) %>%
    mutate(pheno = sub("\\.rds$", "", basename(file)))
})

res_pred <- res %>%
  mutate(nb_var_ldpred2 = sapply(beta_sp, function(beta) sum(beta != 0))) %>%
  select(-node, -time, -postp, -beta, -beta_sp) %>%
  filter(n_keep > 0) %>%
  tidyr::pivot_longer(c(pred, pred_sp), values_to = "pred") %>%
  mutate(sparse = (name == "pred_sp"), name = NULL)

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

res_pred$pcor <- purrr::pmap(res_pred[c("pred", "pheno")], function(pred, pheno) {
  # pred <- res_pred$pred[[29]]
  # pheno <- res_pred$pheno[[29]]
  y <- all_pheno[[print(pheno)]] + 0
  lapply(split(seq_along(pop), pop), function(ind) {
    if (sum(y[ind] != 0, na.rm = TRUE) == 0) return(rep(NA_real_, 3))
    bigstatsr::pcor(pred[ind], y[ind], covar[ind, ])
  })
})

res2 <- res_pred %>%
  select(-pred) %>%
  tidyr::unnest_longer("pcor", indices_to = "pop") %>%
  mutate(pop = factor(pop, levels = POP)) %>%
  tidyr::pivot_wider(names_from = "sparse", values_from = "pcor",
                     names_prefix = "pcor_") %>%
  tidyr::unnest_wider("pcor_FALSE", names_sep = "_") %>%
  tidyr::unnest_wider("pcor_TRUE", names_sep = "_") %>%
  print()

library(ggplot2)

ggplot(res2, aes(pcor_FALSE_1, pcor_TRUE_1)) +
  bigstatsr::theme_bigstatsr() +
  geom_point() +
  # geom_errorbar(aes(xmin = pcor_FALSE_2, xmax = pcor_FALSE_3), color = "green") +
  # geom_errorbar(aes(ymin = pcor_TRUE_2, ymax = pcor_TRUE_3), color = "green") +
  geom_abline(color = "red", linetype = 2) +
  facet_wrap(~ pop, nrow = 2) +
  coord_equal() +
  labs(x = "Not sparse", y = "Sparse")
# ggsave("figures/sparse-ldpred2.pdf", width = 10, height = 5.6)


res3 <- res_pred %>%
  filter(sparse) %>%
  tidyr::unnest_wider("pcor") %>%
  print()

OUTLIERS <- c("less_tanned", "darker_hair", "darker_skin", "red_hair",
              "log_bilirubin", "log_lipoA", "darker_hair0", "darker_skin0")
stopifnot(all(OUTLIERS %in% res$pheno))

slopes <- list("United Kingdom" = 1)

bigstatsr::plot_grid(plotlist = lapply(POP[-1], function(pop) {

  df <- tibble(pheno = res3$pheno, x = res3$`United Kingdom`, y = res3[[pop]]) %>%
    filter(purrr::map_dbl(y, ~ .[3] - .[2]) < 0.2,
           !startsWith(pheno, "LTFH")) %>%
    tidyr::unnest_wider("x", names_sep = "_") %>%
    tidyr::unnest_wider("y", names_sep = "_") %>%
    na.omit() %>%
    mutate(., eps = abs(resid(lm(y_1 ~ x_1 + 0, data = .))))

  robust_slope <- deming::deming(y_1 ~ x_1 + 0,
                                 data = filter(df, !pheno %in% OUTLIERS, y_2 != y_3),
                                 xstd = x_3 - x_2,
                                 ystd = y_3 - y_2)$coef[2]

  cat(pop, ": ", round(100 * (slopes[[pop]] <<- unname(robust_slope^2)), 1), "%\n", sep = "")

  ggplot(df, aes(x_1, y_1, label = pheno)) +
    bigstatsr::theme_bigstatsr(0.85) +
    geom_abline(color = "red", linetype = 2) +
    geom_hline(yintercept = 0, color = "red", linetype = 3) +
    geom_vline(xintercept = 0, color = "red", linetype = 3) +
    geom_abline(slope = robust_slope, color = "blue") +
    geom_errorbar(aes(xmin = x_2, xmax = x_3), width = 0, color = "green") +
    geom_errorbar(aes(ymin = y_2, ymax = y_3), width = 0, color = "green") +
    geom_point() +
    ggrepel::geom_text_repel(data = slice_max(df, eps, n = 5),
                             min.segment.length = 0, seed = 42, box.padding = 0.5) +
    labs(x = "United Kingdom", y = pop) +
    xlim(-0.1, 0.65) + ylim(-0.2, 0.65) +
    ggtitle(paste0(round(100 * unname(robust_slope^2), 1), "%")) +
    coord_equal()
}), nrow = 2, scale = 0.95)
# Poland: 94.1%
# Italy: 86.1%
# Iran: 72.7%
# India: 64.8%
# China: 49.3%
# Caribbean: 25.2%
# Nigeria: 18.4%
# ggsave("figures/ldpred2-ancestry.pdf", width = 12, height = 7)


res_PLR <- purrr::map_dfr(list.files("assess-pred", full.names = TRUE), function(file) {
  readRDS(file) %>%
    mutate(pheno = sub("\\.rds$", "", basename(file))) %>%
    arrange(validation_loss) %>%
    select(-validation_loss) %>%
    slice(1)
}) %>%
  tidyr::pivot_longer(`United Kingdom`:Nigeria,
                      names_to = "pop", values_to = "pcor_PLR") %>%
  print()

res_ldpred2 <- res_pred %>%
  filter(sparse) %>%
  select(-sparse, -pred, -ldsc) %>%
  tidyr::unnest_longer("pcor", values_to = "pcor_ldpred2", indices_to = "pop") %>%
  print()

res_all <-
  inner_join(res_PLR, res_ldpred2, by = c("pheno", "pop")) %>%
  tidyr::unnest_wider("pcor_PLR", names_sep = "_") %>%
  tidyr::unnest_wider("pcor_ldpred2", names_sep = "_") %>%
  tidyr::unnest_wider("h2", names_sep = "_") %>%
  tidyr::unnest_wider("p", names_sep = "_") %>%
  mutate(pop = factor(pop, levels = POP))

p1 <- ggplot(filter(res_all, pcor_PLR_3 - pcor_PLR_2 < 0.3),
             aes(pcor_PLR_1, pcor_ldpred2_1)) +
  bigstatsr::theme_bigstatsr() +
  geom_point() +
  geom_abline(color = "red", linetype = 2) +
  geom_hline(yintercept = 0, color = "blue", linetype = 3) +
  geom_vline(xintercept = 0, color = "blue", linetype = 3) +
  scale_x_continuous(breaks = -5:5 / 5, minor_breaks = -10:10 / 10) +
  facet_wrap(~ pop, nrow = 2) +
  coord_equal() +
  labs(x = "LASSO", y = "LDpred2-auto-sp")

filter(res_all, pcor_PLR_3 < pcor_ldpred2_2 | pcor_PLR_2 > pcor_ldpred2_3) %>%
  select(-h2_ldsc, -h2_2, -h2_3, -p_2, -p_3, -power_adaptive, -power_scale) %>%
  print(n = Inf, width = Inf)
#    nb_var pheno               pop            pcor_PLR_1 pcor_PLR_2 pcor_PLR_3
#    <int> <chr>               <fct>               <dbl>      <dbl>      <dbl>
# 1    690 615                 China              0.0650    0.00739    0.122
# 2     20 hard_falling_asleep United Kingdom    -0.0349   -0.0742     0.00452
# 3 156534 height              United Kingdom     0.634     0.626      0.643
# 4  33210 log_bilirubin       Nigeria            0.546     0.523      0.569
#   n_keep    h2_1      p_1 pcor_ldpred2_1 pcor_ldpred2_2 pcor_ldpred2_3
#    <int>   <dbl>    <dbl>          <dbl>          <dbl>          <dbl>
# 1      1 0.00945 0.0938          -0.0509        -0.108         0.00677
# 2      2 0.0402  0.00244          0.0707         0.0314        0.110
# 3     26 0.546   0.0226           0.613          0.605         0.622
# 4     30 0.361   0.000481         0.475          0.449         0.500

ggplot(filter(res_all, pop == "United Kingdom"),
       aes(pcor_PLR_1, pcor_ldpred2_1, color = p_1)) +
  bigstatsr::theme_bigstatsr() +
  geom_point() +
  geom_abline(color = "red", linetype = 2) +
  geom_hline(yintercept = 0, color = "blue", linetype = 3) +
  geom_vline(xintercept = 0, color = "blue", linetype = 3) +
  scale_x_continuous(breaks = -5:5 / 5, minor_breaks = -10:10 / 10) +
  scale_color_viridis_c(direction = -1, trans = "log10") +
  facet_wrap(~ cut(h2_1, breaks = c(0, 0.05, 0.1, 0.2, 0.3, 0.5, Inf)), scales = "free") +
  # coord_equal() +
  labs(x = "LASSO", y = "LDpred2-auto-sp", color = "p")

res_UK <- filter(res_all, pop == "United Kingdom")
p2 <- ggplot(res_UK, aes(pcor_PLR_1, pcor_ldpred2_1, color = h2_1)) +
  bigstatsr::theme_bigstatsr() +
  geom_abline(color = "red", linetype = 2) +
  geom_hline(yintercept = 0, color = "red", linetype = 3) +
  geom_vline(xintercept = 0, color = "red", linetype = 3) +
  geom_errorbar(aes(xmin = pcor_PLR_2, xmax = pcor_PLR_3)) +
  geom_errorbar(aes(ymin = pcor_ldpred2_2, ymax = pcor_ldpred2_3)) +
  geom_point() +
  scale_x_continuous(breaks = -5:5 / 5, minor_breaks = -10:10 / 10) +
  scale_color_viridis_c(direction = -1, trans = "log10") +
  facet_wrap(~ cut(p_1, breaks = 10^-(6:0))) +
  coord_equal() +
  labs(x = "LASSO", y = "LDpred2-auto-sp", color = "h2")

plot_grid(p1, p2, scale = 0.98, ncol = 1, rel_heights = c(4, 4.5),
          labels = c("A", "B"), label_size = 18)
# ggsave("figures/PLR-ldpred2.pdf", width = 10, height = 12)

M <- length(res$beta_sp[[1]])
ggplot(res_UK, aes(p_1, nb_var_ldpred2 / M, color = h2_1)) +
  theme_bigstatsr() +
  geom_abline(color = "red", linetype = 2) +
  geom_point() +
  scale_color_viridis_c(trans = "log10", direction = -1) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Proportion of causal variants",
       y = "Proportion of non-zero variables",
       color = "SNP h2")
# ggsave("figures/sparsity-ldpred2.pdf", width = 9, height = 6)

ggplot(res_UK, aes(p_1, nb_var / M, color = pcor_PLR_1)) +
  theme_bigstatsr() +
  geom_abline(color = "red", linetype = 2) +
  geom_point() +
  scale_color_viridis_c(trans = "log10", direction = -1) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Proportion of causal variants",
       y = "Proportion of non-zero variables",
       color = "Partial\ncorrelation")
# ggsave("figures/sparsity-PLR.pdf", width = 9, height = 6)
