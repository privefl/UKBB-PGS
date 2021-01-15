library(bigsnpr)
ukbb <- snp_attach("data/UKBB_HM3.rds")
CHR <- as.integer(ukbb$map$chromosome)
POS <- ukbb$map$physical.pos
(NCORES <- parallelly::availableCores() - 1L)
covar <- as.matrix(dplyr::select(ukbb$fam, -eid, -group, -set))
ind_csv <- match(ukbb$fam$eid, readRDS("data/csv_eid.rds"))

pheno <- "log_bilirubin"
pheno <- "log_lipoA"
pheno <- "apoB"

y <- readRDS(paste0("data/ukbb-quant-pheno/", pheno, ".rds"))[ind_csv] + 0
ind_nona <- which(!is.na(y))

gwas <- readRDS(paste0("GWAS/", pheno, ".rds"))
top <- which.max(abs(gwas$score))
chr_top <- CHR[top]
pos_top <- POS[top]

library(bigreadr)
library(dplyr)
snp_id <- fread2(paste0("UKBB/mfi/ukb_mfi_chr", chr_top, "_v3.txt"),
                 select = c(3:6, 8),
                 col.names = c("pos", "a0", "a1", "maf", "info")) %>%
  filter(abs(pos - pos_top) < 500e3, info > 0.3) %>%
  with(paste(chr_top, pos, a0, a1, sep = "_")) %>%
  print()

sample <- fread2("UKBB/ukb58024_imp_chr1_v3_s487296.sample")[-1, ]

system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles   = glue::glue("UKBB/bgen/ukb_imp_chr{chr_top}_v3.bgen"),
    list_snp_id = list(snp_id),
    backingfile = (tmp <- tempfile(tmpdir = "tmp-data")),
    ind_row     = match(ukbb$fam$eid, sample$ID_2),
    ncores      = NCORES
  )
) # < 2 min

zoom <- snp_attach(rds)
G <- zoom$genotypes
POS2 <- zoom$map$physical.pos
mean(in_hm3 <- vctrs::vec_in(zoom$map, ukbb$map))  # < 2%

all_r2_pop <- purrr::map(split(rows_along(G), ukbb$fam$group), function(ind_pop) {
  cat(".")
  ind.train <- intersect(ind_nona, ind_pop)
  gwas <- big_univLinReg(G, y[ind.train], ind.train = ind.train,
                         covar.train = covar[ind.train, ],
                         ncores = `if`(length(ind.train) > 10e3, NCORES, 1))
  N <- as.integer(sub(".*df = ([0-9]+).*", "\\1", body(attr(gwas, "predict"))[2]))
  t <- gwas$score
  list(beta = gwas$estim, beta_se = gwas$std.err,
       r2 = t^2 / (N + t^2), maf = snp_MAF(G, ind.train))
})

POP <- c("United Kingdom", "Poland", "Italy", "Iran",
         "India", "China", "Caribbean", "Nigeria")

res <- all_r2_pop %>%
  purrr::transpose() %>%
  as_tibble() %>%
  mutate(pop = factor(names(all_r2_pop), levels = POP),
         pos = list(POS2), in_hm3 = list(in_hm3)) %>%
  tidyr::unnest() %>%
  na.omit() %>%
  print()

res %>%
  group_by(pop) %>%
  summarize(max_r2 = max(r2), max_r2_hm3 =  max(r2[in_hm3]))
### bilirubin
#   pop            max_r2 max_r2_hm3
# 1 United Kingdom  0.298      0.297
# 2 Poland          0.321      0.320
# 3 Italy           0.313      0.311
# 4 Iran            0.260      0.259
# 5 India           0.294      0.294
# 6 China           0.138      0.138
# 7 Caribbean       0.297      0.295
# 8 Nigeria         0.313      0.313
### lipoA
#   pop            max_r2 max_r2_hm3
# 1 United Kingdom 0.293      0.0524
# 2 Poland         0.228      0.0399
# 3 Italy          0.170      0.0325
# 4 Iran           0.0846     0.0277
# 5 India          0.110      0.0782
# 6 China          0.163      0.163
# 7 Caribbean      0.0494     0.0337
# 8 Nigeria        0.0468     0.0384
### apoB
#   pop            max_r2 max_r2_hm3
# 1 United Kingdom 0.0555    0.0314
# 2 Poland         0.0516    0.0293
# 3 Italy          0.0452    0.0236
# 4 Iran           0.0423    0.0248
# 5 India          0.0227    0.00989
# 6 China          0.0559    0.0487
# 7 Caribbean      0.0848    0.0351
# 8 Nigeria        0.0716    0.0288

library(ggplot2)
ggplot(res) +
  theme_bigstatsr(0.9) +
  geom_vline(xintercept = pos_top / 1e6, linetype = 3) +
  geom_point(aes(pos / 1e6, r2, color = in_hm3, alpha = I(in_hm3 + 0.4))) +
  facet_wrap(~ pop) +
  theme(legend.position = c(0.85, 0.12),
        panel.spacing = unit(1, "lines"),
        plot.margin = margin(r = 1, 0.5, 0.5, 0.5, unit = "lines")) +
  labs(x = paste0("Physical position (Mb, chromosome ", chr_top, ")"),
       y = "Variance explained per variant", color = "In HapMap3?")
# ggsave(paste0("figures/zoom_", pheno, ".png"), width = 14, height = 8)

ggplot(res) +
  theme_bigstatsr(0.9) +
  geom_point(aes(maf, r2)) +
  facet_wrap(~ pop) +
  theme(panel.spacing = unit(1, "lines"),
        plot.margin = margin(r = 1, 0.5, 0.5, 0.5, unit = "lines")) +
  labs(x = "MAF", y = "Variance explained per variant")
# ggsave(paste0("figures/zoom_maf_", pheno, ".png"), width = 14, height = 8)

res %>%
  group_by(pos) %>%
  mutate(r2_max = mean(r2)) %>%
  ungroup() %>%
  slice_max(r2_max, n = 2 * length(POP) + 1) %>%
  print() %>%
  ggplot(aes(beta, r2, color = pop, pch = paste0(chr_top, ":", pos))) +
  theme_bigstatsr(0.9) +
  geom_errorbar(aes(xmin = beta - 2 * beta_se, xmax = beta + 2 * beta_se)) +
  geom_point(size = 2) +
  labs(x = "Effect size (+/- 2SE)", y = "Variance explained",
       color = "Ancestry", pch = "BP (GRCh37)")
# ggsave(paste0("figures/top3_", pheno, ".pdf"), width = 8, height = 5.5)


# cleanup
file.remove(paste0(tmp, c(".bk", ".rds")))
