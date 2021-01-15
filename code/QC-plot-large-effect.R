library(bigsnpr)
ukbb <- snp_attach("data/UKBB_HM3.rds")
G <- ukbb$genotypes
is_train <- (ukbb$fam$set == "train")
ind_csv <- match(ukbb$fam$eid, readRDS("data/csv_eid.rds"))

pheno <- "data/ukbb-quant-pheno/log_bilirubin.rds"
y <- readRDS(pheno)
y.train <- y[ind_csv[is_train]]
n <- sum(!is.na(y.train))

gwas <- readRDS(file.path("GWAS", basename(pheno)))
gwas$lpval <- -predict(gwas)
gwas$id <- rows_along(gwas)
gwas_top <- dplyr::slice_max(gwas, lpval, n = 1000)
sd <- sqrt(big_colstats(G, ind.row = which(is_train), ind.col = gwas_top$id)$var)

se <- gwas_top$std.err
sd_ss <- sd(y.train, na.rm = TRUE) / sqrt(n * se^2)
beta <- gwas_top$estim
sd_ss2 <- sd(y.train, na.rm = TRUE) / sqrt(n * se^2 + beta^2)

library(ggplot2)
theme_set(theme_bigstatsr(0.9))
plot_grid(
  qplot(sd, sd_ss) +
    geom_abline(color = "red", linetype = 2) +
    coord_equal() +
    labs(x = "SD (from genotypes)",
         y = "SD (from summary statistics, previous formula)"),
  qplot(sd, sd_ss2) +
    geom_abline(color = "red", linetype = 2) +
    coord_equal() +
    labs(x = "SD (from genotypes)",
         y = "SD (from summary statistics, updated formula)"),
  scale = 0.95, nrow = 1, labels = c("A", "B"), label_size = 16
)
# ggsave("figures/qc-plot-new-formula.pdf", width = 10.5, height = 6)
