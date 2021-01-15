library(bigsnpr)
library(bigreadr)
library(dplyr)
library(ggplot2)

rel <- fread2("UKBB/ukb58024_rel_s488264.dat")
fam <- snp_attach("data/UKBB_HM3.rds")$fam
is_rel <- fam$eid %in% c(rel$ID1, rel$ID2)

## Ancestry groups
ind_used <- which(fam$set == "test")
ind.train <- intersect(ind_used, which(!is_rel))
table(pop_UKBB <- fam$group[ind.train])
# Caribbean          China          India           Iran
#      2089           1699           5557           1154
#  Italy        Nigeria         Poland United Kingdom
#   5851           3625           3763          14806

ukbb <- snp_attach("data/UKBB_geno.rds")
G <- ukbb$genotypes
CHR <- as.integer(ukbb$map$chromosome)
POS <- ukbb$map$physical.pos
(NCORES <- parallelly::availableCores() - 1L)
obj.svd <- runonce::save_run(
  snp_autoSVD(G, ind.row = ind.train, k = 50, ncores = NCORES,
              infos.chr = CHR, infos.pos = POS),
  file = "tmp-data/svd_UKBB.rds"
) # <2H with 63 cores

plot(obj.svd)
legend <- cowplot::get_legend(
  plot(obj.svd, type = "scores", coeff = 0.7) +
    aes(color = pop_UKBB) +
    labs(color = "Population")
)
plot_grid(plotlist = c(lapply(10:20, function(k) {
  plot(obj.svd, type = "scores", scores = 2 * k - 1:0, coeff = 0.4) +
    aes(color = pop_UKBB) +
    theme(legend.position = "none")
})), legend)
# ggsave("figures/PC-scores-restricted.png", width = 11, height = 7)

PC <- predict(obj.svd)[, 1:32]

slopes <- list(
  "United Kingdom" = 1,
  "Poland" = 0.938,
  "Italy" = 0.856,
  "Iran" = 0.722,
  "India" = 0.647,
  "China" = 0.486,
  "Caribbean" = 0.252,
  "Nigeria" = 0.18
)

centers <- bigutilsr::geometric_median(PC, by_grp = pop_UKBB)
dist_to_UK <- as.matrix(dist(centers))[, "United Kingdom"]

qplot(dist_to_UK[names(slopes)], unlist(slopes), size = I(2)) +
  theme_bigstatsr() +
  geom_smooth(method = "lm") +
  ggrepel::geom_label_repel(aes(label = names(slopes)), min.segment.length = 0) +
  # scale_y_log10() +
  scale_y_continuous(breaks = 0:5 / 5, minor_breaks = 0:10 / 10) +
  labs(x = "PC distance to UK", y = "Relative predictive performance with UK")
# ggsave("figures/ratio-dist-restricted.pdf", width = 8, height = 6)
