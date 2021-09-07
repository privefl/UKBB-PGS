library(bigsnpr)
ukb    <- snp_attach("data/UKBB_geno.rds")
ukb_aj <- snp_attach("data/UKBB_HM3_AJ.rds")
ind_match <- vctrs::vec_match(ukb_aj$map, ukb$map)
nona <- which(!is.na(ind_match))

fam <- snp_attach("data/UKBB_HM3.rds")$fam
ind_test <- which(fam$set == "test")
set.seed(1); ind_row <- unlist(
  by(ind_test, fam$group[ind_test], function(ind) sample(ind, 200)))

ukb$map$genetic.dist <- 0
ukb$fam <- snp_fake(nrow(G), 1)$fam
bed1 <- snp_writeBed(ukb, bedfile = "tmp-data/for_admixture.bed",
                     ind.row = ind_row, ind.col = ind_match[nona])

n_aj <- nrow(ukb_aj$genotypes)
set.seed(1); ind_row_aj <- sample(n_aj, 200)
ukb_aj$map$genetic.dist <- 0
ukb_aj$fam <- snp_fake(n_aj + 1e6, 1)$fam[1:n_aj + 1e6, ]
bed2 <- snp_writeBed(ukb_aj, bedfile = "tmp-data/for_admixture2.bed",
                     ind.row = ind_row_aj, ind.col = nona)

plink <- download_plink("tmp-data")
system(glue::glue("{plink} --memory 30e3",
                  " --bfile tmp-data/for_admixture",
                  " --bmerge tmp-data/for_admixture2",
                  " --out tmp-data/for_admixture_all"))

tgz <- runonce::download_file(
  "http://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz",
  dir = "tmp-data")
untar(tgz, exdir = "tmp-data")

admixture <- "tmp-data/dist/admixture_linux-1.3.0/admixture"
K <- 5
system(glue::glue(
  "{admixture} tmp-data/for_admixture_all.bed {K} -j{nb_cores()}"))

## K = 8
# Fst divergences between estimated populations:
#       Pop0	Pop1	Pop2	Pop3	Pop4	Pop5	Pop6
# Pop0
# Pop1	0.028
# Pop2	0.022	0.026
# Pop3	0.047	0.043	0.033
# Pop4	0.051	0.051	0.044	0.033
# Pop5	0.059	0.062	0.057	0.053	0.061
# Pop6	0.147	0.146	0.142	0.137	0.138	0.128
# Pop7	0.112	0.116	0.111	0.114	0.116	0.077	0.159

## K = 6
# Fst divergences between estimated populations:
#       Pop0	Pop1	Pop2	Pop3	Pop4
# Pop0
# Pop1	0.078
# Pop2	0.049	0.109
# Pop3	0.043	0.108	0.022
# Pop4	0.049	0.112	0.032	0.031
# Pop5	0.125	0.158	0.133	0.134	0.146

## K = 5
# Fst divergences between estimated populations:
#       Pop0	Pop1	Pop2	Pop3
# Pop0
# Pop1	0.085
# Pop2	0.125	0.158
# Pop3	0.038	0.109	0.133
# Pop4	0.037	0.107	0.141	0.022

group <- c(fam$group[sort(ind_row)], rep("Ashkenazi", 200))
Q <- read.table(paste0("for_admixture_all.", K, ".Q"))
main_comp_by_group <- by(Q, group, function(x) which.max(colSums(x)))
val_main_comp_group <- Q[cbind(seq_along(group), main_comp_by_group[group])]

library(dplyr)
tbl <-
  data.frame(grp = group, val = val_main_comp_group) %>%
  bind_cols(Q) %>%
  group_by(grp) %>%
  mutate(id = rank(-val, ties.method = "first"), val = NULL) %>%
  tidyr::pivot_longer(-c(grp, id))

library(ggplot2)
cbPalette <- tail(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                    "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), K)

ggplot(tbl) +
  geom_col(aes(id, value, color = name, fill = name)) +
  theme_bw(15) +
  scale_x_continuous(breaks = 1:10 * 20) +
  scale_color_manual(values = cbPalette) +
  scale_fill_manual(values = cbPalette) +
  facet_wrap(~ grp) +
  labs(x = "Individual # (ordered by main component of group)",
       y = "Ancestry proportion", color = "Ancestry", fill = "Ancestry") +
  theme(legend.position = "none")
# ggsave(paste0("figures/admixture_", K, ".pdf"), width = 12, height = 8)
