library(bigsnpr)
library(ggplot2)
NCORES <- 15

#### Project UKBB individuals onto PCA space from 1000G ####

obj.bed <- bed("../paper-ancestry-matching/data/ukbb.bed")
tarbz2 <- runonce::download_file(
  "https://evolbio.ut.ee/khazar/full_data_panel_for_khazar_paper.tar.bz2",
  dir = "tmp-data")
untar(tarbz2, exdir = "tmp-data")
bed.ref <- bed("tmp-data/khazar3_1-22_maf0.01mind0.035geno0.005.bed")
dim(bed.ref)  # 1774 x 270898

proj <- runonce::save_run(
  bed_projectPCA(bed.ref, obj.bed, match.min.prop = 0.1, k = 25, ncores = NCORES),
  file = "tmp-data/proj-UKBB-to-khazar.rds"
) # < 4 min

PC.ref <- predict(proj$obj.svd.ref)
proj2 <- proj$OADP_proj

info.ref <- readxl::read_excel(
  "~/Downloads/Behar_HumanBiology_2013_KhazarTest_Supplement.xlsx")
fam <- dplyr::left_join(bed.ref$fam, info.ref, by = c("sample.ID" = "Sample ID"))
is_AJ <- fam$Population == "Ashkenazi Jewish"

col <- rep(NA, nrow(fam))
col[fam$Population %in% "Ashkenazi Jewish"] <- "#E69F00"
col[fam$Population %in% c("Algerian Jewish", "Moroccan Jewish",
                          "Libyan Jewish", "Tunisian Jewish")] <- "#009E73"
col[fam$Population %in% c("Iranian Jewish", "Iraqi Jewish")] <- "#0072B2"
col[fam$Population %in% c("Italian Jewish", "Sephardi Jewish")] <- "#CC79A7"

plot_grid(plotlist = lapply(1:12, function(k) {
  k1 <- 2 * k - 1
  k2 <- 2 * k
  qplot(PC.ref[is.na(col), k1], PC.ref[is.na(col), k2],
        color = I("#999999"), alpha = I(0.3)) +
    geom_point(aes(PC.ref[, k1], PC.ref[, k2]),
               color = col, pch = ifelse(is_AJ, 17, 16)) +
    theme_bigstatsr(0.5) +
    labs(x = paste0("PC", k1), y = paste0("PC", k2)) +
    coord_equal() +
    theme(legend.position = "none")
}), nrow = 3)
# ggsave("figures/jews-ref.pdf", width = 9, height = 6)

AJ.center <- bigutilsr::geometric_median(PC.ref[is_AJ, ])

dist_to_AJ0 <- sqrt(bigutilsr:::rowSumsSq(sweep(PC.ref, 2, AJ.center, '-')))
sort(tapply(dist_to_AJ0, fam$Population, min))
# Ashkenazi Jewish     Italian Jewish    Tunisian Jewish    Sephardi Jewish
#         8.544761          12.665727          13.069149          18.292241
#  Moroccan Jewish    Algerian Jewish      Syrian Jewish            Italian
#        24.681102          26.407728          28.070543          29.311714

dist_to_AJ <- sqrt(bigutilsr:::rowSumsSq(sweep(proj2, 2, AJ.center, '-')))
hist(dist_to_AJ)
hist(dist_to_AJ[dist_to_AJ < 20])
THR <- 12.5

info <- readRDS("../paper-ancestry-matching/data/info_UKBB.rds")
PC <- dplyr::select(info, PC1:PC16)
stopifnot(nrow(PC) == length(dist_to_AJ))
fam_groups <- snp_attach("data/UKBB_HM3.rds")$fam
pop <- rep(NA, nrow(PC))
pop[match(fam_groups$eid, info$eid)] <- fam_groups$group
tapply(dist_to_AJ, pop, min)
# Caribbean          China          India           Iran          Italy
#  93.39589      102.64485       46.71575       31.58030       17.32431
#   Nigeria         Poland United Kingdom
# 110.48153       28.95825       24.12444
sum(dist_to_AJ < THR)  # 1780
pop[dist_to_AJ < THR] <- "Ashkenazi Jewish"

source("https://raw.githubusercontent.com/privefl/paper4-bedpca/master/code/plot_grid2.R")
ind <- sample(nrow(proj2), 50e3)
plot_grid2(plotlist = lapply(1:4, function(k) {
  k1 <- 2 * k - 1
  k2 <- 2 * k
  qplot(PC[ind, k1], PC[ind, k2], color = pop[ind],
        alpha = I(ifelse(is.na(pop[ind]), 0.2, 1))) +
    theme_bigstatsr(0.6) +
    labs(x = paste0("PC", k1), y = paste0("PC", k2), color = "Ancestry group")
}), nrow = 2, title_ratio = 0, legend_ratio = 0.18)
# ggsave("figures/PC-nine-groups.png", width = 9, height = 6)


print_table <- function(x) {
  tab <- sort(table(x, exclude = NULL), decreasing = TRUE)
  cat(paste0(names(tab), ": ", tab), sep = " || ")
}
print_table(as.character(info$country[dist_to_AJ < THR]))
# United Kingdom: 1008 || NA: 495 || USA: 97 || South Africa: 57 || Israel: 25
#  || Hungary: 13 || France: 8 || Canada: 7 || Ireland: 7 || Argentina: 5 ||
# Germany: 5 || Australia: 4 || Czech Republic: 4 || Russia: 4 || Caribbean: 3
#  || Romania: 3 || Zimbabwe: 3 || Belarus: 2 || Belgium: 2 || Italy: 2 ||
# New Zealand: 2 || Poland: 2 || Slovakia: 2 || Switzerland: 2 || ... with 1


#### Prepare genotype data ####

# Match indices in BGEN data
sample <- bigreadr::fread2("UKBB/ukb58024_imp_chr1_v3_s487296.sample")[-1, ]
info_AJ <- info[dist_to_AJ < THR, ]
ind.indiv <- match(info_AJ$eid, sample$ID_2)
info_AJ[is.na(ind.indiv), ] <- NA

# Filter for close relatedness
rel <- bigreadr::fread2("UKBB/ukb58024_rel_s488264.dat")
rel2 <- dplyr::filter(rel, ID1 %in% info_AJ$eid, Kinship > 2^-3.5)
info_AJ[info_AJ$eid %in% rel2$ID2, ] <- NA

length(sub <- which(!is.na(info_AJ$eid)))  # 1709 individuals

# Prepare variants to use
map <- snp_attach("data/UKBB_HM3.rds")$map
list_snp_id <- with(map, split(
  paste(chromosome, physical.pos, allele1, allele2, sep = "_"),
  chromosome))


# Prepare bigSNP objects
system.time(
  rds <- snp_readBGEN(
    bgenfiles   = glue::glue("UKBB/bgen/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKBB_HM3_AJ",
    ind_row     = ind.indiv[sub],
    ncores      = NCORES
  )
) # 47 min with 15 cores

ukbb <- snp_attach(rds)
library(dplyr)
df0 <- bigreadr::fread2(
  "UKBB/ukb41181.csv",
  select = c("eid", "34-0.0", "52-0.0", "22001-0.0", "21022-0.0", "189-0.0"),
  col.names = c("eid", "year", "month", "sex", "age", "deprivation_index")
) %>%
  mutate(date = (year - 1900) + (month - 0.5) / 12, year = NULL, month = NULL)


ukbb$fam <- left_join(info_AJ[sub, ], df0)
snp_save(ukbb)

# verif
qplot(PC3, PC4, data = ukbb$fam)
qplot(PC7, PC8, data = ukbb$fam)
