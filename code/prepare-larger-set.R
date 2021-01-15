library(bigreadr)
library(dplyr)
library(data.table)

# UKBB
map_ukbb <- do.call("rbind", lapply(1:22, function(chr) {
  print(chr)
  mfi <- fread2(paste0("UKBB/mfi/ukb_mfi_chr", chr, "_v3.txt"),
                select = c(3:6, 8),
                col.names = c("pos", "a1", "a0", "maf", "info_ukbb"))
  cbind(chr = chr, subset(mfi, maf > 0.01, info_ukbb > 0.6))
}))
str(map_ukbb)

# iPSYCH (2015 + 2012)
map_ipsych <- bigsnpr::snp_attach("../data_ipsych/dosage_ipsych2015.rds")$map %>%
  transmute(chr = CHR, pos = POS, info_2012, freq_2012, info_2015, freq_2015)

merged <- merge(as.data.table(map_ukbb), map_ipsych, all.x = TRUE)
setkey(merged, chr, pos)

wd <- setwd("tmp-data")

#### BRCA ####
# download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/MichailidouK_29059683_GCST004988/oncoarray_bcac_public_release_oct17%20(1).txt.gz",
#               destfile = "sumstats_BRCA.txt.gz")
sumstats_brca <- fread2(
  "sumstats_BRCA.txt.gz", na.strings = "NULL",
  select = c("chr", "position_b37", "bcac_onco2_r2", "bcac_icogs2_r2", "bcac_gwas_all_eaf_controls"),
  col.names = c("chr", "pos", "info_onco_brca", "info_icogs_brca", "freq_brca"))
str(sumstats_brca)

merged <- merge(merged, sumstats_brca, all.x = TRUE)

#### PRCA ####
# download.file("http://practical.icr.ac.uk/blog/wp-content/uploads/uploadedfiles/oncoarray/MetaSummaryData/meta_v3_onco_euro_overall_ChrAll_1_release.zip",
#               destfile = "sumstats_PRCA.zip")
# unzip("sumstats_PRCA.zip"); file.remove("sumstats_PRCA.zip")
sumstats_prca <- fread2(
  "meta_v3_onco_euro_overall_ChrAll_1_release.txt",
  select = c("Chr", "position", "OncoArray_imputation_r2", "Freq1"),
  col.names = c("chr", "pos", "info_onco_prca", "freq_prca")
)
str(sumstats_prca)

merged <- merge(merged, sumstats_prca, all.x = TRUE)

#### CAD ####
# download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/NikpayM_26343387_GCST003116/cad.add.160614.website.txt",
#               destfile = "sumstats_CAD.txt")
sumstats_cad <- fread2("sumstats_CAD.txt",
                       select = c("chr", "bp_hg19", "median_info", "effect_allele_freq"),
                       col.names = c("chr", "pos", "info_cad", "freq_cad"))

merged <- merge(merged, sumstats_cad, all.x = TRUE)

#### T1D ####
# To request a download link, go to
# https://datadryad.org/stash/dataset/doi:10.5061/dryad.ns8q3
# download.file("http://merritt.cdlib.org/cloudcontainer/mrtstore2/35227889.tar.gz",
#               destfile = "tmp-data/sumstats_T1D.tar.gz")
# untar("tmp-data/sumstats_T1D.tar.gz", exdir = "T1D")
sumstats_t1d <- fread2(
  paste0("T1D/meta_chr_", 1:22),
  select = c("chromosome", "position", "info_score.I", "info_score.A",
             "EUR_MAF_1kG", "controls_maf.I", "controls_maf.A"),
  col.names = c("chr", "pos", "info_illu_t1d", "info_affy_t1d",
                "freq_1000G", "freq_illu_t1d", "freq_affy_t1d"))

merged <- merge(merged, sumstats_t1d, all.x = TRUE)

setwd(wd)

## Filter for info score
info_mean <- merged %>%
  select(starts_with("info_")) %>%
  mutate_all(~ ifelse(is.na(.), 0, pmin(., 1))) %>%
  rowMeans(na.rm = TRUE)
hist(info_mean)

merged2 <- filter(merged, info_mean > 0.5) %>%
  select(-starts_with("info_")) %>%
  as_tibble() %>%
  print()
hist(merged2$maf)
table(merged2$chr)

maf_ext <- merged2 %>%
  select(starts_with("freq_")) %>%
  mutate_all(~ pmin(., 1 - .)) %>%
  as.matrix() %>%
  matrixStats::rowMedians(na.rm = TRUE)

ind <- sample(length(maf_ext), 100e3)
library(ggplot2)
theme_set(bigstatsr::theme_bigstatsr())
qplot(merged2$maf[ind], maf_ext[ind], alpha = I(0.2)) +
  labs(x = "MAF in UKBB", y = "Median of external MAFs") +
  geom_abline(color = "red") +
  coord_equal()
# ggsave("figures/compare-MAF.png", width = 6, height = 6)

diff <- abs(merged2$maf - maf_ext)
hist(diff, "FD", xlim = c(0, 0.03))
qplot(y = diff[1:100000]) + labs(x = "Index", y = "Difference between MAFs")
# ggsave("figures/diff-MAF.png", width = 10, height = 6)
diff2 <- bigutilsr::rollmean(diff, 10)
plot(diff2[1:10000], pch = 20, cex = 0.8)
hist(diff2, "FD", xlim = c(0, 0.03))

dups <- vctrs::vec_duplicate_detect(merged2[, c("chr", "pos")])
merged3 <- filter(merged2, diff2 < 0.01, !dups) %>%
  select(-starts_with("freq_"), -maf) %>%
  print()
# saveRDS(merged3, "data/map_UKBB_QC.rds")

list_snp_id <- with(merged3, split(paste(chr, pos, a1, a0, sep = "_"), chr))
lengths(list_snp_id)
#      1      2      3      4      5      6      7      8      9     10     11
# 643770 700731 595289 616012 515429 550013 494053 467463 356363 417506 416805
#     12     13     14     15     16     17     18     19     20     21     22
# 385342 313598 264182 235490 236840 217303 240047 174888 178270 113558 105740


#### Prepare bigSNP objects ####

(NCORES <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")) - 1L)

eid <- bigsnpr::snp_attach("data/UKBB_HM3.rds")$fam$eid
sample <- bigreadr::fread2("UKBB/ukb58024_imp_chr1_v3_s487296.sample")[-1, ]

system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles   = glue::glue("UKBB/bgen/ukb_imp_chr{chr}_v3.bgen",
                             chr = names(list_snp_id)),
    list_snp_id = list_snp_id,
    backingfile = "data/UKBB_large",
    ind_row     = match(eid, sample$ID_2),
    ncores      = NCORES
  )
) # 3.2 H with 15 cores
