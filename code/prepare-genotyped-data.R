#### Prepare variants to use ####

snp_qc <- bigreadr::fread2(runonce::download_file(
  "https://biobank.ctsu.ox.ac.uk/crystal/crystal/auxdata/ukb_snp_qc.txt",
  dir = "tmp-data"
))

af_qc <- bigreadr::fread2(runonce::download_file(
  "http://kunertgraf.com/data/files/snp_maf_comparison.csv",
  dir = "tmp-data"
))


library(dplyr)
snp_qc$qc_all <- rowMeans(select(snp_qc, ends_with("_qc")))
snp_qc2 <- snp_qc %>%
  select(-ends_with("_qc"), -ends_with("_loading")) %>%
  left_join(af_qc, by = c("rs_id" = "rsid"))

plot(as.factor(snp_qc2$oob), snp_qc2$qc_all)

snp_qc_final <- snp_qc %>%
  select(-2, -3, -ends_with("_qc"), -ends_with("_loading")) %>%
  as_tibble() %>%
  filter(array == 2, qc_all == 1) %>%
  anti_join(filter(af_qc, oob == 1), by = c("rs_id" = "rsid")) %>%
  select(2:5)

library(doParallel)
registerDoParallel(cl <- makeCluster(12))
system.time({
  list_snp_id <- foreach(chr = 1:22, .packages = "dplyr") %dopar% {
    paste0("UKBB/mfi/ukb_mfi_chr", chr, "_v3.txt") %>%
      bigreadr::fread2(showProgress = FALSE) %>%
      filter(V6 > 0.01, V8 == 1) %>%   ## MAF > 1% & INFO = 1
      dplyr::semi_join(
        filter(snp_qc_final, chromosome == !!chr),
        by = c("V3" = "position", "V4" = "allele1_ref", "V5" = "allele2_alt")
      ) %>%
      with(paste(chr, V3, V4, V5, sep = "_"))
  }
}) # 80 sec
stopCluster(cl)

sum(lengths(list_snp_id))  # 586,534
lengths(list_snp_id)
#  [1] 46526 47020 39640 37351 35047 40949 32385 30552 25985 29403 29000
# [12] 27976 20264 19012 18572 20691 19026 17675 15755 15404  8836  9465

sample <- bigreadr::fread2("UKBB/ukb58024_imp_chr1_v3_s487296.sample")[-1, ]
fam <- bigsnpr::snp_attach("data/UKBB_HM3.rds")$fam

system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles = glue::glue("UKBB/bgen/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKBB_geno",
    ind_row = match(fam$eid, sample$ID_2),
    ncores = 15
  )
) # 32 min
