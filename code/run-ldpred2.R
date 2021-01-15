library(bigsnpr)
library(dplyr)

# build SFBM (same for all pheno)
corr <- runonce::save_run({

  map_ldref <- readRDS("../paper-ldpred2/ld-ref/map.rds")
  map_ukbb <- snp_attach("data/UKBB_HM3.rds")$map %>%
    transmute(chr = as.integer(chromosome), pos = physical.pos,
              a0 = allele1, a1 = allele2, beta = 1)
  info_snp <- snp_match(map_ukbb, map_ldref)

  for (chr in 1:22) {

    cat(chr, ".. ", sep = "")

    ## indices in 'info_snp'
    ind.chr <- which(info_snp$chr == chr)
    ## indices in 'map_ldref'
    ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
    ## indices in 'corr_chr'
    ind.chr3 <- match(ind.chr2, which(map_ldref$chr == chr))

    ld_file <- paste0("../paper-ldpred2/ld-ref/LD_chr", chr, ".rds")
    corr_chr <- readRDS(ld_file)[ind.chr3, ind.chr3]

    if (chr == 1) {
      corr <- as_SFBM(corr_chr, "tmp-data/LDref")
    } else {
      corr$add_columns(corr_chr, nrow(corr))
    }
  }
  corr
}, file = "tmp-data/LDref.rds")
dim(corr)  # 1040096 x 1040096
rm(corr)

pheno_files <- c(list.files("data/ukbb-quant-pheno", full.names = TRUE),
                 list.files("data/ukbb-binary-pheno", full.names = TRUE))

files <- tibble(
  pheno_file = pheno_files,
  gwas_file = file.path("GWAS", basename(pheno_file)),
  res_file = file.path("ldpred2", basename(pheno_file))
) %>%
  filter(file.exists(gwas_file), !file.exists(res_file)) %>%
  # filter(!startsWith(basename(pheno_file), "LTFH_")) %>%
  # slice_sample(n = 80) %>%
  print(n = Inf)

bigassertr::assert_dir("ldpred2")

NCORES <- 30
library(future.batchtools)
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 2, mem = "250g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

furrr::future_pmap(files, function(pheno_file, gwas_file, res_file) {

  node <- Sys.getenv("SLURMD_NODENAME")

  # Get effective sample size from GWAS
  fam <- snp_attach("data/UKBB_HM3.rds")$fam
  ind_csv <- match(fam$eid, readRDS("data/csv_eid.rds"))
  y <- readRDS(print(pheno_file))[ind_csv]
  ind.train <- which(!is.na(y) & (fam$set == "train"))
  N <- length(ind.train)

  ## Information for the variants provided in the LD reference
  map_ldref <- readRDS("../paper-ldpred2/ld-ref/map.rds")
  map_ukbb <- snp_attach("data/UKBB_HM3.rds")$map %>%
    transmute(chr = as.integer(chromosome), pos = physical.pos,
              a0 = allele1, a1 = allele2)

  sumstats <- readRDS(gwas_file) %>%
    transmute(beta = estim, beta_se = std.err, n_eff = N) %>%
    bind_cols(map_ukbb)

  df_beta <- snp_match(sumstats, map_ldref) %>%
    as_tibble() %>%
    select(-starts_with("pos_")) %>%
    print()

  # Heritability estimation of LD score regression
  # to be used as a starting value in LDpred2-auto
  (ldsc <- with(df_beta, snp_ldsc(ld, ld_size = nrow(map_ldref),
                                  chi2 = (beta / beta_se)^2,
                                  sample_size = n_eff,
                                  ncores = NCORES)))
  h2_est <- ldsc[["h2"]]

  corr <- readRDS("tmp-data/LDref.rds")

  # LDpred2-auto
  time <- system.time(
    multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                                   # burn_in = 200, num_iter = 100,
                                   vec_p_init = seq_log(1e-4, 0.9, 30),
                                   sparse = TRUE, ncores = NCORES)
  )[[3]]

  saveRDS(list(ldsc = ldsc, ldpred2 = multi_auto, time = time, node = node),
          res_file)
})
