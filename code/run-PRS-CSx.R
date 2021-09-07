library(bigsnpr)
library(dplyr)
ukbb <- snp_attach("data/UKBB_HM3.rds")

tmp <- "tmp-data/bim_for_prscsx"
ukbb$map %>%
  transmute(chromosome, marker.ID = rsid, genetic.dist = 0,
            physical.pos, allele1, allele2) %>%
  bigsnpr:::write.table2(paste0(tmp, ".bim"))

info_snp <- ukbb$map %>%
  transmute(SNP = rsid, A1 = allele1, A2 = allele2)

POP <- c("afr", "eur", "sas", "eas")

library(dplyr)
files <- tidyr::expand_grid(
  pheno = c("174.1", "185", "250.2", "401", "411.4",
            "darker_skin", "log_bilirubin", "log_lipoA", "years_of_edu"),
  chr = 1:22
) %>%
  filter(!file.exists(paste0(
    "results-prscsx/", pheno, "_afr_pst_eff_a1_b0.5_phiauto_chr", chr, ".txt"))) %>%
  print(n = Inf)

bigassertr::assert_dir("results-prscsx")

NCORES <- 15
library(future.batchtools)
plan(batchtools_slurm(resources = list(
  t = "12:00:00", c = NCORES + 1, mem = "100g",
  name = basename(rstudioapi::getSourceEditorContext()$path))), workers = Inf)

furrr::future_pmap(files, function(pheno, chr) {

  library(bigstatsr)

  N <- list()
  sumstats <- list()
  for (pop in POP) {
    gwas <- readRDS(paste0("GWAS-multi/", pop, "_", pheno, ".rds"))
    N[pop] <- as.integer(sub(".+, df = ([0-9]+),.+", "\\1", body(attr(gwas, "predict"))[2]))
    sumstats[pop] <- info_snp %>%
      bind_cols(BETA = gwas$estim, P = predict(gwas, log10 = FALSE)) %>%
      na.omit() %>%
      bigreadr::fwrite2(file = tempfile(tmpdir = "tmp-data"), sep = "\t")
  }

  on.exit(file.remove(unlist(sumstats)), add = TRUE)

  prscsx <- "PRScsx/PRScsx.py"
  tmp <- tmp

  system(glue::glue(
    "OMP_NUM_THREADS=", NCORES[[1]],
    " python3 {prscsx}",
    " --ref_dir=PRScsx/ldref",
    " --bim_prefix={tmp}",
    " --pop={paste(POP, collapse = ',')}",
    " --sst_file={paste(sumstats[POP], collapse = ',')}",
    " --n_gwas={paste(N[POP], collapse = ',')}",
    " --chrom={chr}",
    " --out_dir=results-prscsx",
    " --out_name={pheno}"
  ))
})
