
# file.symlink("~/NCRR-PRS/faststorage/UKBB/", ".")

#### Prepare covariates data + ancestry groups ####

library(bigreadr)
library(dplyr)

code_country <- filter(fread2("UKBB/coding89.tsv"), selectable == "Y")

df0 <- fread2(
  "UKBB/ukb41181.csv",
  select = c("eid", "34-0.0", "52-0.0", "22001-0.0", "21022-0.0", "189-0.0",
             "21000-0.0", "20115-0.0", paste0("22009-0.", 1:16)),
  col.names = c("eid", "year", "month", "sex", "age", "deprivation_index",
                "ancestry", "country", paste0("PC", 1:16))
) %>%
  mutate(
    date = (year - 1900) + (month - 0.5) / 12,
    year = NULL, month = NULL,
    country = factor(country, levels = code_country$coding,
                     labels = code_country$meaning)
  )
str(df0)
# saveRDS(df0$eid, "data/csv_eid.rds")
df0$country[with(df0, is.na(country) & ancestry == 1001)] <- "United Kingdom"
df0$country[with(df0, is.na(country) & ancestry == 1002)] <- "Ireland"

# Ancestry grouping (https://doi.org/10.1101/2020.10.06.328203)
PC_UKBB <- as.matrix(select(df0, PC1:PC16))
POP <- c("United Kingdom", "Poland", "Iran", "Italy",
         "India", "China", "Caribbean", "Nigeria")
country2 <- factor(df0$country, levels = POP)
df0$country <- df0$ancestry <- NULL
nona <- complete.cases(PC_UKBB) & !is.na(country2)
all_centers <- bigutilsr::geometric_median(PC_UKBB[nona, ], by_grp = country2[nona])
# saveRDS(all_centers, "data/eight_centers.rds")
all_sq_dist <- apply(all_centers, 1, function(center) {
  bigutilsr:::rowSumsSq(sweep(PC_UKBB, 2, center, '-'))
})
choose_pop <- apply(all_sq_dist, 1, function(x) {
  ind <- which.min(x)
  if (length(ind) == 0) NA else ind
})
min_sq_dist <- all_sq_dist[cbind(seq_along(choose_pop), choose_pop)]
(thr_sq_dist <- 0.002 * (max(dist(all_centers)^2) / 0.16))
hist(log(min_sq_dist)); abline(v = log(thr_sq_dist), col = "red")
df0$group <- ifelse(min_sq_dist > thr_sq_dist, NA, POP[choose_pop])
table(df0$group)
#  Caribbean          China          India           Iran
#       2655           1853           6720           1234
#  Italy        Nigeria         Poland United Kingdom
#   6980           4086           4311         446682
df0[!complete.cases(df0), ] <- NA

# Match indices in BGEN data
sample <- bigreadr::fread2("UKBB/ukb58024_imp_chr1_v3_s487296.sample")[-1, ]
ind.indiv <- match(df0$eid, sample$ID_2)
df0[is.na(ind.indiv), ] <- NA

# Filter for close relatedness
rel <- bigreadr::fread2("UKBB/ukb58024_rel_s488264.dat")
thr_rel <- 2^-3.5
hist(rel$Kinship); abline(v = thr_rel, col = "red")
rel2 <- dplyr::filter(rel, ID1 %in% df0$eid, Kinship > thr_rel)
df0[df0$eid %in% rel2$ID2, ] <- NA

length(sub <- which(!is.na(df0$eid)))  # 437,669 individuals


#### Prepare variants to use ####

library(bigsnpr)

# Provided in https://doi.org/10.1101/2020.04.28.066720
map_ukbb_hm3 <- readRDS(runonce::download_file(
  "https://ndownloader.figshare.com/files/24928055",
  dir = "tmp-data", fname = "map.rds"))
map_ipsych <- snp_attach("../data_ipsych/dosage_ipsych2015.rds")$map %>%
  filter(info_2012 > 0.6, info_2015 > 0.6) %>%
  transmute(chr = CHR, pos = POS, a1, a0 = a2)

info_match <- snp_match(mutate(map_ukbb_hm3, beta = 1), map_ipsych)
# 1,054,330 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 1,040,096 variants have been matched; 0 were flipped and 1,040,096 were reversed.

list_snp_id <- with(info_match, split(paste(chr, pos, a1, a0, sep = "_"), chr))


#### Prepare bigSNP objects ####

(NCORES <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")) - 1L)

system.time(
  rds <- bigsnpr::snp_readBGEN(
    bgenfiles   = glue::glue("UKBB/bgen/ukb_imp_chr{chr}_v3.bgen", chr = 1:22),
    list_snp_id = list_snp_id,
    backingfile = "data/UKBB_HM3",
    ind_row     = ind.indiv[sub],
    ncores      = NCORES
  )
) # 47 min with 15 cores

ukbb <- snp_attach(rds)
ukbb$fam <- df0[sub, ]
set <- ifelse(ukbb$fam$group == "United Kingdom", "train", "test")
set.seed(1); set[sample(which(set == "train"), 20e3)] <- "test"
ukbb$fam$set <- set

snp_save(ukbb)
