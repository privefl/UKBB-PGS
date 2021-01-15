library(bigreadr)
library(dplyr)
csv <- "UKBB/ukb41181.csv"
sex <- fread2(csv, select = "22001-0.0")[[1]]

df_ICD10 <- fread2(csv, colClasses = "character", select = c(
  paste0("40001-", 0:1, ".0"),     # death
  paste0("40002-0.", 0:13),        # death
  paste0("40002-1.", 0:13),        # death
  paste0("40006-", 0:16, ".0"),    # cancer
  paste0("41270-0.", 0:212),       # diagnosis
  paste0("41201-0.", 0:21)         # external cause
))

# Non-cancer illness code, self-reported
# coding609 <- fread2("UKBB/coding609.tsv")
# df_self_reported <- fread2(csv, select = c(paste0("20002-0.", 0:33),
#                                            paste0("20002-1.", 0:33),
#                                            paste0("20002-2.", 0:33),
#                                            paste0("20002-3.", 0:33))) %>%
#   mutate_all(~ as.character(factor(., levels = coding609$coding,
#                                    labels = coding609$meaning)))

# Depression ever diagnosed by a professional
df_dep <- fread2(csv, select = c(paste0("20544-0.", 1:16))) %>%
  mutate_all(~ ifelse(. == 11, "F330", NA)) %>%
  select_if(is.character)

df_ICD9 <- fread2(csv, colClasses = "character", select = c(
  paste0("40013-", 0:14, ".0"),    # cancer
  paste0("41271-0.", 0:46)         # diagnosis
))

id_icd10_count <- bind_cols(df_ICD10, df_dep) %>%
  mutate_all(~ ifelse(. == "", NA, .)) %>%
  mutate(id = row_number()) %>%
  tidyr::pivot_longer(-id, values_to = "code", values_drop_na = TRUE) %>%
  group_by(id, code) %>%
  summarise(count = n(), .groups = "drop") %>%
  transmute(id, vocabulary_id = "ICD10", code, count)
str(id_icd10_count)

id_icd9_count <- df_ICD9 %>%
  mutate_all(~ ifelse(. == "", NA, .)) %>%
  mutate(id = row_number()) %>%
  tidyr::pivot_longer(-id, values_to = "code", values_drop_na = TRUE) %>%
  group_by(id, code) %>%
  summarise(count = n(), .groups = "drop") %>%
  transmute(id, vocabulary_id = "ICD9", code, count)
str(id_icd9_count)

# https://phewascatalog.org/phecodes
map_icd9 <- fread2("tmp-data/phecode_icd9_rolled.csv", colClasses = "character")
phecode_map_icd9 <- transmute(map_icd9, vocabulary_id = "ICD9",
                              code = ICD9, phecode = PheCode)


# https://github.com/PheWAS/PheWAS
library(PheWAS)
phecodes <- createPhenotypes(
  rbind(id_icd10_count, id_icd9_count),
  id.sex = data.frame(id = seq_along(sex), sex = c("F", "M")[sex + 1L]),
  vocabulary.map = mutate_at(rbind(phecode_map_icd10, phecode_map_icd9),
                             "code", ~ sub("\\.", "", .)),
  min.code.count = 1,
  add.phecode.exclusions = TRUE,
  full.population.ids = seq_along(sex)
)
phecodes <- phecodes[order(phecodes$id), -1]
# saveRDS(phecodes[colSums(phecodes, na.rm = TRUE) > 50],
#         "tmp-data/all_phecodes.rds")

table(pheno = phecodes$`185`, sex, exclude = NULL)
#        sex
# pheno        0      1   <NA>
#   FALSE      0 195543  13476
#   TRUE       0  10905    285
#   <NA>  264796  17019    481
table(pheno = phecodes$`174.1`, sex, exclude = NULL)
table(pheno = phecodes$`654.2`, sex, exclude = NULL)
table(pheno = phecodes$`296.2`, sex, exclude = NULL)


codes_kept <- c(
  "153", "165", "172", "174.1", "180", "185", "187.2", "189.2",
  "191.11", "193", "200.1", "208", "211", "218", "241.2", "242",
  "244", "250.1", "250.2", "250.7", "251.1", "252.1", "272", "274.1",
  "275.1", "277.4", "278", "280", "285", "286.12", "290.1", "290.11",
  "296.2", "318", "332", "335", "351", "361", "362.29", "364.5",
  "365", "366", "371", "371.1", "379.3", "383", "401", "411.4",
  "415", "427.2", "428", "433", "433.1", "442.11", "443.9", "451",
  "454", "455", "459.9", "471", "495", "496", "530.1", "531.3",
  "535.6", "540", "550.1", "555.1", "555.2", "557.1", "562.1",
  "564", "565.1", "571.5", "574", "575", "593", "594", "596.1",
  "600", "615", "618", "626.1", "654.2", "681", "695.4", "696.4",
  "697", "702.1", "702.2", "704", "706.2", "709.2", "709.7", "714.1",
  "715.2", "716", "717", "727.4", "728.71", "735.3", "740", "743.1",
  "785", "790.6", "960.2"
)
length(codes_kept) # 106
phe2 <- phecodes[codes_kept]
sort(sapply(phe2, function(.) sum(., na.rm = TRUE)))
# 347 min and 113605 max

bigassertr::assert_dir("data/ukbb-binary-pheno")
purrr::iwalk(phe2, ~ saveRDS(.x, paste0("data/ukbb-binary-pheno/", .y, ".rds")))
