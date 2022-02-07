## Portability of 245 polygenic scores when derived from the UK Biobank and applied to 9 ancestry groups from the same cohort

Preprint: https://doi.org/10.1101/2021.02.05.21251061

AJHG paper: https://doi.org/10.1016/j.ajhg.2021.11.008

### Code to reproduce ancestry groups

```r
PC_UKBB <- bigreadr::fread2(
  "UKBB/ukb41181.csv",  ## REPLACE WITH YOURS
  select = paste0("22009-0.", 1:16))
all_centers <- read.csv(
  "https://raw.githubusercontent.com/privefl/UKBB-PGS/main/pop_centers.csv",
  stringsAsFactors = FALSE)
all_sq_dist <- apply(all_centers[-1], 1, function(one_center) {
  rowSums(sweep(PC_UKBB, 2, one_center, '-')^2)
})
thr_sq_dist <- max(dist(all_centers[-1])^2) * 0.002 / 0.16
group <- apply(all_sq_dist, 1, function(x) {
  grp <- NA
  ind <- which.min(x)
  if (isTRUE(x[ind] < thr_sq_dist)) {
    grp <- all_centers$Ancestry[ind]
    # We used a more stringent cutoff for the Ashkenazi group
    if (grp == "Ashkenazi" && x[ind] > 12.5^2) grp <- NA
  }
  grp
})
table(group, exclude = NULL)
# Ashkenazi      Caribbean          China          India           Iran
#      2500           2655           1853           6720           1234
# Italy        Nigeria         Poland United Kingdom           <NA>
#  6824           4086           4311         446682          25640
```

### Polygenic score effect sizes for 200+ phenotypes

- from the PGS Catalog: https://www.pgscatalog.org/publication/PGP000263/

- from figshare: https://doi.org/10.6084/m9.figshare.14074760.v2 (one large matrix of variants x phenotypes)
