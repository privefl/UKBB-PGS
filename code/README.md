## Structure of the code

### Prepare data

1. `prepare-data`: Prepare covariates + ancestry groups + HapMap3 data

2. `prepare-ashkenazi-jews`: Prepare the ninth ancestry group + figures 1 and S20

2. `prepare-fields` + `prepare-pheno-fields`: Prepare phenotypes from fields of csv

4. `prepare-phecodes`: Prepare phecodes from ICD data

5. `prepare-genotyped-data` + `prepare-larger-set`: prepare additional data with different sets of variants + figures S18 and S19


### Run methods on real data

- `run-PLR-continuous` + `run-PLR-binary`: Run penalized linear/logistic regression (PLR) for all phenotypes using HapMap3 variants

- `run-GWAS` + `run-ldpred2`: Run LDpred2-auto for all phenotypes using HapMap3 variants

- `run-PLR-continuous-geno`: Run PLR using high-quality genotyped variants only

- `run-GWAS-large` + `run-ldpred2-large`: Run LDpred2-auto for a few phenotypes prioritizing 1M variants out of 8M+ common variants

- `run-PLR-continuous-multi` + `run-PLR-binary-multi`: Run PLR for a few phenotypes using a mixture of ancestries in training


### Visualize results

- `assess-pred-AJ` + `assess-pred`: Figures 2 and 3

- `assess-pred-geno`: Figure S2

- `assess-pred-ldpred2`: Figures S3 + S12 to S15

- `recompute-PCA`: Figures S4 and S5

- `plot-zoom-gwas`: Figures 4 and S6 to S10

- `assess-pred-ldpred2-large`: Figure 5 and table S1

- `assess-pred-multi`: Figure S11

- `timings`: Figures S16 and S17

- `QC-plot-large-effect`: Figure S21


### Export results

- `export-effects`
