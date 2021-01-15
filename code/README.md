## Structure of the code

### Prepare data

1. `prepare-data`: Prepare covariates + ancestry groups + HapMap3 data

2. `prepare-fields` + `prepare-pheno-fields`: Prepare phenotypes from fields of csv

4. `prepare-phecodes`: Prepare phecodes from ICD data

5. `prepare-genotyped-data` + `prepare-larger-set`: prepare additional data with different sets of variants + figures S19 and S20


### Run methods on real data

- `run-PLR-continuous` + `run-PLR-binary`: Run penalized linear/logistic regression (PLR) for all phenotypes using HapMap3 variants

- `run-GWAS` + `run-ldpred2`: Run LDpred2-auto for all phenotypes using HapMap3 variants

- `run-PLR-continuous-geno`: Run PLR using high-quality genotyped variants only

- `run-GWAS-large` + `run-ldpred2-large`: Run LDpred2-auto for a few phenotypes prioritizing 1M variants out of 8M+ common variants

- `run-PLR-continuous-multi` + `run-PLR-binary-multi`: Run PLR for a few phenotypes using a mixture of ancestries in training


### Visualize results

- `assess-pred`: Figures 1 and 2

- `assess-pred-geno`: Figure S1

- `assess-pred-ldpred2`: Figures S2 + S13 to S16

- `recompute-PCA`: Figures S3 and S4

- `plot-zoom-gwas`: Figures S5 to S10

- `assess-pred-ldpred2-large`: Figure S11

- `assess-pred-multi`: Figure S12

- `timings`: Figures S17 and S18

- `QC-plot-large-effect`: Figure S21
