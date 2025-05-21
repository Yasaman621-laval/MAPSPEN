# MAPSPEN

**MAPSPEN** (Multi-Ancestry Posterior Shrinkage Penalization) is a novel statistical method for constructing polygenic risk scores (PRS) using GWAS summary statistics across multiple populations. This R package accounts for heterogeneity in genetic architecture by applying adaptive penalization informed by linkage disequilibrium (LD) patterns and cross-population effect distributions.

> Developed and maintained by [Yasaman Tahernezhad](https://github.com/Yasaman621-laval) as part of her PhD thesis.

---

## ✨ Key Features

- Estimates population-specific PRS using only summary statistics
- Adjusts for LD structure via penalized likelihood
- Supports multi-population modeling (EUR, EAS, AFR, etc.)
- Incorporates genetic correlation and effect size variance across ancestries
- Profile likelihood-based hyperparameter tuning
- Benchmarked against PRS-CSX using real and simulated genotype data

---

## 📦 Installation

To install the package from source:

```r
install.packages("MAPSPEN.tar.gz", repos = NULL, type = "source")

 Or clone this repository and install from local:

git clone https://github.com/Yasaman621-laval/MAPSPEN.git
setwd("MAPSPEN")
install.packages(".", repos = NULL, type = "source")

🚀 Usage Example
library(MAPSPEN)

# Load input summary statistics and LD information
# summaryStats: matrix of Z-scores across populations
# LDmat: list of population-specific LD matrices
# Output: posterior mean of SNP effect sizes and predicted PRS
```
---
## 📁 Data Format and Preprocessing

MAPSPEN uses **GWAS summary statistics** and **LD reference panels** as inputs. The expected formats and preprocessing steps are outlined below.

---

### ✅ Required Summary Statistics Format

Each GWAS file (one per population) should contain a single tab-delimited table with the following columns:

| Column | Description                      |
|--------|----------------------------------|
| SNP    | SNP identifier (e.g., rsID)      |
| N      | Sample size                      |
| Z      | Z-score for SNP-trait association |
| A1     | Effect (risk) allele             |
| A2     | Non-effect allele                |

This format matches the standard LDSC and MiXeR input formats. You may optionally preprocess raw GWAS files using `munge_sumstats.py` from [LDSC](https://github.com/bulik/ldsc).

> ⚠️ For case/control studies:  
> Use `neff = 4 / (1/ncase + 1/ncontrol)` for the N column to adjust for imbalanced sample sizes.

---

### 🌐 Public Summary Statistics (Examples)

You can access GWAS summary statistics from these public repositories:

- **Schizophrenia (PGC 2014):**  
  https://www.med.unc.edu/pgc/download-results/  
  Look for: _"Download results: Schizophrenia (49 European-ancestry cohorts)"_

- **Educational Attainment (Lee et al., 2018):**  
  https://www.thessgac.org/data  
  File: `GWAS_EA_excl23andMe.txt`

- **UK Biobank GWAS (Neale Lab):**  
  http://www.nealelab.is/uk-biobank  
  Broad coverage of phenotypes with summary stats by ancestry.

- **IEU OpenGWAS Project:**  
  https://gwas.mrcieu.ac.uk/  
  Search and download thousands of harmonized GWAS summary statistics.

---

### ✅ Reference LD Panel

MAPSPEN requires LD reference data (e.g., from the 1000 Genomes Project Phase 3). For European-based studies, you can download the HapMap3-constrained PLINK files:

```bash
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
tar -xzvf 1000G_Phase3_plinkfiles.tgz
bzip2 -d w_hm3.snplist.bz2
---.
##📚 Citation
If you use this package in your research, please cite:

Tahernezhad, Y. (2025). Development of Novel Statistical Methods for Genetic Architecture Estimation and Polygenic Risk Score Construction Across Traits, Populations, and Sexes. PhD Thesis, Université Laval.

---

Let me know if you’d like a shorter version for CRAN-style, or want me to generate the matching `DESCRIPTION` file or vignette outline.
