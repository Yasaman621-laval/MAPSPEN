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


📚 Citation
If you use this package in your research, please cite:

Tahernezhad, Y. (2025). Development of Novel Statistical Methods for Genetic Architecture Estimation and Polygenic Risk Score Construction Across Traits, Populations, and Sexes. PhD Thesis, Université Laval.


---

Let me know if you’d like a shorter version for CRAN-style, or want me to generate the matching `DESCRIPTION` file or vignette outline.
