# GWAS Pipeline Documentation

## Steps Overview
| Step         | Command                              | Description                  |
|--------------|--------------------------------------|------------------------------|
| Split IID by chip   | splitbychip.ipynb (Python)    | extract IIDs specific to one chip |
|    |                      |          |

# 1. Initial per-chip QC (split) [x]
- why: different chips = different SNPs, genotyping errors etc.

Python script to:
1. Read the metadata and split IIDs by chip.
2. Create `.keep` files for each chip group.

Use `.keep` to filter bed file into chip specific files:
```bash
plink --bfile gwas_data --keep /faststorage/project/populationgenomics/students/amlg/project/data/splitbychip/chipX.keep --make-bed --out /faststorage/project/populationgenomics/students/amlg/project/data/splitbychip/chipX_data
```

## 1.1 Sample-level QC []

### 1.1.1 Missingness []

### 1.1.2 Sex checks []

### 1.1.3 Heterozygosity outliers []

### 1.1.4 Relatedness / Duplicates []

## 1.2 SNP-level QC []

### 1.2.1 Call rate/Missingness []

### 1.2.2 MAF []

### 1.2.3 HWE []


---
## 1. Impute Sex
```bash
plink --bfile gwas_data_indiv_filtered
      --impute-sex 0.8 
      --make-bed --out gwas_data_sex_imputed
      --allow-no-sex
```
| command                 | what it does                                                                      |
|-------------------------|-----------------------------------------------------------------------------------|
| ```--bfile```           | reads in all .bed, .bam & .fam files                                              |
| ```--impute-sex 0.8```  | command to impute sex, F>0.8 likely male, F<0.2 likely female, 0.8 standard value |
| ```--make-bed```        | needed to produce new output files                                                |
| ```--allow-no-sex```    | allows PLINK to process those individuals at all                                  |

- impute sex = infer sex: change sex assignments to imputed values
- get sex of individuals based on X chromosome heterozygosity rates
  - SNPs determine whether an individual is likely male (low heterozygosity on X, since XY) or female (normal heterozygosity on X, since XX)
- outputs a .sexcheck file

```bash
# The terminal output I get
515439 MB RAM detected; reserving 257719 MB for main workspace.
1376653 variants loaded from .bim file.
2009 people (992 males, 916 females, 101 ambiguous) loaded from .fam.
Ambiguous sex IDs written to
../../students/amlg/project/data/gwas_data_sex_imputed.nosex .
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 2009 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Warning: 106246 het. haploid genotypes present (see
../../students/amlg/project/data/gwas_data_sex_imputed.hh ); many commands
treat these as missing.
Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
treat these as missing.
Total genotyping rate is 0.478922.
1376653 variants and 2009 people pass filters and QC.
Note: No phenotypes present.
--impute-sex: 29856 Xchr and 0 Ychr variant(s) scanned, 1889/2009 sexes
imputed. Report written to
../../students/amlg/project/data/gwas_data_sex_imputed.sexcheck .
--make-bed to ../../students/amlg/project/data/gwas_data_sex_imputed.bed +
../../students/amlg/project/data/gwas_data_sex_imputed.bim +
../../students/amlg/project/data/gwas_data_sex_imputed.fam ... done.
```

## 2. Check Sex
---

## How do do .md files

### code blocks
Use triple backticks (```) for code blocks:

```bash
plink --bfile data --mind 0.1 --geno 0.05 --maf 0.01 --make-bed --out qc_filtered
```
You can specify language after the triple backticks (e.g., bash, r, python) for syntax highlighting.

Also possible like this
```bash
plink --bfile data/raw/genotypes \
  --mind 0.1 \
  --geno 0.05 \
  --maf 0.01 \
  --make-bed \
  --out data/clean/genotypes_qc
```

### bullet lists
with - or 1.

### tables
| Step         | Command                              | Description                  |
|--------------|--------------------------------------|------------------------------|
| Initial QC   | plink --mind 0.1                     | Remove individuals with >10% missing data |
| MAF Filter   | plink --maf 0.01                     | Remove rare variants         |

### links
[link name](https://github.com/antoniamlg/GT_GWASproject)

### images
![image name](image.png)

### text formatting
- use (*) or (_) for *italic* or _italic_ text
- use (**) or (__) for __bold__ text