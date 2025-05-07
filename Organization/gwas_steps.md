# GWAS Pipeline Documentation

## Steps Overview
| Step         | Command                              | Description                  |
|--------------|--------------------------------------|------------------------------|
| Impute sex   | xx                     | Remove individuals with >10% missing data |
| MAF Filter   | plink --maf 0.01                     | Remove rare variants         |

## 1. Impute Sex
```bash
plink --bfile gwas_data_indiv_filtered --impute-sex 0.8 --make-bed --out gwas_data_sex_imputed --allow-no-sex
```
| command   | what it does            |
| ```--bfile```   | reads in all .bed, .bam & .fam files  |
| ```--impute-sex 0.8```  | command to impute sex, F>0.8 likely male, F<0.2 likely female, 0.8 standard value |
| ```--make-bed```  | needed to produce new output files  |

- impute sex = infer sex: change sex assignments to imputed values
- get sex of individuals based on X chromosome heterozygosity rates
  - SNPs determine whether an individual is likely male (low heterozygosity on X, since XY) or female (normal heterozygosity on X, since XX)
- outputs a .sexcheck file

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