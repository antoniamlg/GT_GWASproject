# GWAS Pipeline Documentation

- [GWAS Pipeline Documentation](#gwas-pipeline-documentation)
  - [Steps Overview](#steps-overview)
  - [Datanumbers](#datanumbers)
- [0. Filter by phenotype \[X\]](#0-filter-by-phenotype-x)
- [1. Initial per-chip QC (split) \[X\]](#1-initial-per-chip-qc-split-x)
      - [As we do not have that many individuals (now \`\`)](#as-we-do-not-have-that-many-individuals-now-)
  - [1.1 SNP-level QC \[ \]](#11-snp-level-qc--)
    - [1.1.1 Call rate/Missingness \[ \]](#111-call-ratemissingness--)
      - [Results](#results)
    - [1.1.2 Filter SNPs that are missing in 100% of individuals \[ \]](#112-filter-snps-that-are-missing-in-100-of-individuals--)
    - [1.1.3 Merge chips](#113-merge-chips)
    - [1.1.X Replot missingness to see if i actually did something.](#11x-replot-missingness-to-see-if-i-actually-did-something)
    - [1.1.4 Sex-imputation (TRY)](#114-sex-imputation-try)
    - [1.1.5 MAF \& HWE \[ \]](#115-maf--hwe--)
  - [1.2 Sample-level QC \[ \]](#12-sample-level-qc--)
    - [1.2.1 Sex checks \[ \]](#121-sex-checks--)
    - [1.2.2 Missingness \[ \]](#122-missingness--)
    - [1.2.3 Heterozygosity outliers \[ \]](#123-heterozygosity-outliers--)
    - [1.2.4 Relatedness / Duplicates \[ \]](#124-relatedness--duplicates--)
  - [1.3 Population Structure/ Stratification/ Batch Effects](#13-population-structure-stratification-batch-effects)
    - [1.3.1 PCA](#131-pca)
          - [After QC, merge chips again](#after-qc-merge-chips-again)
  - [1. Impute Sex](#1-impute-sex)
  - [2. Check Sex](#2-check-sex)
  - [How do do .md files](#how-do-do-md-files)
    - [code blocks](#code-blocks)
    - [bullet lists](#bullet-lists)
    - [tables](#tables)
    - [links](#links)
    - [images](#images)
    - [text formatting](#text-formatting)

## Steps Overview
| Step         | Command                              | Description                  |
|--------------|--------------------------------------|------------------------------|
| Split IID by chip   | splitbychip.ipynb (Python)    | extract IIDs specific to one chip |
|    |                      |          |

Also check out this link [QC README](https://github.com/kaspermunch/PopulationGenomicsCourse/blob/master/Exercises/GWAS_QC/step_by_step_tutorial.md).

## Datanumbers
| Chip Name             | # Individuals | # Variants          | #Ind. after phenotype filtering  |
|-----------------------|---------------|---------------------|----------------------------------|
| whole dataset         | 2009          | 1376653             | 1071                             |
| HTS_iSelect_HD        | 587           |                     | 275                              |
| Illumina_GSAs         | 248           |                     | 128                              |
| OmniExpress           | 291           |                     | 170                              |
| OmniExpress_plus      | 484           |                     | 266                              |
| unknown_chip          | 399           |                     | 232                              |

# 0. Filter by phenotype [X]
Our height.txt file only contains the phenotype of `1106` individuals. We cannot use the remaining ones without a phenotype for anything so before doing any QC, splitting etc. filter them out.

```bash
# create keep file from individuals
# duplicate column one to get the following format: IID FID HEIGHT
awk '{print $1, $1}' height.txt > /faststorage/project/populationgenomics/students/amlg/project/data/shit_phenotypes_kickedout/phenokeep.txt

# then use plink to filter the b-files
plink --bfile /faststorage/project/populationgenomics/project_data/GWAS/gwas_data --keep phenokeep.txt --make-bed --out pheno_filtered
```

* --keep: 1071 people remaining
*  At least 2 duplicate IDs in --keep file
*  1376653 variants and 1071 people pass filters and QC

# 1. Initial per-chip QC (split) [X]
- why: different chips = different SNPs, genotyping errors etc.

Python script ((find script)[https://github.com/antoniamlg/GT_GWASproject/blob/main/Scripts/splitbychip.ipynb]) to:
1. Read the metadata and split IIDs by chip.
2. Create `.keep` files for each chip group.

Use `.keep` files to filter bed file into chip specific files:
```bash
plink --bfile gwas_data --keep /faststorage/project/populationgenomics/students/amlg/project/data/splitbychip/chipX.keep --make-bed --out /faststorage/project/populationgenomics/students/amlg/project/data/splitbychip/chipX_data
```

This results in 6 different b-file sets. Each of them only containing data specific to the chip. I decided, to not use the all the datapoints with the "chip" label as they do not contain any useful data.

#### As we do not have that many individuals (now ``)
## 1.1 SNP-level QC [ ]

### 1.1.1 Call rate/Missingness [ ]
Run the --missing command again to generate the GWA-data.lmiss with the missing data rate for each SNP. <br>
Use R to make a histogram of the missing data rates (F_MISS).
Run the test-missing command and make a list of all the names of all SNPs where the differential missingness p-value is less than 1e-5. Save the list as fail-diffmiss-qc.txt.

```bash
plink --bfile splitbychip/unknown_chip_data --missing --out ../SNPmiss/unknown_chip_SNPmiss
ls -lht
```

See visualization script at ../Scripts/SNPmiss_QC.ipynb.

A lot of care must be taken when filtering variants, since we can lose potential variants associaated with the phenotypic trait. Therefore, we only filter these variants in which the missigness is associated with the phenotype, i.e case or control, so they could be a source of bias in our study. I compute the Fisher exact test between these two variables by running:

```bash
plink --bfile infile --test-missing --out outfile
```
<font color="pink">=> did not work because we don't have case/control data</font>

[X] Before continuing to the next step, do a plot to see where to set the threshold. Apparently, 0.05 removes almost everything which is annoying.

#### Results

* use `0.9` as a filtering threshold to first get rid of SNPs which are missing in all individuals (data cleaning)
* use `0.05` -> If a SNP has a high F_MISS (e.g., > 0.05), it means that a large proportion of individuals have missing data for that SNP.
* the `unknown_chip` does not have any data below 0.2% missing data rate. So the idea is, to just exclude the whole chip. <font color="pink">=> Does that make any sense though?</font>
=> sooo next step: [ ] do the actual filtering

### 1.1.2 Filter SNPs that are missing in 100% of individuals [ ]

```bash
plink --bfile chipX --geno 1 --make-bed --out chipX...
```
First, I will only filter SNPs that are missing in 100% of individuals in a  chip. Later (add section) I will also filter for higher missingness - but first, do sex imputation.

### 1.1.3 Merge chips
We need to run the sex-imputation on all the datasets, so we at least have a chance to impute sex for the `omics` chip which does not have any information about the sex. <br>
Sooooooooo, apparently you cannot really merge bfiles back together, if they don't have exactly the same SNPs

### 1.1.X Replot missingness to see if i actually did something.

### 1.1.4 Sex-imputation (TRY)
If it works - amazing.
If it doesn't - just continue with QC and don't care too much (-*Bjarke*).

### 1.1.5 MAF & HWE [ ]

after doing all this, remove all low-quality SNPs
```bash
plink --bfile GWA-QC --exclude fail-diffmiss-qc.txt --geno 0.05 --hwe 0.00001 --maf 0.01 --make-bed --out GWA-QC
```

In addition to removing SNPs identified with differential call rates between cases and controls, this command removes SNPs with call rate less than 95% with --geno option and deviation from HWE (p<1e-5) with the --hwe option. It also removes all SNPs with minor allele frequency less than a specified threshold using the --maf option.

## 1.2 Sample-level QC [ ]
At this point still confused what I am doing and how the order matters.

### 1.2.1 Sex checks [ ]

Don't use for *Omics*-chip - it does not have any info on sex

[ ] Is the number of individuals on the Omics chip == ambiguous individuals? Or do we have more individuals somewhere without any sex?

```bash
plink --bfile GWA-data --check-sex --out GWA-QC
```
Then, remove proplematic sex with
```bash
plink --bfile GWA-data --remove wrong_sex.txt --make-bed --out GWA-QC
```

### 1.2.2 Missingness [ ]

```bash
plink --bfile GWA-QC --missing --out GWA-QC
```

### 1.2.3 Heterozygosity outliers [ ]

```bash
plink --bfile GWA-QC --het --out GWA-QC 
```
calculate the observed heterozygosity rate per individual using the formula:

Het = (N(NM) − O(Hom))/N(NM) -> look at the exercises!

* then, do a plot of *F-values* and find outliers manually in Python/R

Make a file with the FID and IID of all individuals that have a genotype missing rate >=0.03 or a heterozygosity rate that is more than 3 s.d. from the mean. Then use plink to remove these individuals from the data set.
```bash
plink --bfile GWA-QC --remove wrong_het_missing_values.txt --make-bed --out GWA-QC
```

### 1.2.4 Relatedness / Duplicates [ ]

* what is this command? flags individuals with pi_hat > 0.185
```bash
plink --bfile chipX_qc4 --genome --min 0.185 --out chipX_related

```

1. prune the data
```bash
plink --bfile GWA-QC --indep-pairwise 500kb 5 0.2 --out GWA-QC
```
2. Calculate IBD between each pair of individuals
```bash
plink --bfile GWA-QC --extract GWA-QC.prune.in --genome --min 0.185 --out GWA-QC
```
3. Remove a member from each of the pairs that are too closely related from the data set. To keep it simple you can just always remove the individual mentioned first. *

4. Remove individuals with wrong IBD.
```bash
plink --bfile  GWA-QC --remove wrong_ibd.txt --make-bed --out GWA-QC
```
## 1.3 Population Structure/ Stratification/ Batch Effects

### 1.3.1 PCA
```bash
plink --bfile chipX_qc4 --pca 10 --out chipX_pca
```

###### After QC, merge chips again ######
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